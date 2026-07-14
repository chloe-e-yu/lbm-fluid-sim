// Phase 1: D2Q9 Engine

// D2Q9 lattice constants

// Direction layout:
    // 6 2 5
    // 3 0 1        (0 = rest particle)
    // 7 4 8 
// opposite pairs -- bounce back
    // 1-3, 2-4, 5-7, 6-8

const NX = 300; // grid width
const NY = 150; // grid height

// x and y components of the 9 discrete velocity (one lattice cell per time step)
const cx = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const cy = [0, 0, 1, 0, -1, 1, 1, -1, -1];

// lattice weights - must be sum to 1 (distribution of a fluid at rest)
// from isotropy conditions: rest(4/9), axis directions(1/9), diagonal directions(1/36)
const w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
// verification
console.log("weights sum:", w.reduce((a, b) => a + b)); 

// opposite direction for bounce back 1-3 2-4 5-7  6-8
const opp = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// Distribution function f + Initialization
    // every lattice point (x,y) stores 9 numbers
    // f[x][y][q] = the amount of particles at point (x,y) in direction q
/* 
 entire simulaton is just updating NX * NY * 9 numbers every time step:
    collision -- mixes the numbers at each point
    streaming -- moves them to neighboring points       
*/

// flat Float32Array instead of nested arrays: contigous memory for cache locality, no pointer chasing
let f = new Float32Array(NX * NY * 9);

// flat index for (x,y,q): each cell stores its 9 directions contiguosly
function idx(x, y, q) {
    return (y * NX + x) * 9 + q;
}

// initialize to equilibrium at res: rho = 1, u = 0 => f_q = w_q
for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
        for (let q = 0; q < 9; q++) {
            f[idx(x, y, q)] = w[q];
        }
    }
}

// verification
    // density at a point = sum of its 9 values
    // total mass = sum of the whole array
    // since every cell starts at density 1 => the toal must equal NX x NY
let totalMass = 0;
for (let i = 0; i < f.length; i++) {
    totalMass += f[i];
}
console.log("Total mass at initialization: " + totalMass + " (expected: " + (NX * NY) + ")");

// macroscopic density and velocity at cell (x,y)
// rho = sum_q f_q; u = (1/rho) * sum_q f_q * c_q
function macroscopic(x,y) {
    let rho = 0, ux = 0, uy = 0;
    for (let q = 0; q < 9; q++) {
        const f_q = f[idx(x, y, q)];
        rho += f_q;
        ux += f_q * cx[q];
        uy += f_q * cy[q];
    }
    ux /= rho;
    uy /= rho;
    return {rho, ux, uy};
}

// equilibrium distribution for firection q given rho, ux, uy
// feq_q = w_q * rho (1 + 3(c_q . u) + 4.5 (c_q . u)^2 - 1.5 (u . u))
function equilibrium(q, rho, ux, uy) {
    const cu = cx[q] * ux + cy[q] * uy; // c_q . u
    const u_sq = ux * ux + uy * uy; // |u|^2
    return w[q] * rho * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * u_sq);
}

// verification
console.log(equilibrium(1,1,0,0), w[1]); //should 0.11111...



// Phase 2: Collision and Streaming
const tau = 0.6; // relaxation time; viscosity nu = (tau - 0.5) / 3
const omega = 1 / tau; // relaxation rate

// BGK collision: relax every f toward local equilibrium
function collide() {
    for (let y = 0; y < NY; y++) {
        for (let x = 0; x < NX; x++) {
            const {rho, ux, uy} = macroscopic(x,y);
            for (let q = 0; q < 9; q++) {
                const feq_q = equilibrium(q, rho, ux, uy);
                const i = idx(x,y,q);
                f[i] = f[i] - omega * (f[i] - feq_q);
                // equivalent: f[i] = (1 - omega)* f[i] + omega * feq_q;
            }
        }
    }
}

// total mass before and after collide() must match
let m0 = 0; for (let i = 0; i < f.length; i++) m0 += f[i];
collide();
let m1 = 0; for (let i = 0; i < f.length; i++) m1 += f[i];
console.log ("mass befor:", m0, "after:", m1);

let fNew = new Float32Array(NX * NY * 9); // temporary array for streaming

// streaming: f values move one cell along their direction
function stream() {
    for (let y = 0; y < NY; y++) {
        for (let x = 0; x < NX; x++) {
            for (let q = 0; q < 9; q++) {
                const xdest = x + cx[q];
                const ydest = y + cy[q];
                // skip if destination is outside the grid (walls come next step)
                if (xdest < 0 || xdest >= NX || ydest < 0 || ydest >= NY) { 
                fNew[idx(x,y,opp[q])] = f[idx(x, y, q)];
                } else {
                    fNew[idx(xdest,ydest,q)] = f[idx(x, y, q)];
                }
            }
        }
    }    
    const temp = f; f = fNew; fNew = temp; // swap arrays    
}
// sealed-box test: with wallas on all sodes, mass must stay constant
f[idx(50, 50, 1)] += 0.5; // disturb the efluid

let mStart = 0; for (let i = 0; i < f.length; i++) mStart += f[i];
for (let n = 0; n < 100; n++) step();
let mEnd = 0; for (let i = 0; i < f.length; i++) mEnd += f[i];
console.log("sealed box, 100 steps: ", mStart, "->", mEnd);

function step() {
    collide();
    stream();
}

