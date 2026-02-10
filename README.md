# N-Body Simulation

### Simulation main file
`./graphics/galsim.c`


### How to compile and rum

```bash
cd graphics/
make clean
make
time ./galsim 3000 ../input_data/ellipse_N_03000.gal 100 0.00001 0
```
### Automated test (compare with all the reference data)
```bash
./test.sh
```
### Program Input

The program `galsim` accepts **five input arguments**:
| Argument   | Description |
|------------|-------------|
| `N`        | Number of stars/particles to simulate (integer) |
| `filename` | The filename of the file to read the initial particle configuration from (string) |
| `nsteps`   | Number of timesteps to simulate (integer) |
| `delta_t`  | The timestep Δt for each simulation step (double/float) |
| `graphics` | 1 to enable graphics visualization, 0 to disable graphics |

### Example Usage

```bash
./galsim 3000 ../input_data/ellipse_N_03000.gal 100 0.00001 0
```
