# Orbiter Propagator (RK4 & Euler) — AI Team 2

Bootcamp homework: a lightweight propagator simulating the motion of an orbiter around Earth, built to mimic what STK does, using our own numerical integration instead of a commercial tool.

## What it is

- **RK4.cpp** — C++ implementation of a 2-body orbital propagator using 4th-order Runge-Kutta (default) or Euler integration. Outputs a CSV of time, position (x, y, z), and velocity (vx, vy, vz).
- **Plotter.py** — Reads the CSV and plots the 3D trajectory around Earth plus position vs. time.

## How to use

1. Compile the C++ propagator:
   ```
   g++ RK4.cpp -o RK4.exe
   ```
2. Run it to generate `output.csv`:
   ```
   ./RK4.exe
   ```
   See all options (method, initial conditions, sample count, time range, output filename):
   ```
   ./RK4.exe --help
   ```
3. Plot the results:
   ```
   python Plotter.py
   ```

## Notes

- Default initial conditions and 1-day simulation window are built in; override any of them via CLI flags.
- `Plotter.py` expects `output.csv` in the same directory.
