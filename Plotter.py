''' 
    Author: Author: Nikola D. Lilov, AI2; 06 August 2025

    Description: Plots the data from output.csv (from RK4.exe).

    Usage: python Plotter.py
'''

import numpy as np
import csv
import matplotlib.pyplot as plt


def parse_orbit_data(filename):
    """
    Parse a .csv file with satellite time, position and velocity (PVT) for plotting
    Expected CSV format: "Time (UTCG),x (km),y (km),z (km),vx (km/sec),vy (km/sec),vz (km/sec)"

    Returns:
    - time[Npoints]: time in [sec] from simulation start
    - pos[3][Npoints]: 3D position in [km]
    - vel[3][Npoints]: 3D velocity in [km/s]
    """
    try:
        with open(filename, "r") as fp:
            data = csv.reader(fp)
            lst = list(data)

        ndata = len(lst) - 1  # Subtract header row
        time = np.zeros(ndata)
        pos = np.zeros((3, ndata))
        vel = np.zeros((3, ndata))

        for i in range(ndata):
            time[i] = float(lst[i + 1][0])
            for j in range(3):
                pos[j][i] = float(lst[i + 1][j + 1])  # x, y, z
                vel[j][i] = float(lst[i + 1][j + 4])  # vx, vy, vz

        return time, pos, vel

    except FileNotFoundError:
        print(f"Error: File {filename} not found")
        return None, None, None
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None, None, None


def main():
    # Define the two CSV files
    # STKfile = "Satellite_PVT_GMAT.csv"
    RK4file = "output.csv"

    # print("Loading STK reference data...")
    # time_stk, pos_stk, vel_stk = parse_orbit_data(STKfile)

    print("Loading RK4 simulation data...")
    time_rk4, pos_rk4, vel_rk4 = parse_orbit_data(RK4file)

    # Check if data was loaded successfully
    # if time_stk is None:
    #    print(f"Failed to load STK file: {STKfile}")
    #    return
    if time_rk4 is None:
        print(f"Failed to load RK4 file: {RK4file}")
        return

    # Print data info
    # print(f"STK data points: {len(time_stk)}")
    print(f"RK4 data points: {len(time_rk4)}")

    # 3D Trajectory Plot
    fig1 = plt.figure(figsize=(12, 9))
    ax = fig1.add_subplot(111, projection='3d')

    # ax.plot3D(pos_stk[0], pos_stk[1], pos_stk[2],
    #          color="red", label="STK Reference", linewidth=2)
    ax.plot3D(pos_rk4[0], pos_rk4[1], pos_rk4[2],
              color="green", label="RK4 Solution", linewidth=2, alpha=0.5)

    # --- Earth ---
    R = 6371.0  # Earth radius in km

    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 25)

    earth_x = R * np.outer(np.cos(u), np.sin(v))
    earth_y = R * np.outer(np.sin(u), np.sin(v))
    earth_z = R * np.outer(np.ones_like(u), np.cos(v))

    ax.plot_surface(
        earth_x,
        earth_y,
        earth_z,
        alpha=0.6
    )

    max_extent = max(
        np.max(np.abs(pos_rk4[0])),
        np.max(np.abs(pos_rk4[1])),
        np.max(np.abs(pos_rk4[2])),
        R
    ) * 1.1  # 10% padding

    ax.set_xlim(-max_extent, max_extent)
    ax.set_ylim(-max_extent, max_extent)
    ax.set_zlim(-max_extent, max_extent)

    ax.set_box_aspect([1, 1, 1])  # cube, not stretched per-axis

    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('Z [km]')
    ax.legend()
    ax.set_title('3D Orbital Trajectory Comparison')
    plt.tight_layout()
    plt.show()

    # 2D Position vs Time Plot
    fig2, ax = plt.subplots(figsize=(14, 8))

    # STK data (reference)
    # ax.plot(time_stk, pos_stk[0], color="darkgreen", linestyle="-",
    #        label="STK pos X", linewidth=2)
    # ax.plot(time_stk, pos_stk[1], color="darkred", linestyle="-",
    #        label="STK pos Y", linewidth=2)
    # ax.plot(time_stk, pos_stk[2], color="darkblue", linestyle="-",
    #        label="STK pos Z", linewidth=2)

    # RK4 data
    ax.plot(time_rk4, pos_rk4[0], color="lightgreen", linestyle="--",
            label="RK4 pos X", linewidth=2, alpha=0.8)
    ax.plot(time_rk4, pos_rk4[1], color="lightcoral", linestyle="--",
            label="RK4 pos Y", linewidth=2, alpha=0.8)
    ax.plot(time_rk4, pos_rk4[2], color="lightblue", linestyle="--",
            label="RK4 pos Z", linewidth=2, alpha=0.8)

    ax.set_xlabel("Time [sec]")
    ax.set_ylabel("Position [km]")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_title('Position Components vs Time')
    plt.tight_layout()
    plt.show()

    # Calculate and plot position error
    # Find common time range for error calculation
    """min_len = min(len(time_stk), len(time_rk4))

    # Ensure we're comparing at the same time points (assuming same sampling)
    pos_err = np.zeros(min_len)
    for i in range(min_len):
        pos_err[i] = np.sqrt(
            (pos_rk4[0][i] - pos_stk[0][i])**2 +
            (pos_rk4[1][i] - pos_stk[1][i])**2 +
            (pos_rk4[2][i] - pos_stk[2][i])**2
        )

    # Position Error Plot
    fig3, ax = plt.subplots(figsize=(12, 6))
    ax.plot(time_stk[:min_len], pos_err, color="red", linewidth=2,
            label="Position Error")
    ax.set_xlabel("Time [sec]")
    ax.set_ylabel("Position Error [km]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('RK4 Position Error vs STK Reference')
    plt.tight_layout()
    plt.show()

    # Print some statistics
    print(f"\nError Statistics:")
    print(f"Maximum error: {np.max(pos_err):.6f} km")
    print(f"Mean error: {np.mean(pos_err):.6f} km")
    print(f"RMS error: {np.sqrt(np.mean(pos_err**2)):.6f} km")
"""


if __name__ == "__main__":
    main()
