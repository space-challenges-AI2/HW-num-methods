/* Author: Nikola D. Lilov, AI2; 06 August 2025
 *
 * Description: Implementation of the IV order Runge-Kutta method, simulating the motion of an orbiter around a planet.
 *
 * Normal usage:   RK4.exe
 * Expanded usage: RK4.exe -h
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

// physics parameters
const double GM = -3.98589156e5; // -<gravitational-constant>*<earth-mass> in km^3/sec^2

class RHSreturn
{
public:
    double x, y, z, vx, vy, vz;
};

// Return type of the program - vectors of coordinates and velocities at a given time
class RK4return
{
public:
    vector<double> x, y, z, vx, vy, vz;
    RK4return(double x0, double y0, double z0, double vx0, double vy0, double vz0, long long samples)
    {
        x.resize(samples);
        y.resize(samples);
        z.resize(samples);
        vx.resize(samples);
        vy.resize(samples);
        vz.resize(samples);
        x[0] = x0;
        y[0] = y0;
        z[0] = z0;
        vx[0] = vx0;
        vy[0] = vy0;
        vz[0] = vz0;
    }
};

// calculate derivatives
RHSreturn RHS(double x, double y, double z, double vx, double vy, double vz)
{
    // ax = -gamma*mass*x/(x^2 + y^2 + z^2)^(3/2)
    double divisor = (x * x + y * y + z * z) * sqrt(x * x + y * y + z * z);

    RHSreturn ans;
    ans.x = vx;
    ans.y = vy;
    ans.z = vz;
    ans.vx = GM * x / divisor;
    ans.vy = GM * y / divisor;
    ans.vz = GM * z / divisor;

    return ans;
}

void Euler(RK4return &ans, const int &samples, const double &dt)
{
    for (long long i = 1; i < samples; i++)
    {
        RHSreturn k1 = RHS(ans.x[i - 1], ans.y[i - 1], ans.z[i - 1], ans.vx[i - 1], ans.vy[i - 1], ans.vz[i - 1]);
        ans.x[i] = ans.x[i - 1] + k1.x * dt;
        ans.y[i] = ans.y[i - 1] + k1.y * dt;
        ans.z[i] = ans.z[i - 1] + k1.z * dt;
        ans.vx[i] = ans.vx[i - 1] + k1.vx * dt;
        ans.vy[i] = ans.vy[i - 1] + k1.vy * dt;
        ans.vz[i] = ans.vz[i - 1] + k1.vz * dt;
    }
}

void RK4(RK4return &ans, const int &samples, const double &dt)
{
    for (long long i = 1; i < samples; i++)
    {
        RHSreturn k1 = RHS(ans.x[i - 1], ans.y[i - 1], ans.z[i - 1], ans.vx[i - 1], ans.vy[i - 1], ans.vz[i - 1]);
        RHSreturn k2 = RHS(ans.x[i - 1] + 0.5 * k1.x * dt,
                           ans.y[i - 1] + 0.5 * k1.y * dt,
                           ans.z[i - 1] + 0.5 * k1.z * dt,
                           ans.vx[i - 1] + 0.5 * k1.vx * dt,
                           ans.vy[i - 1] + 0.5 * k1.vy * dt,
                           ans.vz[i - 1] + 0.5 * k1.vz * dt);
        RHSreturn k3 = RHS(ans.x[i - 1] + 0.5 * k2.x * dt,
                           ans.y[i - 1] + 0.5 * k2.y * dt,
                           ans.z[i - 1] + 0.5 * k2.z * dt,
                           ans.vx[i - 1] + 0.5 * k2.vx * dt,
                           ans.vy[i - 1] + 0.5 * k2.vy * dt,
                           ans.vz[i - 1] + 0.5 * k2.vz * dt);
        RHSreturn k4 = RHS(ans.x[i - 1] + k3.x * dt,
                           ans.y[i - 1] + k3.y * dt,
                           ans.z[i - 1] + k3.z * dt,
                           ans.vx[i - 1] + k3.vx * dt,
                           ans.vy[i - 1] + k3.vy * dt,
                           ans.vz[i - 1] + k3.vz * dt);

        // x[i+1]=x[i]+(k1.x+2*k2.x+2*k3.x+k4.x)*dt/6;
        ans.x[i] = ans.x[i - 1] + (k1.x + 2 * k2.x + 2 * k3.x + k4.x) * dt / 6;
        ans.y[i] = ans.y[i - 1] + (k1.y + 2 * k2.y + 2 * k3.y + k4.y) * dt / 6;
        ans.z[i] = ans.z[i - 1] + (k1.z + 2 * k2.z + 2 * k3.z + k4.z) * dt / 6;
        ans.vx[i] = ans.vx[i - 1] + (k1.vx + 2 * k2.vx + 2 * k3.vx + k4.vx) * dt / 6;
        ans.vy[i] = ans.vy[i - 1] + (k1.vy + 2 * k2.vy + 2 * k3.vy + k4.vy) * dt / 6;
        ans.vz[i] = ans.vz[i - 1] + (k1.vz + 2 * k2.vz + 2 * k3.vz + k4.vz) * dt / 6;
    }
}

// Function to write RK4return data to CSV file (to visualize in Plotter.py)
void writeToCSV(const RK4return &data, const string &filename, const long long &samples, const double &dt)
{
    ofstream file(filename);

    if (!file.is_open())
    {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return;
    }

    file << "Time (UTCG),x (km),y (km),z (km),vx (km/sec),vy (km/sec),vz (km/sec)\n";

    file << fixed << setprecision(6);
    for (long long i = 0; i < samples; i++)
    {
        file << i * dt << ','
             << data.x[i] << ','
             << data.y[i] << ','
             << data.z[i] << ','
             << data.vx[i] << ','
             << data.vy[i] << ','
             << data.vz[i] << '\n';
    }
    file.close();
    cout << "Data written to " << filename << " successfully!\n";
}

void error(char *&name)
{
    cerr << "Normal Usage:  " << name << "\n\nExpanded usage:\n"
         << "--target <filename.csv>  -  Name of the output CSV file (default is output.csv).\n"
         << "--samples <Number_of_samples>  -  Number of samples to be generated (default is 8641).\n"
         << "--method <RK4/Euler>  -  Numerical method to use (default is RK4).\n"
         << "-ts <time_start> -te <time_end>  -  Start and end times of the simulation in [s] (default is 0.00 and 86400.00).\n"
         << "-x <x0> -y <y0> -z <z0> -vx <vx0> -vy <vy0> -vz <vz0>  -  Initial conditions for the orbiter in [km] and [km/s] (default is 7000.00, 0.000001, -0.001608, 0.000002, 1.310359, 7.431412).\n"
         << "--help or -h  -  Brings up the description of the program and this message.\n";
}

int main(int argc, char *argv[])
{
    // orbiter parameters (normalized):
    double time = 0.00, x = 7000.00, y = 0.000001, z = -0.001608, vx = 0.000002, vy = 1.310359, vz = 7.431412;
    // simulation parameters:
    double t_end = 86400.00;
    long long samples = 8641;
    string filename = "output.csv", method = "RK4";

    // CMD arguments
    for (int i = 1; i < argc; i++)
    {
        string arg = argv[i];
        if (arg == "--samples")
            samples = stoll(argv[++i]);
        else if (arg == "--target")
            filename = argv[++i];
        else if (arg == "--method")
            method = argv[++i];
        else if (arg == "-ts")
            time = stod(argv[++i]);
        else if (arg == "-te")
            t_end = stod(argv[++i]);
        else if (arg == "-x")
            x = stod(argv[++i]);
        else if (arg == "-y")
            y = stod(argv[++i]);
        else if (arg == "-z")
            z = stod(argv[++i]);
        else if (arg == "-vx")
            vx = stod(argv[++i]);
        else if (arg == "-vy")
            vy = stod(argv[++i]);
        else if (arg == "-vz")
            vz = stod(argv[++i]);
        else
        {
            if (arg == "--help" || arg == "-h")
                cerr << "\nThis is a program which takes the initial conditions of an orbiter and simulates its motion around Earth, giving the results in a CSV table of the sort <filename.csv>: Time [s], Position (x,y,z) [km], Velocity (vx,vy,vz) [km/s].\nAuthor: Nikola D. Lilov; 06 August 2025.\n\n";
            else
                cerr << "Unknown parameter: " << arg << '\n';
            error(argv[0]);
            return 1;
        }
    }

    double dt = (t_end - time) / (samples - 1);
    RK4return ans(x, y, z, vx, vy, vz, samples);

    if (method == "RK4")
        RK4(ans, samples, dt);
    else if (method == "Euler")
        Euler(ans, samples, dt);
    else
    {
        cout << "Invalid method specified. Use either 'RK4' (default) or 'Euler'.\n";
        return 1;
    }

    writeToCSV(ans, filename, samples, dt);
    return 0;
}