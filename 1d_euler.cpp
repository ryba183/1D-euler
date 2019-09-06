#include <iostream>
#include <complex>
#include <array>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <getopt.h>
#include <fstream>

unsigned int ncells = 1000;
double maxtime = 0.15;
double densityL = 1;
double densityR = 1;
double velocityL = -2;
double velocityR = 2;
double pressureL = 0.4;
double pressureR = 0.4;
std::string outPath = "-";

void PrintHelp() {
    std::cout <<
        "Parameter settings:\n"
        "   -n, --steps           Set number of steps across domain\n"
        "   -t, --time            Set time to simulate until\n"
        "   -d, --density_left    Set density on left side of domain\n"
        "   -D, --density_right   Set density on right side of domain\n"
        "   -v, --velocity_left   Set velocity on left side of domain\n"
        "   -V, --velocity_right  Set velocity on right side of domain\n"
        "   -p, --pressure_left   Set pressure on left side of domain\n"
        "   -P, --pressure_right  Set pressure on right side of domain\n"
        "   -o, --outfile         File to write solution\n"
        "   -h, --help            Show help\n";
    exit(1);
}

// https://gist.github.com/ashwin/d88184923c7161d368a9
void ProcessArgs(int argc, char** argv) {
    const char* const short_opts = "n:t:d:D:v:V:p:P:o:h";
    const option long_opts[] = {
        {"steps", required_argument, nullptr, 'n'},
        {"time", required_argument, nullptr, 't'},
        {"density_left", required_argument, nullptr, 'd'},
        {"density_right", required_argument, nullptr, 'D'},
        {"velocity_left", required_argument, nullptr, 'v'},
        {"velocity_right", required_argument, nullptr, 'V'},
        {"pressure_left", required_argument, nullptr, 'p'},
        {"pressure_right", required_argument, nullptr, 'P'},
        {"outfile", required_argument, nullptr, 'o'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'n':
            ncells = std::stoi(optarg);
            break;
        case 't':
           maxtime = std::stof(optarg);
           break;
        case 'd':
            densityL = std::stof(optarg);
            break;
        case 'D':
            densityR = std::stof(optarg);
            break;
        case 'v':
            velocityL = std::stof(optarg);
            break;
        case 'V':
            velocityR = std::stof(optarg);
            break;
        case 'p':
            pressureL = std::stof(optarg);
            break;
        case 'P':
            pressureR = std::stof(optarg);
            break;
        case 'o':
            outPath = std::string(optarg);
            break;
        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
}

double get_energy(  double pressure, double density, double momentum,
                    double velocity, double gamma = 1.4) {
    double internal_energy = pressure/((gamma - 1)*density);
    return (density * internal_energy) + (0.5*momentum*velocity);
}

double get_velocity(double density, double momentum) {
    return momentum / density;
}

double get_internal_energy(double density, double energy, double velocity) {
    return (energy / density) - (0.5 * velocity * velocity);
}

double get_pressure(
        double density, double internal_energy, double gamma = 1.4) {
    return (gamma - 1) * (density * internal_energy);
}

double get_dx(
        unsigned int ncells, double domain_start = 0, double domain_end = 1) {
    return (domain_end - domain_start)/ncells;
}

double speed_of_sound(std::array<double, 6>& q, double gamma = 1.4) {
    double c_squared = (gamma*q[3])/q[0];             
    return sqrt(c_squared);
}

double get_amax(std::vector<std::array<double, 6>>& q) {
    int size = q.size();
    double c;
    std::vector<double> all_a;
    for (int i = 0; i < size; i++) {
        c = speed_of_sound(q[i]);
        all_a.push_back(std::abs(q[i][4]) + c);
    } 
    return *std::max_element(all_a.begin(), all_a.end());
}

double get_dt(double amax, double dx, double CFL = 0.9) {
    return CFL*(dx/amax);
}

double f_density(std::array<double, 6>& q) {
    return q[1];
}

double f_momentum(std::array<double, 6>& q) {
    return (q[1]*q[4]) + q[3];
}

double f_energy(std::array<double, 6>& q) {
    return (q[2] + q[3])*q[4];
}

std::array<double, 6> get_q_ihalf(std::array<double, 6>& qi , 
                                  std::array<double, 6>& q_i1, 
                                  double dx, double dt) {
                          
    // Declare output vector
    std::array<double, 6> q_out;
    
    double halfdtdx = 0.5 * (dt/dx);
    
    q_out[0] = (0.5 * (qi[0] + q_i1[0])) 
               + (halfdtdx * (f_density(qi) - f_density(q_i1)));
    q_out[1] = (0.5 * (qi[1] + q_i1[1])) 
               + (halfdtdx * (f_momentum(qi) - f_momentum(q_i1)));
    q_out[2] = (0.5 * (qi[2] + q_i1[2])) 
               + (halfdtdx * (f_energy(qi) - f_energy(q_i1)));
    q_out[4] = get_velocity(q_out[0], q_out[1]);
    q_out[5] = get_internal_energy(q_out[0], q_out[2], q_out[4]); 
    q_out[3] = get_pressure(q_out[0], q_out[5]);
    
    return q_out;
}

std::array<double, 6>  get_flux_RI(std::array<double, 6>& q) {
    // Declare output vector
    std::array<double, 6> flux_RI;
  
    flux_RI[0] = f_density(q);
    flux_RI[1] = f_momentum(q);
    flux_RI[2] = f_energy(q);
    
    return flux_RI;
}

std::array<double, 6> get_flux_LF(std::array<double, 6>& qi , 
                                    std::array<double, 6>& q_i1, 
                                    double dx, double dt) {
                          
    // Declare output vector
    std::array<double, 6> flux_LF;
    
    double halfdxdt = 0.5 * (dx/dt);
    
    flux_LF[0] = (0.5 * (f_density(qi) + f_density(q_i1))) 
                 + (halfdxdt * (qi[0] - q_i1[0]));
    flux_LF[1] = (0.5 * (f_momentum(qi) + f_momentum(q_i1))) 
                 + (halfdxdt * (qi[1] - q_i1[1]));
    flux_LF[2] = (0.5 * (f_energy(qi) + f_energy(q_i1))) 
                 + (halfdxdt * (qi[2] - q_i1[2]));
    
    return flux_LF;
}

std::array<double, 6> get_force(std::array<double, 6>& qi, 
                                std::array<double, 6>& q_i1,
                                double dx, double dt) {
                          
    std::array<double, 6>  q_ihalf = get_q_ihalf(qi, q_i1, dx, dt);
    std::array<double, 6>  flux_RI = get_flux_RI(q_ihalf);
    std::array<double, 6>  flux_LF = get_flux_LF(qi, q_i1, dx, dt);
    
    // Declare output vector
    std::array<double, 6> force_flux;
    
    force_flux[0] = 0.5 * (flux_LF[0] + flux_RI[0]);
    force_flux[1] = 0.5 * (flux_LF[1] + flux_RI[1]);
    force_flux[2] = 0.5 * (flux_LF[2] + flux_RI[2]);
    
    return force_flux;
}

std::array<double, 6> get_q_n1(std::array<double, 6>& qi, 
                               std::array<double, 6>& q_iplus1,
                               std::array<double, 6>& q_iless1,
                               double dx, double dt) {
    
    std::array<double, 6> force_plus_half = get_force(qi, q_iplus1, dx, dt);
    std::array<double, 6> force_less_half = get_force(q_iless1, qi, dx, dt);

    // Declare output vector
    std::array<double, 6> q_nplus1;
    
    double dtdx = dt/dx;
    
    q_nplus1[0] = qi[0] + dtdx*(force_less_half[0] - force_plus_half[0]);
    q_nplus1[1] = qi[1] + dtdx*(force_less_half[1] - force_plus_half[1]);
    q_nplus1[2] = qi[2] + dtdx*(force_less_half[2] - force_plus_half[2]);
    q_nplus1[4] = get_velocity(q_nplus1[0], q_nplus1[1]);
    q_nplus1[5] = get_internal_energy(q_nplus1[0], q_nplus1[2], q_nplus1[4]); 
    q_nplus1[3] = get_pressure(q_nplus1[0], q_nplus1[5]);
    
    return q_nplus1;
}

std::array<std::array<double, 6>, 2> get_q_isurround(
    std::vector<std::array<double, 6>>& q, int i) {
    // Returns q[i-1] and q[i+1] accounting for boundary conditions
    
    // Declare output map
    std::array<std::array<double, 6>, 2> q_isurround;

    int size = q.size();
    if (i == 0) {
        q_isurround[0] = q[i];
    } else {
        q_isurround[0] = q[i-1];   
    }
    if (i == size - 1) {
        q_isurround[1] = q[size - 1];
    } else {
        q_isurround[1] = q[i+1];   
    }
    return q_isurround;
}


std::vector<std::array<double, 6>> initialiseData(unsigned int ncells,
    double densityL, double velocityL, double pressureL,
    double densityR, double velocityR, double pressureR) {
    
    // Declare initial data
    std::vector<std::array<double, 6>> q(ncells);

    double energy, density, velocity, pressure, momentum;
    // In odd number of cells, the middle cell is assigned to right domain
    unsigned int midpoint = ncells/2;

    for (unsigned int i = 0; i < ncells; i++) {
        if (i <= midpoint) {
            density = densityL;
            velocity = velocityL;
            pressure = pressureL;
        }
        else {
            density = densityR;
            velocity = velocityR;
            pressure = pressureR;
        }
        momentum = density * velocity;
        energy = get_energy(pressure, density, momentum, velocity);
        q[i][0] = density;
        q[i][1] = momentum;
        q[i][2] = energy;
        q[i][3] = pressure;
        q[i][4] = velocity;
        q[i][5] = get_internal_energy(density, momentum, energy);
    }
    return q;
}
      

int main(int argc, char* argv[]) {
    
    // Extract command line argumenets
    ProcessArgs(argc, argv);
    std::ofstream outFile;
    
    // Set output to file path or stdout (default)
    if (outPath != "-") {
        outFile.open(outPath, std::ios::out);
    } 
    std::ostream& out = (outPath != "-" ? outFile : std::cout);
    
    // Initialise data
    std::vector<std::array<double, 6>> q = initialiseData(
        ncells, densityL, velocityL, pressureL, 
        densityR, velocityR, pressureR);
    
    // Compute delta x across domain space
    double dx = get_dx(ncells);
    
    double t = 0;
    while (t < maxtime) {
        
        // Declare vector for storing q(n+1) for all i
        std::vector<std::array<double, 6>> q_next(ncells);
        
        // Compute dt for q[i]
        double amax = get_amax(q);
        double dt = get_dt(amax, dx);
        
        for (unsigned int i = 0; i < q.size(); i++) {
            // Return q(n, i+1) and q(n, i-1)
            std::array<std::array<double, 6>, 2> q_surround = get_q_isurround(q, i);
            // Compute q(i, n+1)
            q_next[i] = get_q_n1(q[i], q_surround[1], q_surround[0], dx, dt);
        }
        t += dt;
        q = q_next;
    }
    
    int size = q.size();
    for (int i = 0; i < size; i++) {
        out << q[i][0] << "\t" 
            << q[i][4] << "\t" 
            << q[i][3] << "\t" 
            << q[i][5] << std::endl;
    }
    return 0;
}
