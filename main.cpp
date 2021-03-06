#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <fstream>


/* PHYSICAL CONSTANTS */

constexpr double e = 1.60217657e-19;   // Coulombs
constexpr double k_B = 1.3806488e-23;  // Joules per Kelvin


/* USER CONFIGURATION */

// Note: If the simulation halts, try reducing the Runge-Kutta evolution step
// size, listed below. A starting point is 1e-1 but for complex systems this
// may need to be reduced as far as 1e-2, or perhaps further. Larger step size
// means a faster simulation.

/* Semiclassical simulation properties */
constexpr double h = 0.3e-1;  // Runge-Kutta evolution step size
constexpr double convergence_criterion = 1e-5;

/* Temperature */
constexpr double temp = 4.2;  // Kelvin

/* Dot orbital energies */
constexpr double dot_lumo = -5.1*e;  // Joules
constexpr double single_electron_energies[] = {
    dot_lumo + 0.000*e, dot_lumo + 0.000*e,
    dot_lumo + 0.231*e, dot_lumo + 0.231*e,
    dot_lumo + 0.231*e, dot_lumo + 0.231*e,
    dot_lumo + 0.479*e, dot_lumo + 0.479*e,

};  // Joules

/* Electronic properties of leads (s: source; d: drain) */
constexpr double d_fermi_energy = -5.1*e;  // Joules
double source_dos (double energy) {
    // Density of states can be an arbitrary function of energy.
    if (energy > 0.0) {}  // This suppresses an "unused variable" g++ warning.
    return 1;
}

constexpr double s_fermi_energy = -5.1*e;  // Joules
double drain_dos (double energy) {
    if (energy > 0.0) {}
    return 1;
}

/* Tunnel widths (by dot level) in arbitrary units */
constexpr double source_widths[] {
    3.0, 3.0,
    6.0, 6.0,
    6.0, 6.0,
    12.0, 12.0,
};

constexpr double drain_widths[] {
    0.5, 0.5,
    1.0, 1.0,
    1.0, 1.0,
    2.0, 2.0,
};

/* Dot-system capacitances */
constexpr double gate_capacitance =   1e-20;  // Farads
constexpr double source_capacitance = 1.6e-18;  // Farads
constexpr double drain_capacitance =  1.6e-18;  // Farads
constexpr double extra_capacitance =  0.1e-18;  // Farads

/* Voltage-space to be explored */
constexpr double v_g_min = -50;  // Volts
constexpr double v_g_max = 200;  // Volts
constexpr int v_g_steps = 200;  // (y axis resolution)

constexpr double v_sd_min = -1.1;  // Volts
constexpr double v_sd_max = 1.1;  // Volts
constexpr int v_sd_steps = 200;  // (x axis resolution)

/* (end of user configuration) */


/* COMPILE-TIME CALCULATIONS */

constexpr int n_levels = sizeof(single_electron_energies)
                         /sizeof(single_electron_energies[0]);

constexpr int n_source_rates = sizeof(source_widths)
                               /sizeof(source_widths[0]);

constexpr int n_drain_rates = sizeof(drain_widths)
                               /sizeof(drain_widths[0]);

static_assert(n_levels == n_source_rates && n_levels == n_drain_rates,
              "Wrong number of tunnelling rates inputted by user.");

constexpr double total_capacitance = gate_capacitance + source_capacitance
                                     + drain_capacitance + extra_capacitance;


/* TYPES */

typedef struct {
    double gate;
    double sd;  // Source-drain
} v_pair;

typedef unsigned long cfg;
constexpr cfg n_configs = 1 << n_levels;  // 2^n_levels

typedef struct {
    cfg to;    // Row
    cfg from;  // Column
    double value;
} matrix_elem;

typedef double mu_spectrum[n_levels];


/* BINARY REPRESENTATION OF CONFIGURATIONS */

bool occupied (cfg config, int level) {
    return (config >> level) & 1;
}

int sum_occupation (cfg config) {
    int sum = 0;
    for (int level = 0; level < n_levels; ++level)
        sum += (int)occupied(config, level);
    return sum;
}

cfg flipped_occupation (cfg config, int level) {
    return config ^ (1 << level);
}


/* THE PHYSICS */

/* Thermodynamics */

double fermi (double energy, double chem_pot) {
    double exponent = (energy - chem_pot)/(k_B*temp);
    if (exponent > log(std::numeric_limits<double>::max()))
        return 0.0;
    else if (exponent < log(std::numeric_limits<double>::min()))
        return 1.0;
    else
        return 1.0/(exp(exponent) + 1);
}

/* Energy shifts in leads owing to applied bias */

constexpr double s_shift (double energy, double v_sd) {
    return energy - e*v_sd/2.0;
}

constexpr double d_shift (double energy, double v_sd) {
    return energy + e*v_sd/2.0;
}

/* Tunnel rates */

double in_from_source (double mu, double v_sd, int level) {
    return source_widths[level]
           * source_dos(mu - s_shift(0.0, v_sd))
           * fermi(mu, s_shift(s_fermi_energy, v_sd));
}

double in_from_drain (double mu, double v_sd, int level) {
    return drain_widths[level]
           * drain_dos(mu - d_shift(0.0, v_sd))
           * fermi(mu, d_shift(d_fermi_energy, v_sd));
}

double out_to_source (double mu, double v_sd, int level) {
    return source_widths[level]
           * source_dos(mu - s_shift(0.0, v_sd))
           * (1 - fermi(mu, s_shift(s_fermi_energy, v_sd)));
}

double out_to_drain (double mu, double v_sd, int level) {
    return drain_widths[level]
           * drain_dos(mu - d_shift(0.0, v_sd))
           * (1 - fermi(mu, d_shift(d_fermi_energy, v_sd)));
}

double current_through_level (double mu, bool occupied, double v_sd, int level)
{
    if (occupied)
        return -e*(out_to_drain(mu, v_sd, level)
                   - out_to_source(mu, v_sd, level));
    else
        return -e*(in_from_source(mu, v_sd, level)
                   - in_from_drain(mu, v_sd, level));
}

/* Charging energy (constant interaction) */

double chemical_potential (bool occupied, double single_electron_energy,
                           int n_electrons_on_dot, v_pair voltage) {

    return single_electron_energy
           + ((e/total_capacitance)
             *((n_electrons_on_dot + (occupied ? -1 : 1)*0.5)*e
                - gate_capacitance*voltage.gate
                - (source_capacitance - drain_capacitance)*voltage.sd/2.0));
}

/* Rate matrix elements (Beenakker rate equations) */

matrix_elem diag (cfg config, mu_spectrum spectrum, double v_sd) {
    double value = 0.0;
    for (int level = 0; level < n_levels; ++level) {
        auto mu = spectrum[level];
        if (occupied(config, level))
            value -= (out_to_source(mu, v_sd, level)
                      + out_to_drain(mu, v_sd, level));
        else
            value -= (in_from_source(mu, v_sd, level)
                      + in_from_drain(mu, v_sd, level));
    }
    return { config, config, value };
}

matrix_elem offdiag (cfg to, cfg from, int level, double mu, double v_sd) {
    double value;
    if (occupied(from, level))
        value = out_to_source(mu, v_sd, level)
                + out_to_drain(mu, v_sd, level);
    else
        value = in_from_source(mu, v_sd, level)
                + in_from_drain(mu, v_sd, level);
    return { to, from, value };
}


/* MISCELLANEOUS HELPER */

v_pair voltage_pair_from_index (int index) {
    /* 
     * This function is defined such that we iterate through voltage-space
     * like so:
     *
     *        ┌┐┌┐┌┐┌┐┌finish
     *        │││││││││
     *        │││││││││
     *   start┘└┘└┘└┘└┘
     *
     * where y: gate voltage; x: source-drain voltage. In this way, we ensure
     * our weights guesses will mostly be almost correct.
     */
    v_pair voltage;

    if ((int)floor((double)index/(double)v_g_steps) % 2 == 0)
        voltage.gate = v_g_min
                       + ((index % v_g_steps)
                          *(v_g_max - v_g_min)/(double)v_g_steps);
    else
        voltage.gate = v_g_max
                       - (((index % v_g_steps) + 1)
                          *(v_g_max - v_g_min)/(double)v_g_steps);

    voltage.sd = v_sd_min
                   + (floor((double)index/(double)v_g_steps)
                      *(v_sd_max - v_sd_min)/(double)(v_sd_steps));

    return voltage;
}


/* THE SIMULATION */

int main () {

    std::ofstream outfile ("output.csv", std::ofstream::out);
    outfile << n_levels << "\n";  // Needed for visualization later.

    /* 'Guess' that the dot is initially empty; w = { 1, 0, 0, ... }. */
    std::vector<double> guess (n_configs);
    guess[0] = 1;

    /* Iterate through points in voltage-space. */
    for (int voltage_index = 0;
         voltage_index < v_g_steps*v_sd_steps;
         ++voltage_index) {

        /* Choose a point in voltage-space. */
        auto voltage = voltage_pair_from_index(voltage_index);
        outfile << voltage.gate << " ; " << voltage.sd;
        std::cout << "v_g:" << voltage.gate << " ; "
                  << "v_sd:" << voltage.sd;

        /* For each possible config, find all the chemical potentials. */
        std::vector<mu_spectrum> mu;
        mu.reserve(n_configs);
        for (cfg config = 0; config < n_configs; ++config) {
            for (int level = 0; level < n_levels; ++level) {
                mu[config][level] = chemical_potential(
                    occupied(config, level), single_electron_energies[level],
                    sum_occupation(config), voltage);
            }
        }

        /* Generate the rate matrix M. */
        // Use a std::vector because we might decide later to make more matrix
        // elements non-zero, e.g. for radiative decay.
        std::vector<matrix_elem> matrix;
        for (cfg to = 0; to < n_configs; ++to) {
            // TODO: Only store matrix elements greater than some tolerance.
            // For "to" on the following line only, read "away".
            matrix.push_back(diag(to, mu[to], voltage.sd));
            for (int level = 0; level < n_levels; ++level) {
                cfg from = flipped_occupation(to, level);
                matrix.push_back(
                    offdiag(to, from, level, mu[from][level], voltage.sd));
            }
        }

        /* Iterate on dw/dt = Mw until steady state is found. */
        bool converged = false;
        int cycles = 0;
        while (!converged) {
            ++cycles;
            converged = true;  // To be &&'d.
            /* Runge-Kutta (RK4) iteration happens here. */
            double k_1[n_configs] = { 0.0 };
            cfg elem = 0;  // This is a summation variable.
            for (cfg config = 0; config < n_configs; ++config) {
                // We assume that matrix elements are ordered by value of
                // matrix[elem].to (as below).
                while (elem < matrix.size() && matrix[elem].to == config) {
                    k_1[config] += matrix[elem].value*guess[matrix[elem].from];
                    ++elem;
                }
            }
            double k_2[n_configs] = { 0.0 };
            elem = 0;  // Reset summation variable.
            for (cfg config = 0; config < n_configs; ++config) {
                while (elem < matrix.size() && matrix[elem].to == config) {
                    k_2[config] += matrix[elem].value
                                   * (guess[matrix[elem].from]
                                      + (k_1[matrix[elem].from] * h/2.0));
                    ++elem;
                }
            }
            double k_3[n_configs] = { 0.0 };
            elem = 0;
            for (cfg config = 0; config < n_configs; ++config) {
                while (elem < matrix.size() && matrix[elem].to == config) {
                    k_3[config] += matrix[elem].value
                                   * (guess[matrix[elem].from]
                                      + (k_2[matrix[elem].from] * h/2.0));
                    ++elem;
                }
            }
            double k_4[n_configs] = { 0.0 };
            elem = 0;
            for (cfg config = 0; config < n_configs; ++config) {
                while (elem < matrix.size() && matrix[elem].to == config) {
                    k_4[config] += matrix[elem].value
                                   * (guess[matrix[elem].from]
                                      + (k_3[matrix[elem].from] * h));
                    ++elem;
                }
            }
            for (cfg config = 0; config < n_configs; ++config) {
                double increment = (k_1[config]
                                    + 2*k_2[config]
                                    + 2*k_3[config]
                                    + k_4[config])
                                   * h/6.0;
                guess[config] += increment;
                converged = converged && (increment < convergence_criterion);
            }
        }
        // Print the number of Runge-Kutta iterations required for convergence.
        std::cout << " ; " << cycles << " cycles\n";

        /* Find weighted-average current and write it to file. */
        double current = 0.0;
        for (cfg config = 0; config < n_configs; ++config) {
            for (int level = 0; level < n_levels; ++level) {
                current += guess[config]*current_through_level(
                    mu[config][level], occupied(config, level),
                    voltage.sd, level);
            }
        }
        outfile << " ; " << current/e;

        /* Write weights vector to file for later viewing. */
        for (cfg config = 0; config < n_configs; ++config) {
            if (guess[config] > 1e-3)
                outfile << " ; "<< config << " ; " << guess[config];
        }
        
        outfile << "\n";
    }
}
