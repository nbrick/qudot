#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <fstream>


/* PHYSICAL CONSTANTS */

constexpr double e = 1.60217657e-19;   // Coulombs
constexpr double k_B = 1.3806488e-23;  // Joules per Kelvin


/* USER CONFIGURATION */

/* Temperature */
constexpr double temp = 30;  // Kelvin

/* Dot orbital energies */
constexpr double single_electron_energies[] = {
    -0.1*e, -0.1*e, -0.09*e, -0.09*e, -0.07*e, -0.07*e, -0.04*e, -0.04*e
};  // Joules

/* Dot-system capacitances */
constexpr double gate_capacitance =   1e-19;  // Farads
constexpr double source_capacitance = 1e-18;  // Farads
constexpr double drain_capacitance =  3e-18;  // Farads
constexpr double extra_capacitance =  1e-19;  // Farads

/* Tunnel rates (assumed to be independent of electron energy) */
// TODO: Relax this assumption.
constexpr double source_width = 1;  // arb.
constexpr double drain_width =  1;  // arb.

/* Voltage-space to be explored */
constexpr double v_g_min = -5;  // Volts
constexpr double v_g_max = 10;  // Volts
constexpr int v_g_steps = 200;  // (y axis resolution)

constexpr double v_sd_min = 0;  // Volts
constexpr double v_sd_max = 0.08;  // Volts
constexpr int v_sd_steps = 100;  // (x axis resolution)

/* Electronic properties of leads (s: source; d: drain) */
constexpr double source_dos (double energy) {
    if (energy > 0.0) {
        // Do nothing. (Suppress 'unused' warning.)
    }
    return 1;
}
constexpr double d_fermi_energy = 0.0*e;  // Joules
constexpr double drain_dos (double energy) {
    if (energy > 0.0) {
        // Do nothing. (Suppress 'unused' warning.)
    }
    return 1;
}
constexpr double s_fermi_energy = 0.0*e;  // Joules

/* Semiclassical simulation properties */
constexpr double evolution_step_size = 5e-2;
constexpr double convergence_criterion = 1e-10;


/* COMPILE-TIME CALCULATIONS */

constexpr double total_capacitance = gate_capacitance + source_capacitance
                                     + drain_capacitance + extra_capacitance;


/* TYPES */

typedef struct {
    double gate;
    double sd;  // Source-drain
} v_pair;

constexpr int n_levels = sizeof(single_electron_energies)
                         /sizeof(single_electron_energies[0]);

typedef unsigned long cfg;
constexpr cfg n_configs = 1 << n_levels;  // 2^n_levels

typedef struct {
    cfg to;    // Row
    cfg from;  // Column
    double value;
} matrix_elem;

typedef double mu_spectrum[n_levels];
typedef mu_spectrum mu_container[n_configs];
// TODO: typedef double mu_container[occupied][level][occ_num];
// Should do this. Otherwise there are many redundant mu calculations.

typedef double weights[n_configs];
typedef double ddt_weights[n_configs];


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

int flipped_occupation (cfg config, int level) {
    return config ^ (1 << level);
}


/* THERMODYNAMICS */

double fermi (double energy, double chem_pot) {
    double exponent = (energy - chem_pot)/(k_B*temp);
    if (exponent > log(std::numeric_limits<double>::max()))
        return 0.0;
    else if (exponent < log(std::numeric_limits<double>::min()))
        return 1.0;
    else
        return 1.0/(exp(exponent) + 1);
}


/* THE PHYSICS */

/* Energy shifts in leads owing to applied bias */

constexpr double s_shift (double energy, double v_sd) {
    return energy - e*v_sd/2.0;
}

constexpr double d_shift (double energy, double v_sd) {
    return energy + e*v_sd/2.0;
}

/* Tunnel rates */

double in_from_source (double mu, double v_sd) {
    return source_width
           * source_dos(mu - s_shift(0.0, v_sd))
           * fermi(mu, s_shift(s_fermi_energy, v_sd));
}

double in_from_drain (double mu, double v_sd) {
    return drain_width
           * drain_dos(mu - d_shift(0.0, v_sd))
           * fermi(mu, d_shift(d_fermi_energy, v_sd));
}

double out_to_source (double mu, double v_sd) {
    return source_width
           * source_dos(mu - s_shift(0.0, v_sd))
           * (1 - fermi(mu, s_shift(s_fermi_energy, v_sd)));
}

double out_to_drain (double mu, double v_sd) {
    return drain_width
           * drain_dos(mu - d_shift(0.0, v_sd))
           * (1 - fermi(mu, d_shift(d_fermi_energy, v_sd)));
}

double current_through_level (double mu, bool occupied, double v_sd) {
    if (occupied)
        return -e*(out_to_drain(mu, v_sd) - out_to_source(mu, v_sd));
    else
        return -e*(in_from_source(mu, v_sd) - in_from_drain(mu, v_sd));
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

matrix_elem diag (cfg config, mu_spectrum spectrum,
                                        double v_sd) {
    double value = 0;
    for (int level = 0; level < n_levels; ++level) {
        auto mu = spectrum[level];
        if (occupied(config, level))
            value -= (out_to_source(mu, v_sd) + out_to_drain(mu, v_sd));
        else
            value -= (in_from_source(mu, v_sd) + in_from_drain(mu, v_sd));
    }
    return { config, config, value };
}

matrix_elem offdiag (cfg to, cfg from, int level, double mu, double v_sd) {
    double value;
    if (occupied(from, level))
        value = out_to_source(mu, v_sd) + out_to_drain(mu, v_sd);
    else
        value = in_from_source(mu, v_sd) + in_from_drain(mu, v_sd);
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
                          *(double)(v_g_max - v_g_min)/(double)v_g_steps);
    else
        voltage.gate = v_g_max
                       - (((index % v_g_steps) + 1)
                          *(double)(v_g_max - v_g_min)/(double)v_g_steps);

    voltage.sd = v_sd_min
                   + (floor((double)index/(double)v_g_steps)
                      *(double)(v_sd_max - v_sd_min)/(double)(v_sd_steps));

    return voltage;
}


/* THE SIMULATION */

int main () {

    std::ofstream outfile ("output.csv", std::ofstream::out);
    // outfile << n_levels << "\n";  // Needed for visualization later.

    /* 'Guess' that the dot is initially empty; w = { 1, 0, 0, ... }. */
    weights guess = { 1 };

    /* Iterate through points in voltage-space. */
    for (int voltage_index = 0;
         voltage_index < v_g_steps*v_sd_steps;
         ++voltage_index) {

        /* Normalise weights guess to prevent drift (and for ease of input). */
        double sum_weights = 0;
        for (cfg config = 0; config < n_configs; ++config)
            sum_weights += guess[config];
        for (cfg config = 0; config < n_configs; ++config)
            guess[config] /= sum_weights;

        /* Choose a point in voltage-space. */
        auto voltage = voltage_pair_from_index(voltage_index);
        outfile << voltage.gate << " ; " << voltage.sd;
        std::cout << "v_g:" << voltage.gate << " ; "
                  << "v_sd:" << voltage.sd;

        /* For each possible config, find all the chemical potentials. */
        mu_container mu;
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
        // TODO: Store just one row of the matrix in memory.
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

        /* Iterate on dw/dt = M w until steady state is found. */
        // TODO: Use Runge-Kutta.
        ddt_weights ratechange;
        bool converged = false;
        int cycles = 0;
        while (!converged) {
            ++cycles;
            cfg elem = 0;
            converged = true;  // To be &&'d.
            for (cfg config = 0; config < n_configs; ++config) {
                ratechange[config] = 0;
                // We assume that matrix elements are ordered by value of
                // matrix[elem].to.
                while (elem < matrix.size() && matrix[elem].to == config) {
                    ratechange[config] +=
                        matrix[elem].value*guess[matrix[elem].from];
                    ++elem;
                }
                guess[config] += ratechange[config]*evolution_step_size;
                converged = converged
                            && (ratechange[config] < convergence_criterion);
            }
        }
        std::cout << " ; " << cycles << " cycles\n";

        /* Find weighted-average current and write it to file. */
        double current = 0;
        for (cfg config = 0; config < n_configs; ++config) {
            for (int level = 0; level < n_levels; ++level) {
                current += guess[config]*current_through_level(
                    mu[config][level], occupied(config, level), voltage.sd);
            }
        }
        outfile << " ; " << current/e;

        // for (cfg config = 0; config < n_configs; ++config) {
        //     if (guess[config] > 1e-5)
        //         outfile << " ; "<< config << " ; " << guess[config];
        // }

        /* Write weights vector to file for later viewing. */
        for (cfg config = 0; config < n_configs; ++config) {
            outfile << " ; " << guess[config];
        }

        outfile << "\n";
    }
}
