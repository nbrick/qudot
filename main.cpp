#include <iostream>
#include <utility>
#include <cmath>
#include <limits>
#include <vector>


/* PHYSICAL CONSTANTS */

constexpr double e = 1.60217657e-19;  // Coulombs
constexpr double temp = 4.2;  // Kelvin
constexpr double k_B = 1.3806488e-23;  // Joules per Kelvin


/* USER CONFIGURATION */

constexpr double single_electron_energies[] = {
    -0.1*e, 0.1*e, 0.1*e, 0.1*e
};  // Joules

constexpr double gate_capacitance = 1e-19;  // Farads
constexpr double source_capacitance = 1e-18;  // Farads
constexpr double drain_capacitance = 1e-18;  // Farads
constexpr double extra_capacitance = 1e-19;  // Farads

constexpr double source_width = 1;  // arb.
constexpr double drain_width = 1;  // arb.

constexpr double v_g_min = 0;  // Volts
constexpr double v_g_max = 0;  // Volts
constexpr int v_g_steps = 1;  // (y-axis resolution)

constexpr double v_sd_min = -0.1;  // Volts
constexpr double v_sd_max = 1;  // Volts
constexpr int v_sd_steps = 20;  // (x-axis resolution)

constexpr double source_dos (double energy) {
    if (energy > 0) {
        // Do nothing. (Suppress 'unused' warning.)
    }
    return 1;
}
constexpr double d_fermi_energy = 0.0*e;  // Joules

constexpr double drain_dos (double energy) {
    if (energy > 0) {
        // Do nothing. (Suppress 'unused' warning.)
    }
    return 1;
}
constexpr double s_fermi_energy = 0.0*e;  // Joules

constexpr double evolution_step_size = 1e-1;
constexpr double convergence_criterion = 1e-10;


/* THERMODYNAMICS */

double fermi (double energy, double chem_pot) {
    double exponent = (energy - chem_pot)/(k_B*temp);
    if (exponent > log(std::numeric_limits<double>::max()))
        return 0;
    else if (exponent < log(std::numeric_limits<double>::min()))
        return 1;
    else 
        return 1.0/(exp(exponent) + 1);
}

constexpr double s_shift (double energy, double v_sd) {
    return energy - e*v_sd/2.0;
}
constexpr double d_shift (double energy, double v_sd) {
    return energy + e*v_sd/2.0;
}


/* TYPES AND COMPILE-TIME CALCULATIONS */

typedef struct {
    double gate;
    double sd;  // Source-drain
} v_pair;

// typedef double mu_container[occupied][level][occ_num];
// Should do this. Otherwise have many redundant mu calculations.

typedef unsigned long cfg;

constexpr int n_levels = sizeof(single_electron_energies)
                         /sizeof(single_electron_energies[0]);

typedef struct {
    cfg to;  // Row
    cfg from;  // Column
    double value;
} matrix_element;
typedef double mu_spectrum[n_levels];

constexpr cfg n_configs = 1 << n_levels;  // 2^n_levels

typedef std::pair<double, double[n_configs]> iw_list[n_levels];
typedef mu_spectrum mu_container[n_configs];
typedef double weights[n_configs];
typedef double ddt_weights[n_configs];

constexpr double total_capacitance = gate_capacitance + source_capacitance
                                     + drain_capacitance + extra_capacitance;


/* TUNNEL RATES */

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


/* MISCELLANEOUS FUNCTIONS */

v_pair voltage_pair_from_index (int index) {

    v_pair voltage;

    voltage.gate = v_g_min
                   + ((index % v_g_steps)
                     *(double)(v_g_max - v_g_min)/(double)v_g_steps);

    voltage.sd = v_sd_min
                   + (floor((double)index/(double)v_g_steps))
                     *(double)(v_sd_max - v_sd_min)/(double)(v_sd_steps);

    return voltage;
}


/* BINARY REPRESENTATION OF CONFIGURATIONS */

bool occupied (cfg config, int level) {
    return (config >> level) & 1;
}


int sum_occupation (cfg config) {
    int sum = 0;
    for (int level = 0; level < n_levels; ++level) {
        sum += (int)occupied(config, level);
    }
    return sum;
}


int flipped_occupation (cfg config, int level) {
    return config ^ (1 << level);
}


double chemical_potential (bool occupied, double single_electron_energy,
                           int n_electrons_on_dot, v_pair voltage) {

    return single_electron_energy
           + ((e/total_capacitance)
             *((n_electrons_on_dot + (occupied ? -1 : 1)*0.5)*e
                - gate_capacitance*voltage.gate
                - ((source_capacitance - drain_capacitance)*voltage.sd/2)));
}


matrix_element diagonal_matrix_element (cfg config, mu_spectrum spectrum,
                                        double v_sd) {
    double value = 0;
    for (int level = 0; level < n_levels; ++level) {
        auto mu = spectrum[level];
        if (occupied(config, level)) {
            value -= (out_to_source(mu, v_sd) + out_to_drain(mu, v_sd));
        } else {
            value -= (in_from_source(mu, v_sd) + in_from_drain(mu, v_sd));
        }
    }
    return { config, config, value };
}


matrix_element offdiagonal_matrix_element (cfg to, cfg from, int level,
                                           double mu,
                                           double v_sd) {
    double value;
    if (occupied(from, level)) {
        // Tunnelling from dot into leads.
        value = out_to_source(mu, v_sd) + out_to_drain(mu, v_sd);
    } else {
        // Tunnelling from leads into dot.
        value = in_from_source(mu, v_sd) + in_from_drain(mu, v_sd);
    }
    return { to, from, value };
}


double current_through_level (double mu, bool occupied, double v_sd) {
    if (occupied) {
        return -e*(out_to_drain(mu, v_sd) - out_to_source(mu, v_sd));
    } else {
        return -e*(in_from_source(mu, v_sd) - in_from_drain(mu, v_sd));
    }
}


int main () {

    for (int voltage_index = 0;
         voltage_index < v_g_steps*v_sd_steps;
         ++voltage_index) {

        auto voltage = voltage_pair_from_index(voltage_index);
        std::cout << "\nv_g = " << voltage.gate << " ; "
                  << "v_sd = " << voltage.sd << "\n";

        mu_container mu;
        for (cfg config = 0; config < n_configs; ++config) {
            for (int level = 0; level < n_levels; ++level) {
                mu[config][level] = chemical_potential(
                    occupied(config, level), single_electron_energies[level],
                    sum_occupation(config), voltage);
            }
        }

        // Use a std::vector because we might decide later to make more matrix
        // elements non-zero, e.g. for radiative decay.
        std::vector<matrix_element> matrix;

        for (cfg to = 0; to < n_configs; ++to) {
            // "to" should actually be "away" here.
            matrix.push_back(
                diagonal_matrix_element(
                    to, mu[to], voltage.sd));

            for (int level = 0; level < n_levels; ++level) {
                cfg from = flipped_occupation(to, level);
                matrix.push_back(
                    offdiagonal_matrix_element(
                        to, from, level, mu[from][level], voltage.sd));
            }
        }

        weights guess = { 1 };  // Guess that the dot is in an empty config
                                // to start with. The other elements get
                                // initialised to 0.
        ddt_weights ratechange;

        bool converged = false;
        int cycles = 0;
        while (!converged) {
            ++cycles;
            cfg elem = 0;
            converged = true;  // To be &&'d.
            for (cfg config = 0; config < n_configs; ++config) {
                ratechange[config] = 0;
                // We assume that the matrix elements are ordererd by value of
                // matrix[elem].to.
                while (elem < matrix.size() && matrix[elem].to == config) {
                    ratechange[config] +=
                        matrix[elem].value*guess[matrix[elem].from];
                    ++elem;
                }
                guess[config] += ratechange[config]*evolution_step_size;
                converged = converged
                            && (ratechange[config]*evolution_step_size
                                < convergence_criterion);
            }
        }
        std::cout << "\ncycles for convergence: " << cycles << "\n";

        std::cout << "\nweights:\n";
        for (double weight : guess) {
            std::cout << weight << " ";
        }

        double current = 0;
        for (cfg config = 0; config < n_configs; ++config) {
            for (int level = 0; level < n_levels; ++level) {
                current += guess[config]*current_through_level(
                    mu[config][level], occupied(config, level), voltage.sd);
            }
        }
        std::cout << "\ncurrent is " << current/e << "\n\n";
    }

}
