from math import floor, ceil, log2

import numpy as np


def bit(number, index):
    """Return the indexth bit of number."""
    return (number >> index) & 1

def sum_bits(number):
    """Return the sum of digits in the binary representation of number."""
    if number == 0:
        return 0
    sum_ = 0
    for i in range(ceil(log2(number)) + 1):
        sum_ += bit(number, i)
    return sum_


def get_voltage_pair_from_index(n, voltage_area):
    v_sd_range, v_g_range = voltage_area
    v_sd = v_sd_range[n % len(v_sd_range)]
    v_g = v_g_range[floor(n / len(v_sd_range))]
    return v_sd, v_g


def index_of_closest_element(ascending_list, datum):
    """Return index of the list element whose value is closest to datum."""
    old_delta = abs(ascending_list[0] - datum)
    for index, element in enumerate(ascending_list[1:]):
        delta = abs(element - datum)
        if delta > old_delta:
            return index
        old_delta = delta
    return len(ascending_list) - 1


def get_index_from_voltage_pair(v_sd, v_g, voltage_area):
    # It is assumed that v_sd_range, v_g_range are ascending.
    v_sd_range, v_g_range = voltage_area
    small_number = index_of_closest_element(v_sd_range, v_sd)
    big_number = index_of_closest_element(v_g_range, v_g)*len(v_sd_range)
    index = big_number + small_number
    return index


def get_iw_tuple(iw_list, v_sd, v_g, voltage_area):
    index = get_index_from_voltage_pair(v_sd, v_g, voltage_area)
    return iw_list[index]


def get_iw_vs_v_sd(iw_list, v_g, voltage_area):
    v_sd_range, _ = voltage_area
    start = get_index_from_voltage_pair(v_sd_range[0], v_g, voltage_area)
    end = get_index_from_voltage_pair(v_sd_range[-1], v_g, voltage_area) + 1
    return iw_list[start:end]


def get_i_vs_v_sd(iw_list, v_g, voltage_area):
    return np.asarray([iw_tuple[0] for iw_tuple
                       in get_iw_vs_v_sd(iw_list, v_g, voltage_area)])


def get_diff_conductance_vs_v_sd(iw_tuple, v_g, voltage_area):
    return np.gradient(get_i_vs_v_sd(iw_tuple, v_g, voltage_area))


def get_i_vs_v_g(iw_list, v_sd, voltage_area):
    v_sd_range, v_g_range = voltage_area
    start = get_index_from_voltage_pair(v_sd, v_g_range[0], voltage_area)
    end = get_index_from_voltage_pair(v_sd, v_g_range[-1], voltage_area) + 1
    step = len(v_sd_range)
    return np.asarray([ie_tuple[0] for ie_tuple in iw_list[start:end:step]])


def get_diff_conductance_vs_v_g(iw_list, v_sd, voltage_area):
    return np.gradient(get_i_vs_v_g(iw_list, v_sd, voltage_area))


def get_plottable_diff_conductance_in_v_space(iw_list, voltage_area):
    # Informed by http://stackoverflow.com/questions/6323737/
    x, y = voltage_area
    # 0th element of ie_tuple is the current.
    print("here")
    z = np.asarray([get_diff_conductance_vs_v_sd(iw_list, v_g, voltage_area)
                    for v_g in voltage_area[1]])
    nrows, ncols = len(y), len(x)
    grid = z.reshape((nrows, ncols))
    return grid, (x.min(), x.max(), y.max(), y.min())


def mean_occupation(weights):
    return sum(weight*sum_bits(config) for config, weight in weights)


def get_mean_occupation_vs_v_sd(iw_list, v_g, voltage_area):
    return np.asarray([mean_occupation(weights) for current, weights
                       in get_iw_vs_v_sd(iw_list, v_g, voltage_area)])
