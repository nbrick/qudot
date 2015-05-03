from math import ceil, log2
import csv

from termcolor import colored
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.widgets import Slider

import data_handling as dat


with open("output.csv", newline="") as csv_file:
    print("Processing CSV... (This might take a few seconds.)")
    lines = csv.reader(csv_file, delimiter=";", quoting=csv.QUOTE_NONNUMERIC)
    n_levels = int(lines.__next__()[0])
    lines = sorted(sorted(lines, key=lambda x: x[1]), key=lambda x: x[0])
    v_sd_range = np.asarray(sorted(list(set([line[1] for line in lines]))))
    v_g_range = np.asarray(sorted(list(set([line[0] for line in lines]))))
    voltage_area = (v_sd_range, v_g_range)
    iw_list = []
    for line in lines:
        current = line[2]
        weights = []
        for index, configuration in enumerate(line[3::2]):
            weights.append((int(configuration), line[2*index + 4]))
        iw_list.append((current, weights))


# Define some colors in the format required by matplotlib.
red = np.array([[1, 0, 0]])
black = np.array([[0, 0, 0]])
white = np.array([[1, 1, 1]])

whitespace = "".join(["\n" for _ in range(100)])


def i_sd(v_g=0):
    gs = gridspec.GridSpec(1, 1)
    axes = plot.subplot(gs[0])
    data = dat.get_i_vs_v_sd(iw_list, v_g, voltage_area)
    axes.plot(v_sd_range, data)

    axes.set_xlim([v_sd_range[0], v_sd_range[-1]])

    axes.set_title(r"at $V_\mathrm{g} = " + "%.2f" % v_g + "\ \mathrm{V}$")
    axes.set_ylabel(r"$\sigma / \mathrm{arb.}$")
    axes.set_xlabel(r"$V_\mathrm{sd} / \mathrm{V}$")

    plot.show(block=False)


def diff_cond(v_g=0):
    gs = gridspec.GridSpec(1, 1)
    axes = plot.subplot(gs[0])
    data = dat.get_diff_conductance_vs_v_sd(iw_list, v_g, voltage_area)
    axes.plot(v_sd_range, data)

    axes.set_xlim([v_sd_range[0], v_sd_range[-1]])

    axes.set_title(r"at $V_\mathrm{g} = " + "%.2f" % v_g + "\ \mathrm{V}$")
    axes.set_ylabel(r"$\frac{\mathrm{d}\sigma} "
                    + r"{\mathrm{d}V_\mathrm{sd}} / \mathrm{arb.}$")
    axes.set_xlabel(r"$V_\mathrm{sd} / \mathrm{V}$")

    plot.show(block=False)


def mean_occupation(v_g=0):
    gs = gridspec.GridSpec(1, 1)
    axes = plot.subplot(gs[0])
    data = dat.get_mean_occupation_vs_v_sd(iw_list, v_g, voltage_area)
    axes.plot(v_sd_range, data)

    axes.set_xlim([v_sd_range[0], v_sd_range[-1]])

    axes.set_title(r"at $V_\mathrm{g} = " + "%.2f" % v_g + "\ \mathrm{V}$")
    axes.set_ylabel(r"Mean dot occupation number")
    axes.set_xlabel(r"$V_\mathrm{sd} / \mathrm{V}$")

    plot.show(block=False)


def heatmap(axes=None):
    if axes is None:
        show = True
        gs = gridspec.GridSpec(1, 1)
        heatmap_axes = plot.subplot(gs[0])
    else:
        show = False
        heatmap_axes, colorbar_axes = axes

    current_function, extent = dat.get_plottable_diff_conductance_in_v_space(
        iw_list, voltage_area)
    heatmap_ = heatmap_axes.imshow(current_function, extent=extent,
                                   interpolation="nearest", aspect="auto",
                                   cmap=cm.binary)

    heatmap_axes.set_xlim([v_sd_range[0], v_sd_range[-1]])
    heatmap_axes.set_ylim([v_g_range[0], v_g_range[-1]])

    heatmap_axes.set_xlabel(r"$V_\mathrm{sd}/\mathrm{V}$", size=25)
    heatmap_axes.set_ylabel(r"$V_\mathrm{g}/\mathrm{V}$", size=25)
    heatmap_axes.set_title(
        r"$\frac{\partial I}{\partial V_\mathrm{sd}} /\mathrm{arb.\ units}$",
        size=25, y=1.04)


    heatmap_axes.locator_params(axis="x", nbins=5)

    for tick in heatmap_axes.xaxis.get_major_ticks():
                tick.label.set_fontsize(20)

    for tick in heatmap_axes.yaxis.get_major_ticks():
                tick.label.set_fontsize(20)

    plot.subplots_adjust(bottom=0.15)

    if show:
        plot.show(block=False)


def pretty_bin(number, max_number):
    width = ceil(log2(max_number))
    return "".join([colored("*", "white") if dat.bit(number, index)
                    else colored("|", "red")
                    for index in range(width)])


def ui(v_sd=0, v_g=0):
    """Display a graphical user interface for data exploration."""

    # Define the layout of the UI and add "plots" to UI regions.
    gs = gridspec.GridSpec(
        6, 3,
        height_ratios=[10, 2, 10, 2, 1, 1], width_ratios=[20, 1, 2])
    heatmap_axes = plot.subplot(gs[0])
    colorbar_axes = plot.subplot(gs[2])
    line_plot_axes = plot.subplot(gs[6])
    v_sd_slider_axes = plot.subplot(gs[12])
    v_g_slider_axes = plot.subplot(gs[15])

    # Slider UI elements
    # ------------------
    v_sd_slider = Slider(v_sd_slider_axes, r"$V_\mathrm{sd}/\mathrm{V}$",
                         v_sd_range[0], v_sd_range[-1],
                         valinit=v_sd)
    v_g_slider = Slider(
        v_g_slider_axes, r"$V_\mathrm{g}/\mathrm{V}$",
        v_g_range[0], v_g_range[-1],
        valinit=v_g)

    # Use an object to store interactive voltages.
    # This makes it simpler to keep track of slider changes persistently.
    class Voltages:
        sd = v_sd
        g = v_g

    v = Voltages()

    # Conductance heat map
    # --------------------
    heatmap(tuple([heatmap_axes, colorbar_axes]))

    # We will add a small cross marking the current position in v_g-v_sd space
    # to the heatmap plot. To keep track of its position, we use the following
    # object.
    class PointAnnotation():
        point = heatmap_axes.scatter(v_sd, v_g, marker="+")

    annotation = PointAnnotation()

    def replot_point_annotation(v_sd_, v_g_):
        annotation.point.remove()
        annotation.point = heatmap_axes.scatter(v_sd_, v_g_, marker="+")

    # We also add lines annotating values of v_sd and v_g on the two large
    # plots. v_g_line is persistent and has its y_data updated in redraw().
    # There is no persistent equivalent v_sd_line, because the line plot is
    # completely erased and replotted on each redraw() call; we will make a
    # new axvline each time.
    line_plot_axes.axvline(x=v.sd)
    v_g_line = heatmap_axes.axhline(y=v.g)

    # Labels
    # ------
    def label_figures():
        line_plot_axes.set_xlabel(r"$V_\mathrm{sd}/\mathrm{V}$")
        line_plot_axes.set_ylabel(r"Current $I / \mathrm{arb.}$")

    # Slider update actions
    # ---------------------
    def redraw():
        # Move annotations.
        replot_point_annotation(v.sd, v.g)
        v_g_line.set_ydata(v.g)

        # Retrieve from memory and plot current.
        i_vs_v_sd = dat.get_i_vs_v_sd(iw_list, v.g, voltage_area)
        line_plot_axes.clear()
        line_plot_axes.plot(v_sd_range,
                            [current for current in i_vs_v_sd],
                            "black")
        line_plot_axes.set_xlim([v_sd_range[0], v_sd_range[-1]])
        line_plot_axes.axvline(x=v.sd)

        current, weights = dat.get_iw_tuple(iw_list, v.sd, v.g, voltage_area)

        print(whitespace)

        v.sd, v.g = dat.get_voltage_pair_from_index(
            dat.get_index_from_voltage_pair(v.sd, v.g, voltage_area),
            voltage_area)

        for index, weight in sorted([(index, weight)
                                     for index, weight in weights],
                                    key=lambda x: x[1], reverse=False):

                weight_bar = "".join([colored("=", "green") if point < weight
                                      else " "
                                      for point
                                      in np.linspace(0, 1 - 1e-10, 40)])

                print(" " + pretty_bin(index, 2**n_levels),
                      "%.3f" % weight, weight_bar)

        print("\nv_g/V =", "%.3f" % v.g, "; v_sd/V =", "%.3f" % v.sd)
        print("mean occupation =", "%.3f" % dat.mean_occupation(weights))
        print("current =", "%.3f" % (current), "arb. units")

        print("\n>>> ", end="")

        label_figures()

    def update_v_sd(v_sd_):
        v.sd = v_sd_
        redraw()

    def update_v_g(v_g_):
        v.g = v_g_
        redraw()

    # Set sliders to listen for clicks.
    v_sd_slider.on_changed(update_v_sd)
    v_g_slider.on_changed(update_v_g)

    # Draw with default v_g, v_sd values.
    redraw()

    # Launch the matplotlib window.
    plot.show(block=False)

print("Try ui(), or heatmap().")
