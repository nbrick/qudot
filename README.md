# qudot

This program simulates current flowing through a nanocrystal on the single-
electron level, taking into account the capacitance of the nanocrystal and the
single-particle energies of electrons in the nanocrystal. It is assumed that the
nanocrystal has a discrete electron spectrum, i.e. that it constitutes a
quantum dot.

The core program (`main.cpp`) implements rate equations similar to those in
[a paper of C. W. J. Beenakker](http://dx.doi.org/10.1103/PhysRevB.44.1646).
A custom data visualisation tool (based on `matplotlib`) is also included.

For more details, see
[the project report](https://github.com/nbrick/qudot-doc/raw/master/report.pdf).

## Requirements

* GCC
* Python 3
* matplotlib
* termcolor

### Installing the requirements on Debian-based Linuxes (including Ubuntu)

```
sudo apt-get install gcc python3 python3-matplotlib python3-termcolor
```

### Installing elsewhere

This project was developed on Linux but there's no reason it couldn't be run on
Windows, following changes to some of the boilerplate. Please message me at
`nebricks`(whirlpool-symbol)`gmail.com` if you need help, or make a pull request
if you get a Windows version working.

## Usage

* Set the desired simulation parameters in the `USER CONFIGURATION` section of
  `main.cpp`.
* Execute `./run.sh`. This operation could take anywhere from under a second to
  several hours, depending on the parameters set in `main.cpp`. The results of
  the computation are stored in `output.csv` but are not meant to be
  human-readable.
* Execute `./view.sh`. This launches a Python interpreter with the output data
  already loaded. To view the data graphically, type in `ui()` or `heatmap()`.

## License

This project is open-sourced under the MIT License.
