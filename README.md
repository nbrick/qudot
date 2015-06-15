# qudot

## Requirements

* gcc
* Python 3
* matplotlib
* termcolor

### Installing the requirements on Debian-based Linuxes (including Ubuntu)

```
sudo apt-get install gcc python3 python3-matplotlib python3-termcolor
```

## Usage

* Set the desired simulation parameters in the `USER CONFIGURATION` section of
  `main.cpp`.
* Execute `./run.sh` -- this operation could take anywhere from under a second
  to several hours, depending on the parameters set in `main.cpp`. The results
  of the computation are stored in `output.csv`.
* Execute `./view.sh`. This launches a Python interpreter with the output data
  already loaded. To view the data graphicall, type in `ui()` or `heatmap()`.
