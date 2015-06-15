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
* Execute `./run.sh`. This operation could take anywhere from under a second to
  several hours, depending on the parameters set in `main.cpp`. The results of
  the computation are stored in `output.csv` but are not meant to be
  human-readable.
* Execute `./view.sh`. This launches a Python interpreter with the output data
  already loaded. To view the data graphically, type in `ui()` or `heatmap()`.

## License

This project is open-sourced under the MIT License.
