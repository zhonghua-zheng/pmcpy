## Installation

Step 1: create a conda environment

```bash
$ conda create -n pmcpy python=3.8
$ conda activate pmcpy
$ conda install -c conda-forge numpy pandas xarray netcdf4
```

Step 2: install using `pip`

```bash
$ pip install pmcpy
```

(optional) install from source:: 

```bash
$ git clone https://github.com/zhonghua-zheng/pmcpy.git
$ cd pmcpy
$ python setup.py install
```