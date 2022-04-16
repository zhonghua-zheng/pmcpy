Installation
============

Step 1: create an environment::

    $ conda create -n pmcpy python=3.8
    $ conda activate pmcpy
    $ conda install -c conda-forge numpy pandas xarray netcdf4

Step 2: install using pip::

    $ pip install pmcpy

(optional) install from source:: 

    $ git clone https://github.com/zzheng93/pmcpy.git
    $ cd pmcpy
    $ python setup.py install
