pmcpy: A Python package for PartMC post-processing
--------------------------------------------------
|doi| |docs| |GitHub| |binder| |license| |pepy|

.. |doi| image:: https://zenodo.org/badge/409430865.svg
   :target: https://zenodo.org/badge/latestdoi/409430865

.. |docs| image:: https://readthedocs.org/projects/pmcpy/badge/?version=latest
   :target: https://pmcpy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |GitHub| image:: https://img.shields.io/badge/GitHub-pmcpy-brightgreen.svg
   :target: https://github.com/zzheng93/pmcpy

.. |binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/zzheng93/pmcpy/HEAD?filepath=docs%2Fnotebooks

.. |license| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://github.com/zzheng93/pmcpy/blob/master/LICENSE
   
.. |pepy| image:: https://static.pepy.tech/personalized-badge/pmcpy?period=total&units=international_system&left_color=black&right_color=orange&left_text=Downloads
   :target: https://pepy.tech/project/pmcpy

pmcpy is a **Python** package for `PartMC <https://github.com/compdyn/partmc>`_ post-processing.

Author: `Dr. Zhonghua Zheng <https://zzheng93.github.io/>`_

Installation
------------
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

How to use it?
--------------
You may reproduce the jupyter notebook example on `Binder <https://mybinder.org/v2/gh/zzheng93/pmcpy/HEAD?filepath=docs%2Fnotebooks>`_.

Please check `online documentation <https://zhonghuazheng.com/pmcpy/>`_ (or `Read the Docs <https://pmcpy.readthedocs.io/en/latest/>`_ version) for more information.

How to ask for help
-------------------
The `GitHub issue tracker <https://github.com/zzheng93/pmcpy/issues>`_ is the primary place for bug reports. 
