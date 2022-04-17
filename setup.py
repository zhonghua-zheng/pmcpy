from setuptools import setup, find_packages, Extension
import os

classifiers = [
    "Development Status :: 2 - Pre-Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Atmospheric Science"
    ]

with open("README.rst", "r") as fp:
    long_description = fp.read()

setup(
    name="pmcpy",
    version="0.0.2",
    author="Zhonghua Zheng",
    author_email="zhonghua.zheng@outlook.com",
    url="https://github.com/zzheng93/pmcpy",
    description="A Python package for PartMC post-processing",
    long_description=long_description,
    license="MIT",
    classifiers=classifiers,
    install_requires=['numpy', 'pandas', 'netcdf4', 'xarray'],
    packages=find_packages(),
    )
