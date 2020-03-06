import sys
from setuptools import setup
from warnings import warn

if sys.version_info.major != 3:
    raise RuntimeError("Palantir requires Python 3")
if sys.version_info.minor < 6:
    warn("Analysis methods were developed using Python 3.6")

# get version
with open("src/harmony/version.py") as f:
    exec(f.read())
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="harmonyTS",
    version=__version__,  # read in from the exec of version.py; ignore error
    description=(
        "Harmony is a unified framework for data visualization, analysis "
        "and interpretation of scRNA-seq data measured across discrete time points"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dpeerlab/harmony",
    author="Manu Setty",
    author_email="manu.talanki@gmail.com",
    package_dir={"": "src"},
    packages=["harmony"],
    install_requires=[
        "numpy>=1.14.2",
        "pandas>=0.22.0",
        "scipy>=1.0.1",
        "sklearn",
        "fa2",
        "matplotlib>=2.2.2",
        "seaborn>=0.8.1",
        "palantir>=0.2.3"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires='>=3.6',
)
