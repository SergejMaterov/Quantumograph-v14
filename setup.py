from setuptools import setup, find_packages

setup(
    name             = "quantumograph",
    version          = "14.0.0",
    author           = "Sergej Materov",
    author_email     = "sergejmaterov2@gmail.com",
    description      = "Finite quantum graph theory of spacetime",
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = "text/markdown",
    url              = "https://github.com/SergejMaterov/Quantumograph-v14",
    packages         = find_packages(),
    python_requires  = ">=3.10",
    install_requires = [
        "numpy>=1.24",
        "scipy>=1.10",
        "networkx>=3.0",
        "matplotlib>=3.7",
    ],
    extras_require   = {
        "hpc": ["mpi4py>=3.1", "petsc4py>=3.19"],
        "notebook": ["jupyter>=1.0", "ipympl>=0.9"],
    },
    classifiers      = [
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
    ],
)
