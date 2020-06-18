# coding=utf-8

from setuptools import find_packages, setup

dependencies = ["biopython", "pandas", "numpy", ]

setup(
    name="CrisprOpenDB",
    version="0.1",
    packages=find_packages(),

    setup_requires=dependencies,
    install_requires=dependencies,

    author="Edwige Zufferey, Moïra Dion and Pier-Luc Plante",
    maintainer="Edwige Zufferey",
    description="Phage host identification using a CRISPR spacer database",

    license="GL-3",
    keywords="bioinformatics BLAST phage host virome",
    url="https://github.com/edzuf/CrisprOpenDB"
    )