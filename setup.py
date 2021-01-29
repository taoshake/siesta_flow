#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name='siesta_flow',
    version='0.1',
    description="Jiahui's siesta workflow",
    author='Jiahui Jia',
    author_email='jiahui.jia@icn2.cat',
    license='MIT',
    packages=find_packages(),
    install_requires=['numpy','matplotlib','scipy','ase', 'spglib', 'phonopy'],
    scripts=[],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT',
    ])
