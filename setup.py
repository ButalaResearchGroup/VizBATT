"""
VizBATT - Visualization tool for battery materials
"""

from setuptools import setup, find_packages

readme = 'README.md'
long_description = open(readme).read()

setup(
    name='vizbatt',
    version='0.1.0',
    description='A visualization tool for exploring and analyzing battery materials, specifically focusing on Wadsley-Roth materials',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Danielle N. Alverson, Kausturi Parui, Eric Fonseca, Steph J. Meikle, Megan M. Butala',
    author_email='',  # Add appropriate contact email
    url='https://github.com/ButalaResearchGroup/VizBATT',
    packages=find_packages(exclude=['docs', 'tests*']),
    license='MIT',
    python_requires='>=3.9',  # REQUIRED: polyhedral_analysis needs Python 3.9+
    install_requires=[
        'dash==2.14.2',
        'plotly==5.24.1',
        'pandas>=1.5.3',
        'numpy<2.0',  # Critical: polyhedral_analysis requires numpy < 2.0
        'matplotlib>=3.5.0',
        # Use the release version instead of the latest commit with broken dependencies
        'polyhedral-analysis @ git+https://github.com/bjmorgan/polyhedral-analysis.git@0.1',
        'dash-ag-grid==2.4.0',
        'dash-bootstrap-components==1.5.0',
        'dash-core-components==2.0.0',
        'dash-html-components==2.0.0',
        'simpleML',
        'gunicorn==20.1.0',
        # Use the latest available pymatgen version (not the non-existent 2024.7.18)
        'pymatgen>=2020.1.28',  # Compatible with the 0.1 release
        'scipy',
        'vg',
        'monty',
        'tqdm',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
            'coverage',
        ],
        'docs': [
            'sphinx',
            'sphinx-rtd-theme',
        ]
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)