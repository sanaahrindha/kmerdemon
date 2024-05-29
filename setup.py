from setuptools import setup, find_packages
from kmerdemon.__version__ import __version__
setup(
    name='kmerdemon',
    version= __version__,
    author="Natalie Norstad, Michael Brown, Sanaah Arramreddy",
    author_email="nnorstad@ucsd.edu, mibrown@ucsd.edu, sarramre@ucsd.edu",
    description="Bioinformatics tool that estimates genome size and the best k-mer length.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    entry_points={  
        "console_scripts": [
            "kmerdemon=kmerdemon.program:main",
        ],
    },
    install_requires=[
        "argparse",
        "matplotlib"
    ]
)
