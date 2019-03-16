import os
from setuptools import setup, find_packages


setup(
    name = "IRSkyBack",
    version = "0.0.1",
    author = "Nicholas P. Konidaris",
    author_email = "npk@carnegiescience.edu",
    description = ("Swope data reduction"),
    license = "GNU General Public License v3.0",
    keywords = "Python",
    url = "https://github.com/nickkonidaris/IRSkyBack",
    packages=find_packages(exclude=('tests', 'docs', 'sample')),
    long_description="None",
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python',
    ],
)
