from setuptools import setup, find_packages
from plot_nikov import __version__, __authors__

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='plot_nikov',
    version=__version__,
    license='GNU GENERAL PUBLIC LICENSE',
    description='A program to find reaction chains from annotated datasets',
    long_description=readme(),
    packages=find_packages(),
    url='https://github.com/NLPatVCU/ReactionChainLinker',
    author=__authors__,
    author_email='d.y.plotnikov@gmail.com',
    keywords='',
    classifiers=[
        ''
    ],
    install_requires=[
        'anytree'

    ]

)