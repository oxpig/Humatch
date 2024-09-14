from setuptools import setup, find_packages

setup(
    name='Humatch',
    version='1.0.0',
    description='Fast, gene-specific joint humanisation of antibody heavy and light chains.',
    license='BSD 3-clause license',
    maintainer='Lewis Chinery',
    long_description_content_type='text/markdown',
    maintainer_email='lewis.chinery@dtc.ox.ac.uk',
    include_package_data=True,
    package_data={'': ['trained_models/*', 'germline_likeness_lookup_arrays/*']},
    packages=find_packages(include=('Humatch', 'Humatch.*')),
    install_requires=[
        'numpy',
        'pandas',
        'ipykernel',
        'tensorflow',
        'scikit-learn',
        'seaborn',
        'matplotlib',
        'tqdm',
        'sparse',
        'pyyaml',
        'biopython',
        'hmmer',
    ],
)