from setuptools import setup, find_packages

setup(
    name='PARSEbp',
    version='0.1.0',
    packages=find_packages(),
    package_data={
        'PARSEbp': ['bin/*'],
    },
    include_package_data=True,
    install_requires=[
        'numpy',
        'tqdm'
    ],
    entry_points={
        'console_scripts': [
            'PARSEbp = PARSEbp.cli:main',
        ]
    },
    description="Pairwise Agreement-based RNA Scoring with Emphasis on Base Pairings",
    author = ["Sumit Tarafder and Debswapna Bhattacharya"],
    repository = "https://github.com/Bhattacharya-Lab/PARSEbp"
)
