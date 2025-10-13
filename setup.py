from setuptools import setup, find_packages

setup(
    name='PARSEbp',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,   # include data files specified in MANIFEST.in

    install_requires=[
        'numpy',
        'tqdm',
    ],

    package_data={
        'PARSEbp': ['bin/*'],   # include executables inside package
    },

    entry_points={
        'console_scripts': [
            'PARSEbp = PARSEbp.cli:main',   # adds CLI command
        ]
    },

    author="Sumit Tarafder and Debswapna Bhattacharya",
    description="Pairwise Agreement-based RNA Scoring with Emphasis on Base Pairings",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Bhattacharya-Lab/PARSEbp",
    license="MIT",
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
)
