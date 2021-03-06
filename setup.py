from setuptools import setup


setup(
        name = 'pyngs',
        provides = 'pyngs',
        version = "0.4.1",
        author = 'J.A. Marin',
        author_email = 'b32masaj@uco.es',
        url = 'https://github.com/JuantonioMS/pyngs',
        description = 'Simple NGS read quality assessment using Python',
        license = 'MIT',
        packages = ['pystq','pyfilter'],
        install_requires = ['numpy','matplotlib','biopython'],
        entry_points = { 'console_scripts': [ 'pystq = pystq.main:main','pyfilter = pyfilter.main:main' ] },
        classifiers = [
                "Development Status :: 1 - Alpha/test",
                "License :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Bioinformatics",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 2.7",
                "Topic :: Scientific/Engineering :: Bioinformatics"
        ]
)
