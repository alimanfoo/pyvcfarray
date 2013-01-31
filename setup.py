from ast import literal_eval
from distutils.core import setup


def get_version(source='src/vcfarray.py'):
    with open(source) as f:
        for line in f:
            if line.startswith('VERSION'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("VERSION not found")


setup(
    name='vcfarray',
    version=get_version(),
    author='Alistair Miles',
    author_email='alimanfoo@googlemail.com',
    package_dir={'': 'src'},
    py_modules=['vcfarray'],
    url='https://github.com/alimanfoo/pyvcfarray',
    license='MIT License',
    description='Python utility functions for loading data from Variant Call Format (VCF) files into numpy arrays.',
    classifiers=['Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python',
                 'Topic :: Software Development :: Libraries :: Python Modules'
                 ]
)
