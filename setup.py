from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'setuptools',
]
extras_requirements = {
    "predict": [
        'kipoiseq>=0.3.0',
        'mmsplice>=2.1.0'
    ]
}

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'pytest-benchmark']

setup(
    author="Nils Wagner, Aleksandr Neverov",
    author_email='wagnern@in.tum.de',
    classifiers=[
        'Development Status :: 1 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Aberrant splicing prediction across human tissues",
    install_requires=requirements,
    license="GPL-3 license",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='absplice',
    name='absplice2',
    packages=find_packages(include=['absplice2']),
    package_data={
        'absplice': ['absplice/precomputed/*']
    },
    setup_requires=setup_requirements,
    extras_require=extras_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gagneurlab/absplice2',
    version='0.0.1',
    zip_safe=False
)
