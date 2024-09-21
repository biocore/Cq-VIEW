# ----------------------------------------------------------------------------
# Copyright (c) 2018-2022, Amanda Birmingham.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages


setup(
    name='cq-view',
    version="0.1.0",
    license='BSD-3-Clause',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'generate_report=src.generate_reports:generate_report',
            'generate_ucsf_report=src.generate_reports:generate_ucsf_report',
            'generate_ddpcr_report=src.generate_reports:generate_ddpcr_report',
            'intake_new_runs=src.relocate_local_files:intake_new_runs'
        ]}
)
