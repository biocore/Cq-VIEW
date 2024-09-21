from unittest import TestCase
from src.calculate_viral_loads import \
    _calc_single_raw_conc_per_l_from_std_curve

TEST_CONFIG = \
    {'standard_curves': {
        'PPMoV_N1': {
            'slope': -3.713,
            'intercept': 36.5,
            'use_ln': False,
            'reaction_vol_per_l': 1000000},
        'Promega_N1': {
            'slope': -3.493,
            'intercept': 40.6,
            'use_ln': False,
            'reaction_vol_per_l': 1000000},
        'MPX': {'slope': -3.47,
                'intercept': 35.6,
                'use_ln': False,
                'reaction_vol_per_l': 1000000}
    },
        'spikein_curves': {
            'Promega_N1': {
                'slope': -3.3903,
                'intercept': 40.3,
                'use_ln': False,
                'reaction_vol_per_l': 1000000}
        },
        'input_files': {
            'sheet_name': 'Results',
            'rows_to_skip': 46,
            'sample_key': 'Sample Name',
            'target_key': 'Target Name',
            'raw_cq_key': 'CRT',
            'null_value': 'Undetermined',
            'targets': [
                {'target_name': 'PPMoV_N1', 'ladder_prefix': 'PPMoV_'},
                {'target_name': 'Promega_N1', 'ladder_prefix': 'NE_'},
                {'target_name': 'MPX'}
            ],
            'control_prefixes': ['neg', 'pos'],
            'ladders': {'ladder_bottom': 0.006, 'ladder_steps': 8}
        },
        'output_files': {
            'norm_conc_key': 'Mean viral gene copies/L',
            'collection_date_key': 'Sample_Date',
            'site_locations': {
                'PL': 'PointLoma',
                'ENC': 'Encina',
                'SB': 'SouthBay'
            },
            'target_filename_patterns': [
                {'target_name': 'Promega_N1',
                 'filename_pattern': '{}_sewage_qPCR.csv'}
            ]
        }
    }


class ProcessQpcrDataTest(TestCase):
    def test__calc_single_raw_conc_per_l_from_std_curve(self):
        expected_out = 119336018
        real_out = _calc_single_raw_conc_per_l_from_std_curve(
            raw_cq=33.345837, intercept=40.6, slope=-3.493,
            reaction_vol_per_l=1000000, is_natural_log=False)
        self.assertEqual(expected_out, round(real_out, 0))
