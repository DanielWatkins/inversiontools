"""Tests for the inversion tools. Based on the test soundings."""

import pandas as pd
import invtools

data = pd.read_csv('test_data.csv')
params = invtools.default_params()
inv = data.groupby('date').apply(invtools.find_inversions, params)

