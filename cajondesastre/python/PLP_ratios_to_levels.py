# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:44:41 2024

@author: Javier.Burgos-Marmol
"""

import numpy as np
import pandas as pd


class ratios_matrix:
    _kind = 'Ratios Matrix'
    __slots__ = ['matrix', 'input_data', 'total_level', 'levels', 'n_levels', 'groups']

    def __init__(self, input_data):
        self.n_levels = len(input_data)
        self.matrix = np.identity(self.n_levels, dtype=float)
        for m in range(1, self.n_levels):  # COLUMNS
            for n in range(m):  # ROWS
                self.matrix[n][m] = -1
                self.matrix[m][n] = -1
        self.levels = []
        groups = [set()]
        for prop, val in input_data.items():
            ratio_prop = prop.split('_')
            added = False
            if ratio_prop[0] == 'ratio':
                row = int(ratio_prop[1])-1
                column = int(ratio_prop[2])-1
                if isinstance(val, int) or isinstance(val, float):
                    self.matrix[row, column] = val
                elif isinstance(val, str):
                    try:
                        self.matrix[row, column] = float(val)
                    except Exception:
                        err_msg = "Parsed elements in data record must be numerical (int, float)."
                        raise ValueError(err_msg)
                else:
                    err_msg = "Parsed elements in data record must be numerical (int, float)."
                    raise ValueError(err_msg)
                if self.matrix[row, column] > 0:
                    self.matrix[column, row] = 1.0 / self.matrix[row, column]
                elif self.matrix[row, column] == 0:
                    self.matrix[column, row] = np.infty
                if len(groups[0]) == 0:
                    groups[0].add(row)
                    groups[0].add(column)
                else:
                    for s in groups:
                        if row in s and column in s:
                            added = True
                        elif row not in s and column in s:
                            s.add(row)
                            added = True
                        elif column not in s and row in s:
                            s.add(column)
                            added = True
                        else:
                            continue
                    if added is False:
                        groups.append({row, column})
            else:
                if isinstance(val, int) or isinstance(val, float):
                    self.total_level = val
                elif isinstance(val, str):
                    try:
                        self.total_level = float(val)
                    except Exception:
                        err_msg = "Parsed elements in data record must be numerical (int, float)."
                        raise ValueError(err_msg)
                else:
                    err_msg = "Parsed elements in data record must be numerical (int, float)."
                    raise ValueError(err_msg)
                continue
        self.groups = []
        for s1 in range(len(groups)-1):
            for s2 in range(s1, len(groups)):
                if len(groups[s1].intersection(groups[s2])) > 0:
                    groups[s1] = groups[s1].union(groups[s2])
                    groups[s2] = groups[s1]
        for s in groups:
            if s not in self.groups:
                self.groups.append(s)

        if len(self.groups) > 1:
            err_msg = "Two or more independent groups of ingredients given with one total level only. Please, use one group of ingredients per calculation only."
            raise ValueError(err_msg)

        self.input_data = input_data

    def __repr__(self):
        string = (self._kind + " object containing " + str(self.n_levels) +
                  " in " + str(len(self.groups)) + " independent groups.")
        return string

    def __iter__(self):
        return np.nditer(self.matrix)

    def complete_matrix(self):

        for n in range(self.n_levels-1):  # ROWS
            for m in range(n+1, self.n_levels):  # COLUMNS
                if self.matrix[n][m] == -1:
                    for s in self.groups:
                        if n in s and m in s:
                            for k in range(len(s)):
                                if k != n and k != m:
                                    if self.matrix[n][k] >= 0 and self.matrix[k][m] >= 0:
                                        if self.matrix[n][k] != np.infty and self.matrix[k][m] != np.infty:
                                            self.matrix[n][m] = self.matrix[n][k] * self.matrix[k][m]
                                            if self.matrix[n][m] > 0.0:
                                                self.matrix[m][n] = 1.0 / self.matrix[n][m]
                                            else:
                                                self.matrix[m][n] = np.infty
                                        elif self.matrix[n][k] == np.infty and self.matrix[k][m] != np.infty and self.matrix[k][m] > 0:
                                            self.matrix[n][m] == np.infty
                                            self.matrix[m][n] == 0.0
                                        elif self.matrix[n][k] != np.infty and self.matrix[k][m] == np.infty and self.matrix[n][k] > 0:
                                            self.matrix[n][m] = 0.0
                                            self.matrix[m][n] = np.infty
                                        else:
                                            continue
                        else:
                            continue
        return

    def calculate_levels(self):
        """ Calculate levels from ratio matrix and total level.

        :raises ValueError: If too many 0 values.
        :return: levels
        :rtype: list

        """
        dim = self.n_levels
        full_column = [0 for n in range(dim)]
        full_dim = 0
        swap_list = [[], []]
        for m in reversed(range(dim)):  # COLUMNS
            for n in range(dim):  # ROWS
                if n == m or (self.matrix[n][m] >= 0 and self.matrix[n][m] != np.infty):
                    full_column[m] += 1
            if full_column[m] == dim:
                full_dim = m
                break
        if full_dim == 0:
            err_msg = 'Please, check ratio values, too many 0s.'
            raise ValueError(err_msg)

        if full_dim < dim-1 and full_dim > 0:
            self.matrix[[full_dim, dim-1]] = self.matrix[[dim-1, full_dim]]
            self.matrix[:, [ full_dim, dim-1]] = self.matrix[:, [dim-1, full_dim]]
            swap_list[0].append(full_dim)
            swap_list[0].append(dim-1)

        # Calculate levels
        ratios_vector = self.matrix[:,-1][:-1]
        partial_levels = [0 for n in range(dim)]
        for m in reversed(range(dim)):
            if m == dim-1:
                partial_levels[-1] = self.total_level / (1.0 + sum(ratios_vector))
            else:
                partial_levels[m] = ratios_vector[m] * partial_levels[-1]

        #UNDO SWAP IF DONE BEFORE
        if len(swap_list[0]) > 0:
            for m in reversed(range(len(swap_list[0]))):
                self.matrix[[swap_list[0][m], swap_list[1][m]]] = self.matrix[[swap_list[1][m], swap_list[0][m]]]
                self.matrix[:, [swap_list[0][m], swap_list[1][m]]] = self.matrix[:, [swap_list[1][m], swap_list[0][m]]]

        self.levels = partial_levels

        return


def parse_variables(indata):
    varnames = {}
    ingredient_list = indata['Ingredient_names']
    n_ings = len(ingredient_list)
    varnames['Total_level'] = indata['Total_level_property_name']
    ratio = {}
    for key, value in indata.items():
        if "Ratio" in key and "Ingredient" not in key:
            namelist = key.split()
            ratio[value] = [1, 2]
    for key, value in indata.items():
        if "Ratio" in key and "Ingredient" in key:
            namelist = key.split()
            rationame = indata[namelist[0] + " " + namelist[1]]
            ratio[rationame][int(namelist[-1])-1] = value

    for key, value in ratio.items():
        name = "ratio"
        for i in range(len(value)):
            name += "_"
            n = ingredient_list.index(value[i])
            name += str(n+1)
        varnames[name] = key
    for n in range(len(ingredient_list)):
        varnames["level_"+str(n+1)] = ingredient_list[n]

    return varnames


def test_data(data_type=4):
    """Datasets to test the script.

    :param data_type: Choose 1 for "one level in all ratios", 2 for "telescopic", and 3 for non of the above, defaults to 1
    :type data_type: int, optional
    :return: Parsed input data
    :rtype: dict

    """
    # TEST DATA with 4 Ingredients. DO NOT INCLUDE IN PLP
    output = pd.DataFrame()
    if data_type == 1:  # TEST DATA with 4 Ingredients, ONE FIXED
        output['Total_level'] = np.array([68])
        output['ratio_1_4'] = np.array([2])
        output['ratio_2_4'] = np.array([3])
        output['ratio_3_4'] = np.array([0.2])
        varnames={'Total_level': 'Total_level', 'ratio_1_4': 'ratio_1_4',
                  'ratio_2_4': 'ratio_2_4', 'ratio_3_4': 'ratio_3_4',
                  'a_level': 'level_1', 'b_level': 'level_2',
                  'd_level': 'level_3', 'd_level': 'level_4'}
    elif data_type == 2:  # TEST DATA with 4 Ingredients, TELESCOPIC
        output['Total_level'] = np.array([68])
        output['ratio_1_2'] = np.array([2.0/3.0])
        output['ratio_2_3'] = np.array([15])
        output['ratio_3_4'] = np.array([0.2])
        varnames={'Total_level': 'Total_level', 'ratio_1_2': 'ratio_1_2',
                  'ratio_2_3': 'ratio_2_3', 'ratio_3_4': 'ratio_3_4',
                  'a_level': 'level_1', 'b_level': 'level_2',
                  'd_level': 'level_3', 'd_level': 'level_4'}
    elif data_type == 3:  # TEST DATA with 4 Ingredients, NONE OF THE ABOVE
        output['Total_level'] = np.array([68])
        output['ratio_1_2'] = np.array([2.0/3.0])
        output['ratio_1_3'] = np.array([10])
        output['ratio_3_4'] = np.array([0.2])
        varnames={'Total_level': 'Total_level', 'ratio_1_2': 'ratio_1_2',
                  'ratio_1_3': 'ratio_1_3', 'ratio_3_4': 'ratio_3_4',
                  'a_level': 'level_1', 'b_level': 'level_2',
                  'd_level': 'level_3', 'd_level': 'level_4'}
    elif data_type == 4:    # TEST DATA with 4 Ingredients, NONE OF THE ABOVE + PLP names
        output['Total_surf'] = np.array([68]*8)
        output['ad_ratio'] = np.array([2]*8)
        output['bd_ratio'] = np.array([3]*8)
        output['cd_ratio'] = np.array([0.2]*8)
        varnames={'Total_level': 'Total_surf', 'ratio_1_4': 'ad_ratio',
                  'ratio_2_4': 'bd_ratio', 'ratio_3_4': 'cd_ratio',
                  'level_1': 'a_level', 'level_2': 'b_level',
                  'level_3': 'c_level', 'level_4': 'd_level'}
    elif data_type == 5:    # TEST DATA with 5 Ingredients, NONE OF THE ABOVE + PLP names + several 0s
        output['Total_surf'] = np.array([100.0]*8)
        output['ab_ratio'] = np.array([0.0]*8)
        output['bc_ratio'] = np.array([2.0]*8)
        output['bd_ratio'] = np.array([1.0]*8)
        output['ce_ratio'] = np.array([4.0]*8)
        varnames={'Total_level': 'Total_surf', 'ratio_1_2': 'ab_ratio',
                  'ratio_2_3': 'bc_ratio', 'ratio_2_4': 'bd_ratio',
                  'ratio_3_5': 'ce_ratio',
                  'level_1': 'a_level', 'level_2': 'b_level',
                  'level_3': 'c_level', 'level_4': 'd_level',
                  'level_5': 'e_level'}
    elif data_type == 6:    # REAL DATA with 3 Ingredients
        output['Total_Surfactant_level'] = np.array([7.8]*10)
        output['Amphoteric_Anionic_ratio'] = np.array([0.6]*5+[1.1]*5)
        output['Isethionate_Tertiary_ratio'] = np.array([0.5, 3, 5.5, 8, 10.5]*2)
        varnames={'Total_level': 'Total_Surfactant_level',
                  'ratio_1_2': 'Amphoteric_Anionic_ratio',
                  'ratio_2_3': 'Isethionate_Tertiary_ratio',
                  'level_1': 'CAPB', 'level_2': 'SLI',
                  'level_3': 'Glycinate'}
    # elif data_type == 7:    # REAL DATA with 3 Ingredients. REAL PLP NAMES
    else:
        output['Total_Surfactant_level'] = np.array([str(7.8)]*10)
        output['Amphoteric_Anionic_ratio'] = np.array([str(0.6)]*5+[str(1.1)]*5)
        output['Isethionate_Tertiary_ratio'] = np.array([str(0.5), str(3),
                                                         str(5.5), str(8),
                                                         str(10.5)]*2)
        varnames={'Total_level_property_name': 'Total_Surfactant_level',
                  'Ratio_name 1': 'Amphoteric_Anionic_ratio',
                  'Ratio_name 2': 'Isethionate_Tertiary_ratio',
                  'Ratio_name 1 Ingredient 1': 'CAPB',
                  'Ratio_name 1 Ingredient 2': 'SLI',
                  'Ratio_name 2 Ingredient 1': 'SLI',
                  'Ratio_name 2 Ingredient 2': 'Glycinate',
                  'Ingredient_names': ['CAPB', 'SLI', 'Glycinate']}
    return output, varnames


# Test values. DO NOT INCLUDE IN PLP
plp_df, plp_globals = test_data(7)

# Translation from PLP property names to python variables
#plp_inputs = plp_globals
plp_inputs = parse_variables(plp_globals)

# Need to add Pipeline Pilot INPUT PROCESSING
output = {}
for local_name, plp_name in plp_inputs.items():
    if local_name.startswith('level_'):
        output[local_name] = []
for dr in range(len(plp_df)):
    ind_data_record = {}
    for local_name, plp_name in plp_inputs.items():
        if local_name.startswith('level_'):
            continue
        ind_data_record[local_name] = plp_df[plp_name][dr]
    ratios = ratios_matrix(ind_data_record)
    ratios.complete_matrix()
    ratios.calculate_levels()

    for n in range(1, len(ratios.levels)+1):
        key = "level_"+str(n)
        output[key].append(ratios.levels[n-1])

for key in output:
    plp_df[plp_inputs[key]] = np.array(output[key])



















