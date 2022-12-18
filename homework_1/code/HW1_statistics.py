import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest


def check_dge(first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table):
    df_first = pd.read_csv(first_cell_type_expressions_path, index_col=0)
    df_second = pd.read_csv(second_cell_type_expressions_path, index_col=0)
    ci_test_results = check_dge_with_ci(df_first, df_second)
    z_test_results, z_test_p_values = check_dge_with_ztest(df_first, df_second)
    mean_dif = mean_diff(df_first, df_second)
    selected_genes = genes(df_first, df_second)

    results = {
        "genes": selected_genes,
        "ci_test_results": ci_test_results,
        "z_test_results": z_test_results,
        "z_test_p_values": z_test_p_values,
        "mean_diff": mean_dif
    }

    results = pd.DataFrame(results)
    results.to_csv(f"{save_results_table}.csv")

    return


def genes(first_table, second_table):
    columns_first = list(first_table)
    columns_second = list(second_table)
    columns = []
    for i in columns_first:
        if i in columns_second:
            columns.append(i)

    return columns


def mean_diff(first_table, second_table):
    mean_dif = []
    for i in genes(first_table, second_table):
        mean_dif.append(np.mean(first_table[i]) - np.mean(second_table[i]))

    return mean_dif


def check_dge_with_ci(first_table, second_table):
    ci_test_results = []
    for i in genes(first_table, second_table):
        ci_first = st.t.interval(alpha=0.95, df=len(first_table[i]) - 1,
                                 loc=np.mean(first_table[i]),
                                 scale=st.sem(first_table[i]))
        ci_second = st.t.interval(alpha=0.95, df=len(second_table[i]) - 1,
                                  loc=np.mean(second_table[i]),
                                  scale=st.sem(second_table[i]))
        ci_test_results.append(check_intervals_intersect(ci_first, ci_second))

    return ci_test_results


def check_intervals_intersect(first_ci, second_ci):
    if first_ci[0] > second_ci[1] or first_ci[1] < second_ci[0]:
        are_intersect = False
    else:
        are_intersect = True
    return are_intersect  # True or False


def check_dge_with_ztest(first_table, second_table):
    z_test_results = []
    p_values = []
    for i in genes(first_table, second_table):
        z_test = ztest(first_table[i], second_table[i])
        p_values.append(z_test[1])
        if z_test[1] < 0.05:
            z_test_results.append(True)

        else:
            z_test_results.append(False)

    return z_test_results, p_values
