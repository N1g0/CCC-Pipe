import pandas as pd
from typing import Set


def get_start_column(columns: pd.Index) -> int:
    """
    Returns the index of the first column that ends with '_1'.
    This is used to determine where the loci columns start
    in a DataFrame (e.g., PeakCalling data).

    Args:
        columns: pandas Index object containing column names.

    Returns:
        Index (int) of the first column ending with '_1'.

    Raises:
        IndexError: If no column ends with '_1'.
    """
    start_cols = [i for i, item in enumerate(columns) if item.endswith('_1')]
    if not start_cols:
        raise IndexError("No column ending with '_1' found in columns.")
    return start_cols[0]


def get_loci_numbers(columns: pd.Index) -> Set[int]:
    """
    Returns a set of unique locus numbers from the columns of a DataFrame.
    Columns for loci are expected to have names like '1_1', '1_2', '2_1', etc.

    Args:
        columns: pandas Index object containing column names.

    Returns:
        Set of integers representing unique loci numbers.
    """
    start_idx = get_start_column(columns)
    loci_cols = columns[start_idx:]
    return {int(col.split("_")[0]) for col in loci_cols}
