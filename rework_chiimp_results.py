import os
import pandas as pd
from openpyxl import Workbook
from typing import List, Tuple, Dict, Optional

from store import Store
from sequences import Sequence, check_primer, handle_wobble, compare_homo_hetero, follow_up_wobble
from utils import get_loci_numbers
from statistical_ana import do_basic_statistic
from write_to_excel import write_pc_to_excel, write_to_excel
from improve_data import get_imp_seq, set_stutter
from file_handler import get_file_list, get_working_data, get_primer_data, rename_peak_calling
from sorter import Sorter, DatabaseManager
from chiimp import Chiimp
from constants import *

store = Store()


def rework_chimp(db_path: str) -> Tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], pd.DataFrame, pd.DataFrame, List[str], List[Tuple[Sequence, Sequence]]
]:
    """
    Rework CHIIMP output by:
    - Sorting locus-wise output files
    - Checking primer binding and wobble positions
    - Applying stutter filters
    - Generating PeakCalling and FlagsPeakCalling tables
    - Writing results to Excel/CSV
    - Storing all intermediate metadata in SQLite to allow resume-after-crash

    Args:
        db_path: Path to SQLite database for storing processed sample metadata.

    Returns:
        Tuple of:
        - peak_calling_all (DataFrame)
        - flags_peak_calling_all (DataFrame)
        - loci_results (DataFrame)
        - wrong_loci (list of str)
        - id_results (DataFrame)
        - meta_data_id (DataFrame)
        - wrong_id (list of str)
        - wobble_list (list of tuples of Sequence pairs)
    """
    # Initialize database and sorter
    db = DatabaseManager(db_path)
    sorter = Sorter(db)

    chiimp = Chiimp()
    chiimp.handle_assam_folder(sorter.processed_samples_folder, ROOT_DATA)

    os.chdir(ROOT_RESULTS)

    # Load reverse primers and store number of loci
    reverse_primer: pd.DataFrame = pd.read_csv(LOCUS_ATTRS, header=0, dtype='string')
    store.N_LOCI: List[str] = reverse_primer.iloc[:, 0].to_list()

    # Initialize variables
    wobble_list: List[Tuple[Sequence, Sequence]] = []
    peak_calling_all: Optional[pd.DataFrame] = None
    flags_peak_calling_all: Optional[pd.DataFrame] = None

    # Iterate through datasets
    for data_set in LIST_RAW_DATA:
        data_path: str = os.path.join(ROOT_DATA, data_set, FOLDER_PS)

        try:
            folder_list: List[str] = os.listdir(data_path)
            folder_list = sorter.sort_by_natural(folder_list)
        except FileNotFoundError:
            raise FileNotFoundError("Check spelling of <processed-samples> folder name")

        file_list: List[str] = get_file_list(folder_list, data_path)

        # Create column names for peak calling
        column_names: List[str] = ['Sample', 'Run']
        for loci in store.N_LOCI:
            column_names.append(f"{loci}_1")
            column_names.append(f"{loci}_2")

        # Initialize PeakCalling and FlagsPeakCalling tables
        if peak_calling_all is None:
            peak_calling_all = pd.DataFrame(columns=column_names)
        peak_calling: pd.DataFrame = pd.DataFrame(
            [[Sequence.empty() for _ in range(len(column_names))] for _ in range(len(file_list))],
            columns=column_names
        )
        if flags_peak_calling_all is None:
            flags_peak_calling_all = pd.DataFrame(columns=column_names)
        flags_peak_calling: pd.DataFrame = pd.DataFrame(
            [[flagdict.copy() for _ in range(len(column_names))] for _ in range(len(file_list))],
            columns=column_names
        )

        peak_calling['Sample'] = file_list
        peak_calling['Run'] = data_set
        flags_peak_calling['Sample'] = file_list
        flags_peak_calling['Run'] = data_set

        # Iterate folders
        for folder in folder_list:
            primer_seq, wobble_len_min, wobble_len_max, can_wobble = get_primer_data(reverse_primer, folder)
            index_1, index_2 = f"{folder}_1", f"{folder}_2"

            # Iterate files
            for file_name in file_list:
                try:
                    working_data = get_working_data(data_path, folder, file_name)
                    if working_data.empty:
                        continue

                    flag_cell1, flag_cell2 = flagdict.copy(), flagdict.copy()
                    allel_list, prom_list, _ = get_imp_seq(working_data, file_name, folder, primer_seq)
                    seq_list: List[Sequence] = allel_list + prom_list

                    # Check primer binding
                    cell_1: Optional[Sequence] = check_primer(seq_list, flag_cell1)
                    if can_wobble and cell_1:
                        cell_1 = handle_wobble(cell_1, seq_list, wobble_list, primer_seq, flag_cell1)

                    cell_2: Optional[Sequence] = check_primer(seq_list, flag_cell2)
                    if can_wobble and cell_2:
                        cell_2 = handle_wobble(cell_2, seq_list, wobble_list, primer_seq, flag_cell2)

                    # Set stutter flags
                    set_stutter(cell_1, cell_2, flag_cell1, flag_cell2)

                    # Insert into PeakCalling
                    if cell_1 and cell_1.count >= COUNT_LIMIT:
                        flag_cell1['category'] = cell_1.category
                        peak_calling.at[file_list.index(file_name), index_1] = cell_1
                        flags_peak_calling.at[file_list.index(file_name), index_1] = flag_cell1

                    if cell_2 and cell_2.count >= COUNT_LIMIT:
                        flag_cell2['category'] = cell_2.category
                        peak_calling.at[file_list.index(file_name), index_2] = cell_2
                        flags_peak_calling.at[file_list.index(file_name), index_2] = flag_cell2

                    # Record sample into database for resume
                    db.insert_sample(
                        sample_name=file_name,
                        folder=os.path.join(data_path, folder),
                        file_path=os.path.join(data_path, folder, file_name),
                        locus=folder,
                        flags=str(flag_cell1)
                    )

                except FileNotFoundError:
                    continue

        # Write per-dataset Excel
        workbook = Workbook()
        write_pc_to_excel(peak_calling, flags_peak_calling, workbook, is_peak_calling_all=False)
        workbook.save(os.path.join(ROOT_RESULTS, f"PeakCalling_{data_set}.xlsx"))

        # Concatenate into the global DataFrame
        peak_calling_all = pd.concat([peak_calling_all, peak_calling])
        flags_peak_calling_all = pd.concat([flags_peak_calling_all, flags_peak_calling])

    # Final wobble handling
    follow_up_wobble(peak_calling_all, wobble_list, flags_peak_calling_all)

    # Sort by ID and compare homo/hetero
    peak_calling_all = sorter.sort_by_column(peak_calling_all, 'ID')
    flags_peak_calling_all = sorter.sort_by_column(flags_peak_calling_all, 'ID')
    compare_homo_hetero(peak_calling_all, flags_peak_calling_all)

    # Basic statistics
    loci_results, wrong_loci, id_results, meta_data_id, wrong_id = do_basic_statistic(peak_calling_all, flags_peak_calling_all)

    # Write final Excel / CSV
    write_to_excel(peak_calling_all, P_C_XLSX, flags_peak_calling_all, loci_results, id_results, is_peak_calling_all=True)
    peak_calling_all.to_csv(P_C_CSV, index=False, sep=';')
    flags_peak_calling_all.to_csv(F_P_C_CSV, index=False, sep=';')
    peak_calling_all = rename_peak_calling(peak_calling_all)

    write_to_excel(peak_calling_all, P_C_RENAME_XLSX, flags_peak_calling_all, loci_results, id_results, is_peak_calling_all=True)
    peak_calling_all.to_csv(P_C_CSV, index=False, sep=';')

    return peak_calling_all, flags_peak_calling_all, loci_results, wrong_loci, id_results, meta_data_id, wrong_id, wobble_list
