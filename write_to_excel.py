import pandas as pd
import numpy as np
from xlsxwriter import Workbook
from sequences import Sequence
from constants import *
from utils import get_start_column
from typing import Optional


def write_to_excel(
    peak_calling: pd.DataFrame,
    pfad_excel: str,
    flags_peak_calling: pd.DataFrame,
    loci_results: pd.DataFrame,
    id_results: pd.DataFrame,
    is_peak_calling_all: bool
) -> None:
    """
    Writes the PeakCalling DataFrame and FlagsPeakCalling DataFrame to an Excel file,
    including the basic statistics worksheet. Also applies cell coloring based on flags.

    Args:
        peak_calling: DataFrame containing all called sequences.
        pfad_excel: File path to save the Excel workbook.
        flags_peak_calling: DataFrame of flags corresponding to each cell in peak_calling.
        loci_results: DataFrame containing loci statistics.
        id_results: DataFrame containing ID-level statistics.
        is_peak_calling_all: If True, generates additional worksheets for homo/hetero and All sequences.
    """
    workbook = Workbook(pfad_excel)
    write_pc_to_excel(peak_calling, flags_peak_calling, workbook, is_peak_calling_all)
    write_stat_to_excel(workbook, loci_results, id_results, 'Basic_Statistics')
    workbook.close()


def write_stat_to_excel(
    workbook: Workbook,
    loci_results: pd.DataFrame,
    id_results: pd.DataFrame,
    worksheet_name: str
) -> None:
    """
    Writes loci and ID statistics to a new worksheet in the Excel workbook, including conditional formatting.

    Args:
        workbook: xlsxwriter Workbook object.
        loci_results: DataFrame of loci-level statistics.
        id_results: DataFrame of ID-level statistics.
        worksheet_name: Name of the worksheet for statistics.
    """
    header_format = workbook.add_format({'bold': True})
    cell_format = workbook.add_format({'num_format': 10})
    cell_bg = workbook.add_format({'bg_color': BASIC_STAT_COLOR, 'num_format': 10})

    worksheet_stat = workbook.add_worksheet(worksheet_name)

    # Write loci header
    worksheet_stat.write(1, 0, 'Loci')
    worksheet_stat.write_row(1, 1, loci_results.columns, header_format)

    # Write loci data
    count = 2
    for name, data in loci_results.iterrows():
        worksheet_stat.write_row(count, 0, [name] + data.values.tolist(), cell_format)
        count += 1

    # Conditional formatting for loci
    worksheet_stat.conditional_format(2, 1, count, len(loci_results.columns),
                                      {'type': 'cell', 'criteria': '>=',
                                       'value': BASIC_STAT_PERCENTAGE, 'format': cell_bg})

    # Write ID-level statistics
    count += 2
    worksheet_stat.write(count, 0, 'ID/Name')
    id_results_clean = id_results[id_results.columns.drop(list(id_results.filter(regex='nan')))]
    worksheet_stat.write_row(count, 1, id_results_clean.columns, header_format)
    count += 1

    first_row_hh = count
    for name, data in id_results_clean.iterrows():
        worksheet_stat.write_row(count, 0, [name] + data.values.tolist(), cell_format)
        count += 1

    # Freeze first row and format columns
    worksheet_stat.freeze_panes(0, 1)
    worksheet_stat.set_column(0, 0, 55, header_format)

    # Conditional formatting for ID statistics
    worksheet_stat.conditional_format(first_row_hh, 1, count, len(id_results_clean.columns),
                                      {'type': 'cell', 'criteria': '>=', 'value': 0.25, 'format': cell_bg})


def write_pc_to_excel(
    peak_calling: pd.DataFrame,
    flags_peak_calling: pd.DataFrame,
    workbook: Workbook,
    is_peak_calling_all: bool
) -> None:
    """
    Writes the PeakCalling sequences and corresponding flags to Excel worksheets,
    including optional homo/hetero and All sequences sheets.

    Args:
        peak_calling: DataFrame with all sequences.
        flags_peak_calling: DataFrame with flags per cell.
        workbook: xlsxwriter Workbook object.
        is_peak_calling_all: If True, creates additional worksheets for all sequences and homo/hetero.
    """
    header_format = workbook.add_format({'bold': True})

    # Main sequence sheet
    worksheet = workbook.add_worksheet('Sequences')
    worksheet.write_row('A1', peak_calling.columns, header_format)

    # Error/flag sheet
    worksheet_cellcolor = workbook.add_worksheet('Errors_in_CHIIMP')
    worksheet_cellcolor.write_row('A1', peak_calling.columns, header_format)

    start_column: int = get_start_column(peak_calling.columns)

    # Optional additional sheets
    if is_peak_calling_all:
        worksheet_bordercolor = workbook.add_worksheet('Homo_Hetero')
        worksheet_bordercolor.write_row('A1', peak_calling.columns, header_format)

        worksheet_all = workbook.add_worksheet('All')
        worksheet_all.write_row('A1', peak_calling.columns, header_format)

    # Iterate rows and columns
    for index_row in range(len(peak_calling.index)):
        for index_col in range(len(peak_calling.columns)):
            cell_value: str
            if isinstance(peak_calling.iloc[index_row, index_col], Sequence):
                cell_value = peak_calling.iloc[index_row, index_col].seqname
            else:
                cell_value = peak_calling.iloc[index_row, index_col]

            cell_format: dict = {'bg': {}, 'border': {}, 'font': {}}

            # Only apply flags after loci start column
            if index_col >= start_column:
                flags: dict = flags_peak_calling.iloc[index_row, index_col]

                if cell_value == '' or cell_value is np.nan:
                    cell_value = ''
                    cell_format['bg']['bg_color'] = EMPTY_CELL
                else:
                    # Background colors based on flags
                    if flags['is_wobble'] and flags['mod_data'] and flags['is_stutter']:
                        cell_format['bg']['bg_color'] = IMP_DATA_WOBBLE_STUTTER_CELL
                    elif flags['is_wobble'] and flags['is_stutter']:
                        cell_format['bg']['bg_color'] = WOBBLE_STUTTER_CELL
                    elif flags['is_wobble'] and flags['mod_data']:
                        cell_format['bg']['bg_color'] = IMP_DATA_WOBBLE_CELL
                    elif flags['is_stutter']:
                        cell_format['bg']['bg_color'] = STUTTER_CELL
                    elif flags['is_wobble']:
                        cell_format['bg']['bg_color'] = WOBBLE_CELL
                    elif flags['mod_data']:
                        cell_format['bg']['bg_color'] = IMP_DATA_CELL

                    # Font styling
                    if flags['non_primer']:
                        cell_format['font'] = {'font_color': PRIMER_NOT_BOUND_CELL,
                                               'bold': True, 'italic': True, 'font_size': 12}
                    if flags['category'] == 'Allele':
                        cell_format['font']['underline'] = 1
                    if is_peak_calling_all:
                        if flags.get('is_homo'):
                            cell_format['border'] = {'border_color': HOMOZYGOT_CELL, 'border': 5}
                        elif flags.get('is_hetero'):
                            cell_format['border'] = {'border_color': HETEROZYGOT_CELL, 'border': 5}
                        elif flags.get('is_wrong'):
                            cell_format['border'] = {'border_color': NO_HOMO_OR_HETERO, 'border': 5}

            # Write to worksheets
            format_bg = workbook.add_format(cell_format['bg'] | cell_format['font'])
            worksheet.write(index_row + 1, index_col, cell_value)
            worksheet_cellcolor.write(index_row + 1, index_col, cell_value, format_bg)

            if is_peak_calling_all:
                format_border = workbook.add_format(cell_format['border'])
                format_all = workbook.add_format(cell_format['bg'] | cell_format['font'] | cell_format['border'])
                worksheet_bordercolor.write(index_row + 1, index_col, cell_value, format_border)
                worksheet_all.write(index_row + 1, index_col, cell_value, format_all)

    # Add legend for color coding
    add_colour_legend(workbook, is_peak_calling_all)


def add_colour_legend(workbook: Workbook, is_peak_calling_all: bool) -> None:
    """
    Creates a legend worksheet describing the meaning of cell colors and formatting.

    Args:
        workbook: xlsxwriter Workbook object.
        is_peak_calling_all: If True, includes homo/hetero formatting explanations.
    """
    worksheet = workbook.add_worksheet('Colour_Legend')

    # Color legend for different cell types
    color_mapping = [
        (EMPTY_CELL, 'Count too low'),
        (WOBBLE_CELL, 'Wobble'),
        (STUTTER_CELL, 'Stutter'),
        (IMP_DATA_CELL, 'Improved data'),
        (WOBBLE_STUTTER_CELL, 'Wobble and Stutter'),
        (IMP_DATA_WOBBLE_CELL, 'Improved data and Wobble'),
        (IMP_DATA_WOBBLE_STUTTER_CELL, 'Improved data, Wobble and Stutter')
    ]

    row = 1
    for color, desc in color_mapping:
        fmt = workbook.add_format({'bg_color': color})
        worksheet.write(row, 0, '', fmt)
        worksheet.write(row, 1, desc)
        row += 1

    if is_peak_calling_all:
        extra_border = [
            (HOMOZYGOT_CELL, 'Homozygot Cell'),
            (HETEROZYGOT_CELL, 'Heterozygot Cell'),
            (NO_HOMO_OR_HETERO, 'Neither homo- or heterozygot Cell')
        ]
        for color, desc in extra_border:
            fmt = workbook.add_format({'border_color': color, 'border': 5})
            worksheet.write(row, 0, '', fmt)
            worksheet.write(row, 1, desc)
            row += 1

    worksheet.autofit()
