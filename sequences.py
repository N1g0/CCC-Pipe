from __future__ import annotations
import ast
from constants import *
from utils import get_loci_numbers
import pandas as pd


class Sequence:
    """
    Represents a DNA sequence for a sample, along with metadata such as file, folder,
    count, length, category, sequence name, and primer binding status.
    """

    def __init__(self, filename: str, folder: str, seq: str, count: int, length: int,
                 category: str, seqname: str, primer_seq: str | bool) -> None:
        """
        Initializes a Sequence object.

        Args:
            filename: Name of the file containing the sequence.
            folder: Folder name where the file is located.
            seq: DNA sequence string.
            count: Number of times this sequence was observed.
            length: Length of the sequence.
            category: Category label (e.g., allele type).
            seqname: Name identifier of the sequence.
            primer_seq: Either a string with primer sequence to check or a boolean
                        indicating primer binding status.
        """
        self.filename: str = filename
        self.folder: str = folder
        self.seq: str = seq
        self.count: int = count
        self.length: int = length
        self.category: str = category
        self.seqname: str = seqname

        if isinstance(primer_seq, bool):
            self.is_primer_bound: bool = primer_seq
        else:
            self.is_primer_bound: bool = self.__is_primer_bound(primer_seq)

    def __repr__(self) -> str:
        """Returns a readable string representation of the Sequence object."""
        return str([self.filename, self.folder, self.seq, self.count, self.length,
                    self.category, self.seqname, self.is_primer_bound])

    def __eq__(self, other: Sequence) -> bool:
        """
        Checks if two Sequence objects are equal based on the DNA sequence string.

        Args:
            other: Another Sequence object.

        Returns:
            True if sequences are identical, else False.
        """
        return self.seq == other.seq

    def print_imp(self) -> str:
        """Returns a compact string showing the sequence name and category."""
        return self.seqname.ljust(15) + '    ' + self.category.ljust(15)

    def print_all(self) -> str:
        """Returns a detailed string showing count, length, category, and sequence name."""
        return str(self.count).ljust(5) + '  ' + str(self.length).ljust(6) + '  ' + \
               self.category.ljust(9) + ' ' + self.seqname.ljust(10)

    def print_route(self) -> str:
        """Returns a string showing folder number and filename."""
        return 'Folder Nr: ' + str(self.folder).ljust(3) + '  File: ' + str(self.filename).ljust(6)

    def print_route2(self) -> str:
        """Returns a string showing folder name and filename."""
        return 'Folder: ' + str(self.folder) + ', File: ' + str(self.filename)

    def __is_primer_bound(self, primer_seq: str) -> bool:
        """
        Determines if the sequence matches the primer sequence with at most 1 mismatch.

        Args:
            primer_seq: Primer sequence string to compare.

        Returns:
            True if primer binding is valid (<=1 mismatch), else False.
        """
        primer_len: int = len(primer_seq)
        miss_count: int = 0
        for seq_pos, pri_pos in zip(self.seq[-primer_len:], primer_seq):
            if seq_pos != pri_pos:
                miss_count += 1
                if miss_count >= 2:
                    return False
        return True

    @staticmethod
    def seq_from_list(str_list: pd.Series) -> pd.Series:
        """
        Converts a pandas Series of strings representing Sequences into Sequence objects.

        Args:
            str_list: Pandas Series containing string representations of Sequence objects.

        Returns:
            Pandas Series of Sequence objects.
        """
        try:
            return str_list.apply(Sequence.ast_eval)
        except Exception:
            return str_list

    @staticmethod
    def ast_eval(string: str) -> Sequence | str:
        """
        Safely evaluates a string representation of a Sequence object and returns a Sequence.

        Args:
            string: String representation of a Sequence.

        Returns:
            Sequence object if evaluation succeeds, else returns the original string.
        """
        try:
            liste: list = ast.literal_eval(string)
            return Sequence(liste[0], liste[1], liste[2], liste[3], liste[4],
                            liste[5], liste[6], liste[7])
        except Exception:
            return string

    @staticmethod
    def empty() -> Sequence:
        """
        Returns an empty Sequence object with default values.

        Returns:
            Sequence object with empty strings and zeros.
        """
        return Sequence('', '', '', 0, 0, '', '', '')

    @staticmethod
    def is_wobble(all_1: Sequence, all_2: Sequence, primer_seq: str) -> bool:
        """
        Determines if two sequences form a wobble pair, ignoring primer region.

        Args:
            all_1: First Sequence object.
            all_2: Second Sequence object.
            primer_seq: Primer sequence string.

        Returns:
            True if wobble detected, else False.
        """
        if all_1.length == all_2.length:
            mutation_count: int = 0
            len_no_primer: int = len(all_1.seq) - len(primer_seq)
            for point_mutation_1, point_mutation_2 in zip(all_1.seq[:len_no_primer], all_2.seq[:len_no_primer]):
                if point_mutation_1 != point_mutation_2:
                    mutation_count += 1
            if mutation_count == 0:
                for i, wobble_position in enumerate(primer_seq):
                    if wobble_position in WOBBLE_POSITION:
                        wp_all_1: str = all_1.seq[len_no_primer + i]
                        wp_all_2: str = all_2.seq[len_no_primer + i]
                        if wp_all_1 != wp_all_2 and (wp_all_1 in WOBBLE_POSITION[wobble_position] and wp_all_2 in WOBBLE_POSITION[wobble_position]):
                            return True
        return False


def handle_wobble(cell_1: Sequence, seq_list: list[Sequence],
                  wobble_list: list[tuple[Sequence, Sequence]], primer_seq: str,
                  flagdict: dict[str, object]) -> Sequence:
    """
    Checks for wobble pairs in a list of sequences and updates flags/wobble list.

    Args:
        cell_1: Sequence object currently being assigned to a PeakCalling cell.
        seq_list: List of potential Sequence objects for the same cell.
        wobble_list: List of tuples of Sequence pairs already identified as wobble.
        primer_seq: Primer sequence for the locus.
        flagdict: Dictionary storing flags for the current cell.

    Returns:
        Updated Sequence object for the cell.
    """
    for all_2 in seq_list:
        if Sequence.is_wobble(cell_1, all_2, primer_seq):
            flagdict['is_wobble'] = True
            seq_list.remove(all_2)
            if '-blank' in cell_1.seqname:
                cell_1, all_2 = all_2, cell_1
            cell_1.count += all_2.count
            for duplicate in wobble_list:
                if cell_1 == duplicate[0]:
                    return cell_1
                elif Sequence.is_wobble(cell_1, duplicate[0], primer_seq):
                    if '-blank' in duplicate[0].seqname:
                        duplicate[0].seqname = cell_1.seqname
                    else:
                        cell_1.seqname = duplicate[0].seqname
                    return cell_1
            wobble_list.append((cell_1, all_2))
            return cell_1
    return cell_1


def follow_up_wobble(peak_calling_all: pd.DataFrame,
                     wobble_list: list[tuple[Sequence, Sequence]],
                     flag_peak_calling_all: pd.DataFrame) -> None:
    """
    Updates sequence names in PeakCalling DataFrame after wobble detection.

    Args:
        peak_calling_all: Pandas DataFrame of all peak calls.
        wobble_list: List of wobble pairs.
        flag_peak_calling_all: DataFrame of flags for each peak calling cell.
    """
    peak_calling_all.reset_index(drop=True, inplace=True)
    flag_peak_calling_all.reset_index(drop=True, inplace=True)
    for i, _ in peak_calling_all.iterrows():
        for j in range(4, len(peak_calling_all.columns)):
            if flag_peak_calling_all.iloc[i, j]['is_wobble']:
                seq_obj: Sequence = peak_calling_all.iloc[i, j]
                for wobble_1, wobble_2 in wobble_list:
                    if wobble_1 == seq_obj or wobble_2 == seq_obj:
                        seq_obj.seqname = wobble_1.seqname


def check_primer(seq_list: list[Sequence], flagdict_cell: dict[str, object]) -> Sequence | None:
    """
    Finds a sequence with correct primer binding from a list and updates flags.

    Args:
        seq_list: List of Sequence objects for a cell.
        flagdict_cell: Dictionary of flags for the cell.

    Returns:
        The Sequence with primer bound, or None if none found.
    """
    cell: Sequence | None = None
    copy_seq_list: list[Sequence] = seq_list.copy()
    i: int = 0
    while i < len(copy_seq_list):
        if copy_seq_list[i].is_primer_bound:
            cell = copy_seq_list[i]
            break
        else:
            flagdict_cell['non_primer'] = True
            seq_list.remove(copy_seq_list[i])
            i += 1
    if cell is not None:
        seq_list.remove(cell)
    return cell


def compare_homo_hetero(peak_calling_all: pd.DataFrame, flags_peak_calling_all: pd.DataFrame) -> None:
    """
    Compares alleles across samples to determine homozygous, heterozygous, or incorrect calls.
    Updates the flags DataFrame accordingly.

    Args:
        peak_calling_all: DataFrame of peak calling results.
        flags_peak_calling_all: DataFrame of flags for each peak calling cell.
    """
    peak_calling_all = peak_calling_all.reset_index(drop=True)
    start_column: int = peak_calling_all.columns.tolist().index('ID')
    loci_numbers: list[int] = get_loci_numbers(peak_calling_all.columns)
    for i, _ in peak_calling_all.iloc[:-1].iterrows():
        if peak_calling_all.iloc[i, start_column] != 'nan' and peak_calling_all.iloc[i, start_column] == peak_calling_all.iloc[i + 1, start_column]:
            for j in loci_numbers:
                flag: str = ''
                l_1: str = str(j) + '_1'
                l_2: str = str(j) + '_2'
                first_column: pd.Series = peak_calling_all[l_1]
                second_column: pd.Series = peak_calling_all[l_2]
                A: str = first_column[i].seqname
                B: str = second_column[i].seqname
                C: str = first_column[i + 1].seqname
                D: str = second_column[i + 1].seqname
                if A != '' and B != '' and C != '' and D != '':
                    if A == B and A == C and A == D:
                        flag = 'is_homo'
                    elif (A == B and C == D) or (A == C and B == D) or (A == D and B == C):
                        flag = 'is_hetero'
                    else:
                        flag = 'is_wrong'
                    flags_peak_calling_all.iloc[i].loc[l_1][flag] = True
                    flags_peak_calling_all.iloc[i].loc[l_2][flag] = True
                    flags_peak_calling_all.iloc[i + 1].loc[l_1][flag] = True
                    flags_peak_calling_all.iloc[i + 1].loc[l_2][flag] = True
