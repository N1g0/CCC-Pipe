import os
import re
import shutil
import sqlite3
from typing import List, Optional, Any
import pandas as pd


class DatabaseManager:
    """
    SQLite manager to store processed samples, file paths, locus info, and peak calling flags.
    """

    def __init__(self, db_path: str):
        """
        Initialize the database connection and create required tables.

        Args:
            db_path: Path to SQLite database file.
        """
        self.conn = sqlite3.connect(db_path)
        self.create_tables()

    def create_tables(self) -> None:
        """Create tables if they do not exist."""
        cur = self.conn.cursor()
        cur.execute('''
        CREATE TABLE IF NOT EXISTS processed_samples (
            sample_name TEXT,
            folder TEXT,
            file_path TEXT,
            locus TEXT,
            peak_calling_flags TEXT,
            UNIQUE(sample_name, folder, locus)
        )
        ''')
        self.conn.commit()

    def insert_sample(
        self,
        sample_name: str,
        folder: str,
        file_path: str,
        locus: Optional[str] = None,
        flags: Optional[str] = None
    ) -> None:
        """
        Insert a processed sample record into the database.

        Args:
            sample_name: Name of the sample.
            folder: Folder where file is stored.
            file_path: Full file path.
            locus: Optional locus identifier.
            flags: Optional flags (JSON string recommended).
        """
        cur = self.conn.cursor()
        cur.execute('''
        INSERT OR IGNORE INTO processed_samples
        (sample_name, folder, file_path, locus, peak_calling_flags)
        VALUES (?, ?, ?, ?, ?)
        ''', (sample_name, folder, file_path, locus, flags))
        self.conn.commit()

    def get_unprocessed_samples(self, folder: str) -> List[str]:
        """
        Get a list of file paths in a folder that have not been processed yet.

        Args:
            folder: Folder to query.

        Returns:
            List of unprocessed file paths.
        """
        cur = self.conn.cursor()
        cur.execute('''
        SELECT file_path FROM processed_samples
        WHERE folder=? AND peak_calling_flags IS NULL
        ''', (folder,))
        rows = cur.fetchall()
        return [row[0] for row in rows]


class Sorter:
    """
    Provides sorting utilities and folder management for processed sample files.
    """

    def __init__(self, db: DatabaseManager):
        """
        Initialize Sorter with a database manager.

        Args:
            db: DatabaseManager instance for inserting file metadata.
        """
        self.db = db

    def sort_by_natural(self, liste: List[str]) -> List[str]:
        """
        Sorts a list of strings using natural order (e.g., 'file2' before 'file10').

        Args:
            liste: List of strings to sort.

        Returns:
            Naturally sorted list of strings.
        """
        liste.sort(key=Sorter.natural_keys)
        return liste

    @staticmethod
    def atof(text: str) -> Any:
        """Convert a string to float if possible, else return string."""
        try:
            return float(text)
        except ValueError:
            return text

    @staticmethod
    def natural_keys(text: str) -> List[Any]:
        """Generate sorting keys for natural sorting."""
        return [Sorter.atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text)]

    def sort_by_column(self, df: pd.DataFrame, name: str) -> pd.DataFrame:
        """Sorts a DataFrame by a column."""
        df.sort_values(by=[name], inplace=True)
        return df

    def sort_by_numbers(self, liste: List[str]) -> None:
        """Sort a list of strings based on numeric digits in the string."""
        num = lambda a: int(re.search(r'\d+', a).group())
        liste.sort(key=num)

    def processed_samples_folder(self, auswahl_folder: str, folder_list_split: List[str]) -> None:
        """
        Create processed-sample folder structure, copy files, and insert records into DB.

        Args:
            auswahl_folder: Base folder containing raw data.
            folder_list_split: Split folder path to extract Assam folder name.
        """
        assam_path: str = os.path.join(auswahl_folder, folder_list_split[-1])
        p_s_path: str = os.path.join(auswahl_folder, 'processed-samples')
        p_s_start: str = os.path.join(assam_path, 'processed-samples')

        os.makedirs(p_s_path, exist_ok=True)
        file_list: List[str] = os.listdir(p_s_start)

        for file in file_list:
            folder_number: str = file.split('-')[-1].split('.')[0]
            check_folder: str = os.path.join(p_s_path, folder_number)
            os.makedirs(check_folder, exist_ok=True)

            dest_file: str = os.path.join(check_folder, file)
            shutil.copyfile(os.path.join(p_s_start, file), dest_file)

            # Insert file record into database
            self.db.insert_sample(sample_name=file, folder=check_folder, file_path=dest_file)


class Filter:
    """
    Methods to filter peak calling data for CERVUS compatibility
    and detect duplicate names.
    """

    def __init__(self, db: DatabaseManager):
        self.db = db

    def filter_cervus(self, folder: str) -> pd.DataFrame:
        """
        Fetch processed samples from DB and filter for homozygous / heterozygous loci.

        Args:
            folder: Folder name to query processed samples.

        Returns:
            DataFrame formatted for CERVUS.
        """
        cur = self.db.conn.cursor()
        cur.execute('SELECT * FROM processed_samples WHERE folder=?', (folder,))
        records = cur.fetchall()
        df: pd.DataFrame = pd.DataFrame(records, columns=['sample_name', 'folder', 'file_path', 'locus', 'flags'])

        # Filtering logic can be added here using df
        return df

    def double_name(self, df: pd.DataFrame) -> None:
        """
        Check for duplicate sample names in the DataFrame.

        Args:
            df: DataFrame containing sample data.

        Returns:
            None
        """
        for i, row_1 in df.iterrows():
            if not isinstance(row_1['sample_name'], str):
                continue
            for j, row_2 in df.iloc[i + 1:].iterrows():
                if row_1[4:].equals(row_2[4:]):
                    print('Double name in dataset')
                    print(f'{row_1["sample_name"]} row: {i}')
                    print(f'{row_2["sample_name"]} row: {j}')
