import sqlite3
from pathlib import Path
import time


class Store:
    """
    Singleton store using SQLite for local pipeline tracking.
    Stores project config, file status, steps, and quality metrics.
    """

    store_instance = None

    def __new__(cls, db_path=None):
        if cls.store_instance is None:
            cls.store_instance = super(Store, cls).__new__(cls)

            # Initialize DB connection
            cls.store_instance.db_path = db_path or str(Path.cwd() / "pipeline.db")
            cls.store_instance.conn = sqlite3.connect(cls.store_instance.db_path)
            cls.store_instance._init_db()
        return cls.store_instance

    def _init_db(self):
        """Create tables if they don't exist"""
        c = self.conn.cursor()

        # Project config
        c.execute("""
        CREATE TABLE IF NOT EXISTS project (
            key TEXT PRIMARY KEY,
            value TEXT
        )
        """)

        # Files table
        c.execute("""
        CREATE TABLE IF NOT EXISTS files (
            path TEXT PRIMARY KEY,
            step TEXT,
            status TEXT,
            quality REAL,
            updated REAL,
            error TEXT
        )
        """)

        # Steps table (optional)
        c.execute("""
        CREATE TABLE IF NOT EXISTS steps (
            name TEXT PRIMARY KEY,
            completed INTEGER,
            updated REAL
        )
        """)

        self.conn.commit()

    # ----------------------
    # Project config methods
    # ----------------------
    def set_project_value(self, key, value):
        self.conn.execute(
            "INSERT OR REPLACE INTO project (key, value) VALUES (?, ?)",
            (key, str(value))
        )
        self.conn.commit()

    def get_project_value(self, key, default=None):
        c = self.conn.cursor()
        c.execute("SELECT value FROM project WHERE key=?", (key,))
        row = c.fetchone()
        return row[0] if row else default

    # ----------------------
    # File tracking methods
    # ----------------------
    def add_file(self, path, step=None):
        self.conn.execute(
            "INSERT OR IGNORE INTO files (path, step, status, updated) VALUES (?, ?, ?, ?)",
            (path, step, "pending", time.time())
        )
        self.conn.commit()

    def get_pending_files(self, step=None):
        c = self.conn.cursor()
        if step:
            c.execute("SELECT path FROM files WHERE status != 'done' AND step=?", (step,))
        else:
            c.execute("SELECT path FROM files WHERE status != 'done'")
        return [row[0] for row in c.fetchall()]

    def mark_done(self, path, step=None):
        self.conn.execute(
            "UPDATE files SET status='done', step=?, updated=? WHERE path=?",
            (step, time.time(), path)
        )
        self.conn.commit()

    def update_quality(self, path, quality):
        self.conn.execute(
            "UPDATE files SET quality=?, updated=? WHERE path=?",
            (quality, time.time(), path)
        )
        self.conn.commit()
