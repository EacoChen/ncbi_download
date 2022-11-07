# -*- coding: utf-8 -*-
#!/usr/bin/env python3


class DownloadJob(object):
    """Collection of data for a download job."""

    __slots__ = ['full_url', 'filename', 'output']

    def __init__(self, full_url, filename, output):
        """Initialise the download job."""
        self.full_url = full_url
        self.filename = filename
        self.output = output

    def __eq__(self, other):
        """Check for equality."""
        if not isinstance(other, DownloadJob):
            return False

        return self.full_url == other.full_url and \
            self.filename == other.filename and \
            self.output == other.output
