"""
Python library for parsing FastQC results.


Provides methods and classes for parsing a FastQC run results.

Tested with FastQC 0.9.5 output, should work with other versions as well.

Requirements:
    * Python == 2.7.x

Copyright (c) 2013 Wibowo Arindrarto <w.arindrarto@lumc.nl>
Copyright (c) 2013 LUMC Sequencing Analysis Support Core <sasc@lumc.nl>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import os


def load_file(fname, **kwargs):
    """Given a path to a FastQC data file or an open file object pointing to it,
    return a `FastQC` object.

    :param fname: path to FastQC results data (usually `fastqc_data.txt`) or an
                  open handle of a FastQC results data
    :type fname: str or file handle
    :returns: an instance of `FastQC`
    :rtype: `FastQC`

    Keyword arguments are used only if fname is a string. They are passed to the
    file opener.

    """
    if isinstance(fname, basestring):
        with open(fname, **kwargs) as fp:
            return FastQC(fp)
    else:
        return FastQC(fname)


def load_from_dir(dirname, data_fname='fastqc_data.txt', **kwargs):
    """Given a path to a FastQC results directory, return a `FastQC` object.
    
    :param dirname: path to the top directory containing all FastQC results.
    :type dirname: str
    :param data_fname: file name of the FastQC results data
    :type data_fname: str
    :returns: an instance of `FastQC`
    :rtype: `FastQC`

    Keyword arguments are passed to the file opener.

    """
    assert os.path.exists(dirname), "Directory %r does not exist" % dirname
    fqc_path = os.path.join(os.walk(dirname).next()[1][0], data_fname)
    return load_file(fqc_path, **kwargs)


class FastQCModule(object):

    """Class representing a FastQC analysis module."""

    def __init__(self, raw_lines, end_mark='>>END_MODULE'):
        """

        :param raw_lines: list of lines in the module
        :type raw_lines: list of str
        :param end_mark: mark of the end of the module
        :type end_mark: str

        """
        self.raw_lines = raw_lines
        self.end_mark = end_mark
        self._status = None
        self._name = None
        self._data = self._parse()

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,
                '[%r, ...]' % self.raw_lines[0])

    def __str__(self):
        return ''.join(self.raw_lines)

    @property
    def name(self):
        """Name of the module."""
        return self._name

    @property
    def columns(self):
        """Columns in the module."""
        return self._columns

    @property
    def data(self):
        """FastQC data."""
        return self._data

    @property
    def status(self):
        """FastQC run status."""
        return self._status

    def _parse(self):
        """Common parser for a FastQC module."""
        # check that the last line is a proper end mark
        assert self.raw_lines[-1].startswith(self.end_mark)
        # parse name and status from first line
        tokens = self.raw_lines[0].strip().split('\t')
        name = tokens[0][2:]
        self._name = name
        status = tokens[-1]
        assert status in ('pass', 'fail', 'warn'), "Unknown module status: %r" \
            % status
        self._status = status
        # and column names from second line
        columns = self.raw_lines[1][1:].strip().split('\t')
        self._columns = columns
        # the rest of the lines except the last one
        data = []
        for line in self.raw_lines[2:-1]:
            cols = line.strip().split('\t')
            data.append(cols)

        # optional processing for different modules
        if self.name == 'Basic Statistics':
            data = {k: v for k, v in data}

        return data


class FastQC(object):

    """Class representing results from a FastQC run."""

    # module name -- attribute name mapping
    _mod_map = {
        '>>Basic Statistics': 'basic_statistics',
        '>>Per base sequence quality': 'per_base_sequence_quality',
        '>>Per sequence quality scores': 'per_sequence_quality_scores',
        '>>Per base sequence content': 'per_base_sequence_content',
        '>>Per base GC content': 'per_base_gc_content',
        '>>Per sequence GC content': 'per_sequence_gc_content',
        '>>Per base N content': 'per_base_n_content',
        '>>Sequence Length Distribution': 'sequence_length_distribution',
        '>>Sequence Duplication Levels': 'sequence_duplication_levels',
        '>>Overrepresented sequences': 'overrepresented_sequences',
        '>>Kmer content': 'kmer_content',
    }

    def __init__(self, fp):
        """

        :param fp: open file handle pointing to the FastQC data file
        :type fp: file handle

        """
        # get file name
        self.fname = fp.name
        self._modules = {}

        line = fp.readline()
        while True:

            tokens = line.strip().split('\t')
            # break on EOF
            if not line:
                break
            # parse version
            elif line.startswith('##FastQC'):
                self.version = line.strip().split()[1]
            # parse individual modules
            elif tokens[0] in self._mod_map:
                attr = self._mod_map[tokens[0]]
                raw_lines = self._read_module(fp, line, tokens[0])
                self._modules[attr] = FastQCModule(raw_lines)

            line = fp.readline()

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.fname)

    def _filter_by_status(self, status):
        """Filter out modules whose status is different from the given status.

        :param status: module status
        :type status: str
        :returns: a list of FastQC module names with the given status
        :rtype: list of str

        """
        return [x.name for x in self._modules.values() if x.status == status]

    def _read_module(self, fp, line, start_mark):
        """Returns a list of lines in a module.

        :param fp: open file handle pointing to the FastQC data file
        :type fp: file handle
        :param line: first line in the module
        :type line: str
        :param start_mark: string denoting start of the module
        :type start_mark: str
        :returns: a list of lines in the module
        :rtype: list of str

        """
        raw = [line]
        while not line.startswith('>>END_MODULE'):
            line = fp.readline()
            raw.append(line)

            if not line:
                raise ValueError("Unexpected end of file in module %r" % line)

        return raw

    @property
    def modules(self):
        """All modules in the FastQC results."""
        return self._modules

    @property
    def passes(self):
        """All module names that pass QC."""
        return self._filter_by_status('pass')

    @property
    def passes_num(self):
        """How many modules have pass status."""
        return len(self.passes)

    @property
    def warns(self):
        """All module names with warning status."""
        return self._filter_by_status('warn')

    @property
    def warns_num(self):
        """How many modules have warn status."""
        return len(self.warns)

    @property
    def fails(self):
        """All names of failed modules."""
        return self._filter_by_status('fail')

    @property
    def fails_num(self):
        """How many modules failed."""
        return len(self.fails)

    @property
    def basic_statistics(self):
        """Basic statistics module results."""
        return self._modules['basic_statistics']

    @property
    def per_base_sequence_quality(self):
        """Per base sequence quality module results."""
        return self._modules['per_base_sequence_quality']

    @property
    def per_sequence_quality_scores(self):
        """Per sequence quality scores module results."""
        return self._modules['per_sequence_quality_scores']

    @property
    def per_base_sequence_content(self):
        """Per base sequence content module results."""
        return self._modules['per_base_sequence_content']

    @property
    def per_base_gc_content(self):
        """Per base GC content module results."""
        return self._modules['per_base_gc_content']

    @property
    def per_sequence_gc_content(self):
        """Per sequence GC content module results."""
        return self._modules['per_sequence_gc_content']

    @property
    def per_base_n_content(self):
        """Per base N content module results."""
        return self._modules['per_base_n_content']

    @property
    def sequence_length_distribution(self):
        """Per sequence length distribution module results."""
        return self._modules['sequence_length_distribution']

    @property
    def sequence_duplication_levels(self):
        """Sequence duplication module results."""
        return self._modules['sequence_duplication_levels']

    @property
    def overrepresented_sequences(self):
        """Overrepresented sequences module results."""
        return self._modules['overrepresented_sequences']

    @property
    def kmer_content(self):
        """Kmer content module results."""
        return self._modules['kmer_content']
