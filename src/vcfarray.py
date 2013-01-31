VERSION = '0.2'


import vcf
import numpy as np


# VCF types
VCF_TYPE_INTEGER = 'Integer'
VCF_TYPE_FLOAT = 'Float'
VCF_TYPE_CHARACTER = 'Character'
VCF_TYPE_STRING = 'String'
VCF_TYPE_FLAG = 'Flag'


# default mapping from VCF types to numpy dtypes
DEFAULT_DTYPES = {VCF_TYPE_INTEGER: 'i4',
                  VCF_TYPE_FLOAT: 'f4',
                  VCF_TYPE_CHARACTER: 'a1',
                  VCF_TYPE_STRING: 'a20',
                  VCF_TYPE_FLAG: 'b1'}


# default fill values to use in place of missing data
DEFAULT_FILLVALUES = {VCF_TYPE_INTEGER: 0,
                      VCF_TYPE_FLOAT: 0,
                      VCF_TYPE_CHARACTER: '',
                      VCF_TYPE_STRING: '',
                      VCF_TYPE_FLAG: False}


# attributes of PyVCF variant records, including VCF fixed fields
VARIANT_ATTRIBUTES = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                      'is_snp', 'is_indel', 'is_deletion', 'is_transition')


# default numpy dtypes to use for VCF fixed fields and other variant attributes 
DEFAULT_VARIANT_ATTRIBUTE_DTYPES = {'CHROM': 'a20',
                                    'POS': 'i4',
                                    'ID': 'a20',
                                    'REF': 'a20',
                                    'ALT': 'a20',
                                    'QUAL': 'f4',
                                    'FILTER': 'a20',
                                    'is_snp': 'b1', 
                                    'is_indel': 'b1', 
                                    'is_deletion': 'b1', 
                                    'is_transition': 'b1'}


def fromvcfinfo(filename, fields=None, types=None, arities=None, fillvalues=None, converters=None):
    """
    Load a numpy structured array from data in the fixed and INFO
    fields of a Variant Call Format (VCF) file.
    
    If just the filename is given, then all fixed and INFO fields will
    be loaded.
    
    To only extract specific fields, provide a sequence of strings as
    the `fields` keyword argument.
    
    To override the default numpy dtype for one or more fields,
    provide a dict as the `types` keyword argument mapping field names
    to numpy dtype specifications.
    
    To override the default arity (i.e., number of expected values)
    for one or more fields, provide a dict as the `arities` keyword
    argument mapping field names to integers.
    
    To override the default fill value to use when a value is missing
    or None, provide a dict as the `fillvalues` keyword argument
    mapping field names to fill values.
    
    To inject functions for preprocessing values from specific
    fields before loading into the array, provide a dict as the
    `converters` keyword argument mapping field names to functions.

    """

    # set up VCF reader
    vcf_reader = vcf.Reader(filename=filename)

    # determine fields to use
    if fields is None:
        # start with variant attributes
        fields = list(VARIANT_ATTRIBUTES)
        # add in all INFO fields declared in VCF header
        fields.extend(vcf_reader.infos.keys())
    else:
        # check all requested fields are available
        for f in fields:
            assert f in VARIANT_ATTRIBUTES or f in vcf_reader.infos.keys(), 'bad field name: %s' % f
        
    # determine a numpy dtype to use for each field
    if types is None:
        types = dict()
    for f in fields:
        if f not in types:
            if f in VARIANT_ATTRIBUTES:
                t = DEFAULT_VARIANT_ATTRIBUTE_DTYPES[f]
            else:
                vcf_t = vcf_reader.infos[f].type
                t = DEFAULT_DTYPES[vcf_t]
            types[f] = t
    
    # determine expected arity for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f in VARIANT_ATTRIBUTES:
            arities[f] = 1 # expect only one value
        elif f not in arities:
            vcf_n = vcf_reader.infos[f].num
            if vcf_n > 1:
                n = vcf_n
            else:
                n = 1 # fall back to expecting one value
            arities[f] = n
    
    # decide what fill values to use for each INFO field if value is missing
    if fillvalues is None:
        fillvalues = dict()
    for f in fields:
        if f not in fillvalues:
            if f in VARIANT_ATTRIBUTES:
                v = '' # should only apply to ID
            else:
                vcf_t = vcf_reader.infos[f].type
                v = DEFAULT_FILLVALUES[vcf_t]
            fillvalues[f] = v
            
    # pad out converters
    if converters is None:
        converters = dict()
    for f in fields:
        if f not in converters:
            converters[f] = None

    # construct a numpy dtype for structured array
    dtype = list()
    for f in fields:
        t = types[f]
        n = arities[f]
        if n == 1:
            dtype.append((f, t))
        else:
            dtype.append((f, t, (n,)))
    
    # set up an iterator over the VCF records
    it = _itervcfinfo(vcf_reader, fields, arities, fillvalues, converters)

    # build an array from the iterator
    a = np.fromiter(it, dtype=dtype)

    return a


def _itervcfinfo(vcf_reader, fields, arities, fillvalues, converters):
    for rec in vcf_reader:
        yield tuple(_mkival(rec, f, arities[f], fillvalues[f], converters[f]) for f in fields)


def _mkival(rec, f, num, fill, conv):
    if f in VARIANT_ATTRIBUTES:
        val = getattr(rec, f)
        if conv is not None: # user-provided value converter
            val = conv(val)
        # special case ALT and FILTER because of variable length
        elif f in {'ALT', 'FILTER'}:
            val = ','.join(map(str, val))
        elif val is None:
            val = fill
    elif f in rec.INFO:
        val = rec.INFO[f]
        if conv is not None: # user-provided value converter, overrides everything else
            val = conv(val)
        elif num > 1:
            if isinstance(val, basestring) and ',' in val:
                val = val.split(',')
            val = tuple(val[:num]) # try to pick off as many values as requested
        elif isinstance(val, (list, tuple)) and len(val) > 0:
            val = val[0] # fall back to picking off first value
        elif isinstance(val, (list, tuple)) and len(val) == 0:
            # edge case
            val = fill
        else:
            pass # leave val as-is
    else:
        val = fill
    return val    


def fromvcfcalldata(filename, samples=None, fields=None, types=None, arities=None, fillvalues=None, converters=None):   
    """
    Load a numpy structured array from data in the sample columns of a
    Variant Call Format (VCF) file.
    
    If just the filename is given, then all FORMAT fields for all
    samples will be loaded.
   
    To only extract data for specific samples, provide a sequence of
    strings as the `samples` keyword argument.
    
    To only extract specific fields, provide a sequence of strings as
    the `fields` keyword argument.
    
    To override the default numpy dtype for one or more fields,
    provide a dict as the `types` keyword argument mapping field names
    to numpy dtype specifications.
    
    To override the default arity (i.e., number of expected values)
    for one or more fields, provide a dict as the `arities` keyword
    argument mapping field names to integers.
    
    To override the default fill value to use when a value is missing
    or None, provide a dict as the `fillvalues` keyword argument
    mapping field names to fill values.
    
    To inject functions for preprocessing values from specific
    fields before loading into the array, provide a dict as the
    `converters` keyword argument mapping field names to functions.

    """
    
    # set up VCF reader
    vcf_reader = vcf.Reader(filename=filename)
    
    # determine samples to extract calldata for
    if samples is None:
        samples = vcf_reader.samples
    else:
        # check all requested samples are available
        for s in samples:
            assert s in vcf_reader.samples, 'bad sample: %s' % s

    # determine fields to use
    if fields is None:
        # use all format fields
        fields = list(vcf_reader.formats.keys())
    else:
        # check all requested fields are available
        for f in fields:
            assert f in vcf_reader.formats.keys(), 'bad field name: %s' % f
        
    # determine a numpy dtype to use for each field
    if types is None:
        types = dict()
    for f in fields:
        if f not in types:
            vcf_t = vcf_reader.formats[f].type
            t = DEFAULT_DTYPES[vcf_t]
            types[f] = t
    
    # determine expected arity for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f not in arities:
            vcf_n = vcf_reader.formats[f].num
            if vcf_n > 1:
                n = vcf_n
            else:
                n = 1 # fall back to expecting one value
            arities[f] = n
    
    # decide what fill values to use for each field if value is missing
    if fillvalues is None:
        fillvalues = dict()
    for f in fields:
        if f not in fillvalues:
            vcf_t = vcf_reader.formats[f].type
            v = DEFAULT_FILLVALUES[vcf_t]
            fillvalues[f] = v
            
    # pad out converters
    if converters is None:
        converters = dict()
    for f in fields:
        if f not in converters:
            converters[f] = None
            
    # construct a numpy dtype for calldata cells
    cell_dtype = list()
    for f in fields:
        t = types[f]
        n = arities[f]
        if n == 1:
            cell_dtype.append((f, t))
        else:
            cell_dtype.append((f, t, (n,)))
            
    # construct a numpy dtype for structured array
    dtype = [(s, cell_dtype) for s in samples]
    
    # set up iterator
    it = _itervcfcalldata(vcf_reader, samples, fields, arities, fillvalues, converters)

    # build an array from the iterator
    a = np.fromiter(it, dtype=dtype)

    return a


def _itervcfcalldata(vcf_reader, samples, fields, arities, fillvalues, converters):
    for rec in vcf_reader:
        yield tuple(_mkcvals(rec.genotype(s), fields, arities, fillvalues, converters) for s in samples)
    
    
def _mkcvals(call, fields, arities, fillvalues, converters):
    return tuple(_mkcval(call, f, arities[f], fillvalues[f], converters[f]) for f in fields)


def _mkcval(call, f, num, fill, conv):
    val = call[f]
    if conv is not None: # user-provided value converter, overrides everything else
        val = conv(val)
    elif num > 1:
        if isinstance(val, basestring) and ',' in val:
            val = val.split(',')
        val = tuple(val[:num]) # try to pick off as many values as requested
    elif isinstance(val, (list, tuple)) and len(val) > 0:
        val = val[0] # fall back to picking off first value
    elif isinstance(val, (list, tuple)) and len(val) == 0:
        # edge case
        val = fill
    elif val is None:
        val = fill
    return val


def view2d(a):
    """
    Utility function to view a structured 1D array where all fields have a uniform dtype 
    (e.g., an array constructed by :func:fromvcfcalldata) as a 2D array.
    
    """
    
    rows = a.size
    cols = len(a.dtype)
    dtype = a.dtype[0]
    b = a.view(dtype).reshape(rows, cols)
    return b
    
