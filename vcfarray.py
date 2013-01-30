VERSION = 0.1


import vcf
import numpy as np


VCF_TYPE_INTEGER = 'Integer'
VCF_TYPE_FLOAT = 'Float'
VCF_TYPE_CHARACTER = 'Character'
VCF_TYPE_STRING = 'String'
VCF_TYPE_FLAG = 'Flag'


VCF_FIXED_FIELDS = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER')


# default dtypes for fixed fields
DEFAULT_FIXED_FIELD_DTYPES = {'CHROM': 'a20',
                              'POS': 'i4',
                              'ID': 'a20',
                              'REF': 'a20',
                              'ALT': 'a20',
                              'QUAL': 'f4',
                              'FILTER': 'a20'}


# default mapping from VCF types to numpy dtypes
DEFAULT_DTYPES = {VCF_TYPE_INTEGER: 'i4',
                  VCF_TYPE_FLOAT: 'f4',
                  VCF_TYPE_CHARACTER: 'a1',
                  VCF_TYPE_STRING: 'a20',
                  VCF_TYPE_FLAG: 'b1'}


# default fill values
DEFAULT_FILLVALUES = {VCF_TYPE_INTEGER: 0,
                      VCF_TYPE_FLOAT: 0,
                      VCF_TYPE_CHARACTER: '',
                      VCF_TYPE_STRING: '',
                      VCF_TYPE_FLAG: False}


def fromvcfinfo(filename, dtype=None, fillvalues=None):

    # set up VCF reader
    vcf_reader = vcf.Reader(filename=filename)

    # determine dtype for structured array
    if dtype is not None:
        # user has provided dtype, check all field names are either fixed fields or declared INFO fields
        for dt in dtype:
            fieldname = dt[0]
            if fieldname not in VCF_FIXED_FIELDS: 
                assert fieldname in vcf_reader.infos, 'field name not declared in VCF INFO header: %s' % fieldname
    else:
        # add all fixed fields
        dtype = [('CHROM', DEFAULT_FIXED_FIELD_DTYPES['CHROM']), 
                 ('POS', DEFAULT_FIXED_FIELD_DTYPES['POS']), 
                 ('ID', DEFAULT_FIXED_FIELD_DTYPES['ID']),  
                 ('REF', DEFAULT_FIXED_FIELD_DTYPES['REF']),  
                 ('ALT', DEFAULT_FIXED_FIELD_DTYPES['ALT']), 
                 ('QUAL', DEFAULT_FIXED_FIELD_DTYPES['QUAL']), 
                 ('FILTER', DEFAULT_FIXED_FIELD_DTYPES['FILTER'])]
        # add all INFO fields declared in VCF header
        for info in vcf_reader.infos.values():
            # map VCF type to numpy type
            t = DEFAULT_DTYPES[info.type]
            # how many values should we expect?
            if info.type == VCF_TYPE_FLAG:
                # expect one value (ignore num)
                dtype.append((info.id, t))
            elif info.num == 1:
                # expect only one value
                dtype.append((info.id, t))
            elif info.num > 1:
                # expect multiple values
                dtype.append((info.id, t, (info.num,)))
            else:
                # fall back to one value
                dtype.append((info.id, t))

    # decide what fill values to use for each INFO field if value is missing
    if fillvalues is None:
        fillvalues = dict()
    for dt in dtype:
        fieldname = dt[0]
        # use default fill value for type if not specified by user
        if fieldname not in fillvalues and fieldname not in VCF_FIXED_FIELDS:
            fillvalues[fieldname] = DEFAULT_FILLVALUES[vcf_reader.infos[fieldname].type]

    # set up an iterator over the VCF records
    it = itervcfinfo(vcf_reader, dtype, fillvalues)

    # build an array from the iterator
    a = np.fromiter(it, dtype=dtype)

    return a


def itervcfinfo(vcf_reader, dtype, fillvalues):
    infos = vcf_reader.infos
    for rec in vcf_reader:
        out = tuple(mkval(dt, rec, fillvalues, infos) for dt in dtype)
        yield out


def mkval(dt, rec, fillvalues, infos):
    # what's the field name?
    f = dt[0]
    if f in VCF_FIXED_FIELDS:
        val = getattr(rec, f)
        if f in {'ALT', 'FILTER'}:
            val = ','.join(str(v) for v in val)
        elif val is None:
            val = ''
    elif f in rec.INFO:
        val = rec.INFO[f]
        # work around PyVCF handling of multiple string values
        if infos[f].type == VCF_TYPE_STRING and ',' in val:
            val = val.split(',')
        # how many values are expected in the dtype?
        if len(dt) == 2:
            # expect one value
            if isinstance(val, (tuple, list)):
                if len(val) > 0:
                    # fall back to first value
                    val = val[0]
                else:
                    val = fillvalues[f]
            else:
                pass # leave val as-is
        elif len(dt) == 3:
            # expect multiple values
            n = dt[2][0]
            val = tuple(val[:n]) # attempt to pick out as many values as requested
    else:
        val = fillvalues[f]
    return val    
    
