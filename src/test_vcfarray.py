"""
Unit tests for the vcfarray package.

"""


import numpy as np
import vcfarray
from nose.tools import eq_ 


def test_fromvcfinfo_default():

    # default array
    a = vcfarray.fromvcfinfo('fixture/example-4.1.vcf').view(np.recarray)
    
    # check fixed fields
    eq_('20', a.CHROM[0])
    eq_(14370, a.POS[0])
    eq_('rs6054257', a.ID[0])
    eq_('G', a.REF[0])
    eq_('A', a.ALT[0])
    eq_(29, a.QUAL[0])

    # check FILTER field
    eq_(True, a.FILTER.PASS[0])
    eq_(False, a.FILTER.PASS[1])
    eq_(False, a.FILTER.q10[0])
    eq_(True, a.FILTER.q10[1])
    
    # check other attributes
    eq_(True, a.is_snp[0])
    eq_(False, a.is_indel[0])
    eq_(False, a.is_deletion[0])
    eq_(True, a.is_transition[0])
    eq_(3, a.num_called[0])
    eq_(0, a.num_unknown[0])
    eq_(1., a.call_rate[0])
    
    # check special attributes
    eq_(2, a.num_alleles[0])
    eq_(3, a.num_alleles[2])
    
    # check INFO fields
    eq_(3, a.NS[0])
    eq_(14, a.DP[0])
    eq_(.5, a.AF[0])
    eq_(True, a.DB[0])
    eq_(True, a.H2[0])

    

def test_fromvcfinfo_override_fields():

    # default array
    a = vcfarray.fromvcfinfo('fixture/example-4.1.vcf', fields=['CHROM', 'POS', 'DP']).view(np.recarray)
    
    eq_('20', a.CHROM[0])
    eq_(14370, a.POS[0])
    eq_(14, a.DP[0])

    for f in ('ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'NS', 'AF', 'DB', 'H2'):
        assert not hasattr(a, f), 'unexpected attribute on array: %s' % f
    
    
def test_fromvcfcalldata_default():
    
    # default array
    a = vcfarray.fromvcfcalldata('fixture/example-4.1.vcf')
    
    # check FORMAT fields
    eq_('0|0', a['NA00001']['GT'][0])
    eq_(48, a['NA00001']['GQ'][0])
    eq_(1, a['NA00001']['DP'][0])
    eq_(51, a['NA00001']['HQ'][0][0])
    eq_(51, a['NA00001']['HQ'][0][1])
    
    # check other attributes
    eq_(True, a['NA00001']['called'][0])
    eq_(0, a['NA00001']['gt_type'][0])
    eq_(False, a['NA00001']['is_het'][0])
    eq_(False, a['NA00001']['is_variant'][0])
    
    
def test_fromvcfcalldata_view2d():

    # default array
    a = vcfarray.fromvcfcalldata('fixture/example-4.1.vcf')
    a = vcfarray.view2d(a).view(np.recarray)

    # check FORMAT fields
    eq_('0|0', a.GT[0,0])
    eq_(48, a.GQ[0,0])
    eq_(1, a.DP[0,0])
    eq_(51, a.HQ[0,0][0])
    eq_(51, a.HQ[0,0][1])
    
    # check other attributes
    eq_(True, a.called[0,0])
    eq_(0, a.gt_type[0,0])
    eq_(False, a.is_het[0,0])
    eq_(False, a.is_variant[0,0])



    

