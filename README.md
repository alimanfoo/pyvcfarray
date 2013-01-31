pyvcfarray
==========

Python utility functions for loading data from Variant Call Format
(VCF) files into numpy arrays.

Installation
------------

```
$ pip install vcfarray
```

Usage
-----

```python
import vcfarray
fn = '/path/to/my/variants.vcf'
info = vcfarray.fromvcfinfo(fn)
calldata = vcfarray.fromvcfcalldata(fn)
calldata2d = vcfarray.view2d(calldata)
```

