#!/usr/bin/env python

import sys
import time
import pstats
import cProfile
import timeit


sys.path.append('src')
import vcfarray


def profile():
    a = vcfarray.fromvcfinfo(sys.argv[1])


prof_fn = 'profile.prof'
cmd = 'profile()'
cProfile.runctx(cmd, globals(), locals(), prof_fn)
s = pstats.Stats(prof_fn)
s.strip_dirs().sort_stats('time').print_stats()
print timeit.repeat(cmd, number=20, repeat=3, setup='from __main__ import profile')


