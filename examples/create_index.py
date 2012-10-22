#!/usr/bin/python
"""Example script that creates the example gem-index"""

import gem
import os

if __name__ == "__main__":
    basedir = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    target = "%s/results/chr21.gem" % basedir
    source = "%s/data/chr21.fa" % basedir
    print("Creating index for %s" % (source))
    print("GEM index is created in %s" % (target))

    if not os.path.exists(source):
        print("Unable to find source file %s. "
              "Please make sure you ran the bootstrap "
              "script to download and prepare demo data" % (source))
        exit(1)
    index = gem.index(source, target)
    print("Created index %s " % index)


