#!/usr/bin/python
"""Example that reads in a GTF annotation and prints the number of entries"""

import gem
import os
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Please specify the gtf file to read"
        exit(1)
    junctions = gem.junctions.from_gtf(sys.argv[1])
    print len(set(junctions)), "Possible Junction Sites"  # we create a set because from_gtf returns a generator 

