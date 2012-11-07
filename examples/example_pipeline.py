#!/usr/bin/env python
import sys
import os
import gem
import time


## global configuration
# number of threads to use for each run
THREADS = 4
# remove intermediate files
REMOVE_FILES = False

## make sure we have access to the demo data
base_dir = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
result_dir = "%s/results" % base_dir
data_dir = "%s/data" % base_dir

if not os.path.exists(data_dir):
    print("Unable to find teh demo data."
          "Please make sure you ran the bootstrap "
          "script to download and prepare demo data")
    exit(1)


def create_index():
    """Create the example index"""
    print "Preparing index"
    return gem.index("%s/chr21.fa" % data_dir, "%s/chr21.gem" % result_dir, threads=THREADS)


if __name__ == "__main__":
    ## define the input
    #
    # the gem index
    index = create_index()
    # the input reads (multiplexed fastq)
    reads = "%s/data_chr21.fastq" % data_dir
    # the annotation used for junction site detection
    annotation = "%s/chr21.gtf" % data_dir

    # define some output file name
    #
    # global name for the resulting data-set
    name = "%s/chr21_mapping" % result_dir
    # initial mapping output
    initial_out = "%s_initial.map" % name
    # output from the denovo junction detection step
    denovo_out = "%s_denovo.map" % name
    # denovo junction sites
    junctions_out = "%s.junctions" % name
    # initial splitmap output
    initial_split_out = "%s_initial_split.map" % name
    # second round mapping with 3" end trimmed by 20 bases
    trim_20_out = "%s_trim_20.map" % name
    # second round split mapping with 3" end trimmed by 20 bases
    trim_20_split_out = "%s_trim_20_split.map" % name
    # first end second mapping step merged
    merge_out = "%s_merged.map" % name
    # paired mappings
    paired_out = "%s_paired.map" % name
    # scored mappings
    final_out = "%s.map" % name
    # sam/bam output
    sam_out = "%s.bam" % name

    ## Create initial mapping
    # we deal with a single file with interleaved paired reads here, but
    # creating input from two files for the read pairs is straight forward using the interleave filter
    #
    # input_1 = gem.files.open(input_file)
    # input_2 = gem.files.open(input_file2)
    # input = gem.filter.interleave([input_1, input_2])
    print "Running initial mapping"
    input = gem.files.open(reads)
    initial_mapping = gem.mapper(input, index, initial_out, mismatches=0.07, delta=1, threads=THREADS)

    ## junction sites
    # before we can do the split mapping, we have to load
    # the junction sites from a gtf annotation and
    # run the denovo-junction detection. This will also give
    # us a mapping that preserves short indels detected during the
    # extraction run.
    print "Loading GTF junctions from %s" % annotation
    junctions = gem.junctions.from_gtf(annotation)

    # now the denovo run. This returns a tuple : (mapping, junctions)
    # and here we use the merge_with parameter to merge the denovo junctions with
    # the previously loaded gtf junctions.
    #
    # Also note that we pass only unmapped reads from the initial mapping to the junction
    # extraction. This is done using the "unmapped filter"
    print "Getting de-novo junctions"
    (denovo_mapping, junctions) = gem.extract_junctions(
        gem.filter.unmapped(initial_mapping), # only unmapped reads from the initial mapping
        index,
        denovo_out,
        mismatches=0.04,
        threads=THREADS,
        merge_with=set(junctions)
    )

    ## we filter the junnctions now by their distance and
    # write all junctions with a distance <= 500000 to a file
    print "Writing junctions file"
    gem.junctions.write_junctions(gem.junctions.filter_by_distance(junctions, 500000), junctions_out, index)


    ## Initial split map run with junction sites
    # Here we take all unmapped reads after the denovo junction detection and
    # pass them to the split mapper along with the junctions
    print "Running initial split-map"
    initial_split_mapping = gem.splitmapper(
        gem.filter.unmapped(denovo_mapping),
        index,
        initial_split_out,
        junctions_file=junctions_out,
        mismatches=0.06,
        threads=THREADS)


    ## Phase two, hard trim 20 bases from the right side
    # and try mapping again
    print "Running trim 20 mapping"
    trim_20_mapping = gem.mapper(
        gem.filter.unmapped(initial_split_mapping),
        index,
        trim_20_out,
        trim=(0, 20),  # this is the additional parameter to control the trimming
        threads=THREADS)

    ## second round split mapping
    print "Running trim 20 split-map"
    trim_20_split_mapping = gem.splitmapper(
        gem.filter.unmapped(trim_20_mapping),
        index,
        trim_20_split_out,
        trim=(0, 20),  # this is the additional parameter to control the trimming
        junctions_file=junctions_out,
        threads=THREADS)


    ## We are done with the basic mapping pipeline,
    # now we merge everything into the final file
    # Note that we call the .clone() method on most
    # of the objects we pass to the merger. This is because the
    # mapping steps return iterators over the results and
    # we eventually iterated through them already once. The
    # .clone() makes sure we new iterator that starts at the beginning
    print "Merging results"
    merged = gem.merger(
        initial_mapping.clone(), ## this is the main target and must contain ALL reads
        [denovo_mapping.clone(),
         initial_split_mapping.clone(),
         trim_20_mapping.clone(),
         trim_20_split_mapping]
    ).merge(merge_out) # write to file


    ## remove files
    if REMOVE_FILES:
        print "Removing intermediate files"
        os.remove(initial_out)
        os.remove(initial_split_out)
        os.remove(denovo_out)
        os.remove(junctions_out)
        os.remove(trim_20_out)
        os.remove(trim_20_split_out)

    ## pair align the mappings
    print "Running pair aligner"
    paired_mapping = gem.pairalign(merged, index, paired_out, max_insert_size=100000, threads=THREADS)

    if REMOVE_FILES:
        os.remove(merge_out)

    ## validate and score the alignment
    print "Validating and scoring alignment"
    #scored = gem.validate_and_score(paired_mapping, index, final_out, threads=8)
    #scored = gem.score(paired_mapping, index, final_out, threads=8)
    scored = gem.validate(paired_mapping, index, final_out, threads=THREADS)
    if REMOVE_FILES:
        os.remove(paired_out)

    ## convert the result to sam and then to bam
    print "Converting to sam"
    sam = gem.gem2sam(scored, index, threads=4)
    bam = gem.sam2bam(sam, sam_out, sorted=True)

    print "Done :)"
