#!/bin/bash

function download(){
    file=$1
    if [ ! -f "${file%.*}" ]; then
        echo "Getting $file"
        curl -L -O https://github.com/downloads/gemtools/gemtools-examples/${file}
        gunzip ${file}
    fi
}

if [ ! -d "data" ]; then
    mkdir -p data
fi

cd data
download "data_chr21.fastq.gz"
download "chr21.fa.gz"
download "chr21.gtf.gz"

cd ..
mkdir -p results
