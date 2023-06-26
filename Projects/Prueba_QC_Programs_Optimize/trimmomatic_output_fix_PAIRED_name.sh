#!/bin/bash
for f in *.fastq; do mv "$f" "$(echo "$f" | sed s/kneaddata.trimmed/kneaddata_trimmed/)"; done
