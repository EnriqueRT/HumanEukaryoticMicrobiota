#!/bin/bash
for f in *.fastq; do mv "$f" "$(echo "$f" | sed s/.denovo_duplicates_marked.trimmed./_denovo_duplicates_marked_trimmed_/)"; done
