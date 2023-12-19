#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
READ1=$1
READ2=$2
BASE=$3
RULE=$4
MISMATCH=$5

echo "Using ${BASE}"

# Split files
zcat ../$1 | parallel --pipe -N1000000 "gzip > SPLIT_{#}_R1.fq.gz"
zcat ../$2 | parallel --pipe -N1000000 "gzip > SPLIT_{#}_R2.fq.gz"

parallel -j30 --link python3 -W ignore /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/wuyuxuan/dat/Hairpin-Seq/couplet/bin/run_couplet.py  --fq1 {1}  --fq2 {2} --mismatch-threshold ${MISMATCH} --end-gap-penalty -2 --gap-penalty -2 --match-award 1 --mismatch-penalty -1 --phred=prob --rule ${RULE} ::: SPLIT_*_R1.fq.gz ::: SPLIT_*_R2.fq.gz

cat SPLIT_*_resolved.fq.gz > "${BASE}_resolved.fq.gz"
cat SPLIT_*_discarded_R1.fq.gz > "${BASE}_discarded_R1.fq.gz"
cat SPLIT_*_discarded_R2.fq.gz > "${BASE}_discarded_R2.fq.gz"

python3 /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/wuyuxuan/dat/Hairpin-Seq/couplet/bin/postprocess_stats.py --input-stats-files SPLIT_*couplet.yaml --output-prefix ${BASE}

# rm SPLIT_*
