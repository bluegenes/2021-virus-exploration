#conda activate virus-smash

# dna
python count-hashes.py --siglist /group/ctbrowngrp/virus-references/pigeon/dna-input/pigeon1.0.signatures.txt --output-csv pigeon1.0.stats2.csv --fastafile /group/ctbrowngrp/virus-references/pigeon/PIGEONv1.0.fa.gz -s 200 -s 1000 -s 2000 -s 10000

# test protein
#python count-hashes.py --siglist /home/ntpierce/2021-virus-exploration/output.protein-pigeon-test/compare/pigeon1.0.prodigal.siglist.txt --output-csv pigeon1.0.proteintest.stats.csv --length-csv output.protein-pigeon-test/fastasplit/pigeon1.0.lengths.txt -s 200 -s 1000 -s 2000 -s 10000
