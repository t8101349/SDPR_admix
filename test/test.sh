
# train the PRS model
./SDPR_admix -vcf test/chr22_train.vcf.gz -msp test/chr22_train.msp.tsv -pheno test/train.pheno -covar test/covar.txt -iter 100 -burn 0  -out res.txt
