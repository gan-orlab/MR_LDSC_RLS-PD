
## Linkage disequilibrium score regression (LDSC) between RLS and PD
#Open R and install the required packages if necessary
require(data.table)
#Change name of the file
data1 <-fread("RLS_GWAS_withSE.txt")
#See the names of the columns
head(data1)
#Calculate Zscore
data1$Zscore <- data1$Effect/data1$StdErr
#Select the data
data_selected<-select(data1,snpid,A1, A2, Zscore, N, Pvalue)
#save the file
write.table(data_selected, file="RLS_GWAS_fiveColumns.txt", row.names=FALSE, quote=FALSE, sep = "\t")
#Perform LDSC using python package (Bulik-Sullivan B et al., 2015)
munge_sumstats.py --out RLS --sumstats RLS_GWAS_fiveColumns.txt --merge-alleles w_hm3.snplist --N-col N --snp snpid --a1 A1 --a2 A2 --p Pvalue --signed-sumstats Zscore,0
ldsc.py --rg  PD_with_UKBB.sumstats.gz --out PD_RLS --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/
