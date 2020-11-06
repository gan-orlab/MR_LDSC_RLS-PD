## Mendelian randomization when RLS is exposure and PD is outcome
# Select SNPs with p value <= 5e-08 and create a new file named “RLS_meta.GCyes_prevalence_Ncases_Ncontrols_SampleSize_significantPvalues.txt”

awk -F"\t" '{if $10<=5e-08) print $0}' RLS_meta.GCyes_prevalence_Ncases_Ncontrols_SampleSize.txt > RLS_meta.GCyes_prevalence_Ncases_Ncontrols_SampleSize_significantPvalues.txt

#Open R, install packages if necessary

install.package("TwoSampleMR")
require(TwoSampleMR)

#Read exposure data and perform clumping with standard parameters 

exp_data <- read_exposure_data("RLS_meta.GCyes_prevalence_Ncases_Ncontrols_SampleSize_significantPvalues.txt", sep="\t", snp_col = "SNP", beta_col = "Effect",  eaf_col = "Freq1", se_col="StdErr", effect_allele_col = "Allele1", other_allele_col= "Allele2", pval_col = "Pvalue", samplesize_col = "SampleSize", ncase_col= "Ncases", ncontrol_col="Ncontrols",  clump=TRUE)
#Save and quit
write.csv(exp_data,file="RLS_exposure",quote=FALSE)
q()

#Open R, install packages if necessary
install.package("ggplot2")
require(ggplot2)	
install.package("dplyr")
library(dplyr)
install.package("MRPRESSO")
library(MRPRESSO)
#Read exposure and outcome data
exp_data <- read_exposure_data("RLS_exposure", sep=",", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", ncase_col= "ncase.exposure", ncontrol_col="ncontrol.exposure",  clump=FALSE)
out_data <- read_outcome_data(snps = exp_data$SNP,filename = "META_no23_yesUKBB_SampleSize.txt", sep="\t",  snp_col = "SNP", beta_col = "b", se_col = "se", eaf_col = "freq", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p", ncase_col = "N_cases", ncontrol_col = "N_controls",samplesize_col="SampleSize")
out_data$r.outcome<- get_r_from_lor(out_data$beta.outcome, out_data$eaf.outcome, out_data$ncase.outcome, out_data$ncontrol.outcome, 0.01,  model = "logit")
#Harmonization of data
dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
#To calculate r2 and perform Steiger filtering, we added the columns. Prevalence should also be specified in each exposure and outcome file or manually added

dat$units.outcome<-"log odds"
dat$units.exposure<-"log odds"
dat1<-subset(dat, dat$eaf.exposure!="NA")
dat1$r.exposure<- get_r_from_lor(dat1$beta.exposure, dat1$eaf.exposure, dat1$ncase.exposure, dat1$ncontrol.exposure, 0.1,  model = "logit")
dat1$prevalence.exposure<- 0.1
dat1$prevalence.outcome<- 0.01
#Steiger filtering was performed to exclude SNPs that explain more variance in the outcome than in the exposure
steiger <- steiger_filtering(dat1)
#MR PRESSO was performed to detect horizontal pleiotropy
sig<-subset(steiger, steiger$steiger_dir==TRUE)
presso <-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data= sig, NbDistribution = 1000,  SignifThreshold = 0.05) 
capture.output(print(presso), file = "presso.txt")
#Perform MR
mr(sig)
#F-statitistics and R2 were calculated for each exposure
R2<-mean(sig$r.exposure)
capture.output(print(R2), file="r2.txt")
n<-mean(sig$samplesize.exposure)
k <- nrow(subset(sig, sig$ambiguous == FALSE))
F<-(R2*(n-1-k))/((1-R2)*k)
capture.output(print(F), file = "f.txt")
#Output from the `mr_report` function will generate a report containing tables and graphs summarising the results.
mr_report(sig)
res_single <- mr_singlesnp(sig)
p5 <- mr_forest_plot(res_single)
p5[[1]]
ggsave(p5[[1]], file= "plot.jpg", width=7, height=12)
