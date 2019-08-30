setwd("/mnt/volume1/T2D/analysis")

#################################
library("PerformanceAnalytics") #paired-wise correlation
phenotype = read.table("/mnt/volume1/T2D/cleanData/feature_max_Tab_YuanAug25.txt",header = T,sep = '\t')
chart.Correlation(phenotype[,c("BMI","age","bp_high","triglyceride", "ldl", "hdl", "a1c","gulcose", "Fasting")],histogram = T)

##################################

plink.lmiss = read.table("plink.lmiss",header = T)
plink.lmiss[1:5,]
lmiss_hist = hist(plink.lmiss$N_MISS, main = "Histogram of missing counts at all loci",xlab = "No. Missing")
text(lmiss_hist$mids,lmiss_hist$counts+15000,labels = lmiss_hist$counts)

plink.imiss = read.table("plink.imiss",header = T)
plink.imiss[1:5,1:5]
imiss_hist = hist(plink.imiss$N_MISS,main = "Histogram of missing counts at all individuals",xlab = "No. Missing")
text(imiss_hist$mids,imiss_hist$counts+55,labels = imiss_hist$counts)


plink.frq.counts = read.table("plink.frq.counts",header = T,stringsAsFactors = F)
plink.frq.counts[1:5,]
plink.frq.counts$C1+plink.frq.counts$C2
plink.frq.counts_with_geno = plink.frq.counts[which((plink.frq.counts$C1+plink.frq.counts$C2)>0),]
fre.count_hist = hist(plink.frq.counts_with_geno$C1/(plink.frq.counts_with_geno$C1+plink.frq.counts_with_geno$C2),
                 main = "Minor allele frequency distribution",xlab = NA,
                 breaks = seq(0,0.5,0.05))
text(fre.count_hist$mids,fre.count_hist$counts+8000,labels = fre.count_hist$counts,cex = 1)



###################################
library("qqman")

#####
#BMI cofactor
BMI_cofactor_a1c = read.table("/mnt/volume1/T2D/analysis/Aug25/BMI_cofactor/plink.a1c.assoc.linear",header = T)
manhattan(BMI_cofactor_a1c[which(BMI_cofactor_a1c$P<0.05),],main= "A1C - BMI as cofactor")
dim(BMI_cofactor_a1c[which(BMI_cofactor_a1c$P<0.0000005),])

BMI_cofactor_triglyceride = read.table("/mnt/volume1/T2D/analysis/Aug25/BMI_cofactor/plink.triglyceride.assoc.linear",header = T,stringsAsFactors = F) 
manhattan(BMI_cofactor_triglyceride[which(BMI_cofactor_triglyceride$P<0.0000005),])
dim(BMI_cofactor_triglyceride[which(BMI_cofactor_triglyceride$P<0.0000005),])

BMI_cofactor_glucose = read.table("/mnt/volume1/T2D/analysis/Aug25/BMI_cofactor/plink.gulcose.assoc.linear",header = T)
manhattan(BMI_cofactor_glucose[which(BMI_cofactor_glucose$P<0.0000005),])
dim(BMI_cofactor_glucose[which(BMI_cofactor_glucose$P<0.0000005),])

BMI_cofactor_Fasting = read.table("/mnt/volume1/T2D/analysis/Aug25/BMI_cofactor/plink.Fasting.assoc.linear",header = T)
manhattan(BMI_cofactor_Fasting[which(BMI_cofactor_Fasting$P<0.05),],main="Fasting Glu - BMI as cofactor")
dim(BMI_cofactor_Fasting[which(BMI_cofactor_Fasting$P<0.0000005),])

BMI_cofactor_OrdinaryGlucose = read.table("/mnt/volume1/T2D/analysis/Aug25/BMI_cofactor/plink.OrdinaryGlucose.assoc.linear",header = T)
manhattan(BMI_cofactor_OrdinaryGlucose[which(BMI_cofactor_OrdinaryGlucose$P<0.0000005),])
dim(BMI_cofactor_OrdinaryGlucose[which(BMI_cofactor_OrdinaryGlucose$P<0.0000005),])

BMI_cofactor_bp_high = read.table("/mnt/volume1/T2D/analysis/Aug25/BMI_cofactor/plink.bp_high.assoc.linear",header = T)
manhattan(BMI_cofactor_bp_high[which(BMI_cofactor_bp_high$P<0.0000005),])
dim(BMI_cofactor_bp_high[which(BMI_cofactor_bp_high$P<0.0000005),])

#No_factor
No_cofactor_a1c = read.table("/mnt/volume1/T2D/analysis/Aug25/No_cofactor/plink.a1c.qassoc",header = T)
manhattan(No_cofactor_a1c[which(No_cofactor_a1c$P<0.01),],main="A1C - no cofactor")
dim(No_cofactor_a1c[which(No_cofactor_a1c$P<0.0000005),])

No_cofactor_triglyceride = read.table("/mnt/volume1/T2D/analysis/Aug25/No_cofactor/plink.triglyceride.qassoc",header = T) 
manhattan(No_cofactor_triglyceride[which(No_cofactor_triglyceride$P<0.0000005),])
dim(No_cofactor_triglyceride[which(No_cofactor_triglyceride$P<0.0000005),])

No_cofactor_glucose = read.table("/mnt/volume1/T2D/analysis/Aug25/No_cofactor/plink.gulcose.qassoc",header = T)
manhattan(No_cofactor_glucose[which(No_cofactor_glucose$P<0.0000005),])
dim(No_cofactor_glucose[which(No_cofactor_glucose$P<0.0000005),])

No_cofactor_Fasting = read.table("/mnt/volume1/T2D/analysis/Aug25/No_cofactor/plink.Fasting.qassoc",header = T)
manhattan(No_cofactor_Fasting[which(No_cofactor_Fasting$P<0.05),],main="Fasting Glucose - no cofactor")
dim(No_cofactor_Fasting[which(No_cofactor_Fasting$P<0.0000005),])

No_cofactor_OrdinaryGlucose = read.table("/mnt/volume1/T2D/analysis/Aug25/No_cofactor/plink.OrdinaryGlucose.qassoc",header = T)
manhattan(No_cofactor_OrdinaryGlucose[which(No_cofactor_OrdinaryGlucose$P<0.0000005),])
dim(No_cofactor_OrdinaryGlucose[which(No_cofactor_OrdinaryGlucose$P<0.0000005),])

No_cofactor_bp_high = read.table("/mnt/volume1/T2D/analysis/Aug25/No_cofactor/plink.bp_high.qassoc",header = T)
manhattan(No_cofactor_bp_high[which(No_cofactor_bp_high$P<0.0000005),])

No_cofactor_BMI =read.table("/mnt/volume1/T2D/analysis/Aug20/plink.BMI.qassoc",header = T)
manhattan(No_cofactor_BMI[which(No_cofactor_BMI$P <0.001),],main="BMI")






intersect_check = function(sig_set){
  print(intersect(Mahajan2018$Index.variant,gsub("GSA-","",sig_set$SNP)))
}

######

#####################################
Mahajan2018 = read.csv("/mnt/volume1/T2D/analysis/Mahajan_2018_Nature_genetics.csv",header = T,stringsAsFactors = F)
Mahajan2018[1:5,]$Index.variant

