
setwd("E:/reeditedsurvival")

load("LmqcmData.RData")
library(survival)
library(survminer)
###debug(utils:::unpackPkgZip)#######use this before using "install.packages" command if "unable to move temporary..." error comes
###install.packages("survivalROC")
library(survivalROC)
##devtools::install_github("matthieugomez/statar")
library(statar)##for x-tile

##########follow: "http://bioconnector.org/workshops/r-survival.html"##############

############Squamous sample##############
ADENO_samplename=read.csv("ADENO_SAMPLEID_su.csv",header=F,sep=",") ###updated##with cpgnames (rows) and samples (columns)
SCC_samplename=read.csv("SCC_SAMPLEID_su.csv",header=F,sep=",") ###updated###with cpgnames (rows) and samples (columns)


ADENO_samplename$V1<-gsub(".", "-", ADENO_samplename$V1[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") 
SCC_samplename$V1<-gsub(".", "-", SCC_samplename$V1[], fixed = TRUE)


ADENO_samplename$V1 = substr(ADENO_samplename$V1,1,nchar(ADENO_samplename$V1)-3)##discard last three characters e.g., "-01", "-11", "-06", etc.
SCC_samplename$V1 = substr(SCC_samplename$V1,1,nchar(SCC_samplename$V1)-3)##discard last three characters e.g., "-01", "-11", "-06", etc.

allsubtypes<-rbind(ADENO_samplename,SCC_samplename)

############read the clinical information##############
clinicalinfo=read.csv("CESC_clinicalMatrix.csv",header=T,sep=",") ###with cpgnames (rows) and samples (columns)##nrow=185

clinicalinfo_filt<-clinicalinfo[which(clinicalinfo$bcr_patient_barcode %in% allsubtypes$V1),]##nrow=178##filtered out the samples that are not in the set of all the experimental subtypes 
clinicalinfo_filt$clinical_stage
clinicalinfo_filt$clinical_stage<-gsub("A", "", clinicalinfo_filt$clinical_stage[], fixed = TRUE)##replace "dot"("A") by "hypen" ("") since here e.g., "Stage IIA" is replaced by "Stage II"
clinicalinfo_filt$clinical_stage<-gsub("B", "", clinicalinfo_filt$clinical_stage[], fixed = TRUE)##replace "dot"("B") by "hypen" ("") since here e.g., "Stage IIB" is replaced by "Stage II"
clinicalinfo_filt$clinical_stage<-gsub("1", "", clinicalinfo_filt$clinical_stage[], fixed = TRUE)
clinicalinfo_filt$clinical_stage<-gsub("2", "", clinicalinfo_filt$clinical_stage[], fixed = TRUE)
clinicalinfo_filt<-clinicalinfo_filt[which(nchar(clinicalinfo_filt$clinical_stage[])!=0),]##nrow=268##eliminate the rows that has no stage information
nrow(clinicalinfo_filt)
clinicalinfo_filt<-clinicalinfo_filt[order(clinicalinfo_filt$clinical_stage),]###ordering the whole matrix according to the alphabetical order of its stages (i.e., c("Stage I","Stage II","Stage III","Stage IV"))



length(clinicalinfo_filt$clinical_stage)##268

sum(!is.na(match(clinicalinfo_filt$clinical_stage, "Stage I")))##41####check how many among 268 samples belonging to "Stage I"
sum(!is.na(match(clinicalinfo_filt$clinical_stage, "Stage II")))##53####check how many among 268 samples belonging to "Stage II"
sum(!is.na(match(clinicalinfo_filt$clinical_stage, "Stage III")))#44####check how many among 268 samples belonging to "Stage III"
sum(!is.na(match(clinicalinfo_filt$clinical_stage, "Stage IV")))##17####check how many among 268 samples belonging to "Stage IV"


clinicalinfo_filt$bcr_patient_barcode
clinicalinfo_filt$vital_status
clinicalinfo_filt$OS.time
unique(clinicalinfo_filt$clinical_stage)
clinicalinfo_filt$days_to_last_followup ##updated
  
for(i in 1:length(clinicalinfo_filt$vital_status)){
  print(i)
  print(clinicalinfo_filt$vital_status[i])
  if(match(as.character(clinicalinfo_filt$vital_status[i]),"LIVING",nomatch = 0)==TRUE) ###nrow=116=number of alive
  {
    clinicalinfo_filt$recordedStatus[i]=1 #living will be denoted by 1
    clinicalinfo_filt$final_time[i]<-clinicalinfo_filt$days_to_last_followup[i]
  }
  if(match(as.character(clinicalinfo_filt$vital_status[i]),"DECEASED",nomatch = 0)==TRUE)  ###nrow=59=number of dead
  {
    clinicalinfo_filt$recordedStatus[i]=0 #dead will be denoted by 0
    clinicalinfo_filt$final_time[i]<-clinicalinfo_filt$OS.time[i]
  }
}
write.csv(clinicalinfo_filt,"clinical_su.csv")
mySurv<-Surv(time=clinicalinfo_filt$final_time, event = clinicalinfo_filt$recordedStatus)
class(mySurv)
head(mySurv)
mySurv


myfit<-survfit(mySurv~1 ,data=clinicalinfo_filt)##single curve by previously fitted Cox model for all (filtered) patients that are matched with experimental samples 
myfit##median=636
plot(myfit)
plot(myfit,conf.int="both")
ggsurvplot(myfit)




#############Overall survival analysis of different pathological stages by Kaplan-Meier Curves#######
myfit<-survfit(mySurv~clinicalinfo_filt$clinical_stage)
myfit##median for stage I samples= 232, median for stage II samples= 444, and median for stage III samples= 394
plot(myfit,col=c("red","blue","green","brown"),mark=3,xlab = "Overall survival time (days to death)", ylab = "survival probability")
legend("topright",c("Stage I","Stage II","Stage III","Stage IV"),col=c("red","blue","green","brown"),lty=1)
title("Overall survival analysis of different pathological stages by Kaplan-Meier Curves")
tp<-clinicalinfo_filt$clinical_stage
Stage<-tp##updated
#lt<-survfit(mySurv~tp)
lt<-survfit(mySurv~Stage)##included new
ggsurvplot(lt,main="Overall survival of pathological stages by Kaplan-Meier Curves",data = clinicalinfo_filt)

ggsurvplot(lt,data = clinicalinfo_filt,conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c( "Stage I","Stage II",  "Stage III", "Stage IV"), legend.title="Pathological Stage")##pval is here from lor-rank test

ggsurvplot(lt, data = clinicalinfo_filt,  main = "Survival curve",
           font.main = c(16, "bold", "darkblue"),
           font.x = c(14, "bold.italic", "red"),
           font.y = c(14, "bold.italic", "darkred"),
           font.tickslab = c(12, "plain", "darkgreen"),legend.title="Pathological Stage")

# Customized survival curves

# Customized survival curves
require("survival")
time=clinicalinfo_filt$final_time; status = clinicalinfo_filt$recordedStatus
fit<- survfit(Surv(time, status) ~ clinicalinfo_filt$clinical_stage, data = clinicalinfo_filt)
unique( clinicalinfo_filt$clinical_stage)
# Fit survival curves
#%%%%%%%%%%%%%%%%%%%%%%%
library(survival)
fit_clinical_stage<- survfit(Surv(time, status) ~( clinicalinfo_filt$clinical_stage), data = clinicalinfo_filt)
fit_all <- survfit(Surv(time, status) ~ 1, data = clinicalinfo_filt)

# Visualize
#%%%%%%%%%%%%%%%%%%%%%%
library(survminer)
plot_stage <- ggsurvplot(fit_clinical_stage)

sum_all <- surv_summary(fit_all)
plot_stage$plot+geom_step(data = sum_all,
                          mapping = aes(x = time, y = surv),
                          color = "black")




survdiff(formula=mySurv~clinicalinfo_filt$clinical_stage,rho=0)##Chisq= 1.3  on 2 degrees of freedom, p= 0.521###log-rank test (more weight on higher values of time)##compare the hazards accross groups
survdiff(formula=mySurv~clinicalinfo_filt$clinical_stage,rho=1)##Chisq= 0.5  on 2 degrees of freedom, p= 0.789###generalized wilcoxon test (more weight on lower values of time)
survdiff(formula=mySurv~clinicalinfo_filt$clinical_stage,rho=1.5)##Chisq= 0.3  on 2 degrees of freedom, p= 0.867 ###Tarone-Ware test (in between log-rank and wilcoxon)


ADENO_samplename=read.csv("ADENO_SAMPLEID_su.csv",header=F,sep=",") ###with cpgnames (rows) and samples (columns)
SCC_samplename=read.csv("SCC_SAMPLEID_su.csv",header=F,sep=",")
ADENO_samplename$V1<-gsub(".", "-", ADENO_samplename$V1[], fixed = TRUE)##replace "dot"(" .") by "hypen" ("-") since all the mirnas of mirna target prediction tools are in the format of "-" it is not supported in write.table
SCC_samplename$V1<-gsub(".", "-", SCC_samplename$V1[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") since all the mirnas of mirna target prediction tools are in the format of "-" it is not supported in write.table


ADENO_samplename$V1 = substr(ADENO_samplename$V1,1,nchar(ADENO_samplename$V1)-3)##discard last three characters e.g., "-01", "-11", "-06", etc.
SCC_samplename$V1 = substr(SCC_samplename$V1,1,nchar(SCC_samplename$V1)-3)##discard last three characters e.g., "-01", "-11", "-06", etc.



allsubtypes<-rbind(ADENO_samplename,SCC_samplename)
allsubtypes.cls<-c(rep("ADENO",nrow(ADENO_samplename)),rep("SCC",nrow(allsubtypes)-nrow(ADENO_samplename)))###classinformation rest (2) vs squamous (1) 


allsubtypes.all<-cbind(allsubtypes,allsubtypes.cls)  ##nrow=275
colnames(allsubtypes.all)[1]<-"samplename"



TCGACESC_mirdat<-read.delim("HiSeqV2.csv")
TCGACESC_mirdat1<-TCGACESC_mirdat[,2:ncol(TCGACESC_mirdat)]
rownames(TCGACESC_mirdat1)<-TCGACESC_mirdat$sample
colnames(TCGACESC_mirdat1)<-gsub(".", "-", colnames(TCGACESC_mirdat1)[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") since all the matched samples in "allsubtypes" are in the format of "-" 
colnames(TCGACESC_mirdat1)<-substr(colnames(TCGACESC_mirdat1)[],1,nchar(colnames(TCGACESC_mirdat1)[])-3)##discard last three characters
TCGACESC_mirdat2<-TCGACESC_mirdat1[match(genename2,rownames(TCGACESC_mirdat1)),]##select the data of the 9 dysregulated mirnas only (for Squamous subtype)#ncol=183
TCGACESC_mirdat3<-TCGACESC_mirdat2

TCGACESC_mirdat3fd<-as.data.frame(TCGACESC_mirdat3,stringsAsFactors = TRUE)#ncol=313
ncol(TCGACESC_mirdat3fd)##313

TCGACESC_mirdat3f<-TCGACESC_mirdat3fd[,which(colnames(TCGACESC_mirdat3fd) %in% allsubtypes.all$samplename)]

class.TCGACESC_mirdat3f<-allsubtypes.all[match(colnames(TCGACESC_mirdat3f), allsubtypes.all$samplename),]#nrow=313 #in the order of the samples of "TCGAPAAD_mirdat3f" matrix
nrow(class.TCGACESC_mirdat3f)##275

class.TCGACESC_mirdat3f.ord<-(class.TCGACESC_mirdat3f)
class.TCGACESC_mirdat3f.ord.cls1<-as.data.frame(class.TCGACESC_mirdat3f.ord$allsubtypes.cls) #select the class information of "class.TCGAPAAD_mirdat3f.ord"
class.TCGACESC_mirdat3f.ord.cls=na.omit(class.TCGACESC_mirdat3f.ord.cls1)
nrow(class.TCGACESC_mirdat3f.ord.cls)


###########row-wise (gene-wise) standardization of each TCGA data#######################################

TCGACESC_mirdat3f_scaleddat1<- t(apply(TCGACESC_mirdat3f, 1, function(x)(x-min(x))/(max(x)-min(x))))##where the second argument 1 tells apply to work with rows.##min-max normalization
select <- !is.na(class.TCGACESC_mirdat3f.ord.cls1)
sum(select)
TCGACESC_mirdat3f_scaleddat <- TCGACESC_mirdat3f_scaleddat1[, select]
ncol(TCGACESC_mirdat3f_scaleddat)





TCGAPAAD_mirdat2<-TCGACESC_mirdat3f_scaleddat
TCGAPAAD_mirdat3<-TCGAPAAD_mirdat2[,match(clinicalinfo_filt$bcr_patient_barcode,colnames(TCGAPAAD_mirdat2))]##select the data of the matched samples only (i.e., "clinicalinfo_filt$bcr_patient_barcode"=175)
TCGAPAAD_mirdat4<-as.data.frame(TCGAPAAD_mirdat3,stringsAsFactors = TRUE)
ncol(TCGAPAAD_mirdat4)##268

###################Overall survival analysis of each dysregulated mirnas by Kaplan-Meier Curves#####
####find the z-score of each mirna, and find which sample (patient) is up-regulated (higher expression) and which sample (patient) id down-regulated (lower expression)###########
zscore_DEMIR_squavsrest=matrix(nrow=nrow(TCGAPAAD_mirdat4),ncol=ncol(TCGAPAAD_mirdat4))
zscore_DEMIR_squavsrest[]<-999###put a random value for initialization as it will be updated in the next loop
for(i in 1:nrow(TCGAPAAD_mirdat4))
{
  avg=mean(as.numeric(TCGAPAAD_mirdat4[i,]))
  std=sd(as.numeric(TCGAPAAD_mirdat4[i,]))
  for(j in 1:ncol(TCGAPAAD_mirdat4))
  {
    zscore_DEMIR_squavsrest[i,j]<-(as.numeric(TCGAPAAD_mirdat4[i,j])-avg)/std
  }
}
zscore_DEMIR_squavsrest1<-as.data.frame(zscore_DEMIR_squavsrest)
rownames(zscore_DEMIR_squavsrest1)<-rownames(TCGAPAAD_mirdat4)
colnames(zscore_DEMIR_squavsrest1)<-colnames(TCGAPAAD_mirdat4)

####for all the dysregulated mirnas###
exppattern_DEMIR_squavsrest<-matrix( nrow=nrow(TCGAPAAD_mirdat4), ncol=ncol(TCGAPAAD_mirdat4))

#exppattern_DEMIR_squavsrest<-t(apply(zscore_DEMIR_squavsrest1, 1, function(x) ifelse(x > 0,"high expression","low expression")))##using mean
####exppattern_DEMIR_squavsrest<-t(apply(zscore_DEMIR_squavsrest1, 1, function(x) ifelse(x > 1.96,"Over expression",ifelse(x < -1.96,"under expression","Non-differential expression"))))##need to be changed
# for(i in 1:nrow(zscore_DEMIR_squavsrest1)) ##using median
# {
#   med_mir<-median(as.numeric(zscore_DEMIR_squavsrest1[i,]))
#   for(j in 1:ncol(zscore_DEMIR_squavsrest1))
#   {
#      exppattern_DEMIR_squavsrest[i,j]<-ifelse(zscore_DEMIR_squavsrest1[i,j]>med_mir,"high expression","low expression")
#   }
# }


require(statar)
for(i in 1:nrow(zscore_DEMIR_squavsrest1)) ##using x-tile r
{
  gr34=xtile(zscore_DEMIR_squavsrest1[i,], n = 2)
  for(j in 1:ncol(zscore_DEMIR_squavsrest1))
  {
    exppattern_DEMIR_squavsrest[i,j]<-ifelse(gr34[j]==1,"low expression","high expression")
  }
}
rownames(exppattern_DEMIR_squavsrest)<- rownames(zscore_DEMIR_squavsrest1)
zscore_DEMIR_squavsrest1<-as.data.frame(zscore_DEMIR_squavsrest1)



#mySurv<-Surv(time=clinicalinfo_filt$days_to_death, event = clinicalinfo_filt$recordedStatus)

# 1. Determine the optimal cutpoint of variables for single gene
zscore_DEMIR_squavsrest1_trans<-t(zscore_DEMIR_squavsrest1)
days_to_death<-clinicalinfo_filt$OS.time
recordedStatus<-clinicalinfo_filt$recordedStatus
mysurv_vec_forcutofffinding<-cbind(days_to_death,recordedStatus,zscore_DEMIR_squavsrest1_trans)
names(mysurv_vec_forcutofffinding)[1]<-"days_to_death"
names(mysurv_vec_forcutofffinding)[2]<-"recordedStatus"

mysurv_vec_forcutofffinding<-data.frame(mysurv_vec_forcutofffinding)
mysurv_vec_forcutofffinding<-mysurv_vec_forcutofffinding[which(!is.na(mysurv_vec_forcutofffinding$days_to_death)),]
res.cut <- surv_cutpoint(mysurv_vec_forcutofffinding, time = "days_to_death", event = "recordedStatus",
                         variables = colnames(mysurv_vec_forcutofffinding)[3:ncol(mysurv_vec_forcutofffinding)])

summary(res.cut)
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
genename2



####for 1st mirna#####
library("survival")
myfit <- survfit(Surv(days_to_death, recordedStatus) ~FGF9, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[1],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)


####for 2nd mirna#####
##tt<-res.cut[,(2+2)]
myfit <- survfit(Surv(days_to_death, recordedStatus) ~FGF18, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[2],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 3rd mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~PPP1R9A, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[3],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 4th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~ERBB4, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[4],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 5th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~DCDC2, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[5],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 6th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~TOX3, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[6],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 7th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~ARMC3, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[7],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 8th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~DNALI1, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[8],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

####for 9th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~RGL3, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[9],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)
####for 10th mirna#####
myfit <- survfit(Surv(days_to_death, recordedStatus) ~ENPP3, data = res.cat)
ggsurvplot(myfit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title=rownames(zscore_DEMIR_squavsrest1)[10],xlab = "Overall survival time (days to death)",
           risk.table.height=.20)

######for integrated signature (10 mrna together)##########

####find multiple-genetic score(MGS)=sum(xi)-sum(xj)########
##see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0045894 ####

MGS.10mir<-vector("integer", length=nrow(zscore_DEMIR_squavsrest1_trans))
for(i in 1: nrow(zscore_DEMIR_squavsrest1_trans))
{
  set.seed(300)###to keep the randomly picked training samples same everytime for the specified seed value
  cl <- kmeans(zscore_DEMIR_squavsrest1_trans[i,], 2)##here second parameter "2" is the number of cluster
  ##plot(zscore_DEMIR_squavsrest1_trans[i,], col = cl$cluster)
  #points(cl$centers, col = 1:2, pch = 8, cex = 2)
  gene.clust<-cl$cluster
  
  MGS.10mir[i]<- abs(mean(zscore_DEMIR_squavsrest1_trans[i,which(gene.clust==1)])-mean(zscore_DEMIR_squavsrest1_trans[i,which(gene.clust==2)])) 
}


# 2. Determine the optimal cutpoint of variables for multiple-gene signature

mysurv_vec_forcutofffinding.integrated<-cbind(days_to_death,recordedStatus,MGS.10mir)
mysurv_vec_forcutofffinding.integrated<-data.frame(mysurv_vec_forcutofffinding.integrated)
names(mysurv_vec_forcutofffinding.integrated)[3]<-"Ten_mRNA_Signature"

mysurv_vec_forcutofffinding.integrated<-mysurv_vec_forcutofffinding.integrated[which(!is.na(mysurv_vec_forcutofffinding.integrated$days_to_death)),]
res.cut.integrated <- surv_cutpoint(mysurv_vec_forcutofffinding.integrated, time = "days_to_death", event = "recordedStatus",
                                    variables = colnames(mysurv_vec_forcutofffinding.integrated)[3])

summary(res.cut.integrated)
# 3. Categorize variables
res.cat.integrated <- surv_categorize(res.cut.integrated)
head(res.cat.integrated)

myfit.integrated <- survfit(Surv(days_to_death, recordedStatus)~ Ten_mRNA_Signature , data = res.cat.integrated)
ggsurvplot(myfit.integrated, data = res.cat.integrated, risk.table = TRUE, conf.int = TRUE, pval=TRUE, legend.title="Ten-mRNA Signature",xlab = "Overall survival time (days to death)",
           risk.table.height=.20)





##############survivalROC###############
############for 1st mirna##########
####\references{Heagerty, P.J., Lumley, T., Pepe, M. S. (2000)
###Time-dependent ROC Curves for Censored Survival Data and a Diagnostic
###Marker \emph{Biometrics}, \bold{56}, 337 -- 344} 
###\author{Patrick J. Heagerty }
#####https://www.rdocumentation.org/packages/survivalROC/versions/1.0.1/topics/survivalROC ###
#####https://cran.r-project.org/web/packages/survivalROC/survivalROC.pdf ####


mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$FGF9))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[1],", method =  Kaplan-Meier survival estimate \n Year = 1"))




############for 2nd mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$FGF18))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[2],", method =  Kaplan-Meier survival estimate \n Year = 1"))



############for3rd mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$PPP1R9A))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[3],", method =  Kaplan-Meier survival estimate \n Year = 1"))


############for 4th mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$ERBB4))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[4],", method =  Kaplan-Meier survival estimate \n Year = 1"))


############for 5th mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$DCDC2))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[5],", method =  Kaplan-Meier survival estimate \n Year = 1"))


############for 6th mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$TOX3))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[6],", method =  Kaplan-Meier survival estimate \n Year = 1"))


############for 7th mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$ARMC3))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,3)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[7],", method =  Kaplan-Meier survival estimate \n Year = 1"))

abline(0,1)
############for 8th mirna##########
mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$DNALI1))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM
abline(0,1)
##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,2)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[8],", method =  Kaplan-Meier survival estimate \n Year = 1"))

abline(0,1)

#############for Integrated 9 mirna-signature ########

mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$RGL3))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM
abline(0,1)
##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, cut.values=NULL, method="KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,3)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[9],", method =  Kaplan-Meier survival estimate \n Year = 1"))
abline(0,1)
#############for Integrated 10th mirna-signature ########

mayo1<-as.data.frame(cbind(res.cat$days_to_death,res.cat$recordedStatus,mysurv_vec_forcutofffinding$ENPP3))
colnames(mayo1)<-c("event_time","censor","marker_score")
nrow(mayo1)##258

nobs <- NROW(mayo1)
cutoff <- 365
##METHOD = KM

##mayo11<-mayo1[which(!is.na(mayo1$event_time)),]
mayo11<-mayo1
NROW(mayo11)##59

Mayo4.21= survivalROC(Stime=mayo11$event_time,
                      status=mayo11$censor,
                      marker = mayo11$marker_score,
                      predict.time = cutoff, span = 0.25*nobs^(-0.20),method = "KM")

plot(Mayo4.21$FP, Mayo4.21$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.21$AUC,3)),
     ylab="TP",main=paste("For",rownames(zscore_DEMIR_squavsrest1)[10],", method =  Kaplan-Meier survival estimate \n Year = 1"))


abline(0,1)

