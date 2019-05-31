setwd("C:/Users/dell/Desktop/pamr_svm -RF")

library(pamr)###for pamr
library(caret)
library(rpart)
library(pracma)#for isempty
library(AUC)
library(e1071)
library(DrugClust)###for createFolds (newly added)
library(ROCR)

############Squamous sample##############
ADENO_samplename=read.csv("ADENO_SAMPLEID.csv",header=F,sep=",") ###with cpgnames (rows) and samples (columns)
SCC_samplename=read.csv("SCC_SAMPLEID.csv",header=F,sep=",")
ADENO_samplename$V1<-gsub(".", "-", ADENO_samplename$V1[], fixed = TRUE)##replace "dot"(" .") by "hypen" ("-") since all the mirnas of mirna target prediction tools are in the format of "-" it is not supported in write.table
SCC_samplename$V1<-gsub(".", "-", SCC_samplename$V1[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") since all the mirnas of mirna target prediction tools are in the format of "-" it is not supported in write.table


ADENO_samplename$V1 = substr(ADENO_samplename$V1,1,nchar(ADENO_samplename$V1)-3)##discard last three characters e.g., "-01", "-11", "-06", etc.
SCC_samplename$V1 = substr(SCC_samplename$V1,1,nchar(SCC_samplename$V1)-3)##discard last three characters e.g., "-01", "-11", "-06", etc.



allsubtypes<-rbind(SCC_samplename,ADENO_samplename)
allsubtypes.cls<-c(rep("SCC",nrow(SCC_samplename)),rep("Adeno",nrow(allsubtypes)-nrow(SCC_samplename)))###classinformation rest (2) vs squamous (1) 


allsubtypes.all<-cbind(allsubtypes,allsubtypes.cls)  ##nrow=275
colnames(allsubtypes.all)[1]<-"samplename"


############select the dysregulated (upregulated and downregulated mirnas) for squamous subtype##########
URMIRinfo_adeno_vs_scc=read.csv("cervical.ADENO.vs.SCC.URmiRNA.csv",header=T,sep=",")  ###with cpgnames (rows) and samples (columns)
DRMIRinfo_adeno_vs_scc=read.csv("cervical.ADENO.vs.SCC.DRmiRNA.csv",header=T,sep=",") ###with cpgnames (rows) and samples (columns)

URMIRinfo_adeno_vs_scc$X<-gsub(".", "-", URMIRinfo_adeno_vs_scc$X[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") since all the mirnas of "CESC_mirna.csv" are in the format of "-" 
DRMIRinfo_adeno_vs_scc$X<-gsub(".", "-", DRMIRinfo_adeno_vs_scc$X[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") since all the mirnas of "CESC_mirna.csv" are in the format of "-" 
AllDEMIRinfo_adeno.vs.scc=rbind(URMIRinfo_adeno_vs_scc,DRMIRinfo_adeno_vs_scc)

TCGACESC_mirdat<-read.delim("genomic_miRNA.csv",header=TRUE,sep=",")###TCGA CESC mirna data
TCGACESC_mirdat1<-TCGACESC_mirdat[,2:ncol(TCGACESC_mirdat)]
rownames(TCGACESC_mirdat1)<-TCGACESC_mirdat$sample
colnames(TCGACESC_mirdat1)<-gsub(".", "-", colnames(TCGACESC_mirdat1)[], fixed = TRUE)##replace "dot"(".") by "hypen" ("-") since all the matched samples in "allsubtypes" are in the format of "-" 
colnames(TCGACESC_mirdat1)<-substr(colnames(TCGACESC_mirdat1)[],1,nchar(colnames(TCGACESC_mirdat1)[])-3)##discard last three characters e.g., "-01", "-11", "-06", etc. to make match with all matched subtypes (i.e., "allsubtypes$V1")
TCGACESC_mirdat2<-TCGACESC_mirdat1[match(AllDEMIRinfo_adeno.vs.scc$X,rownames(TCGACESC_mirdat1)),]##select the data of the 9 dysregulated mirnas only (for Squamous subtype)#ncol=183
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
#############SVM,PAM,Random Forest##########################################
library(doSNOW)
library(MASS,quietly = TRUE)
library(caTools)
set.seed(123)
library(caret)
DataFrame<-data.frame(x=t(TCGACESC_mirdat3f_scaleddat),y=as.factor(class.TCGACESC_mirdat3f.ord.cls$`class.TCGACESC_mirdat3f.ord$allsubtypes.cls`))

SampleID=createFolds(DataFrame$y,k=5,returnTrain = FALSE,list = TRUE)
for(i in 1:length(SampleID))
  {
  Test_ID=SampleID[[i]]
  
       # Train_id=unlist(setdiff(SampleID,Test_ID));
        
        Test_MAt=DataFrame[Test_ID,];
        train_MAT=DataFrame[-Test_ID,];
        ControlParameters<-trainControl(method="cv",number=10,savePredictions = TRUE,classProbs = TRUE)
        modelRandom<-train(y~.,data = train_MAT,method="rf",trControl=ControlParameters)
        modelRandom$modelInfo
        
        print(modelRandom$results)
        pred<-predict(modelRandom,Test_MAt,type = 'prob')
        x<-ifelse(Test_MAt$y =="SCC", 1,0)
        labels = as.factor(x)
        predicted_label<-predict(modelRandom,Test_MAt)
        tab<-table(predicstions=predicted_label,actual=Test_MAt$y)
        print(tab)
        predictions <- prediction(pred$SCC, labels)
        perf <- performance(predictions, measure = "tpr", x.measure = "fpr")
        plot(perf, col=rainbow(10))
        h = roc(labels = labels,predictions = pred$SCC)
        print(paste("AUC:",auc(h)))
        
}

library(MLmetrics)
prob <- predict(modelRandom,Test_MAt) # Prediction
#F1_Score(y_pred = y1, y_true =x , positive = "0")
pred1 <- ifelse(modelRandom$pred$pred=="Adeno",0,1)
x<-ifelse(modelRandom$pred$obs =="Adeno", 0,1)
F1_Score(y_pred = pred1, y_true = x, positive = "1")

View(Test_MAt)

















