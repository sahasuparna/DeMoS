DEGstatNew_sm<-function(datamatrix,pval_cutoff,upperfoldchange_cutoff,lowerfoldchange_cutoff,pval_adj_used,ctl.subtype,exp.subtype,control.size,tissue)
     {
  
          
          print("Yes, we are in the significant gene-finding extraction precedure");
          #grep("\\?",rownames(datamatrix))
         #datamatrix1<- datamatrix[-grep("\\?",rownames(datamatrix)),]
         nrow(datamatrix)
          myCPM<-cpm(datamatrix)
          thresh_myCPM <- myCPM > 0.5 ##Which values in myCPM are greater than 0.5?
          table(rowSums(thresh_myCPM)) ##Summary of how many TRUEs there are in each row
          keeprna <- rowSums(thresh_myCPM) >= 2 ##keep genes that have at least 2 TRUES in each row of thresh
          RNA_Data<-datamatrix
         # grep("\\?",rownames(RNA_Data1))
          #RNA_Data<- RNA_Data1[-grep("\\?",rownames(RNA_Data1)),]
         # nrow(RNA_Data)#19685
          RNA_Data.keeprna <- RNA_Data[keeprna,] ##Subset the rows of countdata to keep the more highly expressed genes
          summary(keeprna)
          dim(RNA_Data.keeprna)#19603 genes 275 samples
      
          plot(myCPM[,1],RNA_Data[,1],ylab="tcga_rnadata[,1]")
          
          # Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
          plot(myCPM[,1],RNA_Data[,1],ylim=c(0,5),xlim=c(0,5),ylab="tcga_rnadata[,1]")
          # Add a vertical line at 0.5 CPM
          abline(v=0.5)
          
          
          #######Convert counts to DGEList object for rna########
          y_r <- DGEList(RNA_Data.keeprna)
          names(y_r)### See what slots are stored in y
          # Library size information is stored in the samples slot
          y_r$samples
          nrow(y_r$samples)#275
          
          
         
          ########*************for ADENO (control) vs cervical squamous cell carcinoma (SCC or experimental or rest) (rna)**********##########################
          ###Voom transform the data
          sample_clslabel_cntl.rest_r <- rep(c(ctl.subtype,exp.subtype),c(control.size, ncol(RNA_Data)- control.size))
          sample_clslabel_cntl.rest_r
          design_cntlvsrest_r<-model.matrix(~ sample_clslabel_cntl.rest_r)
          head(design_cntlvsrest_r,n=15)
          #png("mean-variance-trend.png", 490, 350)
          setEPS()
          postscript(paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DEG.mean_variance_trend.eps",sep=""))
          ###postscript(file=sprintf("imagesDEG/mean_variance_trend%svs%s.eps",ctl.subtype,exp.subtype), onefile=TRUE, horizontal=FALSE)
          #postscript(file=sprintf("imagesDEG/mean-variance-trend.Adenovs%s.eps",subtype), onefile=TRUE, horizontal=FALSE)
          par(mar = rep(5, 4))
          v_cntlvsrest_r <- voom(y_r,design_cntlvsrest_r,plot = TRUE)
          dev.off()
          
          v_cntlvsrest_r
          names(v_cntlvsrest_r)
          v_cntlvsrest_r$E
          logcounts_cntlvsrest_r <- cpm(y_r,log=TRUE)
          
          setEPS()
          postscript(paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DEG.boxplot.normalization.eps",sep=""))
          #postscript(file=sprintf("imagesDEG/boxplot.Adenovs%s.eps",subtype), onefile=TRUE, horizontal=FALSE)
          par(mfrow=c(1,2),mar = c(5, 5, 4, 2) + 0.1)
          #par(mar = rep(2, 4))
          par(cex.axis=0.3)##resizing the boxplot label
          par(cex.lab=0.7)
          par(cex.main=0.8)
          boxplot(logcounts_cntlvsrest_r[,1:25], xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")##for first 25 samples since the total number of samples are very high (=275)
          
          ## Let's add a blue horizontal line that corresponds to the median logCPM
          abline(h=median(logcounts_cntlvsrest_r),col="blue")
          par(cex.axis=0.3)##resizing the boxplot label
          par(cex.lab=0.7)
          par(cex.main=0.8)
          boxplot(v_cntlvsrest_r$E[,1:25], xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
          abline(h=median(v_cntlvsrest_r$E),col="blue")
          dev.off()
          
          ############Testing for differential expression for rna###
          # Fit the linear model
          fit_cntlvsrest_r <- lmFit(v_cntlvsrest_r,design_cntlvsrest_r)
          names(fit_cntlvsrest_r)
          fit_ebayes_cntlvsrest_r <- eBayes(fit_cntlvsrest_r)
          dim(fit_ebayes_cntlvsrest_r)
          options(digits=3)
          top2_cntlvsrest_r <- topTable(fit_ebayes_cntlvsrest_r,number=Inf,coef=2,adjust=pval_adj_used,sort.by="P")###use this if not use "makeContrasts" function
          sum(top2_cntlvsrest_r$adj.P.Val<pval_cutoff)
          
          
          DEGlist_cntlvsrest_r<-top2_cntlvsrest_r[which(top2_cntlvsrest_r$adj.P.Val<pval_cutoff),]
          DEGlist_srt_cntlvsrest_r<-DEGlist_cntlvsrest_r[order(DEGlist_cntlvsrest_r$adj.P.Val, decreasing = FALSE),]
          fURGlist_srt_cntlvsrest_r<-DEGlist_srt_cntlvsrest_r[which(DEGlist_srt_cntlvsrest_r$logFC>upperfoldchange_cutoff),]##FC filtering for upregulation
          fDRGlist_srt_cntlvsrest_r<-DEGlist_srt_cntlvsrest_r[which(DEGlist_srt_cntlvsrest_r$logFC<lowerfoldchange_cutoff),]##FC filtering for downregulation
          fDEGlist_pvalFCfilt_cntlvsrest_r<-DEGlist_srt_cntlvsrest_r[which(DEGlist_srt_cntlvsrest_r$logFC>upperfoldchange_cutoff | DEGlist_srt_cntlvsrest_r$logFC<(lowerfoldchange_cutoff)),]###final DEGlist for ADENO vs SCC
         
          
         
          n1=nrow(fDEGlist_pvalFCfilt_cntlvsrest_r)
          n2=nrow(fURGlist_srt_cntlvsrest_r);
          n3= nrow(fDRGlist_srt_cntlvsrest_r);
          print("The total no of significant genes:");print(n1);
          print("The total no of Upregulated genes:"); print(n2);
          print("The total no of Downregulated genes:");print(n3);
          fURGlist_srt_cntlvsrest_r$Gene_Name<-rownames(fURGlist_srt_cntlvsrest_r)
          fDRGlist_srt_cntlvsrest_r$Gene_Name<-rownames(fDRGlist_srt_cntlvsrest_r)
          fDEGlist_pvalFCfilt_cntlvsrest_r$GeneName<-rownames(fDEGlist_pvalFCfilt_cntlvsrest_r)
          a<-list(UG_Name= fURGlist_srt_cntlvsrest_r$Gene_Name,DG_Name=fDRGlist_srt_cntlvsrest_r$Gene_Name)
          print("the Upregulated gene Names are as follows....")
          print(a$UG_Name)
          print("the Downregulated gene Names are as follows....")
          print(a$DG_Name)
          
          write.table(fURGlist_srt_cntlvsrest_r, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".URG.csv",sep=""),sep=",",col.names = TRUE,row.names = FALSE)
          #write.table(fURGlist_srt_cntlvsrest_r, paste("UG_ADENOvs",subtype, ".csv",sep=""), sep=",", row.names=FALSE, col.names=TRUE,quote=FALSE)
          write.table(fDRGlist_srt_cntlvsrest_r, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DRG.csv",sep=""),sep=",",col.names = TRUE,row.names = FALSE)
          #write.table(fDRGlist_srt_cntlvsrest_r,paste("DG_ADENOvs",subtype, ".csv",sep=""), sep=",", row.names=FALSE, col.names=TRUE,quote=FALSE)
          write.table(fDEGlist_pvalFCfilt_cntlvsrest_r, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DEG.csv",sep=""),sep=",",col.names =TRUE,row.names = TRUE)
          #write.table(fDEGlist_pvalFCfilt_cntlvsrest_r,paste("UG+DG_ADENOvs",subtype, ".csv",sep=""), sep=",", row.names=FALSE, col.names=TRUE,quote=FALSE)
         
          top2_cntlvsrest_r$ID<-rownames(top2_cntlvsrest_r)
          
          
          ###volcanoplot####
          setEPS()
          postscript(file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".volcanoplot.eps",sep=""))
          par(mfrow=c(1,1),mar = c(5, 5, 4, 2) + 0.1)
          yaxis.max<-max(-log10(top2_cntlvsrest_r$adj.P.Val))+3
          with(top2_cntlvsrest_r, plot(logFC, -log10(adj.P.Val), pch=10, main="Volcano plot", xlim=c(-5,3),ylim=c(0,yaxis.max)))
          with(subset(top2_cntlvsrest_r, adj.P.Val<pval_cutoff & logFC>upperfoldchange_cutoff), points(logFC, -log10(adj.P.Val), pch=yaxis.max,col="red"))
          with(subset(top2_cntlvsrest_r, adj.P.Val<pval_cutoff & logFC<lowerfoldchange_cutoff), points(logFC, -log10(adj.P.Val), pch=yaxis.max,col="green"))
          #with(subset(top2_cntlvsrest_r, adj.P.Val<pval_cutoff & abs(logFC)>upperfoldchange_cutoff), textxy(logFC, -log10(adj.P.Val), labs=ID, cex=.87))###for genename printing on valcanoplot
          dev.off()
          
          
          ####heatmap for fDEG (fURG+fDRG)##########
          heatmap1<-datamatrix[which(rownames(datamatrix) %in% rownames(fDEGlist_pvalFCfilt_cntlvsrest_r)),]
          heatmap_matrix=heatmap1
          row1<-rownames(heatmap1)
          heatmap_matrix<-cbind(V2=row1,heatmap_matrix)
          write.table(heatmap_matrix, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DEG.dataval.csv",sep=" "),sep=",",col.names=TRUE,row.names = FALSE)
          write.table(heatmap_matrix, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".input_WGCNA_dataval.csv",sep=""),sep=",",col.names = TRUE,row.names = TRUE)
          gene_list<-rownames(heatmap1)
          a<-as.data.frame(gene_list)
          gene_number<-rownames(a)
          write.table(gene_list, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DEG.GeneNames.csv",sep=""),sep=",",col.names = FALSE,row.names = FALSE)
          write.table(gene_number, file=paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".DEG.GeneNumbers.csv",sep=""),sep=",",col.names = FALSE,row.names = FALSE)

          setEPS()
          postscript(paste0(tissue,".",ctl.subtype,".vs.",exp.subtype,".heatmap.eps",sep=""))
          expval_fDRGlist_srt_cntlvsrest_r<-heatmap1
          DG<-expval_fDRGlist_srt_cntlvsrest_r
          D_nd1<-data.matrix(DG)#row=gene, col=sample
          D_nd <- scale(t(D_nd1)) ##row=sample, ##col=gene#scaling the data
          D_nd_t <- t(D_nd)#row=gene, col=sample#transpose the data
          round(colMeans(D_nd),1)# Check means are zero and std devs are 1
          apply(D_nd,2,sd)# Check means are zero and std devs are 1
          my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)# Set colours for heatmap, 25 increments
          par(cex.main=0.75) # Shrink title fonts on plot# Plot heatmap with heatmap.2
          d1 <- dist(D_nd,method = "euclidean", diag = FALSE, upper = TRUE)# Calculate distance between experiments (samples) between rows
          round(d1,3)
          d2 <- dist(D_nd_t,method = "euclidean", diag = FALSE, upper = TRUE)# Calculate distance between proteins (genes) between rows
          round(d2,3)
          c1 <- hclust(d1, method = "ward.D2", members = NULL)# Clustering distance between experiments using Ward linkage
          c2 <- hclust(d2, method = "ward.D2", members = NULL)# Clustering distance between proteins using Ward linkage
          # Check clustering by plotting dendrograms
          par(mfrow=c(2,1),cex=0.5) # Make 2 rows, 1 col plot frame and shrink labels
          #plot(c1); plot(c2) # Plot both cluster dendrograms
          library(d3heatmap)
          library(gplots)
          heatmap.2(D_nd_t,                     # Tidy, normalised data
                    Colv=as.dendrogram(c1),     # Experiments clusters in cols
                    Rowv=as.dendrogram(c2),     # Protein clusters in rows
                    trace="none",               # Turn of trace lines from heat map
                    scale = c("row"),
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.5,cexCol=0.75,     # Amend row and column label fonts
                    margins=c(6,8))
          
          dev.off()
          
}       
          
          