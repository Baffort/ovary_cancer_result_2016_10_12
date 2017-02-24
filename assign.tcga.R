#P_value<-read.table("RRBS-Step1.update3.txt",header=T,sep="\t")
P_value<-read.table("RRBS-Step1.update3.txt",header=T,sep="\t",stringsAsFactors = FALSE)
P_value[,"qvalue"]<-p.adjust(P_value[,"P_value"],method="fdr",length(P_value[,"P_value"]))
#P_value<-P_value[P_value$P_value<=0.05,]
tcga_cor<-read.table("Correlate_Methylation_vs_mRNA_OV-TP_matrix.txt",header=T,sep="\t",stringsAsFactors = FALSE)
#for (g in 1:nrow(P_value)){
#        P_value[g,"TCGA_Corr_Coeff"]<-mean(tcga_cor[tcga_cor$Gene==P_value[g,"gene"],"Corr_Coeff"])
#        P_value[g,"TCGA_Pval"]<-mean(tcga_cor[tcga_cor$Gene==P_value[g,"gene"],"Pval"])
#        P_value[g,"TCGA_Qval"]<-mean(tcga_cor[tcga_cor$Gene==P_value[g,"gene"],"Qval"])
#        P_value[g,"TCGA_Expr_Mean"]<-mean(tcga_cor[tcga_cor$Gene==g,"Expr_Mean"])
#        P_value[g,"TCGA_Meth_Mean"]<-mean(tcga_cor[tcga_cor$Gene==g,"Meth_Mean"])
#}
a<-lapply(1:nrow(P_value),function (g) c(mean(tcga_cor[tcga_cor$Gene==P_value[g,"gene"],"Corr_Coeff"]),mean(tcga_cor[tcga_cor$Gene==P_value[g,"gene"],"Pval"]),mean(tcga_cor[tcga_cor$Gene==P_value[g,"gene"],"Qval"])))
a<-data.frame(do.call(rbind,a))
colnames(a)<-c("TCGA_Corr_Coeff","TCGA_Pval","TCGA_Qval")
P_value<-cbind(P_value,a)
write.table(P_value,"RRBS-Step1.update3.tcga.txt",quote=FALSE,sep="\t",row.names=FALSE)
