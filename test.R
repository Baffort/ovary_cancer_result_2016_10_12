library(data.table)
DMR<-read.table("RRBS-Step1.update2.unique.txt",header=F,sep="\t",stringsAsFactors = FALSE)
exp<-read.table("all.sample.exp.txt",header=T,sep="\t")
rownames(exp)<-exp[,"gene"]
exp[,"gene"]<-NULL
MeanExp<-function(x){
	Hyper_index<-which(x>0)
	Hypo_index<-which(x<=0)
	Hyper_exp_sample_number <- length(Hyper_index)
	Hyper_exp_mean <- mean(x[Hyper_index])
	Hypo_exp_sample_number <- length(Hypo_index)
        Hypo_exp_mean <- mean(x[Hypo_index])
	c(Hyper_exp_sample_number,Hyper_exp_mean,Hypo_exp_sample_number,Hypo_exp_mean)
}
MeanExpData<-data.frame(t(apply(exp,1,MeanExp)))
colnames(MeanExpData)<-c("Hyper_exp_sample_number","Hyper_exp_mean","Hypo_exp_sample_number","Hypo_exp_mean")
exp<-cbind(exp,MeanExpData)
mean_exp<-data.frame(do.call(rbind,lapply(1:nrow(DMR),function (x) MeanExpData[DMR[x,4],])))
DMR<-cbind(DMR,mean_exp)
write.table(DMR,"RRBS-Step1.update2.unique.addMeanExp.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
write.table(exp,"Expression-Step2.addMeanExp.txt",quote=FALSE,sep="\t")
