library(ggplot2)
library(data.table)
setwd("/Users/baifeng/Desktop/BGIproject/ovarian_DMR/ovary_cancer_result_2016_10_12")
ov.sample <- data.frame(fread("sample.all.list"))
ov.clean.data <- read.table("Clean_Rate_stat.xls",header = T,sep="\t")
ov.clean.data <- ov.clean.data[ov.clean.data$Sample%in%c(ov.sample[,1],ov.sample[,2]),]
ov.clean.data[ov.clean.data$Sample%in%ov.sample[,1],"Type"] <- "Tumor"
ov.clean.data[ov.clean.data$Sample%in%ov.sample[,2],"Type"] <- "Normal"

ov.bisulfite.rate <- read.table("BS_conversion_rate.xls",header = T,sep="\t")
rownames(ov.bisulfite.rate) <- ov.bisulfite.rate[,1]
ov.mapping.rate <- read.table("Mapping_rate.xls",header = T,sep="\t")
rownames(ov.mapping.rate) <- ov.mapping.rate[,1]
ov.cpg.coverage <- read.table("cpg.coverage.xls",header=T,sep="\t")
colnames(ov.cpg.coverage) <- c("Sample","Theoretical","1X","4X","10X")
rownames(ov.cpg.coverage) <- ov.cpg.coverage[,1]

ov.sig.rrbs <- data.frame(fread("RRBS-Step1.update3.tcga.txt"))
ov.sig.rrbs.cgi <- unique(ov.sig.rrbs[,c("seqnames","start","end","UCSC_REFGENE_GROUP","RELATION_TO_UCSC_CPG_ISLAND","mean.meth.diff")])
ov.sig.rrbs.cgi[ov.sig.rrbs.cgi$mean.meth.diff < 0,"METH_TYPE"] <- "Hypo"
ov.sig.rrbs.cgi[ov.sig.rrbs.cgi$mean.meth.diff > 0,"METH_TYPE"] <- "Hyper"
draw_pie <- function (data) {
  library(scales)
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  ggplot(data, aes(x="",fill=RELATION_TO_UCSC_CPG_ISLAND))+geom_bar(width = 1,position = "fill")+coord_polar("y", start=0)+blank_theme+theme(axis.text.x=element_blank())+facet_grid(METH_TYPE~UCSC_REFGENE_GROUP)# +geom_text(aes(y = 1/5 + c(0, cumsum(value)[-length(value)]),label = percent(value/100)), size=5)
}
draw_pie(data = ov.sig.rrbs.cgi)

heatmap.rrbs <- function(data=ov.sig.rrbs){
  hyper.rrbs <- ov.sig.rrbs[,c("gene","sample","Hyper_meth_sample_number","Hyper.group.mean.meth.diff")]
  hypo.rrbs <- ov.sig.rrbs[,c("gene","sample","Hypo_meth_sample_number","Hypo.group.mean.meth.diff")]
  hyper.rrbs.T10.gene <- head(unique(hyper.rrbs[do.call(order,-hyper.rrbs["Hyper_meth_sample_number"]),"gene"]),10)
  hyper.rrbs.T10 <- hyper.rrbs[hyper.rrbs$gene%in%hyper.rrbs.T10.gene,-3]
  hypo.rrbs.T10.gene <- head(unique(hypo.rrbs[do.call(order,-hypo.rrbs["Hypo_meth_sample_number"]),"gene"]),10)
  hypo.rrbs.T10 <- hypo.rrbs[hypo.rrbs$gene%in%hypo.rrbs.T10.gene,-3]
  colnames(hyper.rrbs.T10) <- c("gene","sample","Meth.Diff")
  hyper.rrbs.T10 <- hyper.rrbs.T10[do.call(order,-hyper.rrbs.T10["Meth.Diff"]),]
  hyper.rrbs.T10 <- hyper.rrbs.T10[!duplicated(hyper.rrbs.T10[,c("gene","sample")]),]
  colnames(hypo.rrbs.T10) <- c("gene","sample","Meth.Diff")
  hypo.rrbs.T10 <- hypo.rrbs.T10[do.call(order,hypo.rrbs.T10["Meth.Diff"]),]
  hypo.rrbs.T10 <- hypo.rrbs.T10[!duplicated(hypo.rrbs.T10[,c("gene","sample")]),]
  rrbs.T20 <- rbind(hyper.rrbs.T10,hypo.rrbs.T10)
  rrbs.T20 <- rrbs.T20[!is.na(rrbs.T20$Meth.Diff),]
  rrbs.T20 <- dcast(rrbs.T20[,c("gene","sample","Meth.Diff")],gene~sample)
  rownames(rrbs.T20) <- rrbs.T20[,1]
  rrbs.T20[,1] <- NULL
  pheatmap(rrbs.T20)
}


sig.cpg.coverage <- function(rrbs) {
  cpg.len.bed <- read.table("hg19_cpgisland_all.bed",header=F,sep="\t")
  cpg.len <- sum(cpg.len.bed$V5)
  cpg.coverage <- c()
  for (i in unique(ov.rrbs$sample)){
    rrbs.cpg.length <- sum(ov.rrbs[ov.rrbs$sample==i & ov.rrbs$RELATION_TO_UCSC_CPG_ISLAND=="Island","length"])
    cpg.coverage <- rbind(cpg.coverage,c(i,rrbs.cpg.length,rrbs.cpg.length*100/cpg.len))
  }
  cpg.coverage <- data.frame(cpg.coverage)
  colnames(cpg.coverage) <- c("sample","cpg.length","cpg.coverage")
  cpg.coverage
}
ov.sig.cpg.coverage <- cpg.coverage(rrbs = ov.rrbs)
ov.clean.data[,"bisulfite.rate"] <- ov.bisulfite.rate[as.character(ov.clean.data[,"Sample"]),2]
ov.clean.data[,"mapping.rate"] <- ov.mapping.rate[as.character(ov.clean.data[,"Sample"]),"Mapping_rate..."]
ov.clean.data <- cbind(ov.clean.data,ov.cpg.coverage[as.character(ov.clean.data[,"Sample"]),-1])
test <- melt(ov.clean.data[,c("Sample","bisulfite.rate","mapping.rate","Type","Theoretical","1X","4X","10X")])
colnames(test) <- c("Sample","Type","variable","Fraction")
ggplot(ov.clean.data,aes(x=Type,y=Clean_Bases,color=Type))+geom_boxplot()+geom_jitter()+labs(title="Distribution of Data Output")
ggplot(test,aes(x=variable,y=Fraction,fill=Type))+geom_boxplot()+labs(title="Distribution of Base Quality in Normal and Tumor Samples")

Hypo_HiEx <- data.frame(fread("Hypo_HiEx.list"))
Hypo_HiEx <- Hypo_HiEx[,c("gene","Hypo.group.mean.meth.diff","Hyper_exp_mean")]
colnames(Hypo_HiEx) <- c("gene","Meth","Exp")
Hypo_HiEx <- Hypo_HiEx[do.call(order,Hypo_HiEx["Meth"]),]
Hypo_HiEx <- Hypo_HiEx[!duplicated(Hypo_HiEx[,c("gene","Exp")]),]
Hyper_LoEx <- data.frame(fread("Hyper_LoEx.list"))
Hyper_LoEx <- Hyper_LoEx[,c("gene","Hyper.group.mean.meth.diff","Hypo_exp_mean")]
colnames(Hyper_LoEx) <- c("gene","Meth","Exp")
Hyper_LoEx <- Hyper_LoEx[do.call(order,-Hyper_LoEx["Meth"]),]
Hyper_LoEx <- Hyper_LoEx[!duplicated(Hyper_LoEx[,c("gene","Exp")]),]

meth.exp <- rbind(Hyper_LoEx,Hypo_HiEx)
gene.name <- meth.exp$gene
plot(Meth~Exp,meth.exp)+abline(lm(Meth~Exp,meth.exp))+mtext("Adjusted R-squared = 0.8026   p-value: < 2.2e-16")
meth.exp <- melt(meth.exp)
#meth.exp <- meth.exp[do.call(order,meth.exp("value")),]
#meth.exp <- meth.exp[!duplicated(meth.exp[,c("gene","variable")]),]
meth.exp[meth.exp$variable=="Meth","value"] <- meth.exp[meth.exp$variable=="Meth","value"]
for (i in 1:length(gene.name)){
  meth.exp[meth.exp$gene==gene.name[i],"pos"] <- i
}
ggplot(meth.exp,aes(x=pos,y=value,fill=variable))+geom_bar(stat = "identity")+scale_x_discrete(limits = gene.name)+theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=8,face="bold"))+labs(title = "Genes with correlation between expression and methylation")

ov.exp <- load("Expression/Ovarian_Cancer_Gene_Expression.RData")
ov.exp <- cbind(normal,tumor)
outcome <- c(rep(1,ncol(normal)),rep(2,ncol(tumor)))
ov.exp.dat <- list(n=ov.exp,y=outcome,pair=TRUE,type="twoclass") 
ov.exp.poisson <- function(data){
  library(PoissonSeq)
  ov.exp.ps <- PS.Main(ov.exp)
  ov.exp.ps[ov.exp.ps$pval < 0.01,"significance"] <- TRUE
  ov.exp.ps[ov.exp.ps$pval >= 0.01,"significance"] <- FALSE
  ov.exp.ps.cut <- ov.exp.ps[ov.exp.ps$pval < 0.01,]
  ov.exp.ps.cut.down <- ov.exp.ps.cut[ov.exp.ps.cut$log.fc < -1,]
  ov.exp.ps.cut.down[,"fc"] <- 1/(2**ov.exp.ps.cut.down$log.fc)
  ov.exp.ps.cut.up <- ov.exp.ps.cut[ov.exp.ps.cut$log.fc > 1,]
  ov.exp.ps.cut.up[,"fc"] <- 2**ov.exp.ps.cut.up$log.fc
  ov.exp.ps.cut.fc <- rbind(ov.exp.ps.cut.up,ov.exp.ps.cut.down)
  interval <- c(2,2.5,3)
  ov.exp.ps.cut.fc.stat <- c()
  for (i in interval){
    ov.exp.ps.cut.fc.stat <- rbind(ov.exp.ps.cut.fc.stat,c(i,nrow(ov.exp.ps.cut.fc[ov.exp.ps.cut.fc$fc>i & ov.exp.ps.cut.fc$log.fc<0,]),nrow(ov.exp.ps.cut.fc[ov.exp.ps.cut.fc$fc>i & ov.exp.ps.cut.fc$log.fc>0,])))
  }
  ov.exp.ps.cut.fc.stat <- data.frame(ov.exp.ps.cut.fc.stat)
  colnames(ov.exp.ps.cut.fc.stat) <- c("FoldChange","Down","Up")
  ov.exp.ps.cut.fc.stat
  
  ggplot(ov.exp.ps[-log10(ov.exp.ps$fdr) < 0.2,],aes(x=log.fc,y=-log10(pval),colour=significance))+geom_point(alpha=1,size=1.5)+scale_color_manual(values = c("black","red"))+xlab('log2 (fold change)')+ylab('-log10 (pval)')+theme(legend.position = 'top',axis.text = element_text(color='black'),panel.background = element_rect(fill = 'transparent'),panel.grid = element_line(color = 'grey'),panel.border = element_rect(fill = 'transparent',color = 'black'),axis.title = element_text(size = 15))
  
  log2Ratio.2fc <- log2Ratio[rownames(ov.exp.ps.cut.fc),]
  pheatmap(log2Ratio.2fc,show_rownames=F)
  
}





Hyper_HiEx <- data.frame(fread("Hyper_HiEx.region"))
Hyper_LoEx <- data.frame(fread("Hyper_LoEx.region"))
Hypo_HiEx <- data.frame(fread("Hypo_HiEx.region"))
Hypo_LoEx <- data.frame(fread("Hypo_LoEx.region"))
select.region <- rbind(Hyper_HiEx,Hyper_LoEx,Hypo_HiEx,Hypo_LoEx)
ggplot(select.region,aes(x=mean.meth.diff,y=-log10(DMR.qvalue)))+geom_point()












library(TCGAbiolinks)
setwd("/Users/baifeng/Desktop/BGIproject/ovarian_DMR/ovary_cancer_result_2016_10_12/TCGAmethylation")
query.met <- GDCquery(project = "TCGA-OV", 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)

ov.tcga.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "ov.TCGA.DNAmet.rda",
                      summarizedExperiment = TRUE)

query.exp <- GDCquery(project = "TCGA-OV", 
                     legacy = TRUE,
                     data.category = "Gene expression",
                     data.type = "Gene expression quantification",
                     platform = "Illumina HiSeq", 
                     file.type = "results")
GDCdownload(query.exp)

ov.tcga.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "ov.TCGA.Exp.rda")

# na.omit
ov.tcga.met <- subset(ov.tcga.met,subset = (rowSums(is.na(assay(acc.met))) == 0))

# Volcano plot
ov.tcga.met <- TCGAanalyze_DMR(ov.tcga.met, groupCol = "subtype_MethyLevel",
                           group1 = "CIMP-high",
                           group2="CIMP-low",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "CIMP-highvsCIMP-low_metvolcano.png")

