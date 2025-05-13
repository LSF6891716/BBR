rm(list=ls())
library(dplyr)
library(org.Hs.eg.db)
library(limma) 
library(GEOquery)
library(affy)
library(ggplot2)
library(ggpubr)
library(limma)
library(ggplot2)
library(ggpubr)
keytypes(org.Hs.eg.db)
mycounts<-read.csv("mycounts.csv",header=T)
x<-mycounts$X  
get<-c("SYMBOL","GENENAME")
result<-select(org.Hs.eg.db,
               keys = x, 
               columns=get, 
               keytype="ENSEMBL")
result<-na.omit(result) 
table(duplicated(result$SYMBOL))
result_1<-mapIds(org.Hs.eg.db,
                 keys = x, 
                 column="SYMBOL", 
                 keytype="ENSEMBL")
result_1<-as.data.frame(result_1)
result_1$ENSEMBL<-rownames(result_1)
colnames(result_1)<-c("SYMBOL","ENSEMBL")
options(stringsAsFactors = F) 
rm(list=ls()) 
gse<-getGEO("GSE41446",destdir = ".",  
            getGPL = T,     
            AnnotGPL = T)   
mode(gse) 
gse[[1]]  
exp<-exprs(gse[[1]])  
cli<-pData(gse[[1]])  
group<-c(rep("control",3),rep("hht",3))
GPL<-fData(gse[[1]])  
gpl<-GPL[,c(1,3)]  
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,
                                     function(x)unlist(strsplit(x,"///"))[1]),
                              stringsAsFactors = F)[,1]
exp<-as.data.frame(exp) 
exp$ID<-rownames(exp)  
exp_symbol<-merge(exp,gpl,by="ID")  
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$`Gene symbol`))
exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`Gene symbol`)
table(duplicated(rownames(exp_unique)))
save(gse,exp,cli,GPL,gpl,exp_symbol,exp_unique,group,file="GSE41446.Rdata")
rawdata = ReadAffy()  
eset.rma = rma(rawdata)   
dim(eset.rma)
data = t(data.frame(eset.rma))
save(data,file="raw.Rdata")
g<-ifelse(cli$source_name_ch1=="Kasumi-1_control","control","hht")
options(stringsAsFactors = F)
rm(list=ls())
load("GSE41446.Rdata")
boxplot(exp_unique,cex=0.2,las=2,col="red") 
design<-model.matrix(~0+factor(group))  
colnames(design)<-levels(factor(group))  
rownames(design)<-colnames(exp_unique) 
contrast.matrix<-makeContrasts(hht-control,levels=design) 
fit<-lmFit(exp_unique,design) 
fit2<-contrasts.fit(fit,contrast.matrix) 
fit2<-eBayes(fit2)
options(digits = 4)
DEG<-topTable(fit2,coef=1,n=Inf) 
DEG$group<-ifelse(DEG$P.Value>0.05,"no_change",
                  ifelse(DEG$logFC>1.5,"up",
                         ifelse(DEG$logFC< -1.5,"down","no_change")))
table(DEG$group)
DEG$gene<-rownames(DEG) 
save(DEG,file="DEG.Rdata")
data<-data.frame(gene=exp_unique["PTGS2",],group=group)
ggplot(data=data,aes(x=group,y=gene,fill=group))+geom_boxplot()+
  stat_compare_means()+geom_point(position="jitter")
pp<-function(g){
  library(ggplot2)
  library(ggpubr)
  frame<-data.frame(gene=g,group=group)
  ggplot(data=frame,aes(x=group,y=gene,fill=group))+geom_boxplot()+
    stat_compare_means()+geom_point(position="jitter")
}
pp(exp_unique[2,]) 
pp(exp_unique["PTGS2",]) 
options(stringsAsFactors = F)
rm(list=ls())
load("GSE41446.Rdata")
load("DEG.Rdata")
DEG$p<- -log10(DEG$P.Value)
ggplot(data=DEG,aes(x=logFC,y=p,color=group))+geom_point()+
  theme_bw()  
library(ggrepel)
ggplot(data=DEG,aes(x=logFC,y=p,color=group))+
  geom_point(data=DEG[DEG$P.Value<0.01&abs(DEG$logFC)>2,],size=3)+  
  geom_point(data=DEG[DEG$P.Value>0.01|abs(DEG$logFC)<2,],size=1)+#没变的
  theme_bw()+#白色主题
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+#设置颜色
  ylab("-log10(P.Value)")+ #y标题改名
  xlab("log2FoldChange")+#x标题改名
  geom_text_repel(
    data=DEG[DEG$P.Val<0.0001&abs(DEG$logFC)>3.8,],
    aes(label=gene),
    size=3,
    color="black")
library(pheatmap)  
pheatmap(exp_unique[1:50,],show_rownames = F,show_colnames = F)
x<-DEG$logFC 
names(x)<-rownames(DEG) 
x[1:8]
upgene<-names(tail(sort(x),10))
downgene<-names(head(sort(x),10))
top20<-c(upgene,downgene) 
pheatmap(exp_unique[top20,],
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         scale = "row")
group<-data.frame(group=group) 
rownames(group)<-colnames(exp_unique) 
pheatmap(exp_unique[top20,],
         show_rownames = T,
         show_colnames = F,
         annotation_col = group,
         cluster_rows = F,
         cluster_cols = F,
         scale = "row",
         col=c("powderblue" ,"pink","deeppink4","red2"),
         main = "Top20 gene Pheatmap")
GO富集分析
rm(list=ls())
load("GSE41446.Rdata")
load("DEG.Rdata")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(AnnotationDbi)
up<-DEG[DEG$group=="up",]
down<-DEG[DEG$group=="down",]
upgene<-up$gene
downgene<-down$gene
id_down<-mapIds(org.Hs.eg.db,
                keys=downgene,
                keytype="SYMBOL",
                column="ENTREZID") 
id_down<-data.frame(id=id_down)
id_down<-na.omit(id_down)
library(clusterProfiler)
ALL_down<-enrichGO(id_down$id,
                   keyType = "ENTREZID",
                   ont="ALL",
                   org.Hs.eg.db,
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2)
all_down<-as.data.frame(ALL_down)
library(ggplot2)
dotplot(ALL_down, x = "GeneRatio",color="p.adjust",title='enrichment'
        ,showCategory=Inf,font.size=12)
dotplot(ALL_down, x = "GeneRatio",color="p.adjust",title='enrichment'
        ,showCategory=Inf,font.size=12)+ 
  scale_color_continuous(low="blue",high="red")+
  scale_size(range = c(2,12))
library(DOSE)
dotplot(ALL_down, split="ONTOLOGY")+ facet_grid(.~ONTOLOGY,scale="free")
dotplot(ALL_down, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dotplot(ALL_down,x="count",
        showCategory=Inf,
        font.size=12,
        title="Enrichment GO",
        split="ONTOLOGY")+
  scale_color_continuous(low="blue",high="red")+
  scale_size(range=c(2,8))+
  facet_grid(ONTOLOGY~.,scales = "free")
dotplot(ALL_down, split="ONTOLOGY")+ facet_grid(.~ONTOLOGY,scale="free")
dotplot(ALL_down, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
library(clusterProfiler)
x <-mutate(all_down,
           FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
FoldEnrichment
library(ggplot2)
ggplot(x,aes(x = FoldEnrichment,y = Description))+
  geom_point(aes(color = p.adjust,
                 size = Count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("Fold Enrichment")+
  theme_bw()+ylim(rev(x$Description))+ 
  guides(color = guide_colorbar(reverse = TRUE))
ev = function(x){
  eval(parse(text = x))
}
all_down$generatio = round(sapply(all_down$GeneRatio,ev),3)
ggplot(all_down,aes(x = generatio,y = Description))+
  geom_point(aes(color = p.adjust,
                 size = Count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("GeneRatio")+
  theme_bw()+ylim(rev(x$Description))+
  guides(color = guide_colorbar(reverse = TRUE))
barplot(ALL_down, showCategory=10,title="EnrichmentGO")
rm(list=ls())
gene<-read.table("gene.txt",header=F)
gene<-gene$V1
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(AnnotationDbi)
id<-mapIds(org.Hs.eg.db,
           keys = gene,
           keytype="SYMBOL",
           column="ENTREZID")
id<-data.frame(id=id)
id<-na.omit(id)
#KEGG enrichment
library(clusterProfiler)
KEGG<-enrichKEGG(id$id,
                 organism = "hsa") 
kegg<-as.data.frame(KEGG)

dotplot(KEGG,showCategory=5,title="KEGG ENrichment")
