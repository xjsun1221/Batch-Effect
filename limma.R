rm(list = ls())
load("normlize.Rdata")
#处理批次效应
library(limma)
#?removeBatchEffect()
y <- exp
batch <- c(rep("A",12),rep("B",5))
y2 <- removeBatchEffect(y, batch)
par(mfrow=c(1,2))
boxplot(as.data.frame(y),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")
exp = y2
#exp = log2(exp+1)
#PCA
{
  dat=as.data.frame(t(exp))
  library(FactoMineR)#画主成分分析图需要加载这两个包
  library(factoextra) 
  # pca的统一操作走起
  dat.pca <- PCA(dat, graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               #palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  dir.create(gse)
  ggsave(paste0(gse,"/PCA.png"))
}


#热图 
cg=names(tail(sort(apply(exp,1,sd)),1000))
if(!require(pheatmap))install.packages("pheatmap")
library(pheatmap)
n=exp[cg,]

#加注释分组

annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row")

exp[1:4,1:4]
boxplot(exp[1,]~group_list) 

#差异分析，用limma包来做
#需要表达矩阵和group_list，其他都不要动
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)

#差异基因排名
deg=topTable(fit,coef=2,number = Inf)
head(deg)
boxplot(exp[rownames(deg)[1],]~group_list)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
#tibble::rownames_to_column()
head(deg)


if(F){
  #2.加symbol列，火山图要用
  #id转换，查找芯片平台对应的包
  eSet[[1]]@annotation
  #http://www.bio-info-trainee.com/1399.html
  #hgu133a
  if(!require(hgu133a.db))BiocManager::install("hgu133a.db")
  library(hgu133a.db)
  ls("package:hgu133a.db")
  ids <- toTable(hgu133aSYMBOL)
  head(ids)
  
  #merge
  deg <- inner_join(deg,ids,by="probe_id")
  head(deg)
}

if(T){
  ids = data.table::fread("GSE83521_family.soft",header = T)[,c("ID","circRNA")]
  colnames(ids)[1] = "probe_id"
  deg <- inner_join(deg,ids,by="probe_id")
  head(deg)
}
#3.加change列：上调或下调，火山图要用

logFC_t=1 #不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
change=ifelse(deg$P.Value>0.05,'stable', 
              ifelse( deg$logFC >logFC_t,'up', 
                      ifelse( deg$logFC < -logFC_t,'down','stable') )
)
deg <- mutate(deg,change)
head(deg)
table(deg$change)
deg <- mutate(deg,v = -log10(P.Value))

head(dat)
p <- ggplot(data = deg, 
            aes(x = logFC, 
                y = v)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
#ggsave(paste0(gse,"/volcano.png"))


#x=deg$logFC 
#names(x)=deg$probe_id 
#cg=c(names(head(sort(x),30)),
#     names(tail(sort(x),30)))
cg = deg$circRNA[deg$change != "stable"]
library(pheatmap)
exp2 = exp[deg$probe_id,]
rownames(exp2) = deg$circRNA
n=exp2[cg,]

annotation_col=data.frame(group=group_list,
                          gse = c(rep("GSE83521",12),rep("GSE89143",5)))
rownames(annotation_col)=colnames(n) 

pheatmap(n,show_colnames =F,
         show_rownames = T,
         scale = "row",
         #cluster_cols = F, 
         annotation_col=annotation_col) 

