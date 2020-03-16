rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
gse = "GSE83521"
eSet1 <- getGEO("GSE83521", 
                destdir = '.', 
                getGPL = F)
eSet2 <- getGEO("GSE89143", 
                destdir = '.', 
                getGPL = F)
#(1)提取表达矩阵exp
exp1 <- exprs(eSet1[[1]])
exp1[1:4,1:4]
exp2 <- exprs(eSet2[[1]])
exp2[1:4,1:4]
exp2 = log2(exp2+1)
table(rownames(exp1) %in% rownames(exp2))
length(intersect(rownames(exp1),rownames(exp2)))
exp1 <- exp1[intersect(rownames(exp1),rownames(exp2)),]

exp2 <- exp2[intersect(rownames(exp1),rownames(exp2)),]
if(F){
  exp1 = limma::normalizeBetweenArrays(exp1)
  exp2 = limma::normalizeBetweenArrays(exp2)
}else{
  exp2 = exp2[,-3]
}


exp = cbind(exp1,exp2)
boxplot(exp)

#exp = log2(exp+1)
#(2)提取临床信息
pd1 <- pData(eSet1[[1]])
pd2 <- pData(eSet2[[1]])

#(3)提取芯片平台编号
gpl <- eSet2[[1]]@annotation

if(!identical(rownames(pd1),colnames(exp1))) exp1 = exp1[,match(rownames(pd1),colnames(exp1))]
if(!identical(rownames(pd2),colnames(exp2))) exp2 = exp2[,match(rownames(pd2),colnames(exp2))]

group_list1 = ifelse(str_detect(pd1$title,"Tumour"),"Tumour","Normal")
group_list2 = ifelse(str_detect(pd2$source_name_ch1,"Paracancerous"),"Normal","Tumour")[-3]

group_list = c(group_list1,group_list2)
table(group_list)
group_list = factor(group_list,levels = c("Normal","Tumour"))
save(gse,group_list,exp,gpl,file = "normlize.Rdata")

