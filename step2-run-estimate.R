rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  

if(!require("estimate")){
  library(utils)
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
library(estimate)
library(stringr)

# 参考：http://www.bio-info-trainee.com/6602.html 

###### step1: 设置 estimate 的包和函数  ###### 
estimate_RNAseq <- function(RNAseq_logCPM,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(RNAseq_logCPM,file = input.f,sep = '\t',quote = F)
  
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  scores=data.frame(scores)
  return(scores)
} 

dir.create('estimate_results')


###### step2: 一个个癌症内部运行estimate  ###### 
fs=list.files('Rdata/',pattern = 'htseq_counts')
fs 

# 首先，针对全部的TCGA数据库癌症的表达量矩阵批量运行  estimateScore {estimate} 函数
# StromalScorenumeric scalar specifying the presence of stromal cells in tumor tissue
# ImmuneScorenumeric scalar specifying the level of infiltrating immune cells in tumor tissue
# ESTIMATEScorenumeric scalar specifying tumor cellularity
# TumorPuritynumeric scalar specifying ESTIMATE-based tumor purity with value in range[0,1]

lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro)
  
  load(file =  file.path('Rdata/',x)) 
  
  dat=log2(edgeR::cpm(pd_mat)+1)
  scores=estimate_RNAseq(dat,pro)
  
  head(scores)
  scores$group = ifelse(
    as.numeric(substring(rownames(scores),14,15)) < 10,
    'tumor','normal'
  )
  table(scores$group ) 
  
  save(scores,file = file.path('estimate_results',
                               paste0('estimate_RNAseq-score-for-',pro,'.Rdata'))) 
  
  })
 
# 然后批量载入各自 estimate 结果 进行绘图
fs=list.files('estimate_results/',
              pattern = 'estimate_RNAseq-score-forTCGA')
fs
#lapply(head(fs), function(x){ 
lapply( fs , function(x){ 
  # x=fs[1]
  pro=gsub('.Rdata','',
           gsub('estimate_RNAseq-score-forTCGA-','',x))
  print(pro)
  
  load(file =  file.path('estimate_results/',x))  
  
  #可视化
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(ggsci)
  library(ggstatsplot) 
  b1 = ggplot(dat = scores,aes(group,StromalScore))+
    geom_boxplot(aes(fill = group))+ 
    stat_compare_means()+ ggtitle('StromalScore') +
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  b2 = ggplot(dat = scores,aes(group,ImmuneScore ))+
    geom_boxplot(aes(fill = group))+ 
    stat_compare_means()+ggtitle('ImmuneScore') +
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  
  b3 = ggplot(dat = scores,aes(group,ESTIMATEScore ))+
    geom_boxplot(aes(fill = group))+ 
    stat_compare_means()+ggtitle('ESTIMATEScore') +
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  
  b1+b2+b3
  ggsave(file.path('estimate_results',
                   paste0('estimate_RNAseq-score-boxplot-for-',pro,'.pdf'))
  )
  
  
})


###### step3: 接下来提取相应的基因列表绘制热图 ######  
gmtFile=system.file("extdata", "SI_geneset.gmt",
                    package="estimate")
gmtFile
library(GSVA)
library(limma)
library(GSEABase)
library(data.table)
geneSet=getGmt(gmtFile,
               geneIdType=SymbolIdentifier())

geneSet
StromalGenes=geneSet[[1]]@geneIds
ImmuneGenes=geneSet[[2]]@geneIds 

fs=list.files('Rdata/',pattern = 'htseq_counts')
fs

lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro) 
  load(file =  file.path('Rdata/',x))  
  dat=log2(edgeR::cpm(pd_mat)+1)
  StromalGenes=StromalGenes[StromalGenes %in% rownames(dat)]
  ImmuneGenes=ImmuneGenes[ImmuneGenes %in% rownames(dat)]
  library(pheatmap)
  group = ifelse(
    as.numeric(substring(colnames(dat),14,15)) < 10,
    'tumor','normal'
  )
  ac=data.frame(group)
  rownames(ac)=colnames(dat)
  
  
  pheatmap(dat[StromalGenes,],show_colnames = F,
           annotation_col = ac,width =  ncol(dat)/10 ,height = 15,
           filename =file.path('estimate_results',
                               paste0('StromalGenes-pheatmap-for-',pro,'.pdf')) )  
  pheatmap(dat[ImmuneGenes,],show_colnames = F,
           annotation_col = ac,
           filename =file.path('estimate_results',
                               paste0('ImmuneGenes-pheatmap-for-',pro,'.pdf')) )  
  n=rbind(dat[ImmuneGenes,],
          dat[StromalGenes,])
  ar=data.frame(geneGroup=c(
    rep('ImmuneGenes',length(ImmuneGenes)),
    rep('StromalGenes',length(StromalGenes))
  ))
  rownames(ar)=rownames(n)
  pheatmap(n,show_colnames = F,show_rownames = F,
           annotation_col = ac,annotation_row = ar,
           filename =file.path('estimate_results',
                               paste0('combine-pheatmap-for-',pro,'.pdf')) ) 
})




###### step4: 针对前面的seurat对象运行 estimate  ######   

library(Seurat)
load(file = 'sce.Rdata')
gp=substring(colnames(sce),14,15)
table(gp)
sce@meta.data$gp=gp
pd_mat=sce@assays$RNA@counts
dat=log2(edgeR::cpm(pd_mat)+1) 
pro='estimate_for_seurat'

scores=estimate_RNAseq(dat,pro)
head(scores) 
scores$group= sce@meta.data$gp  
scores$type= sce@meta.data$group  

save(scores,file = file.path('estimate_results',
                             paste0('estimate_RNAseq-score-for-',
                                    pro,'.Rdata'))) 

#  载入针对seurat对象的 estimate 结果 
pro='estimate_for_seurat'
load(file = file.path('estimate_results',
                      paste0('estimate_RNAseq-score-for-',
                             pro,'.Rdata')))
head(scores)
#可视化
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsci)
library(ggstatsplot)

# figure 1 
{
  
  cl=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17')
  b1 = ggplot(dat = scores,aes(group,StromalScore))+
    geom_boxplot(aes(fill = group))+
    scale_fill_manual(values = cl)+
    stat_compare_means()+ ggtitle('StromalScore') +
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  b1
  b2 = ggplot(dat = scores,aes(group,ImmuneScore ))+
    geom_boxplot(aes(fill = group))+
    scale_fill_manual(values = cl)+
    stat_compare_means()+ggtitle('ImmuneScore') +
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15)) 
  b3 = ggplot(dat = scores,aes(group,ESTIMATEScore ))+
    geom_boxplot(aes(fill = group))+
    scale_fill_manual(values = cl)+
    stat_compare_means()+ggtitle('ESTIMATEScore') +
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  
  b1/b2/b3
  ggsave(file.path('estimate_results',
                   paste0('estimate_RNAseq-score-boxplot-for-',pro,'.pdf')),
         height = 15,
  )
  
  
}

# figure 2 
cg_scores=scores[as.numeric(scores$group )< 11,]
{
  
  b1 = ggplot(dat = cg_scores,aes(type,StromalScore))+
    geom_boxplot(aes(fill = type))+ 
    stat_compare_means()+ ggtitle('StromalScore') +
    theme_ggstatsplot()+
    theme(legend.position='none') +
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,
                                 hjust = 1,angle = 90,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  b1
  b2 = ggplot(dat = cg_scores,aes(type,ImmuneScore ))+
    geom_boxplot(aes(fill = type))+ 
    stat_compare_means()+ggtitle('ImmuneScore') +
    theme_ggstatsplot()+
    theme(legend.position='none') +
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,
                                 hjust = 1,angle = 90,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  b3 = ggplot(dat = cg_scores,aes(type,ESTIMATEScore ))+
    geom_boxplot(aes(fill = type))+ 
    stat_compare_means()+ggtitle('ESTIMATEScore') +
    theme_ggstatsplot()+
    theme(legend.position='none') +
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,
                                 hjust = 1,angle = 90,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  
  b1/b2/b3
  
  ggsave(file.path('estimate_results',
                   paste0('cancerType-in-tumors-',pro,'.pdf') ),
         height = 15
  )
}


# figure 3
cg_scores=scores[scores$group=='11',]
{
  
  b1 = ggplot(dat = cg_scores,aes(type,StromalScore))+
    geom_boxplot(aes(fill = type))+ 
    stat_compare_means()+ ggtitle('StromalScore') +
    theme_ggstatsplot()+
    theme(legend.position='none') +
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,
                                 hjust = 1,angle = 90,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  b1
  b2 = ggplot(dat = cg_scores,aes(type,ImmuneScore ))+
    geom_boxplot(aes(fill = type))+ 
    stat_compare_means()+ggtitle('ImmuneScore') +
    theme_ggstatsplot()+
    theme(legend.position='none') +
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,
                                 hjust = 1,angle = 90,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  b3 = ggplot(dat = cg_scores,aes(type,ESTIMATEScore ))+
    geom_boxplot(aes(fill = type))+ 
    stat_compare_means()+ggtitle('ESTIMATEScore') +
    theme_ggstatsplot()+
    theme(legend.position='none') +
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,
                                 hjust = 1,angle = 90,face = 'bold'),
          plot.title = element_text(hjust = 0.5,size =15))
  
  b1/b2/b3
  
  ggsave(file.path('estimate_results',
                   paste0('cancerType-in-normal-',pro,'.pdf')),
         height = 15
  )
}
 

###### step5: 对比不同 estimate 结果 ######   
# todo

rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  
# 首先载入针对seurat对象的 estimate 结果 
pro='estimate_for_seurat'
load(file = file.path('estimate_results',
                      paste0('estimate_RNAseq-score-for-',
                             pro,'.Rdata')))
head(scores)
all_scores = scores

# 然后批量载入各自 estimate 结果 
fs=list.files('estimate_results/',
              pattern = 'estimate_RNAseq-score-forTCGA')
fs

lapply(fs, function(x){ 
  # x=fs[1]
  pro=gsub('.Rdata','',
           gsub('estimate_RNAseq-score-forTCGA-','',x))
  print(pro)
  
  load(file =  file.path('estimate_results/',x))  
  df=cbind(scores[,1:3],
           all_scores[rownames(scores) ,1:3])
  kp=apply(df, 2,sd) > 0
  df=df[,kp]
  m=cor(df)
  # pheatmap::pheatmap(  cor(df))
  pheatmap::pheatmap(  cor(df),
                       filename = file.path('estimate_results',
                                            paste0('cor-by-merge-and-single-',pro,'.pdf')))
  
})


 
###### step6: 使用ssGSEA代替estimate ######   
# todo
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table) 

library(Seurat)
load(file = 'sce.Rdata')
gp=substring(colnames(sce),14,15)
table(gp)
sce@meta.data$gp=gp
pd_mat=sce@assays$RNA@counts
dat=log2(edgeR::cpm(pd_mat)+1) 

pro='ssGSEA_for_seurat'
gmtFile=system.file("extdata", "SI_geneset.gmt",
                    package="estimate")
gmtFile
library(GSVA)
library(limma)
library(GSEABase)
library(data.table)
geneSet=getGmt(gmtFile,
               geneIdType=SymbolIdentifier())

geneSet
StromalGenes=geneSet[[1]]@geneIds
ImmuneGenes=geneSet[[2]]@geneIds 
# 参考：http://www.bio-info-trainee.com/6602.html
geneSet
as.matrix(dat)[1:4,1:4]
ssgseaScore=gsva(as.matrix(dat), geneSet, 
                 method='ssgsea', kcdf='Gaussian', 
                 abs.ranking=TRUE)
ssgseaScore=as.data.frame(t(ssgseaScore))
head(ssgseaScore)
save(ssgseaScore,file = file.path('estimate_results',
                             paste0('ssgseaScore-for-',
                                    pro,'.Rdata'))) 

pro='ssGSEA_for_seurat'
load(file = file.path('estimate_results',
                      paste0('ssgseaScore-for-',
                             pro,'.Rdata')))
head(ssgseaScore)
pro='estimate_for_seurat'
load(file = file.path('estimate_results',
                      paste0('estimate_RNAseq-score-for-',
                             pro,'.Rdata')))
head(scores)
rownames(ssgseaScore)=gsub('-','.', rownames(ssgseaScore))
head(ssgseaScore)

identical(rownames(scores),
          rownames(ssgseaScore) ) 
cor(scores$StromalScore,ssgseaScore$StromalSignature)
cor(scores$ImmuneScore,ssgseaScore$ImmuneSignature)

plot(scores$StromalScore,
     ssgseaScore$StromalSignature)

df=cbind(ssgseaScore,scores)
head(df)
library(ggpubr)
library(papatchwork)
p1=ggscatter(df, x = "StromalSignature", y = "StromalScore",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson",  label.sep = "\n")
)
p1
p2=ggscatter(df, x = "ImmuneSignature", y = "ImmuneScore",
             color = "black", shape = 21, size = 3, # Points color, shape and size
             add = "reg.line",  # Add regressin line
             add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE, # Add confidence interval
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             cor.coeff.args = list(method = "pearson",  label.sep = "\n")
)
p2
p1+p2




###### step7: 不同癌症内部按照estimate的两个打分值高低分组看蛋白编码基因表达量差异  ######   

rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  
library(Seurat)
load(file = 'sce.Rdata')
gp=substring(colnames(sce),14,15)
table(gp)
sce@meta.data$gp=gp
pd_mat=sce@assays$RNA@counts
pd_mat[1:4,1:4]
dat=log2(edgeR::cpm(pd_mat)+1)  
# 首先载入针对seurat对象的 estimate 结果 
pro='estimate_for_seurat'
load(file = file.path('estimate_results',
                      paste0('estimate_RNAseq-score-for-',
                             pro,'.Rdata'))) 
cg_scores=scores[as.numeric(scores$group )< 11,]
head(cg_scores)


suppressMessages(library(limma))
suppressMessages(library(edgeR))
tp=unique(cg_scores$type) 
stromal_DEG_limma_voom_list <- lapply(tp, function(i){
  # i=tp[1]
  tp_scores=cg_scores[cg_scores$type==i,]
  group_list=ifelse(tp_scores$StromalScore > median(tp_scores$StromalScore ),'high','low')
  exprSet=pd_mat[,match(gsub('[.]','-',rownames(tp_scores)),
                      colnames(pd_mat)) ]
  
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design 
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design) 
  con='high-low'
  cont.matrix=makeContrasts(contrasts=c(con),
                            levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef=con, n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
  return(DEG_limma_voom)
})

Immune_DEG_limma_voom_list <- lapply(tp, function(i){
  # i=tp[1]
  tp_scores=cg_scores[cg_scores$type==i,]
  group_list=ifelse(tp_scores$ImmuneScore > median(tp_scores$ImmuneScore ),'high','low')
  exprSet=pd_mat[,match(gsub('[.]','-',rownames(tp_scores)),
                        colnames(pd_mat)) ]
  
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design 
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design) 
  con='high-low'
  cont.matrix=makeContrasts(contrasts=c(con),
                            levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef=con, n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
  return(DEG_limma_voom)
})
names(stromal_DEG_limma_voom_list)=tp
names(Immune_DEG_limma_voom_list)=tp
save(stromal_DEG_limma_voom_list,Immune_DEG_limma_voom_list,
     file='DEG_limma_voom_list_by_estimate_in_all_tcga_cancers.Rdata')


load(   file='DEG_limma_voom_list_by_estimate_in_all_tcga_cancers.Rdata' )
head(stromal_DEG_limma_voom_list[[1]])

# 批量可视化上下调基因数量：
# 这里仅仅是展现 stromal_DEG_limma_voom_list 的操作

logFC_cut = 1
padj_cut = 0.05
 
stat = lapply( stromal_DEG_limma_voom_list , function(DEG){ 
  DEG_sig = subset(DEG, abs(logFC)> logFC_cut & 
                     adj.P.Val< padj_cut)
  tb=table(DEG_sig$logFC>0)
  up_total = as.numeric(tb["TRUE"])
  down_total = as.numeric(tb["FALSE"])
  return(c(up_total, down_total))
})
stat = as.data.frame(do.call(rbind,stat))
colnames(stat) = c("UP","DOWN")
rownames(stat) =  names(stromal_DEG_limma_voom_list)
stat$Type = rownames(stat)
stat
# stat_reshaped=reshape2::melt(stat)
# head(stat_reshaped)
stat_up = subset(stat, select = -DOWN)
stat_up$group = "UP" 
stat_down = subset(stat, select = -UP)
stat_down$label = stat_down$DOWN
stat_down$DOWN = -1*stat_down$DOWN
stat_down$group = "DOWN" 
p_multi_DEG=ggplot() + 
  geom_bar(data = stat_up, aes(x=Type, y=UP, fill=group),stat = "identity",position = 'dodge') +
  geom_text(data = stat_up, aes(x=Type,  y=UP, label=UP, vjust=-0.25)) +
  geom_bar(data = stat_down, aes(x=Type, y=DOWN, fill=group),stat = "identity",position = 'dodge')+
  geom_text(data = stat_down, aes(x=Type,  y=DOWN, label=label, vjust= 1.25)) +
  scale_fill_manual(values=c("#0072B5","#BC3C28")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8, size = 12)) + 
  xlab(label = "") + ylab(label = "DEG numbers") + 
  ggtitle("logFC cut-off : 1\nPadj cut-off : 0.05")
p_multi_DEG
ggsave(filename = "multi_DEG_stromal_DEG_limma_voom_list.pdf", 
       plot = p_multi_DEG, width = 8, height = 10)


# 基因的上下调排序，进行批量ssGSEA分析
# 这里仅仅是展现 stromal_DEG_limma_voom_list 的操作
{
  
  lapply(stromal_DEG_limma_voom_list,dim)
  allgenes=rownames(stromal_DEG_limma_voom_list[[1]])
  stromal_logFC_df = as.data.frame( do.call(cbind,
                                            lapply(stromal_DEG_limma_voom_list,function(deg){
                                              deg[allgenes,1]
                                            })))
  rownames(stromal_stromal_logFC_df)=allgenes
  colnames(stromal_logFC_df)=names(stromal_DEG_limma_voom_list)
  stromal_logFC_df[1:5,1:5]
  # install.packages("msigdbr")
  library(msigdbr)
  h_gene_sets = msigdbr(species = "human", category = "H")
  
  msigdbr_list = split(x = h_gene_sets$gene_symbol,
                       f = h_gene_sets$gs_name)
  msigdbr_list 
  
  library(GSVA)
  library(limma)
  library(GSEABase)
  library(data.table) 
  msigdbr_list
  # 参考：http://www.bio-info-trainee.com/6602.html
  gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
    GeneSet(unique(geneIds), geneIdType=EntrezIdentifier(),
            collectionType=KEGGCollection(keggId),
            setName=keggId)
  }, msigdbr_list, names(msigdbr_list)
  ))
  gsc 
  stromal_ssgseaScore=gsva(as.matrix(stromal_logFC_df), gsc, 
                           method='ssgsea', kcdf='Gaussian', 
                           abs.ranking=TRUE)
  rownames(stromal_ssgseaScore)=gsub('HALLMARK_','',rownames(stromal_ssgseaScore))
  pheatmap(stromal_ssgseaScore) 
  
  
}


{
  
  lapply(Immune_DEG_limma_voom_list,dim)
  allgenes=rownames(Immune_DEG_limma_voom_list[[1]])
  Immune_logFC_df = as.data.frame( do.call(cbind,
                                           lapply(Immune_DEG_limma_voom_list,function(deg){
                                             deg[allgenes,1]
                                           })))
  rownames(Immune_logFC_df)=allgenes
  colnames(Immune_logFC_df)=names(Immune_DEG_limma_voom_list)
  Immune_logFC_df[1:5,1:5]
  # install.packages("msigdbr")
  library(msigdbr)
  h_gene_sets = msigdbr(species = "human", category = "H")
  
  msigdbr_list = split(x = h_gene_sets$gene_symbol,
                       f = h_gene_sets$gs_name)
  msigdbr_list 
  
  library(GSVA)
  library(limma)
  library(GSEABase)
  library(data.table) 
  msigdbr_list
  # 参考：http://www.bio-info-trainee.com/6602.html
  gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
    GeneSet(unique(geneIds), geneIdType=EntrezIdentifier(),
            collectionType=KEGGCollection(keggId),
            setName=keggId)
  }, msigdbr_list, names(msigdbr_list)
  ))
  gsc 
  Immune_ssgseaScore=gsva(as.matrix(Immune_logFC_df), gsc, 
                          method='ssgsea', kcdf='Gaussian', 
                          abs.ranking=TRUE)
  rownames(Immune_ssgseaScore)=gsub('HALLMARK_','',rownames(Immune_ssgseaScore))
  pheatmap(Immune_ssgseaScore) 
  
  
}




