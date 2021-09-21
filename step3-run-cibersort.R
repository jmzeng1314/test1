rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  
library(preprocessCore)
library(parallel)
library(e1071) 
library(dplyr)
library(tidyr)
library(tidyverse)
# 表达矩阵文件:行名变成一列
# LM22.txt: 22种免疫细胞的基因表达特征数据
dir.create('CIBERSORT_results') 

###### step1: 一个个癌症内部运行 CIBERSORT  ###### 
fs=list.files('Rdata/',pattern = 'htseq_counts')
fs
lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro)
  
  load(file =  file.path('Rdata/',x)) 
  
  dat=log2(edgeR::cpm(pd_mat)+1)
  # 只能说假装是 tpm，去搞基因长度信息很麻烦
  exp_tpm=dat
  exp_tpm[1:4,1:4]  
  source("CIBERSORT.R")  
  load('lm22.rda')
  ciber=CIBERSORT(dat) 
  cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
    rename("Patients" = "Mixture") %>%
    select(-c("P.value","Correlation","RMSE"))
  # 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Patiens”。
  #并赋值给cibersort_raw。 
  
  save(cibersort_raw, file = file.path('CIBERSORT_results',
                                     paste0('cibersort_raw-for-',pro,'.Rdata'))) 
  
})



###### step2: 批量可视化（没有分组的），细胞比例的热图，柱状图，箱线图，相关性图  ######   

fs=list.files('CIBERSORT_results/',
              pattern = 'cibersort_raw')
fs

lapply(fs, function(x){ 
  # x=fs[1]
  pro=gsub('.Rdata','',
           gsub('cibersort_raw-for-TCGA-','',x))
  print(pro)
  
  load(file =  file.path('CIBERSORT_results/',x)) 
  
  cibersort_tidy <- cibersort_raw %>%
    remove_rownames() %>%
    column_to_rownames("Patients")
  
  # 将cibersort_raw第一列变为列名后赋值给cibersort_tidy。
  
  flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) <
                  dim(cibersort_tidy)[1]/2)
  # 筛选出0值太多的一些细胞。
  
  cibersort_tidy <- cibersort_tidy[,which(flag)] %>%
    as.matrix() %>%
    t()
  # 留下在大部分样本中有所表达的细胞。 
  bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
  # breaks用来定义数值和颜色的对应关系。
  
  # Step4:将CIBERSORT_Result进行可视化
  
  #1）热图
  library(pheatmap)
  library(RColorBrewer)
  pheatmap(
    cibersort_tidy,
    breaks = bk,
    cluster_cols = T,
    scale = "row",
    cluster_row = T,
    border_color = NA,
    show_colnames = F,
    show_rownames = T,
    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
              colorRampPalette(colors = c("white","red"))(length(bk)/2)
    ),
    filename = file.path('CIBERSORT_results',
                         paste0('pheatmap-for-',pro,'.pdf'))
      )
  
  #调整参数让热图更加美观。
  
  #柱状图可视化细胞占比预测
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  cibersort_barplot <- cibersort_raw %>%
    gather(key = Cell_type,value = Proportion,2:23)
  #使用RColorBrewer包配置需要的色彩方案，使用gather函数中的key-value对应关系重建细胞名称和比例的对应关系并赋值给cibersort_barplot
  
  #cibersort_barplot$Patient1 <- factor(cibersort_barplot$Patient,
  #                                   levels = str_sort(unique(cibersort_barplot$Patient),
  #                                                      numeric = T))
  
  ggplot(cibersort_barplot,aes(Patients,Proportion,fill = Cell_type)) +
    geom_bar(position = "stack",stat = "identity") +
    labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + theme_bw() +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_y_continuous(expand = c(0.01,0)) +
    scale_fill_manual(values = mypalette(23))
  ggsave(filename  = file.path('CIBERSORT_results',
                              paste0('barplot-for-',pro,'.pdf')))
  #调整参数让柱状图更加美观。
  
  #直观箱线图
  ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) +
    geom_boxplot(outlier.shape = 21,color = "black") + theme_bw() +
    labs(x = "", y = "Estimated Proportion") +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = mypalette(23)) 
  ggsave(filename  = file.path('CIBERSORT_results',
                               paste0('boxplot-for-',pro,'.pdf')))
  
  library(corrplot)
  M = cor( t(cibersort_tidy)  )
  p.mat <- cor.mtest( t(cibersort_tidy)   )$p
  library(paletteer)
  my_color = rev(paletteer_d("RColorBrewer::RdYlBu")[-1])
  my_color = colorRampPalette(my_color)(10)
  pdf(file.path('CIBERSORT_results',
                                paste0('corrplot-for-',pro,'.pdf')) 
      )
  corrplot(M, type="upper", 
           order="hclust", 
           col = my_color,
           p.mat = p.mat, 
           sig.level = 0.01, 
           insig = "blank",
           tl.col = "black",
           tl.srt=45)
  dev.off()
  
   
  save(cibersort_tidy, file = file.path('CIBERSORT_results',
                                       paste0('cibersort_tidy-for-',pro,'.Rdata'))) 
  
  
})




###### step3: 加上分组的可视化  ######   
# todo




###### step4: 直接对seurat对象走 CIBERSORT ######   
# todo

rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  
library(preprocessCore)
library(parallel)
library(e1071) 
library(dplyr)
library(tidyr)
library(tidyverse)
library(Seurat)
load(file = 'sce.Rdata')
gp=substring(colnames(sce),14,15)
table(gp)
sce@meta.data$gp=gp
pd_mat=sce@assays$RNA@counts
dat=log2(edgeR::cpm(pd_mat)+1) 
pro='CIBERSORT_for_seurat' 
# 只能说假装是 tpm，去搞基因长度信息很麻烦
exp_tpm=dat
exp_tpm[1:4,1:4]  
source("CIBERSORT.R")  
load('lm22.rda')
ciber=CIBERSORT(dat) 

cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
  rename("Patients" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE"))

dim(cibersort_raw)
sce
cibersort_raw$group= sce@meta.data$gp  
cibersort_raw$type= sce@meta.data$group  

pro='seurat'
save(cibersort_raw, 
     file = file.path('CIBERSORT_results',
    paste0('cibersort_raw-for-',pro,'.Rdata'))) 
 

pro='seurat'
load(  file = file.path('CIBERSORT_results',
                        paste0('cibersort_raw-for-',pro,'.Rdata')))
library(reshape2) 
cgDat=as.data.frame(reshape2::melt( cibersort_raw,
                                    id=c('Patients','group','type') ))
head(cgDat) 

library(ggpubr)
library(ggstatsplot)
ggboxplot(cgDat, "variable", "value", color = "type" ) +
  ylab('expression value  ')+  
  stat_compare_means(aes(group =  type ),  
                     label = "p.signif",hide.ns = T) +
  theme_ggstatsplot()+
  theme(legend.position='none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
 

## 下面是画PCA的必须操作，需要看说明书。
cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
  rename("Patients" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE"))
cibersort_tidy <- cibersort_raw %>%
  remove_rownames() %>%
  column_to_rownames("Patients")

dat=as.data.frame(cibersort_tidy) 
dim(dat)
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
# The variable group_list (index = 54676) is removed
# before PCA analysis
dat.pca <- PCA(dat , 
               graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             #col.ind =  sce$gp, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('cibersort_all_tcga_PCA_by_position.pdf',
       height = 4,width = 8)  
 


###### step5: 对比不同 CIBERSORT ######   
# todo

cibersort_all <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
  rename("Patients" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE")) 

# 首先区分癌症各自运行时候有免疫细胞比例推断

fs=list.files('CIBERSORT_results/',
              pattern = 'cibersort_raw-for-TCGA')
fs

lapply(fs, function(x){ 
  # x=fs[1]
  pro=gsub('.Rdata','',
           gsub('cibersort_raw-for-TCGA-','',x))
  print(pro)
  
  load(file =  file.path('CIBERSORT_results/',x)) 
  pos=match(cibersort_raw$Patients,cibersort_all$Patients)
  df=cbind(cibersort_raw[,-1],
           cibersort_all[pos,-1])
  kp=apply(df, 2,sd) > 0
  df=df[,kp]
  pheatmap::pheatmap(  cor(df),
                       filename = file.path('CIBERSORT_results',
                                            paste0('cor-by-merge-and-single-',pro,'.pdf')))
   
})

# 其次，pan-cancer官网自带一个免疫细胞比例

cibersort_all <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
  rename("Patients" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE")) 
cibersort_all[1:4,1:4]
library(data.table)
TCGA.Kallisto.fullIDs.cibersort = fread('TCGA.Kallisto.fullIDs.cibersort.relative.tsv',data.table = F)
TCGA.Kallisto.fullIDs.cibersort[1:4,1:4]

identical(colnames(cibersort_all)[2:22],
          colnames(TCGA.Kallisto.fullIDs.cibersort)[3:23])
pid=gsub('[.]','-',substring(TCGA.Kallisto.fullIDs.cibersort$SampleID,1,16))
head(pid)
table(cibersort_all$Patients %in% pid)
pos=match(cibersort_all$Patients,pid) 
TCGA.Kallisto.fullIDs.cibersort=TCGA.Kallisto.fullIDs.cibersort[pos,]
pheatmap::pheatmap(cor(
  cbind(cibersort_all[2:22],
        TCGA.Kallisto.fullIDs.cibersort[,3:23])
))



