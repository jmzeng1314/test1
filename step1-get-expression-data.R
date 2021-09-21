rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  
dir.create('Rdata')

###### step1: Counts文件批量转为Rdata  ###### 

# dataset: gene expression RNAseq - HTSeq - Counts 
# unitlog2(count+1)
fs=list.files('ucsc_xena',pattern = 'htseq_counts.tsv.gz')
fs
lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('tsv.gz','',x)
  print(pro)
  a=fread(file.path('ucsc_xena',x),data.table = F)
  head(a[ ,1:4])
  tail(a[ ,1:4])
  mat=a[1:60483,] 
  rownames(mat)=mat$Ensembl_ID
  mat[1:4,1:4]
  mat=mat[,-1]
  
  exprSet=floor(2^mat - 1 )
  exprSet[1:4,1:4] 
  
  n=floor(ncol(exprSet)/20)
  n
  keep_feature <- rowSums (exprSet > 2) > n
  table(keep_feature)
  mat <- exprSet[keep_feature, ]
  mat=mat[, colSums(mat) > 1000000]
  mat[1:4,1:4]
  dim(mat) 
  colnames(mat)
  
  b=read.table('ucsc_xena/gencode.v22.annotation.gene.probeMap',header = T)
  head(b)
  d=read.table('ucsc_xena/ensembID2type.txt')
  head(d)
  b=merge(b,d,by.x='id',by.y='V1')
  head(b)
  
  length(unique(b[match(rownames(mat),b$id),2]))
  # 可以看到，多个ensembl的ID对应同一个基因的symbol，这个是需要处理掉的。
  # 下面的代码略微有一点点复杂
  dat=mat
  ids <-  b[match(rownames(dat),
                  b$id),1:2] #取需要的列 
  head(ids)
  colnames(ids)=c('probe_id','symbol')  
  ids=ids[ids$symbol != '',]
  ids=ids[ids$probe_id %in%  rownames(dat),]
  dat[1:4,1:4]   
  dat=dat[ids$probe_id,] 
  
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
  
  mat=dat
  mat[1:4,1:4]
  print(fivenum(mat['GAPDH',]))
  print(fivenum(mat['ACTB',]))
  
  tp=b[match(rownames(mat),b$gene),7]
  print(  tail(sort(table(tp))) )
  
  pd_mat=mat[tp=='protein_coding',]
  non_mat=mat[tp !='protein_coding',]
  
  save(pd_mat,non_mat,
       file = file.path('Rdata', paste0(pro,'Rdata')))
  
}) 

###### step2: 合并全部癌症表达量矩阵，自己走PCA,tSNE,DBSCAN  ###### 

fs=list.files('Rdata/',pattern = 'htseq_counts')
fs
# 首先加载每个癌症的蛋白编码基因表达量矩阵
# 成为一个 list 
dat_list = lapply(fs, function(x){
 # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro)
  
  load(file =  file.path('Rdata/',x)) 
  
  dat=log2(edgeR::cpm(pd_mat)+1)
  return(dat)
})
tmp=table(unlist(lapply(dat_list, rownames)))
cg=names(tmp)[tmp==33]
length(cg)
# 在33个癌症都存在的蛋白编码基因是一万六千个
all_tcga_matrix=do.call(cbind,
            lapply(dat_list, function(x){
              x[cg,]
            })) 
dim(all_tcga_matrix)
group=rep(gsub('TCGA-','',gsub('.htseq_counts..Rdata','',fs)),
      unlist(lapply(dat_list, ncol)))
table(group)

# 为了加快运算速度，我们取top1000基因即可，按照sd派系
cg=names(tail(sort(apply(all_tcga_matrix,1,sd)),1000))  

# 首先查看pca情况 
library("FactoMineR")
library("factoextra")  
dat.pca <- PCA(t(all_tcga_matrix[cg,]), graph = FALSE)  
fviz_pca_ind(dat.pca,
             repel =T,
             geom.ind = "point", 
             col.ind = group,
             # addEllipses = TRUE,
             ggtheme = theme_minimal(),
             legend.title = "Groups"
) 
plot(dat.pca$ind$coord[,1:2],
     col=rainbow(33)[as.numeric(as.factor(group))] )

library(Rtsne)
tsne_out <-Rtsne(t(all_tcga_matrix[cg,]),
                 initial_dims=30,verbose=FALSE,
                 check_duplicates=FALSE,
                 perplexity=27, dims=2,max_iter=5000)
plot(tsne_out$Y,
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2",
     col=rainbow(33)[as.numeric(as.factor(group))])

# 比较kmeans和dbscan的差异

library(gplots)
kmeans_cluster_tsne=kmeans(tsne_out$Y,centers = 33)$clust
library(dbscan)
dbscan_cluster_tsne=dbscan(tsne_out$Y,eps=3.1)$cluster 
balloonplot(table(dbscan_cluster_tsne,
                  kmeans_cluster_tsne))
balloonplot(table(dbscan_cluster_tsne,
                  group))
plot(tsne_out$Y,col=kmeans_cluster_tsne,
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
plot(tsne_out$Y,col=dbscan_cluster_tsne,
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
library(RColorBrewer)
plot(tsne_out$Y,col=kmeans_cluster_tsne,
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
plot(tsne_out$Y,col=dbscan_cluster_tsne,
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")


###### step3: 合并全部癌症表达量矩阵,走seurat ###### 

fs=list.files('Rdata/',pattern = 'htseq_counts')
fs
# 首先加载每个癌症的蛋白编码基因表达量矩阵
# 成为一个 list 
dat_list = lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro)
  
  load(file =  file.path('Rdata/',x)) 
  # 因为要走单细胞 seurat流程，所以无需处理表达量矩阵
  # 直接最原始的counts即可
  # dat=log2(edgeR::cpm(pd_mat)+1)
  dat=pd_mat
  return(dat)
})
tmp=table(unlist(lapply(dat_list, rownames)))
cg=names(tmp)[tmp==33]
length(cg)
# 在33个癌症都存在的蛋白编码基因是一万六千个
all_tcga_matrix=do.call(cbind,
                        lapply(dat_list, function(x){
                          x[cg,]
                        })) 
dim(all_tcga_matrix)
all_tcga_matrix[1:4,1:2]
group=rep(gsub('TCGA-','',gsub('.htseq_counts..Rdata','',fs)),
          unlist(lapply(dat_list, ncol)))
table(group)

library(Seurat)
sce=CreateSeuratObject(counts = all_tcga_matrix) 
sce
group=as.data.frame(group)
rownames(group)=colnames(sce)
sce=AddMetaData(sce,group)
table(sce$group)
# 这里并不是真正的单细胞表达量矩阵，所以无需过滤
sce.filt=sce
dim(sce.filt)   
sce.filt <- NormalizeData(sce.filt, normalization.method =  "LogNormalize", 
                              scale.factor = 10000)
GetAssay(sce.filt,assay = "RNA")
sce.filt = FindVariableFeatures(sce.filt)
sce.filt = ScaleData(sce.filt, vars.to.regress = c("nFeature_RNA" ))
sce.filt = RunPCA(sce.filt, npcs = 20)
sce.filt = RunTSNE(sce.filt, npcs = 20)
sce.filt = RunUMAP(sce.filt, dims = 1:10) 
sce=sce.filt

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.2)
table(sce@meta.data$RNA_snn_res.0.2) 

library(patchwork)
p1=DimPlot(sce,reduction = "tsne",label=T,repel = T,
           group.by  ='group')
p2=DimPlot(sce,reduction = "umap",label=T,repel = T,
           group.by ='group')
p1+p2


library(patchwork)
p1=DimPlot(sce,reduction = "tsne",label=T,repel = T )
p2=DimPlot(sce,reduction = "umap",label=T,repel = T )
p1+p2


library(gplots)
balloonplot(table(sce$seurat_clusters,
                  sce$group))


# 接下来对seurat默认的分群进行差异分析
pro='basic_seurat'
table(Idents(sce))  
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers) 
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=3)
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),width = 18,height = 21)
p <- DotPlot(sce, features = unique(top10$gene),
             assay='RNA'  )  + coord_flip()

p
ggsave(paste0(pro,'-DotPlot_check_top10_markers_by_clusters.pdf'),width = 18,height = 21)

library(dplyr) 
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(sce,top3$gene,size=3)
ggsave(paste0(pro,'-DoHeatmap_check_top3_markers_by_clusters.pdf'),width = 18,height = 11)
p <- DotPlot(sce, features = unique(top3$gene),
             assay='RNA'  )  + coord_flip()

p
ggsave(paste0(pro,'-DotPlot_check_top3_markers_by_clusters.pdf'),width = 18,height = 11)

save(sce.markers,file=paste0(pro,file = '-sce.markers.Rdata'))
save(sce,file = 'sce.Rdata')


genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'MKI67','TOP2A',
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
 
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'   )  + coord_flip()

p

# 
# ggsave(plot=p, filename="check_marker_by_celltype.pdf")
# table(sce@meta.data$celltype,sce@meta.data$seurat_clusters)
# 
# 
# DimPlot(sce, reduction = "umap", group.by = "celltype",label = T)  
# ggsave('umap_by_celltype.pdf')
# 
# library(patchwork)
# th=theme(axis.text.x = element_text(angle = 45, 
#                                     vjust = 0.5, hjust=0.5))
# p_all_markers=DotPlot(sce, features = genes_to_check,
#                       assay='RNA' ,group.by = 'celltype' )  + coord_flip()+ th
# p_umap=DimPlot(sce, reduction = "umap", group.by = "celltype",label = T) 
# p_all_markers+p_umap
# ggsave('markers_umap_by_celltype.pdf',width = 13)
# phe=sce@meta.data
# save(phe,file = 'phe_by_markers.Rdata')
# 
# 
# 

###### step4: 不同癌症的差异难道大于其与正常对照差异吗 ？ ###### 

library(Seurat)
load(file = 'sce.Rdata')
gp=substring(colnames(sce),14,15)
table(gp)
sce@meta.data$gp=gp

library(gplots)
balloonplot(table(sce@meta.data$gp,sce@meta.data$group))

library(patchwork)
p1=DimPlot(sce,reduction = "tsne",label=T,repel = T,
           group.by  ='gp')
p2=DimPlot(sce,reduction = "umap",label=T,repel = T,
           group.by ='gp')
p1+p2
cl=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17')
p1=DimPlot(sce,reduction = "tsne",label=T,repel = T,
           cols = cl,
           group.by  ='gp')
p2=DimPlot(sce,reduction = "umap",label=T,repel = T,
           cols = cl,
           group.by ='gp')
p1+p2
library(gplots)
balloonplot(table(sce$seurat_clusters,
                  sce$gp))

balloonplot(table(sce$seurat_clusters,
                  sce$group))


# 尝试使用harmony 抹去癌症差异 
# 参考；https://mp.weixin.qq.com/s/5ivvtVFwCutUpJkqRgXYow
table( sce$group)
library(harmony)
seuratObj <- RunHarmony(sce, "group")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=T ) +NoLegend()
library(cowplot)
p1.compare=plot_grid(ncol = 3,
                     DimPlot(sce, reduction = "pca", group.by = "group")+NoAxes() +NoLegend()+ggtitle("PCA"),
                     DimPlot(sce, reduction = "tsne", group.by = "group")+NoAxes() +NoLegend()+ggtitle("tSNE"),
                     DimPlot(sce, reduction = "umap", group.by = "group")+NoAxes() +NoLegend()+ggtitle("UMAP"),
                     DimPlot(seuratObj, reduction = "pca", group.by = "group")+NoAxes() +NoLegend()+ggtitle("PCA"),
                     DimPlot(seuratObj, reduction = "tsne", group.by = "group")+NoAxes() +NoLegend()+ggtitle("tSNE"),
                     DimPlot(seuratObj, reduction = "umap", group.by = "group")+NoAxes() +NoLegend()+ggtitle("UMAP")
)
p1.compare
ggsave(plot=p1.compare,filename="before_and_after_inter_cancerType.pdf")


sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
table(sce@meta.data$seurat_clusters)
DimPlot(sce,reduction = "umap",label=T) 
ggsave(filename = 'harmony_umap_recluster_by_0.5.pdf') 

p1=DimPlot(sce,reduction = "umap",label=T,
        group.by = 'group') 
p1
ggsave(filename = 'harmony_umap_sce_recluster_by_cancerType.pdf')

cl=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17')
p2=DimPlot(sce,reduction = "umap",label=T,
        cols = cl,
        group.by = 'gp') 
p2
ggsave(filename = 'harmony_umap_sce_recluster_by_position.pdf')
library(patchwork)
p1+p2

# 接下来对seurat默认的分群进行差异分析

pro='harmony_seurat-0.5'
table(Idents(sce))  
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers) 
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=3)
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),width = 18,height = 21)
p <- DotPlot(sce, features = unique(top10$gene),
             assay='RNA'  )  + coord_flip()

p
ggsave(paste0(pro,'-DotPlot_check_top10_markers_by_clusters.pdf'),width = 18,height = 21)

library(dplyr) 
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(sce,top3$gene,size=3)
ggsave(paste0(pro,'-DoHeatmap_check_top3_markers_by_clusters.pdf'),width = 18,height = 11)
p <- DotPlot(sce, features = unique(top3$gene),
             assay='RNA'  )  + coord_flip()

p
ggsave(paste0(pro,'-DotPlot_check_top3_markers_by_clusters.pdf'),width = 18,height = 11)

save(sce.markers,file=paste0(pro,file = '-sce.markers.Rdata'))


load(file='seurat_results/harmony_seurat-0.5-sce.markers.Rdata')
table(sce.markers$cluster)
library(org.Hs.eg.db)
library(clusterProfiler)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)
gcSample # entrez id , compareCluster 
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
p=dotplot(xx) 
p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
p
ggsave('seurat_results/harmony_seurat-0.5-compareCluster-enrichKEGG.pdf',width = 15)

library(gplots)
balloonplot(table(sce$seurat_clusters,
                  sce$gp))

balloonplot(table(sce$seurat_clusters,
                  sce$group))

# 可以看到 0，3，6，10是核糖体富集
# 

###### step5:   ###### 








