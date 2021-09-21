rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  

###### step1: 批量读取前面的Rdata  ###### 

fs=list.files('Rdata/',
              pattern = 'htseq_counts')
fs
# 首先加载每个癌症的表达量矩阵（提取非蛋白编码基因部分）
# 成为一个 list 
dat_list = lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro)
  
  load(file =  file.path('Rdata/',x)) 
  
  dat=log2(edgeR::cpm(non_mat)+1)
  dat[1:4,1:4]
  return(dat)
})

###### step2: 区分基因  ######

tmp=table(unlist(lapply(dat_list, rownames)))
all_genes=names(tmp) #[tmp==33]
length(all_genes)
ubiquitous_genes=names(tmp)[tmp==33]
length(ubiquitous_genes)
specific_genes =names(tmp)[tmp==1]
length(specific_genes)
intermediately_genes=setdiff(setdiff(all_genes,specific_genes),ubiquitous_genes)

# 在33个癌症都存在的非编码基因不到九千个
# 但是33个癌症总共涉及到31455个非编码基因
# 独特存在于33种癌症的仅仅是一个里面的是 2236个基因。


all_tcga_matrix=do.call(cbind,
                        lapply(dat_list, function(x){
                          y=matrix(0,length(all_genes),ncol(x))
                          rownames(y)=all_genes
                          y[match(rownames(x),rownames(y)),] = x
                          y
                        })) 

dim(all_tcga_matrix)
# 上面的代码难度有点高，但实际上没有必要

###### step3: 计算各个基因在各个癌症的均值  ######

# 但是突然间发现，下面的代码难度更高了：
sm = do.call(cbind,
        lapply(dat_list, function(x){
          m=rowMeans(x)
          ubiquitous_m=rep(0,length(ubiquitous_genes))
          specific_m=rep(0,length(specific_genes))
          intermediately_m=rep(0,length(intermediately_genes))
          
          cg=intersect(names(m),ubiquitous_genes)
          ubiquitous_m[match(cg,ubiquitous_genes)]=m[match(cg,names(m))]
          cg=intersect(names(m),specific_genes)
          specific_m[match(cg,specific_genes)]=m[match(cg,names(m))]
          cg=intersect(names(m),intermediately_genes)
          intermediately_m[match(cg,intermediately_genes)]=m[match(cg,names(m))]
           
          return(c(ubiquitous_m,intermediately_m,specific_m))
          })) 
colnames(sm)=gsub('TCGA-','',gsub('.htseq_counts..Rdata','',fs))
rownames(sm)=c(ubiquitous_genes,
               intermediately_genes,
               specific_genes) 

library(pheatmap)
pheatmap(sm)


###### step4: 美化  ###### 


sm_bak=sm
sm=sm_bak
sm[sm>6]=6

df1=sm[ubiquitous_genes,]
df1=df1[order(rowMeans(df1),decreasing = T),]
p1=pheatmap(t(df1), 
            cluster_rows = F,
            cluster_cols = F,
            show_rownames = T,
            show_colnames = F)
p1


df2=sm[intermediately_genes,]
n_genes = apply(df2, 1, function(x) sum(x>0)) 
tmp_list= split(names(n_genes),n_genes)
lapply(tmp_list, length)
ordd=do.call(c,
        lapply(tmp_list, function(cg){
         # cg=tmp_list[[1]]
          ord=do.call(c,lapply(dat_list, function(x){
            cg[cg %in% rownames(x)]
          }))
          length(ord)
          return(unique(ord))
        }) )
length(ordd)
df2=df2[order( n_genes ,ordd,rowMeans(df2),decreasing = T),]
p2=pheatmap(t(df2), 
            cluster_rows = F,
            cluster_cols = F,
            show_rownames = T,
            show_colnames = F)
p2


df3=sm[specific_genes,] 
ord=do.call(c,lapply(dat_list, function(x){
  specific_genes[specific_genes %in% rownames(x)]
}))

p3=pheatmap(t(df3[ord,]),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F)
p3

library(cowplot) 
plot_grid(p1$gtable,
          p2$gtable,
          p3$gtable,nrow = 1,
          labels=c('ubiquitous_genes',
                   'intermediately_genes',
                   'specific_genes') )




###### step5: 如果是蛋白编码基因矩阵  ######
fs=list.files('Rdata/',
              pattern = 'htseq_counts')
fs
# 首先加载每个癌症的表达量矩阵（提取非蛋白编码基因部分）
# 成为一个 list 
dat_list = lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('.htseq_counts..Rdata','',x)
  print(pro)
  
  load(file =  file.path('Rdata/',x)) 
  
  dat=log2(edgeR::cpm(pd_mat)+1)
  dat[1:4,1:4]
  return(dat)
})


tmp=table(unlist(lapply(dat_list, rownames)))
all_genes=names(tmp) #[tmp==33]
length(all_genes)
ubiquitous_genes=names(tmp)[tmp==33]
length(ubiquitous_genes)
specific_genes =names(tmp)[tmp==1]
length(specific_genes)
intermediately_genes=setdiff(setdiff(all_genes,specific_genes),ubiquitous_genes)
length(intermediately_genes)

 
