
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing

install.packages("extrafont")
install.packages("Cairo")
install.packages("showtext")
install.packages("rlang")
devtools::install_github("dreamRs/essquisse")
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(scCancer)
library(extrafont)
library(Cairo)
library(showtext)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(scales)
library(cowplot)


#monocle3
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
library(devtools)
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
packageVersion("monocle3")



dataPath <- "/bcc_scRNA_counts.txt"     # The path of cell ranger processed data
savePath <- "/bcc_scRNA_counts.txt"  # A path to save the results
statPath <-"/bcc_scRNA_counts.txt"
sampleName <- "BCC"          # The sample name
authorName <- "JiangYQ"           # The author name to mark the report
memory.limit(300000)


sampleName = "BCC"
bool.filter.cell = T
bool.filter.gene = T
anno.filter = c("mitochondrial","ribosome","dissociation")
nCell.min = 3
bgPercent.max = 1
bool.rmContamination = F
contamination.fraction = NULL
vars.add.meta = c("mito.percent","ribo.percent","diss.percent")
vars.to.regress = c("nCount_RNA","mito.percent","ribo.percent")
pc.use = 30
resolution = 0.8
clusterStashName = "default"
show.features = NULL
bool.add.features = T
bool.runDiffExpr = T
n.markers = 5
species = "human"
genome = "hg19"
hg.mm.mix = F
bool.runDoublet = T
doublet.method = "bcds"
bool.runCellClassify = T
ct.templates = NULL
coor.names = c("tSNE_1","tSNE_2")
bool.runMalignancy = T
cnv.ref.data = NULL
cnv.referAdjMat = NULL
cutoff = 0.1
p.value.cutoff = 0.5
bool.intraTumor = T
bool.runCellCycle = T
bool.runStemness = T
bool.runGeneSets = T
geneSets = NULL
geneSet.method = "average"
bool.runExprProgram = T
nmf.rank = 50
bool.runInteraction = T
genReport = T

results <- as.list(environment())
checkAnnoArguments(results)

if(is.null(savePath)){
  savePath <- sataPath
}

if(species == "mouse" & genome == "hg19"){
  genome <- "mm10"
}

if(!dir.exists(file.path(savePath, "figures/"))){
  dir.create(file.path(savePath, "figures/"), recursive = T)
}

suppressWarnings( dataPath <- normalizePath(dataPath, "/") )
suppressWarnings( statPath <- normalizePath(statPath, "/") )
suppressWarnings( savePath <- normalizePath(savePath, "/") )
results[["dataPath"]] <- dataPath
results[["statPath"]] <- statPath
results[["savePath"]] <- savePath

#input data
rt=read.table("bcc_scRNA_counts.txt",sep="\t",header=T,check.names=F)


#create seurat
expr <- CreateSeuratObject(counts = rt,project = "seurat")
rm(rt)
gc()
expr[["percent.mt"]] <- PercentageFeatureSet(object = expr, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6)           
VlnPlot(object = expr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#SCTransform
expr <- SCTransform(expr, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = TRUE)

runSeurat <- function(expr,
                      savePath,
                      pc.use = 30, resolution = 0.8,
                      clusterStashName = "default",
                      bool.runDiffExpr = T,
                      comb.method = NULL){

    if(!dir.exists(file.path(savePath))){
        dir.create(file.path(savePath), recursive = T)
    }

    if(!is.null(comb.method)){
        if(comb.method == "Harmony"){
            reduction.type = "harmony"
        }else if(comb.method == "LIGER"){
            reduction.type = "inmf"
            pc.use <- min(pc.use, ncol(expr@reductions$inmf@cell.embeddings))
        }else{
            message("[", Sys.time(), "] -----: PCA")
            expr <- RunPCA(expr, verbose = F)
            reduction.type <- "pca"
        }
    }else{
        message("[", Sys.time(), "] -----: PCA")
        expr <- RunPCA(expr, verbose = F)
        reduction.type <- "pca"
    }

    message("[", Sys.time(), "] -----: clustering")
    expr <- FindNeighbors(expr, reduction = reduction.type, dims = 1:pc.use, verbose = F)
    expr <- FindClusters(expr, resolution = resolution, verbose = F)
    expr[[clusterStashName]] <- as.numeric(Idents(object = expr))

    if(is.null(comb.method)){
        message("[", Sys.time(), "] -----: tSNE")
        expr <- RunTSNE(object = expr, dims = 1:pc.use, reduction = reduction.type)
    }else{
        if(comb.method != "LIGER"){
            message("[", Sys.time(), "] -----: tSNE")
            expr <- RunTSNE(object = expr, dims = 1:pc.use, reduction = reduction.type)
        }
    }

    message("[", Sys.time(), "] -----: UMAP")
    suppressWarnings(
        tryCatch(expr <- RunUMAP(expr, dims = 1:pc.use, reduction = reduction.type, verbose = F),
                 error = function(err) {
                     cat("Error in 'RunUMAP': Please use 'pip install umap-learn' to install UMAP firstly.\n")}
        )
    )

    if(bool.runDiffExpr){
        message("[", Sys.time(), "] -----: differential expression analysis")
        if(!dir.exists(file.path(savePath, "diff.expr.genes"))){
            dir.create(file.path(savePath, "diff.expr.genes"), recursive = T)
        }
        diff.expr.genes <- FindAllMarkers(expr, only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25,
                                          verbose = F)
        # write.table(diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")],
        #             file = file.path(savePath, "diff.expr.genes.txt"), quote = F, sep = "\t", row.names = F)

        diff.expr.genes <- diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")]
        diff.expr.genes$cluster <- as.numeric(diff.expr.genes$cluster)
        nCluster <- length(unique(diff.expr.genes$cluster))
        for(ci in 1:nCluster){
            cur.diff.genes <- subset(diff.expr.genes, cluster == ci)
            cur.diff.genes <- cur.diff.genes[order(cur.diff.genes$avg_logFC, decreasing = T), ]
            write.csv(cur.diff.genes,
                      file = file.path(savePath, "diff.expr.genes", paste0("cluster", ci ,".csv")),
                      quote = F, row.names = F)
        }
        # write.xlsx(subset(diff.expr.genes, cluster == 0),
        #            file = file.path(savePath, "diff.expr.genes.xlsx"),
        #            sheetName = "cluster0", row.names = F)
        # nCluster <- length(unique(diff.expr.genes$cluster))
        # for(ci in 1:(nCluster - 1)){
        #     write.xlsx(subset(diff.expr.genes, cluster == ci),
        #                file = file.path(savePath, "diff.expr.genes.xlsx"),
        #                sheetName = paste0("cluster", ci), append = T, row.names = F)
        # }
    }else{
        diff.expr.genes <- NULL
    }

    cell.annotation <- data.frame(barcodes = colnames(x = expr), stringsAsFactors = F)
    cell.annotation <- cbind(cell.annotation, expr@reductions$tsne@cell.embeddings)
    if("umap" %in% names(expr@reductions)){
        cell.annotation <- cbind(cell.annotation, expr@reductions$umap@cell.embeddings)
    }
    # cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]],
    #                                   levels = 0:(length(unique(expr@meta.data[[clusterStashName]]))-1))
    cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]])

    return(list(expr = expr,
                diff.expr.genes = diff.expr.genes,
                cell.annotation = cell.annotation))
}



t.results <- runSeurat(
  expr = expr,
  savePath = savePath,
  pc.use = pc.use,
  resolution = 0.8,
  clusterStashName = clusterStashName,
  bool.runDiffExpr = F
)
expr = t.results$expr
cell.annotation = t.results$cell.annotation
results[["diff.expr.genes"]] = t.results$diff.expr.genes
rm(t.results)
gc()

SampleTag <- vector("logical",length=ncol(expr))
names(SampleTag) <- colnames(expr)

SampleTag[expr@meta.data[["Treatment.Patients"]]=="su001"]<-"Responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su002"]<-"Responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su003"]<-"Responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su004"]<-"Responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su009"]<-"Responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su012"]<-"Responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su005"]<-"Non-responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su006"]<-"Non-responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su007"]<-"Non-responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su008"]<-"Non-responders"
SampleTag[expr@meta.data[["Treatment.Patients"]]=="su010"]<-"Non-responders"

expr[["efficacy"]] <- SampleTag

#other metadata
metadata<-read.table("bcc_all_metadata.txt",sep="\t",header=T,check.names=F)
metadata$cluster<-as.character(metadata$cluster)
metadata$cluster1<-metadata$cluster
metadata$cluster[metadata$cluster1=="B_cells_1"]<-"B_cells"
metadata$cluster[metadata$cluster1=="B_cells_2"]<-"B_cells"
cellSets<-unique(metadata$cluster)
cell.id<-colnames(expr)
cell.id<-as.data.frame(cell.id)
metadata1<-dplyr::left_join(cell.id,metadata,by="cell.id")

expr[["Cell.Type"]]<-metadata1$cluster
expr[["Treatment.Status"]]<-metadata1$treatment
expr[["Treatment.Patients"]]<-metadata1$patient

premetadata<-metadata[metadata$treatment=="pre",]
postmetadata<-metadata[metadata$treatment=="post",]


preeffectivemetadata<-premetadata[premetadata$patient %in% c("su001","su002","su003","su004","su009","su012"),]
prenoneffectivemetadata<-premetadata[premetadata$patient %in% c("su005","su006","su007","su008","su010"),]
posteffectivemetadata<-postmetadata[postmetadata$patient %in% c("su001","su002","su003","su004","su009","su012"),]
postnoneffectivemetadata<-postmetadata[postmetadata$patient %in% c("su005","su006","su007","su008","su010"),]



#cell type frequency
preeffectivetable<-with(preeffectivemetadata,table(cluster))
preeffectivetable<-as.data.frame(preeffectivetable)
fix(preeffectivetable)

prenoneffectivetable<-with(prenoneffectivemetadata,table(cluster))
prenoneffectivetable<-as.data.frame(prenoneffectivetable)
fix(prenoneffectivetable)

posteffectivetable<-with(posteffectivemetadata,table(cluster))
posteffectivetable<-as.data.frame(posteffectivetable)
fix(posteffectivetable)

postnoneffectivetable<-with(postnoneffectivemetadata,table(cluster))
postnoneffectivetable<-as.data.frame(postnoneffectivetable)
fix(postnoneffectivetable)

celltypefrequency<-cbind(preeffectivetable,prenoneffectivetable,posteffectivetable,postnoneffectivetable)

write.csv(celltypefrequency,file="celltypefrequency.csv")


#subgroup by efficacy and treatment status
Idents(expr)<-"Treatment.Status"
pretreat<-SubsetData(expr,ident.use ="pre" )
posttreat<-SubsetData(expr,ident.use ="post" )

#pretreat effective
Idents(pretreat)<-"Treatment.Patients"
pretreateffective<-SubsetData(pretreat,ident.use =c("su001","su002","su003","su004","su009","su012"))

#posttreat effective
Idents(posttreat)<-"Treatment.Patients"
posttreateffective<-SubsetData(posttreat,ident.use =c("su001","su002","su003","su004","su009","su012"))

  
#pretreat non-effective
  Idents(pretreat)<-"Treatment.Patients"
  pretreatnoneffective<-SubsetData(pretreat,ident.use =c("su005","su006","su007","su008","su010") )
  
  #posttreat non-effective
  Idents(posttreat)<-"Treatment.Patients"
  posttreatnoneffective<-SubsetData(posttreat,ident.use =c("su005","su006","su007","su008","su010") )
  

  
  
  

#Fig1 UMAP
  #pretreateffective
  pretreateffective<-AddMetaData(pretreateffective,pretreateffective@reductions[["umap"]]@cell.embeddings,col.name = colnames(pretreateffective@reductions[["umap"]]@cell.embeddings))
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  class_avg <- pretreateffective@meta.data %>%
    group_by(Cell.Type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  CairoPDF(file="pretreateffective.pdf",width=12.5,height=9) 
  showtext_begin()
  ggplot(pretreateffective@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=Cell.Type))+
    scale_color_manual(values = allcolour)+
    geom_text(aes(label = Cell.Type), data = class_avg)+
    theme(text=element_text(family="A",size=18)) +
    theme(panel.background = element_rect(fill='white', colour='black'), 
          panel.grid=element_blank(), axis.title = element_text(color='black',
                                                                family="A",size=18),axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'), 
          axis.ticks.margin = unit(1,"lines"),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=18),
          axis.title.y=element_text(colour='black', size=18),
          axis.text=element_text(colour='black',size=18),
          legend.title=element_blank(),
          legend.text=element_text(family="A", size=18),
          legend.key=element_blank())+
    theme(plot.title = element_text(size=22,colour = "black",family = "A",face = "bold"))  + 
    guides(colour = guide_legend(override.aes = list(size=6)))
  showtext_end()
  dev.off()
  
  #posttreateffective
  posttreateffective<-AddMetaData(posttreateffective,posttreateffective@reductions[["umap"]]@cell.embeddings,col.name = colnames(posttreateffective@reductions[["umap"]]@cell.embeddings))
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  class_avg <- posttreateffective@meta.data %>%
    group_by(Cell.Type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  CairoPDF(file="posttreateffective.pdf",width=12.5,height=9) 
  showtext_begin()
  ggplot(posttreateffective@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=Cell.Type))+
    scale_color_manual(values = allcolour)+
    geom_text(aes(label = Cell.Type), data = class_avg)+
    theme(text=element_text(family="A",size=18)) +
    theme(panel.background = element_rect(fill='white', colour='black'), 
          panel.grid=element_blank(), axis.title = element_text(color='black',
                                                                family="A",size=18),axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'), 
          axis.ticks.margin = unit(1,"lines"),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=18),
          axis.title.y=element_text(colour='black', size=18),
          axis.text=element_text(colour='black',size=18),
          legend.title=element_blank(),
          legend.text=element_text(family="A", size=18),
          legend.key=element_blank())+
    theme(plot.title = element_text(size=22,colour = "black",family = "A",face = "bold"))  + 
    guides(colour = guide_legend(override.aes = list(size=6)))
  showtext_end()
  dev.off()


  #pretreatnoneffective
  pretreatnoneffective<-AddMetaData(pretreatnoneffective,pretreatnoneffective@reductions[["umap"]]@cell.embeddings,col.name = colnames(pretreatnoneffective@reductions[["umap"]]@cell.embeddings))
  allcolour=c("#DC143C","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  class_avg <- pretreatnoneffective@meta.data %>%
    group_by(Cell.Type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  CairoPDF(file="pretreatnoneffective.pdf",width=12.5,height=9) 
  showtext_begin()
  ggplot(pretreatnoneffective@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=Cell.Type))+
    scale_color_manual(values = allcolour)+
    geom_text(aes(label = Cell.Type), data = class_avg)+
    theme(text=element_text(family="A",size=18)) +
    theme(panel.background = element_rect(fill='white', colour='black'), 
          panel.grid=element_blank(), axis.title = element_text(color='black',
                                                                family="A",size=18),axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'), 
          axis.ticks.margin = unit(1,"lines"),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=18),
          axis.title.y=element_text(colour='black', size=18),
          axis.text=element_text(colour='black',size=18),
          legend.title=element_blank(),
          legend.text=element_text(family="A", size=18),
          legend.key=element_blank())+
    theme(plot.title = element_text(size=22,colour = "black",family = "A",face = "bold"))  + 
    guides(colour = guide_legend(override.aes = list(size=6)))
  showtext_end()
  dev.off()  
  

  #posttreatnoneffective
  posttreatnoneffective<-AddMetaData(posttreatnoneffective,posttreatnoneffective@reductions[["umap"]]@cell.embeddings,col.name = colnames(posttreatnoneffective@reductions[["umap"]]@cell.embeddings))
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  class_avg <- posttreatnoneffective@meta.data %>%
    group_by(Cell.Type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  CairoPDF(file="posttreatnoneffective.pdf",width=12.5,height=9) 
  showtext_begin()
  ggplot(posttreatnoneffective@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=Cell.Type))+
    scale_color_manual(values = allcolour)+
    geom_text(aes(label = Cell.Type), data = class_avg)+
    theme(text=element_text(family="A",size=18)) +
    theme(panel.background = element_rect(fill='white', colour='black'), 
          panel.grid=element_blank(), axis.title = element_text(color='black',
                                                                family="A",size=18),axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'), 
          axis.ticks.margin = unit(1,"lines"),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=18),
          axis.title.y=element_text(colour='black', size=18),
          axis.text=element_text(colour='black',size=18),
          legend.title=element_blank(),
          legend.text=element_text(family="A", size=18),
          legend.key=element_blank())+
    theme(plot.title = element_text(size=22,colour = "black",family = "A",face = "bold"))  + 
    guides(colour = guide_legend(override.aes = list(size=6)))
  showtext_end()
  dev.off()    
  
#TMB Fig 1J
TMB <- read_csv("C:/bcc_scRNA_counts.txt/TMB.csv")
ggplot(TMB) +
  aes(x = Efficacy, y = TMB, fill = Efficacy) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Spectral") +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(vars(TreatmentStatus))





#PDL1 Fig 1K
VlnPlot(PretreatSeurat,features = "CD274",pt.size = 1,group.by = "Cell.Type",split.by = "efficacy")

VlnPlot(PretreatSeurat,features = "CD274",pt.size = 1,group.by = "efficacy")

#celltype specific DEG between responders and non-responders demo
Idents(PretreatSeurat)<-"Cell.Type"
PretreatSeurat$celltype.stim <- paste(Idents(PretreatSeurat), PretreatSeurat$efficacy, sep = "_")
PretreatSeurat$celltype <- Idents(PretreatSeurat)
Idents(PretreatSeurat) <- "celltype.stim"
b.interferon.response <- FindMarkers(PretreatSeurat, ident.1 = "Tumor_2_Responders", ident.2 = "Tumor_2_Non-responders", verbose = FALSE,min.pct = 0,subset.ident = PretreatSeuratCD274)
b.interferon.response$gene<-rownames(b.interferon.response)

DefaultAssay(PretreatSeurat)<-"RNA"
PretreatSeuratCD274 <- WhichCells(object = PretreatSeurat, expression = CD274 > 0)

PretreatSeuratCD274seurat<-PretreatSeurat[,colnames(PretreatSeurat)%in% PretreatSeuratCD274]
b.interferon.response <- FindMarkers(PretreatSeuratCD274seurat, ident.1 = "Tumor_2_Responders", ident.2 = "Tumor_2_Non-responders", verbose = FALSE,min.pct = 0)
b.interferon.response$gene<-rownames(b.interferon.response)

VlnPlot(PretreatSeuratCD274seurat,features = "CD274",pt.size = 1,group.by = "Cell.Type",split.by = "efficacy")



#cell.annotation
cell.annotation<-cbind(cell.annotation,Cell.Type)

cell.annotation$Cell.Type<-cell.annotation$`expr@meta.data[["default"]]`
cell.annotation$Cell.Type<-as.factor(cell.annotation$Cell.Type)

#getcelltypecolor
getCellTypeColor <- function(cell.types){
  cell.colors <- c(
    "T.cells.CD4" = "#07a2a4",
    "T.cells.CD8" = "#9a7fd1",
    "B.cells" = "#588dd5",
    "NK.cells" = "#f5994e",
    "Myeloid.cells" = "#c05050",
    "Endothelial" = "#59678c",
    "Fibroblast" = "#c9ab00",
    "Epithelial" = "#7eb00a",
    "Unknown" = "grey")
  cti = 1
  new.types <- setdiff(cell.types, names(cell.colors))
  for(ct in new.types){
    cell.colors[ct] <- getDefaultColors(n = length(new.types), type = 3)[cti]
    cti = cti + 1
  }
  return(cell.colors)
}





#cellphoneDB process demo
cellphonedb method statistical_analysis \
       --project-name tissue \
       --iterations 1000 \
       --threshold 0.2 \
       --result-precision 3 \
       --output-path /ObsPath/Cells_Communications \
       --verbose \
       --threads 16 \
       /ObsPath/Cell_type.txt \
       /ObsPath/Expression.txt 


#fdr adjustment for cellphoneDB
pvaluesPreeffective<-read.table(file = "pvaluesPreeffective.txt",sep = "\t",header = T)
pvaluesPreeffectiveReshape<-melt(pvaluesPreeffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                                                      "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin"))
pvaluesPreeffectiveReshape<-pvaluesPreeffectiveReshape[pvaluesPreeffectiveReshape$value<0.05,]
pValue=pvaluesPreeffectiveReshape$value
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
pvaluesPreeffectiveReshape=cbind(pvaluesPreeffectiveReshape,fdr=fdr)
pvaluesPreeffectiveReshape<-pvaluesPreeffectiveReshape[pvaluesPreeffectiveReshape$fdr<0.05,]


pvaluesPreNonEffective<-read.table(file = "pvaluesPreNonEffective.txt",sep = "\t",header = T)
pvaluesPreNonEffectiveReshape<-melt(pvaluesPreNonEffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                          "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin"))
pvaluesPreNonEffectiveReshape<-pvaluesPreNonEffectiveReshape[pvaluesPreNonEffectiveReshape$value<0.05,]
pValue=pvaluesPreNonEffectiveReshape$value
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
pvaluesPreNonEffectiveReshape=cbind(pvaluesPreNonEffectiveReshape,fdr=fdr)
pvaluesPreNonEffectiveReshape<-pvaluesPreNonEffectiveReshape[pvaluesPreNonEffectiveReshape$fdr<0.05,]


pvaluesPostNonEffective<-read.table(file = "pvaluesPostNonEffective.txt",sep = "\t",header = T)
pvaluesPostNonEffectiveReshape<-melt(pvaluesPostNonEffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                          "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin"))
pvaluesPostNonEffectiveReshape<-pvaluesPostNonEffectiveReshape[pvaluesPostNonEffectiveReshape$value<0.05,]
pValue=pvaluesPostNonEffectiveReshape$value
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
pvaluesPostNonEffectiveReshape=cbind(pvaluesPostNonEffectiveReshape,fdr=fdr)
pvaluesPostNonEffectiveReshape<-pvaluesPostNonEffectiveReshape[pvaluesPostNonEffectiveReshape$fdr<0.05,]


PostEffective<-read.table(file = "PostEffective.txt",sep = "\t",header = T)
PostEffectiveReshape<-melt(PostEffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                          "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin"))
PostEffectiveReshape<-PostEffectiveReshape[PostEffectiveReshape$value<0.05,]
pValue=PostEffectiveReshape$value
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
PostEffectiveReshape=cbind(PostEffectiveReshape,fdr=fdr)
PostEffectiveReshape<-PostEffectiveReshape[PostEffectiveReshape$fdr<0.05,]


#cellphoneDB data reshape
#PreEffective
significant_meansPreeffective<-read.table("significant_meansPreeffective.txt",sep="\t",header = T)
significant_meansPreeffective<-significant_meansPreeffective[!duplicated(significant_meansPreeffective$id_cp_interaction),]
write.csv(significant_meansPreeffective,file="significant_meansPreeffectiveUnique.csv")

significant_meansPreeffectiveReshape<-melt(significant_meansPreeffective,id=c("id_cp_interaction","interacting_pair","partner_a",
"partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"))
significant_meansPreeffectiveReshape<-significant_meansPreeffectiveReshape[complete.cases(significant_meansPreeffectiveReshape),]
write.csv(significant_meansPreeffectiveReshape,file = "significant_meansPreeffectiveReshape.csv")

#PreNonEffective
significant_meansPreNonEffective<-read.table("significant_meansPreNonEffective.txt",sep="\t",header = T)
significant_meansPreNonEffective<-significant_meansPreNonEffective[!duplicated(significant_meansPreNonEffective$id_cp_interaction),]
write.csv(significant_meansPreNonEffective,file="significant_meansPreNonEffectiveUnique.csv")

significant_meansPreNonEffectiveReshape<-melt(significant_meansPreNonEffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                                              "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"))
significant_meansPreNonEffectiveReshape<-significant_meansPreNonEffectiveReshape[complete.cases(significant_meansPreNonEffectiveReshape),]
write.csv(significant_meansPreNonEffectiveReshape,file = "significant_meansPreNonEffectiveReshape.csv",row.names = F)


#PostEffective
significant_meansPostEffective<-read.table("significant_meansPostEffective.txt",sep="\t",header = T)
significant_meansPostEffective<-significant_meansPostEffective[!duplicated(significant_meansPostEffective$id_cp_interaction),]
write.csv(significant_meansPostEffective,file="significant_meansPostEffectiveUnique.csv",row.names = F)

significant_meansPostEffectiveReshape<-melt(significant_meansPostEffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                                                    "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"))
significant_meansPostEffectiveReshape<-significant_meansPostEffectiveReshape[complete.cases(significant_meansPostEffectiveReshape),]
write.csv(significant_meansPostEffectiveReshape,file = "significant_meansPostEffectiveReshape.csv",row.names = F)


#PostNonEffective
significant_meansPostNonEffective<-read.table("significant_meansPostNonEffective.txt",sep="\t",header = T)
significant_meansPostNonEffective<-significant_meansPostNonEffective[!duplicated(significant_meansPostNonEffective$id_cp_interaction),]
write.csv(significant_meansPostNonEffective,file="significant_meansPostNonEffectiveUnique.csv",row.names = F)

significant_meansPostNonEffectiveReshape<-melt(significant_meansPostNonEffective,id=c("id_cp_interaction","interacting_pair","partner_a",
                                                                                "partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"))
significant_meansPostNonEffectiveReshape<-significant_meansPostNonEffectiveReshape[complete.cases(significant_meansPostNonEffectiveReshape),]
write.csv(significant_meansPostNonEffectiveReshape,file = "significant_meansPostNonEffectiveReshape.csv",row.names = F)

#cellphoneDB Summary Score Map
SummaryScore<-read.csv("SummaryScore.csv",header = T)


#Fig 3A-D demo
#Pre Post Res Non-Res Status
#Pre Res
CairoPDF(file="PreRes.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                                    color = log2(PreEffectiveScore+1),size=PreEffectiveNum)) +
  geom_point() +
  coord_fixed() +
  scale_color_gradient2(low = "#ffffff",high = "#ea1a59") +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Score+1)", size = "Number of Ligand-Receptor pairs")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = 0, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Responsers (Pretreatment Status)")
dev.off() 

#Pre NonRes
CairoPDF(file="PreNonRes.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                         color = log2(PreNonEffectiveScore+1),size=PreNonEffectiveNum)) +
  geom_point() +
  coord_fixed() +
  scale_color_gradient2(low = "#ffffff",high = "#ea1a59") +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Score+1)", size = "Number of Ligand-Receptor pairs")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = 0, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Non-Responsers (Pretreatment Status)")

dev.off() 


#Post Res
CairoPDF(file="PostRes.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                         color = log2(PostEffectiveScore+1),size=PostEffectiveNum)) +
  geom_point() +
  coord_fixed() +
  scale_color_gradient2(low = "#ffffff",high = "#ea1a59") +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Score+1)", size = "Number of Ligand-Receptor pairs")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = 0, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Responsers (Posttreatment Status)")

dev.off() 

#Post NonRes
CairoPDF(file="PostNonRes.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                         color = log2(PostNonEffectiveScore+1),size=PostNonEffectiveNum)) +
  geom_point() +
  coord_fixed() +
  scale_color_gradient2(low = "#ffffff",high = "#ea1a59") +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Score+1)", size = "Number of Ligand-Receptor pairs")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = 0, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Non-Responsers (Posttreatment Status)")

dev.off() 

#Fig2 C demo
#PreEffectiveVsNonEffective
#SummaryScore[,3:18]<-lapply(SummaryScore[,3:18], function(x) as.numeric(as.character(x)))
CairoPDF(file="PreEffectiveVsNonEffective.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                         size =PreEffectiveVsNonEffectiveNumRatio , color = log2(PreEffectiveVsNonEffectiveScoreRatio))) +
  geom_point() +
  coord_fixed() +
  scale_size(limits=c(0,5),range = c(1,5))+
  scale_color_gradient2(low = "#2e68b7",mid = "#ffffff",high = "#ea1a59",midpoint = 0) +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Ratio of Score)", size = "Ratio of Pair Counts")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = +0.5, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Pre-Treatment Patients (Responsers Vs. Non-Responsers)")

dev.off()   


Fig 3C/D demo
#Effective PostVsPre
#SummaryScore[,3:18]<-lapply(SummaryScore[,3:18], function(x) as.numeric(as.character(x)))
CairoPDF(file="EffectivePostVsPre.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                            size =EffectivePostVsPreNumRatio , color = log2(EffectivePostVsPreScoreRatio))) +
  geom_point() +
  coord_fixed() +
  scale_size(limits=c(0,5),range = c(1,5))+
  scale_color_gradient2(low = "#2e68b7",mid = "#ffffff",high = "#ea1a59",midpoint = 0) +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Ratio of Score)", size = "Ratio of Pair Counts")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = +0.5, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Responsers (Post Vs. Pre)")

dev.off()   


#NonEffective PostVsPre
#SummaryScore[,3:18]<-lapply(SummaryScore[,3:18], function(x) as.numeric(as.character(x)))
CairoPDF(file="NonEffectivePostVsPre.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                         size =NonEffectivePostVsPreNumRatio , color = log2(NonEffectivePostVsPreScoreRatio))) +
  geom_point() +
  coord_fixed() +
  scale_size(limits=c(0,5),range = c(1,5))+
  scale_color_gradient2(low = "#2e68b7",mid = "#ffffff",high = "#ea1a59",midpoint = 0) +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Ratio of Score)", size = "Ratio of Pair Counts")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = +0.5, vjust = 3)) +labs(title = "Cell-Cell Crosstalk in Non-Responsers (Post Vs. Pre)")

dev.off()   

Fig3 G demo
#ResVsNonResofPostVsPre
#SummaryScore[,3:18]<-lapply(SummaryScore[,3:18], function(x) as.numeric(as.character(x)))
CairoPDF(file="ResVsNonResofPostVsPre.pdf",width=8,height=7) 

ggplot(SummaryScore, aes(x = Ligand.Cellset, y = Receptor.Cellset,
                         size =ResVsNonResNumRatio , color = log2(ResVsNonResScoreRatio))) +
  geom_point() +
  coord_fixed() +
  scale_size(limits=c(0,5),range = c(1,5))+
  scale_color_gradient2(low = "#2e68b7",mid = "#ffffff",high = "#ea1a59",midpoint = 0) +
  theme_light() +
  theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        axis.title.y = element_text(size = 13, vjust = 1.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black")) + theme(axis.line = element_line(size = 0.8, 
                                                                                       linetype = "solid"), axis.text = element_text(family = "Times", 
                                                                                                                                     size = 7, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                             colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                                                                           colour = "black"), legend.text = element_text(family = "Times"), 
                                                              legend.title = element_text(family = "Times"), 
                                                              panel.background = element_rect(fill = NA), 
                                                              legend.key = element_rect(fill = NA), 
                                                              legend.background = element_rect(fill = NA)) +labs(colour = "Log2(Ratio Change Comparison)", size = "Pair Counts Change Comparison")+
  theme(axis.title = element_text(family = "Times"), plot.title = element_text(family = "Times", hjust = +0.5, vjust = 3)) +labs(title = "Change Comparison between Responsers and Non-Responsers")

dev.off()   


#Fig2D demo
#significant_meansReshapePreEffectiveVsNonEffective
significant_meansReshapePreEffectiveVsNonEffective<-read.csv("significant_meansReshapePreEffectiveVsNonEffective.csv",header = T)
significant_meansReshapePreEffectiveVsNonEffective$interacting_cell<-gsub("\\.","-",significant_meansReshapePreEffectiveVsNonEffective$interacting_cell)
#all
CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffective.pdf",height=60,width=65) 

ggplot(significant_meansReshapePreEffectiveVsNonEffective,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                       axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                       axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                       axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                  colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                       legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                            fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))

dev.off()

#Obvious
significant_meansReshapePreEffectiveVsNonEffective<-significant_meansReshapePreEffectiveVsNonEffective[significant_meansReshapePreEffectiveVsNonEffective$FoldChange==">2 or <0.5",]
CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffectiveObvious.pdf",height=35,width=35) 

ggplot(significant_meansReshapePreEffectiveVsNonEffective,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values = 23)+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Comparison of Cell-Cell Crosstalk between Responsers and Non-Responsers (Pre-Treatment Status)")+
  theme(axis.title = element_text(size = 30), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20), 
    legend.key = element_rect(size = 5), 
    legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
    vjust = 2))
dev.off()


#significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq<-significant_meansReshapePreEffectiveVsNonEffective[(significant_meansReshapePreEffectiveVsNonEffective$gene_aName %in% significant_meansReshapePreEffectiveVsNonEffective$anti.PD1RNASEQDEG.preEffectiveVsNonEffective.GSE78220)|(significant_meansReshapePreEffectiveVsNonEffective$gene_bName %in% significant_meansReshapePreEffectiveVsNonEffective$anti.PD1RNASEQDEG.preEffectiveVsNonEffective.GSE78220),]

CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq.pdf",height=20,width=25) 

ggplot(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 25)) +labs(title = "Intersecting Ligand-Receptor Pairs with Hugo et al. Database")+
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()


#significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious<-significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq[significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$FoldChange==">2 or <0.5",]
CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious.pdf",height=9,width=9.5) 

ggplot(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values = 23)+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 13), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 16)) +labs(title = "Intersecting Ligand-Receptor Pairs with Hugo et al. Database (Pairs with Foldchange >2 or <0.5)")+
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))
dev.off()



#significant_meansReshapeEffectivePostVsPre
significant_meansReshapeEffectivePostVsPre<-read.csv("significant_meansReshapeEffectivePostVsPre.csv",header = T)
significant_meansReshapeEffectivePostVsPre$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectivePostVsPre$interacting_cell)
#all
CairoPDF(file="significant_meansReshapeEffectivePostVsPre.pdf",height=60,width=65) 

ggplot(significant_meansReshapeEffectivePostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectivePostVsPre$Ratio)),0,max(log(significant_meansReshapeEffectivePostVsPre$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Cell-Cell Crosstalk of Post- Vs. Pre-Treatment Status(Responsers)")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()

#Obvious
significant_meansReshapeEffectivePostVsPre<-significant_meansReshapeEffectivePostVsPre[significant_meansReshapeEffectivePostVsPre$FoldChange==">2 or <0.5",]
CairoPDF(file="significant_meansReshapeEffectivePostVsPreObvious.pdf",height=28,width=30) 

ggplot(significant_meansReshapeEffectivePostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values = 23)+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectivePostVsPre$Ratio)),0,max(log(significant_meansReshapeEffectivePostVsPre$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Cell-Cell Crosstalk of Post- Vs. Pre-Treatment Status (Responsers)")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))
dev.off()

#significant_meansReshapeNonEffectivePostVsPre
significant_meansReshapeNonEffectivePostVsPre<-read.csv("significant_meansReshapeNonEffectivePostVsPre.csv",header = T)
significant_meansReshapeNonEffectivePostVsPre$interacting_cell<-gsub("\\.","-",significant_meansReshapeNonEffectivePostVsPre$interacting_cell)
#all
CairoPDF(file="significant_meansReshapeNonEffectivePostVsPre.pdf",height=75,width=70) 

ggplot(significant_meansReshapeNonEffectivePostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeNonEffectivePostVsPre$Ratio)),0,max(log(significant_meansReshapeNonEffectivePostVsPre$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Cell-Cell Crosstalk of Post- Vs. Pre-Treatment Status (Non-Responsers)")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()

#Obvious
significant_meansReshapeNonEffectivePostVsPre<-significant_meansReshapeNonEffectivePostVsPre[significant_meansReshapeNonEffectivePostVsPre$FoldChange==">2 or <0.5",]
CairoPDF(file="significant_meansReshapeNonEffectivePostVsPreObvious.pdf",height=25,width=32) 

ggplot(significant_meansReshapeNonEffectivePostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values = 23)+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeNonEffectivePostVsPre$Ratio)),0,max(log(significant_meansReshapeNonEffectivePostVsPre$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Cell-Cell Crosstalk of Post- Vs. Pre-Treatment Status (Non-Responsers)")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()


#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll
#Ligand Scene(all)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell)
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.pdf",height=60,width=60) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))

dev.off()

#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll[(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_aName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061|significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_bName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061),]
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq2.pdf",height=25,width=60) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio),shape=FoldChange.EffectiveCombineNonEffective.))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  scale_shape_manual(values=c(23,22))+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 40)) +labs(title = "Change Comparison between Responsers and Non-Responsers Intersecting Riaz et al. Database")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()

#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll[(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_aName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061|significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_bName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061),]
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$FoldChange==">2 or <0.5",]
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious.pdf",height=20,width=40) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",size=5,shape=23)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                         axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                         axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                         axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                    colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                         legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                              fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 14), legend.text = element_text(size = 18,family = "Times"),legend.title = element_text(size = 18,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 40)) +labs(title = "Change Comparison between Responsers and Non-Responsers Intersecting Riaz et al. Database")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()



#Ligand Scene(Obvious)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$interacting_cell)
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreObvious.pdf",height=40,width=43) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",size=5,shape=23)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                         axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                         axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                         axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                    colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                         legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                              fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 50)) +labs(title = "Change Comparison between Responsers and Non-Responsers")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()

#Receptor Scene(Obvious)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$interacting_cell)
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreObviousReceptorScene.pdf",height=40,width=43) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=Receptor.LigandCellset,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",size=5,shape=23)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                 axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                 axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                 axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                            colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                 legend.title = element_text(family = "Times")) +labs(x = "Receptor-Ligand Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                      fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Change Comparison between Responsers and Non-Responsers")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()

#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreIntersectingAntiPd1RnaSeq
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreIntersectingAntiPd1RnaSeq<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$gene_aName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$anti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061)|(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$gene_bName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$anti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061),]
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreIntersectingAntiPd1RnaSeq.pdf",height=22,width=36) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreIntersectingAntiPd1RnaSeq,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",size=5,shape=23)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                 axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                 axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                 axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                            colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                 legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                      fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 40)) +labs(title = "Change Comparison between Responsers and Non-Responsers Intersecting Riaz et al. Database")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()



#Fig4/5 demo
#Sole Map of Cell Types

significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll[(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_aName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061|significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_bName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061),]
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$FoldChange==">2 or <0.5",]

#Fig4/5 demo
#CD8_T_cells/NK
CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$LigandCellset %in% c("CD8_act_T_cells","CD8_ex_T_cells","CD8_mem_T_cells","NK_cells"),]

PDF(file="CD8_T_cellsNK_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre2.pdf",height=12,width=32) 

CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1<-CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveRatio),]
p1<-ggplot(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio)),0,max(log2(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                 axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                 axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                 axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                            colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                 legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                      fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))
CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2<-CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$NonEffectiveRatio),]
p2<-ggplot(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(NonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio)),0,max(log2(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                       axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                  colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                            fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Non-Responders (Post Vs. Pre)") + theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                           vjust = 2))



p3<-ggplot(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log2(CD8_T_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                         axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                         axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                         axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                    colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                         legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                              fill = "Log2(Ratio Change Comparison)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5, vjust = 2)) +labs(title = "Change Comparison between Responders and Non-Responders")
plot_grid(p1, p2, p3, align = "w", ncol = 3, rel_widths =  c(0.9, 0.9, 1))


dev.off()


#EndothelialCD4TregTprolifBcellPlasma
CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$LigandCellset %in% c("Endothelial","CD4_T_cells","Tregs","Tcell_prolif","B_cells","Plasma_cells"),]

PDF(file="EndothelialCD4TregTprolifNK_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre3.pdf",height=12,width=32) 

CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1<-CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveRatio),]
p1<-ggplot(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio)),0,max(log2(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                 axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                 axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                 axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                            colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                 legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                      fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))
CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2<-CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$NonEffectiveRatio),]
p2<-ggplot(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(NonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio)),0,max(log2(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                       axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                  colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                            fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Non-Responders (Post Vs. Pre)") + theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                           vjust = 2))



p3<-ggplot(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log2(CD4TregTprolifsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                         axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                         axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                         axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                    colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                         legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                              fill = "Log2(Ratio Change Comparison)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5, vjust = 2)) +labs(title = "Change Comparison between Responders and Non-Responders")
plot_grid(p1, p2, p3, align = "w", ncol = 3, rel_widths =  c(0.9, 0.9, 1))


dev.off()


#CAFsMyofibroblasts
CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$LigandCellset %in% c("CAFs","Myofibroblasts"),]

PDF(file="CAFsMyofibroblastssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre3.pdf",height=12,width=32) 

CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1<-CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveRatio),]
p1<-ggplot(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio)),0,max(log2(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                       axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                  colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                       legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                            fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))
CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2<-CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$NonEffectiveRatio),]
p2<-ggplot(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(NonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio)),0,max(log2(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                             axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                             axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                             axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                        colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                             legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                  fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Non-Responders (Post Vs. Pre)") + theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                           vjust = 2))



p3<-ggplot(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log2(CAFsMyofibroblastsEndothelialsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                               axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                               axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                               axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                          colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                               legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                    fill = "Log2(Ratio Change Comparison)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5, vjust = 2)) +labs(title = "Change Comparison between Responders and Non-Responders")
plot_grid(p1, p2, p3, align = "w", ncol = 3, rel_widths =  c(0.9, 0.9, 1))


dev.off()



#MacrophagesDCspDCs
MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$LigandCellset %in% c("Macrophages","DCs","pDCs"),]

PDF(file="MacrophagesDCspDCsNK_cellssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre3.pdf",height=12,width=32) 

MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1<-MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveRatio),]
p1<-ggplot(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio)),0,max(log2(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                     axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                     axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                     axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                     legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                          fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))
MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2<-MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$NonEffectiveRatio),]
p2<-ggplot(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(NonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio)),0,max(log2(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                           axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                           axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                           axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                      colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                           legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                                fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Non-Responders (Post Vs. Pre)") + theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                           vjust = 2))



p3<-ggplot(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log2(MacrophagesDCspDCssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                                        colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                                             legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                                                  fill = "Log2(Ratio Change Comparison)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5, vjust = 2)) +labs(title = "Change Comparison between Responders and Non-Responders")
plot_grid(p1, p2, p3, align = "w", ncol = 3, rel_widths =  c(0.9, 0.9, 1))


dev.off()


#Tumors
Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$LigandCellset %in% c("Tumor_1","Tumor_2"),]

PDF(file="Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre3.pdf",height=12,width=32) 

Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1<-Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveRatio),]
p1<-ggplot(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio)),0,max(log2(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                     axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                     axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                     axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                     legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                          fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))
Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2<-Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$NonEffectiveRatio),]
p2<-ggplot(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(NonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio)),0,max(log2(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                           axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                           axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                           axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                      colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                           legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                                fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Non-Responders (Post Vs. Pre)") + theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                           vjust = 2))



p3<-ggplot(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log2(Tumorssignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                                        colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                                             legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                                                  fill = "Log2(Ratio Change Comparison)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5, vjust = 2)) +labs(title = "Change Comparison between Responders and Non-Responders")
plot_grid(p1, p2, p3, align = "w", ncol = 3, rel_widths =  c(0.9, 0.9, 1))


dev.off()

#CD4T
CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqObvious$LigandCellset %in% c("CD4_T_cells"),]

PDF(file="CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre3.pdf",height=12,width=32) 

CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1<-CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveRatio),]
p1<-ggplot(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio)),0,max(log2(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP1$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                       axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                       axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                       axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                  colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                       legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                            fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))
CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2<-CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre[complete.cases(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$NonEffectiveRatio),]
p2<-ggplot(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(NonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio)),0,max(log2(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreP2$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                             axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                             axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                             axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                        colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                             legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                  fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Non-Responders (Post Vs. Pre)") + theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                           vjust = 2))



p3<-ggplot(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=23,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio)),0,max(log2(CD4Tsignificant_meansReshapeEffectiveVsNonEffectiveOfPostVsPre$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                               axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                               axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                               axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                          colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                               legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                    fill = "Log2(Ratio Change Comparison)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10, 
                                                                                         family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5, vjust = 2)) +labs(title = "Change Comparison between Responders and Non-Responders")
plot_grid(p1, p2, p3, align = "w", ncol = 3, rel_widths =  c(0.9, 0.9, 1))


dev.off()


#Fig S2 demo
significant_meansReshapePreEffectiveVsNonEffective<-read.csv("significant_meansReshapePreEffectiveVsNonEffective.csv",header = T)
significant_meansReshapePreEffectiveVsNonEffective$interacting_cell<-gsub("\\.","-",significant_meansReshapePreEffectiveVsNonEffective$interacting_cell)
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq<-significant_meansReshapePreEffectiveVsNonEffective[(significant_meansReshapePreEffectiveVsNonEffective$gene_aName %in% significant_meansReshapePreEffectiveVsNonEffective$anti.PD1RNASEQDEG.preEffectiveVsNonEffective.GSE78220)|(significant_meansReshapePreEffectiveVsNonEffective$gene_bName %in% significant_meansReshapePreEffectiveVsNonEffective$anti.PD1RNASEQDEG.preEffectiveVsNonEffective.GSE78220),]
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious<-significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq[significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$FoldChange==">2 or <0.5",]
ggplot(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values = 23)+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqObvious$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                             axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                             axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                             axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                        colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                             legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                  fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 13), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 16)) +labs(title = "Intersecting Ligand-Receptor Pairs with Hugo et al. Database (Pairs with Foldchange >2 or <0.5)")+
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

#Wilcoxon rank-sum (Mann-Whitney U) tests demo
#Fig 2E demo
data<-read.csv("significant_meansReshapePreEffectiveVsNonEffective.csv",header = T)
fdrFilter=0.05                                                  
logFCfilter=0.5                                                     

outTab=data.frame()
for(i in data$interacting_pair){
  pair1<-data[data$interacting_pair==i,]
  pair1<-pair1[,"PreNonEffectiveScore"]
  pair1<-na.omit(pair1)
  pair1<-as.vector(pair1)
  pair2<-data[data$interacting_pair==i,]
  pair2<-pair2[,"PreEffectiveScore"]
  pair2<-na.omit(pair2)
  pair2<-as.vector(pair2)
  grade<-c(rep(1,length(pair1)),rep(2,length(pair2)))
  grade<-as.data.frame(grade)
  rt<-c(pair1,pair2)
  rt<-as.data.frame(rt)
  rt<-cbind(rt,grade)
  colnames(rt)[1]<-"interactionscore"
  if ((length(pair1)>0) & (length(pair2)>0)){
    wilcoxTest<-wilcox.test(interactionscore ~ grade, data=rt)
    conNum=length(pair1)                                                   
    treatNum=length(pair2)
    conGeneMeans=mean(pair1)
    treatGeneMeans=mean(pair2)
    logFC=log2(treatGeneMeans+0.001)-log2(conGeneMeans+0.001)
    pvalue=wilcoxTest$p.value
    conMed=median(pair1)
    treatMed=median(pair2)
    diffMed=treatMed-conMed
    if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
      outTab=rbind(outTab,cbind(interactionpair=i,NonResponderMean=conGeneMeans,ResponderMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
    }
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
PreResVsNonResoutDiff<-outDiff[!duplicated(outDiff$interactionpair),]
write.csv(PreResVsNonResoutDiff,file = "PreResVsNonResDEinteraction.csv")

#significant_meansReshapePreEffectiveVsNonEffective
significant_meansReshapePreEffectiveVsNonEffective<-read.csv("significant_meansReshapePreEffectiveVsNonEffective.csv",header = T)
significant_meansReshapePreEffectiveVsNonEffective$interacting_cell<-gsub("\\.","-",significant_meansReshapePreEffectiveVsNonEffective$interacting_cell)
antiPD1RNASEQDEGpreEffectiveVsNonEffectiveGSE78220<-read.csv("antiPD1RNASEQDEGpreEffectiveVsNonEffectiveGSE78220.csv",header = T)
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq<-significant_meansReshapePreEffectiveVsNonEffective[(significant_meansReshapePreEffectiveVsNonEffective$gene_aName %in% antiPD1RNASEQDEGpreEffectiveVsNonEffectiveGSE78220$antiPD1RNASEQDEGpreEffectiveVsNonEffectiveGSE78220)|(significant_meansReshapePreEffectiveVsNonEffective$gene_bName %in% antiPD1RNASEQDEGpreEffectiveVsNonEffectiveGSE78220$antiPD1RNASEQDEGpreEffectiveVsNonEffectiveGSE78220),]
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq<-significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq[significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$interacting_pair %in% PreResVsNonResoutDiff$interactionpair,]

CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1RnaseqWilcoxon.pdf",height=10,width=10.7) 

ggplot(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio)))+
  geom_point(colour="black",size=8.4,shape=22)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                               axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                               axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                               axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                          colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                               legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                    fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 13)) +labs(title = "Intersecting Ligand-Receptor Pairs (p<0.05) with Hugo et al. Database")+
  theme(axis.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 12))+
  theme(axis.ticks = element_line(linetype = "blank"), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    panel.background = element_rect(fill = NA), 
    legend.key = element_rect(fill = NA), 
    legend.background = element_rect(fill = NA))
dev.off()

#heatmap
outDiffreshape<-PreResVsNonResoutDiff[,c("interactionpair","NonResponderMean","ResponderMean")]
outDiffreshape<-outDiffreshape[outDiffreshape$interactionpair %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$interacting_pair,]
outDiffreshape<-melt(outDiffreshape,id="interactionpair")
outDiffreshape$value<-as.character(outDiffreshape$value)
outDiffreshape$value<-as.numeric(outDiffreshape$value)
outDiffreshape<-outDiffreshape[order(outDiffreshape$value),]
CairoPDF(file="IntersectingAntiPd1RnaSeqWilcoxonHeatmapPreResVsNonRes.pdf",height=4,width=8) 

ggplot(outDiffreshape,aes(x=variable,y=interactionpair,fill=log(value)))+
  geom_point(colour="black",size=8.4,shape=22)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(outDiffreshape$value)),0,max(log(outDiffreshape$value)))))+ 
  theme(axis.ticks = element_line(linetype = "blank"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(family = "Times", 
                                  size = 9), axis.text = element_text(family = "Times", 
                                                                      size = 9, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                              size = 9, colour = "black", angle = 90), 
        plot.title = element_text(family = "Times", 
                                  size = 10), legend.text = element_text(size = 9, 
                                                                         family = "Times"), legend.title = element_text(size = 9, 
                                                                                                                        family = "Times"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA)) +labs(x = "Pre-treatment Responders vs. Pre-treatment Non-responders", 
                                                           y = "Cell-Cell Interaction Pairs", fill = "Log2 (Score)")
dev.off()




##predicted model Fig 8 demo
#"set1roc.csv" and "set2roc.csv" for ¡°Ligand-receptor Pairs Related to Response Before Treatment¡± model
#"set1prepostcalculated.csv" for ¡°Ligand-receptor Pairs Related to Response on Treatment¡± model

#set1
finalgenes<-read.csv("final genes.csv",header = T)

rt=read.table("Melanoma_PRJEB23709_reRun_fpkm.txt",sep="\t",header=T,check.names=F)

rt<-rt[rt$gene_symbol %in% finalgenes$gene,]
rt<-t(rt)
colnames(rt)<-rt[1,]
rt<-as.data.frame(rt)
rt<-rt[-1,]

metadata<-read.table("Melanoma_PRJEB23709_reRun_clinicaldata.txt",sep="\t",header=T,check.names=F)
metadata<-metadata[metadata$treatment=="Pembrolizumab"|metadata$treatment=="Nivolumab",]
metadatapre<-metadata[metadata$Treatment=="PRE",]

rt<-rt[row.names(rt) %in% metadatapre$Sample_ID,]
rt$Sample_ID<-row.names(rt)
rt<-dplyr::left_join(rt,metadatapre,by="Sample_ID")

rt$efficacy[rt$response=="CR"|rt$response=="PR"|rt$response=="SD"]<-"responser"
rt$efficacy[rt$response=="PD"]<-"nonresponser"
fix(rt)

write.csv(rt,file="set1roc.csv")


#set1 pre vs post

rt=read.table("Melanoma_PRJEB23709_LR.txt",sep="\t",header=T,check.names=F)

rt<-rt[rt$gene_symbol %in% finalgenes$gene,]
rt<-t(rt)
colnames(rt)<-rt[1,]
rt<-as.data.frame(rt)
rt<-rt[-1,]

metadataset1prepost<-read.csv("prevspostID.csv",header = T)
rt$SampleID<-row.names(rt)
set1prepost<-dplyr::left_join(metadataset1prepost,rt,by="SampleID")
write.csv(set1prepost,file = "set1prepost.csv")




#set2
rt2<-read.table("GSE145996_Melanoma_Immune_LR.txt",sep="\t",header=T,check.names=F)
rt2<-rt2[,-2]
rt2<-rt2[rt2$`Gene Name` %in% finalgenes$gene,]

rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp=rt2[,2:ncol(rt2)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data<-t(data)
data<-as.data.frame(data)
data$sampleID<-row.names(data)
data$sampleID<-str_sub(data$sampleID,-7,-1)
data$sampleID<-gsub("_","",data$sampleID)
metadata<-read.csv("sample response.csv",sep=",",header=T,check.names=F)
metadata<-metadata[,-1]
colnames(metadata)[1]<-"sampleID"
rt2<-dplyr::left_join(data,metadata,by="sampleID")
fix(rt2)
write.csv(rt2,"set3 roc.csv")




install.packages("xgboost")
require(xgboost)


##model demo
model1Label<-read.table(file = "model1_SRP070710_reRun_clinicaldata.txt",sep="\t",header=T,check.names=F)
model1Label<-model1Label[model1Label$Treatment=="PRE",]
model1Label$efficacy[model1Label$response=="PRCR"]<-"1"
model1Label$efficacy[model1Label$response=="PD"]<-"0"

model1LR<-read.table(file = "model1LR.txt",sep="\t",header=T,check.names=F)
model1LR=as.matrix(model1LR)
rownames(model1LR)=model1LR[,1]
exp=model1LR[,2:ncol(model1LR)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
model1LR<-data
model1LR<-model1LR[,colnames(model1LR) %in% model1Label$Sample_ID]
model1LR<-model1LR[rownames(model1LR) %in% finalgenes$gene,]
rm(data,exp)
gc()
model1LR<-t(model1LR)
model1LR<-as(as.matrix(model1LR),"dgCMatrix")
model1ID<-rownames(model1LR)
model1ID<-as.data.frame(model1ID)
colnames(model1ID)[1]<-"Sample_ID"
model1ID<-dplyr::left_join(model1ID,model1Label,by="Sample_ID")
model1ID$efficacy<-as.double(model1ID$efficacy)
trainmodel1<-list(model1LR,model1ID$efficacy)
names(trainmodel1) <- c("data", "label")



model1vLR<-read.csv(file="set1roc.csv",header = T)
model1vLR=as.matrix(model1vLR)
rownames(model1vLR)=model1vLR[,1]
exp=model1vLR[,2:ncol(model1vLR)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
model1vLR<-data
model1vLR<-model1vLR[,colnames(model1vLR) %in% model1Label$Sample_ID]
model1vLR<-model1vLR[rownames(model1vLR) %in% finalgenes$gene,]
rm(data,exp)
gc()
model1vLR<-t(model1vLR)
model1vLR<-as(as.matrix(model1vLR),"dgCMatrix")
model1ID<-rownames(model1vLR)
model1ID<-as.data.frame(model1ID)
colnames(model1ID)[1]<-"Sample_ID"

testmodel1<-model1vLR
names(trainmodel1) <- c("data", "label")

library(xgboost)

dtrain = xgb.DMatrix(data=train$data, label=train$label)
dtest = xgb.DMatrix(data=test$data, label=test$label)
watchlist = list(train=dtrain, test=dtest)

bst = xgb.train(data=dtrain, max_depth=2, eta=1, nthread=2, nrounds=2, watchlist=watchlist, objective='binary:logistic')

pred = predict(bst, test$data)
prediction = as.numeric(pred > 0.5)
print(head(prediction, 5))
err = mean(as.numeric(prediction != test$label))

rocl<-cbind(model1vLR,pred)

#plot

rocobj3<-roc(rocl$efficacy,rocl$pred)
auc(rocobj3)
plot(rocobj3)
g <- ggroc(rocobj3)
g + theme_minimal() + ggtitle("ROC curve") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
  
  
  

#Fig6 demo
install.packages("DMwR")
library(DMwR)
data<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.csv",header = T)

fdrFilter=0.05                                                  
logFCfilter=1                                                   

outTab=data.frame()
for(i in data$interacting_pair){
  pair1<-data[data$interacting_pair==i,]
  pair1<-pair1[,"NonEffectiveRatio"]
  pair1<-na.omit(pair1)
  pair1<-as.vector(pair1)
  pair2<-data[data$interacting_pair==i,]
  pair2<-pair2[,"EffectiveRatio"]
  pair2<-na.omit(pair2)
  pair2<-as.vector(pair2)
  grade<-c(rep(1,length(pair1)),rep(2,length(pair2)))
  grade<-as.data.frame(grade)
  rt<-c(pair1,pair2)
  rt<-as.data.frame(rt)
  rt<-cbind(rt,grade)
  colnames(rt)[1]<-"interactionscore"
  if ((length(pair1)>0) & (length(pair2)>0)){
  wilcoxTest<-wilcox.test(interactionscore ~ grade, data=rt)
  conNum=length(pair1)                                                   
  treatNum=length(pair2)
  conGeneMeans=mean(pair1)
  treatGeneMeans=mean(pair2)
  logFC=log2(treatGeneMeans+0.001)-log2(conGeneMeans+0.001)
  pvalue=wilcoxTest$p.value
  conMed=median(pair1)
  treatMed=median(pair2)
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(interactionpair=i,NonResponderChangeMean=conGeneMeans,ResponderChangeMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
ResVsNonResChangesoutDiff<-outDiff[!duplicated(outDiff$interactionpair),]
write.csv(ResVsNonResChangesoutDiff,file = "ResVsNonResChangesDEinteraction.csv")



#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll
#Ligand Scene(all)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell)

#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq
IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061<-read.csv("IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll[(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_aName %in% IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061|significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_bName %in% IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061),]
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq[significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$interacting_pair %in% ResVsNonResChangesoutDiff$interactionpair,]
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq<-na.zero(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq)
#responders
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqResponderWilcoxon.pdf",height=10,width=17) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveRatio)))+
  geom_point(colour="black",size=5,shape=22)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                         axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                         axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                         axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                                    colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                         legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                              fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Change from Pre- to Post-treatment in Responders")+
  theme(axis.title = element_text(size = 25), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()

#non-responders
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqNonresponderWilcoxon.pdf",height=10,width=17) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq,aes(x=interacting_cell,y=interacting_pair,fill=log(NonEffectiveRatio)))+
  geom_point(colour="black",size=5,shape=22)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$NonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$NonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                             axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                             axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                             axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                        colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                             legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                  fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Change from Pre- to Post-treatment in Non-Responders")+
  theme(axis.title = element_text(size = 25), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()


#Relative differences
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqRelativeDifferenceWilcoxon.pdf",height=10,width=17) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",size=5,shape=22)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Relative Differences of Changes between Responders and Non-Responders")+
  theme(axis.title = element_text(size = 25), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()


#heatmap
outDiffreshape<-ResVsNonResChangesoutDiff[,c("interactionpair","NonResponderChangeMean","ResponderChangeMean")]
outDiffreshape<-outDiffreshape[outDiffreshape$interactionpair %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$interacting_pair,]
outDiffreshape<-melt(outDiffreshape,id="interactionpair")
outDiffreshape$value<-as.character(outDiffreshape$value)
outDiffreshape$value<-as.numeric(outDiffreshape$value)
outDiffreshape<-outDiffreshape[order(outDiffreshape$value),]
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeqWilcoxonHeatmap.pdf",height=4,width=8) 

ggplot(outDiffreshape,aes(x=variable,y=interactionpair,fill=log(value)))+
  geom_point(colour="black",size=8.4,shape=22)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(outDiffreshape$value)),0,max(log(outDiffreshape$value)))))+ 
  theme(axis.ticks = element_line(linetype = "blank"), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(family = "Times", 
        size = 9), axis.text = element_text(family = "Times", 
        size = 9, colour = "black"), axis.text.x = element_text(family = "Times", 
        size = 9, colour = "black", angle = 90), 
    plot.title = element_text(family = "Times", 
        size = 10), legend.text = element_text(size = 9, 
        family = "Times"), legend.title = element_text(size = 9, 
        family = "Times"), panel.background = element_rect(fill = NA), 
    legend.key = element_rect(fill = NA), 
    legend.background = element_rect(fill = NA)) +labs(x = "Responders vs. Non-responders", 
    y = "Cell-Cell Interaction Pairs", fill = "Log2 (Changes)")
dev.off()






#update
expr<-UpdateSeuratObject(object = expr)

##monocle3 Fig7 demo

data <- as(as.matrix(expr@assays$RNA@counts+0.001), 'sparseMatrix')
pd <-  expr@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
colnames(pd)

#responders
Idents(expr)<-"Treatment.Patients"
respondersSeurat<-SubsetData(expr,ident.use =c("su001","su002","su003","su004","su009","su012"))
data <- as(as.matrix(respondersSeurat@assays$RNA@counts+0.001), 'sparseMatrix')
pd <-  respondersSeurat@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
colnames(pd)

#non-responders
nonrespondersSeurat<-SubsetData(expr,ident.use =c("su005","su006","su007","su008","su010"))
data <- as(as.matrix(nonrespondersSeurat@assays$RNA@counts+0.001), 'sparseMatrix')
pd <-  nonrespondersSeurat@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
colnames(pd)

#Construct monocle cds
cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)

#Pre-process the data

cds = preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = c("UMAP"))
plot_cells(cds)
plot_cells(cds, color_cells_by="Cell.Type",label_cell_groups = TRUE)
#batch effect detection
plot_cells(cds, color_cells_by="Treatment.Patients", label_cell_groups=FALSE)


#cluster
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "Cell.Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
plot_cells(cds,
           color_cells_by = "Treatment.Status",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

get_earliest_principal_node <- function(cds, time_bin="pre"){

  cell_ids <- which(colData(cds)[, "Treatment.Status"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

ciliated_genes<-read.csv(file="ligand-receptor list.csv",header = T)
ciliated_genes<-ciliated_genes$Gene
ciliated_genes<-as.character(ciliated_genes)
plot_cells(cds,
           genes=c("SELPLG","WNT5A"),
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE)

cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
endothelialcds_subset<-cds_subset[,colData(cds_subset)$Cell.Type =="CAFs"]
gene_fits = fit_models(cds_subset, model_formula_str = "~Treatment.Status+Cell.Type")
fit_coefs = coefficient_table(gene_fits)

emb_time_terms = fit_coefs %>% filter(term != "(Intercept)")

emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

#responder change
responderTimeRelatedGenes<-emb_time_terms[,c("gene_short_name", "term", "q_value", "estimate","normalized_effect")]
responderTimeRelatedGenes<-responderTimeRelatedGenes[responderTimeRelatedGenes$q_value<0.05,]
write.csv(responderTimeRelatedGenes,file = "responderTimeRelatedLigandReceptorGenes.csv")

#nonresponder change
nonresponderTimeRelatedGenes<-emb_time_terms[,c("gene_short_name", "term", "q_value", "estimate","normalized_effect")]
nonresponderTimeRelatedGenes<-nonresponderTimeRelatedGenes[nonresponderTimeRelatedGenes$q_value<0.05,]
write.csv(nonresponderTimeRelatedGenes,file = "nonresponderTimeRelatedLigandReceptorGenes.csv")

plot_genes_violin(cds_subset, group_cells_by="Treatment.Status", ncol=4) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

#monocle3 using ligand-receptor genes 
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
cds_subset = preprocess_cds(cds_subset, num_dim = 100)

plot_pc_variance_explained(cds_subset)

cds_subset <- reduce_dimension(cds_subset, preprocess_method = "PCA",reduction_method = c("UMAP"))
plot_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="Cell.Type",label_cell_groups=TRUE,label_groups_by_cluster = TRUE,labels_per_group = 2,show_trajectory_graph = TRUE,
           cell_size = 0.8,group_label_size=4,graph_label_size=1.5)+scale_color_manual(values = allcolour)

#cluster
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset,
           color_cells_by = "Cell.Type",
           label_cell_groups = TRUE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,cell_size = 1.5,group_label_size=5,show_trajectory_graph = TRUE)
plot_cells(cds_subset,
           color_cells_by = "Treatment.Status",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=2)

get_earliest_principal_node <- function(cds_subset, time_bin="pre"){

  cell_ids <- which(colData(cds_subset)[, "Treatment.Status"] == time_bin)
  
  closest_vertex <-
    cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_subset), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds_subset)[["UMAP"]])$name[as.numeric(names
                                                                     (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds_subset = order_cells(cds_subset, root_pr_nodes=get_earliest_principal_node(cds_subset))
plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=2)



#Finding genes that change as a function of pseudotime
plot_cells(cds_subset,
           color_cells_by = "Cell.Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ciliated_cds_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

plot_cells(cds_subset, genes=c("SELE","SELPLG"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)


AFD_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% c("SELE"),
                              colData(cds_subset)$Cell.Type %in% c("Endothelial")]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0)

endothelial_cds <- cds[,grepl("Endothelial", colData(cds_subset)$Cell.Type, ignore.case=TRUE)]
plot_cells(endothelial_cds, color_cells_by="Treatment.Status")
pr_graph_test_res <- graph_test(endothelial_cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))





#CD8 ex
CD8cds_subset <- cds_subset[,grepl("CD8_ex", colData(cds_subset)$Cell.Type, ignore.case=TRUE)]
CD8cds_subset <- cds_subset[,filtercelltype]
#CD8
CD8cds_subset <- choose_cells(CD8cds_subset)
CD8cds_subset <- cluster_cells(CD8cds_subset)
CD8cds_subset <- learn_graph(CD8cds_subset)
plot_cells(CD8cds_subset,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
plot_cells(CD8cds_subset,
           color_cells_by = "Treatment.Status",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3,cell_size = 1.5)

get_earliest_principal_node <- function(CD8cds_subset, time_bin="pre"){

  cell_ids <- which(colData(CD8cds_subset)[, "Treatment.Status"] == time_bin)
  
  closest_vertex <-
    CD8cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(CD8cds_subset), ])
  root_pr_nodes <-
    igraph::V(principal_graph(CD8cds_subset)[["UMAP"]])$name[as.numeric(names
                                                                        (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

CD8cds_subset = order_cells(CD8cds_subset, root_pr_nodes=get_earliest_principal_node(CD8cds_subset))


CD8cds_subset <- order_cells(CD8cds_subset)

plot_cells(CD8cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,cell_size = 1.5)

plot_cells(CD8cds_subset,
           color_cells_by = "Cell.Type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,cell_size = 1.5)

plot_cells(CD8cds_subset,
           color_cells_by = "Treatment.Status",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,cell_size = 1.5,alpha = 0.5)+ scale_color_manual(values=c("orange","darkblue"))

plot_cells(CD8cds_subset,
           color_cells_by = "efficacy",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,cell_size = 1,alpha = 0.3)+ scale_color_manual(values=c("orange","darkblue"))

plot_cells(CD8cds_subset,
           color_cells_by = "exCluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,cell_size = 1.3,alpha = 0.3)




#co-regulated modules
subset_pr_test_res <- graph_test(CD8cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.01))
subset_pr_test_res<-subset_pr_test_res[pr_deg_ids,]
write.csv(subset_pr_test_res,file = "subset_pr_test_res.csv")

gene_module_df <- find_gene_modules(CD8cds_subset[pr_deg_ids,], resolution=0.2)
agg_mat <- aggregate_gene_expression(CD8cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

cell_group_df <- tibble::tibble(cell=row.names(colData(CD8cds_subset)), 
                                cell_group=partitions(CD8cds_subset))
agg_mat <- aggregate_gene_expression(CD8cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(CD8cds_subset,norm_method = "log",
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

write.csv(gene_module_df,file = "gene_module_df.csv")

#co-regulated modules pre vs post
subset_pr_test_res <- graph_test(CD8cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.01))
subset_pr_test_res<-subset_pr_test_res[pr_deg_ids,]
write.csv(subset_pr_test_res,file = "subset_pr_test_res.csv")

gene_module_df <- find_gene_modules(CD8cds_subset[pr_deg_ids,], resolution=0.01)
agg_mat <- aggregate_gene_expression(CD8cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

cell_group_df <- tibble::tibble(cell=row.names(colData(CD8cds_subset)), 
                                cell_group=colData(CD8cds_subset)$Treatment.Status)
agg_mat <- aggregate_gene_expression(CD8cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(CD8cds_subset,norm_method = "log",
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

write.csv(gene_module_df,file = "gene_module_dfprepost.csv")

#pre vs post 
gene_fits <- fit_models(CD8cds_subset, model_formula_str = "~Treatment.Status")
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "Treatment.Statuspre")
emb_time_terms <-emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(emb_time_terms,file = "emb_time_terms(prevspost).csv")

plot_cells(CD8cds_subset,genes = c("ADGRE5"),color_cells_by = 
             +scale_shape_manual(values = c(18, 16, 17))) 

ggplot(data=markers_exprs2,aes(x=markers_exprs2$UMAP1,y=markers_exprs2$UMAP2,size=log(value+1)))+
  geom_point(aes(shape = markers_exprs2$treatmentstatus,color=log(markers_exprs2$value)))+
  scale_color_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log(markers_exprs2$value+1)),0,max(log(markers_exprs2$value+1)))))+
  scale_shape_manual(values = c(16, 17))+
  
  plot_cells(CD8cds_subset,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,min_expr = 0)

cd8exprepostcds_subset = CD8cds_subset[rowData(CD8cds_subset)$gene_short_name %in% emb_time_terms$gene_short_name,]
plot_genes_violin(cd8exprepostcds_subset,group_cells_by="Treatment.Status", ncol=6) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


#cd8ex pre
CD8exprecds_subset <- CD8cds_subset[,grepl("pre", colData(CD8cds_subset)$Treatment.Status, ignore.case=TRUE)]
p1 <- plot_cells(CD8exprecds_subset,
           color_cells_by = "exCluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,cell_size = 1.3,alpha = 0.3)

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = 'black',
  fill = '#D4C8C8'
)


CD8expreClusterProportion<-cbind(CD8exprecds_subset@colData@listData[["exCluster"]],as.character(CD8exprecds_subset@colData@listData[["Treatment.Status"]]))
CD8expreClusterProportion<-as.data.frame(CD8expreClusterProportion)
colnames(CD8expreClusterProportion)[1]<-c("cluster")


#cd8ex post
CD8expostcds_subset <- CD8cds_subset[,grepl("post", colData(CD8cds_subset)$Treatment.Status, ignore.case=TRUE)]
p2 <- plot_cells(CD8expostcds_subset,
                 color_cells_by = "exCluster",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=TRUE,
                 graph_label_size=5,cell_size = 1.3,alpha = 0.3)

ggExtra::ggMarginal(
  p = p2,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = 'black',
  fill = '#D4C8C8'
)

CD8expostClusterProportion<-cbind(CD8expostcds_subset@colData@listData[["exCluster"]],as.character(CD8expostcds_subset@colData@listData[["Treatment.Status"]]))
CD8expostClusterProportion<-as.data.frame(CD8expostClusterProportion)
colnames(CD8expostClusterProportion)[1]<-c("cluster")

cd8exClusterPrePostProportion<-rbind(CD8expreClusterProportion,CD8expostClusterProportion)
ggplot(cd8exClusterPrePostProportion) +
  aes(x = V2, fill = cluster) +
  geom_bar() +
  scale_fill_hue() +
  theme_classic()

#cd8ex cluster0
CD8exCluster0cds_subset <- CD8cds_subset[,grepl("0", colData(CD8cds_subset)$exCluster, ignore.case=TRUE)]
gene_fitsExCluster0 <- fit_models(CD8exCluster0cds_subset, model_formula_str = "~Treatment.Status")
fit_coefsExCluster0 <- coefficient_table(gene_fitsExCluster0)
emb_time_termsExCluster0 <- fit_coefsExCluster0 %>% filter(term == "Treatment.Statuspre")
emb_time_termsExCluster0 <-emb_time_termsExCluster0 %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

CD8exCluster0cds_subsetMatrix<-SingleCellExperiment::counts(CD8exCluster0cds_subset)
CD8exCluster0cds_subsetMatrix_exprs = matrix(CD8exCluster0cds_subsetMatrix, nrow=nrow(CD8exCluster0cds_subset))
colnames(CD8exCluster0cds_subsetMatrix_exprs) = colnames(SingleCellExperiment::counts(CD8exCluster0cds_subset))
row.names(CD8exCluster0cds_subsetMatrix_exprs) = row.names(CD8exCluster0cds_subset)


#cd8ex cluster1
CD8exCluster1cds_subset <- CD8cds_subset[,grepl("1", colData(CD8cds_subset)$exCluster, ignore.case=TRUE)]
gene_fitsExCluster1 <- fit_models(CD8exCluster1cds_subset, model_formula_str = "~Treatment.Status")
fit_coefsExCluster1 <- coefficient_table(gene_fitsExCluster1)
emb_time_termsExCluster1 <- fit_coefsExCluster1 %>% filter(term == "Treatment.Statuspre")
emb_time_termsExCluster1 <-emb_time_termsExCluster1 %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

CD8exCluster1cds_subsetMatrix<-SingleCellExperiment::counts(CD8exCluster1cds_subset)
CD8exCluster1cds_subsetMatrix_exprs = matrix(CD8exCluster1cds_subsetMatrix, nrow=nrow(CD8exCluster1cds_subset))
colnames(CD8exCluster1cds_subsetMatrix_exprs) = colnames(SingleCellExperiment::counts(CD8exCluster1cds_subset))
row.names(CD8exCluster1cds_subsetMatrix_exprs) = row.names(CD8exCluster1cds_subset)

#cd8ex cluster2
CD8exCluster2cds_subset <- CD8cds_subset[,grepl("2", colData(CD8cds_subset)$exCluster, ignore.case=TRUE)]
gene_fitsExCluster2 <- fit_models(CD8exCluster2cds_subset, model_formula_str = "~Treatment.Status")
fit_coefsExCluster2 <- coefficient_table(gene_fitsExCluster2)
emb_time_termsExCluster2 <- fit_coefsExCluster2 %>% filter(term == "Treatment.Statuspre")
emb_time_termsExCluster2 <-emb_time_termsExCluster2 %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

CD8exCluster2cds_subsetMatrix<-SingleCellExperiment::counts(CD8exCluster2cds_subset)
CD8exCluster2cds_subsetMatrix_exprs = matrix(CD8exCluster2cds_subsetMatrix, nrow=nrow(CD8exCluster2cds_subset))
colnames(CD8exCluster2cds_subsetMatrix_exprs) = colnames(SingleCellExperiment::counts(CD8exCluster2cds_subset))
row.names(CD8exCluster2cds_subsetMatrix_exprs) = row.names(CD8exCluster2cds_subset)


#cd8ex cluster3
CD8exCluster3cds_subset <- CD8cds_subset[,grepl("3", colData(CD8cds_subset)$exCluster, ignore.case=TRUE)]
gene_fitsExCluster3 <- fit_models(CD8exCluster3cds_subset, model_formula_str = "~Treatment.Status")
fit_coefsExCluster3 <- coefficient_table(gene_fitsExCluster3)
emb_time_termsExCluster3 <- fit_coefsExCluster3 %>% filter(term == "Treatment.Statuspre")
emb_time_termsExCluster3 <-emb_time_termsExCluster3 %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

CD8exCluster3cds_subsetMatrix<-SingleCellExperiment::counts(CD8exCluster3cds_subset)
CD8exCluster3cds_subsetMatrix_exprs = matrix(CD8exCluster3cds_subsetMatrix, nrow=nrow(CD8exCluster3cds_subset))
colnames(CD8exCluster3cds_subsetMatrix_exprs) = colnames(SingleCellExperiment::counts(CD8exCluster3cds_subset))
row.names(CD8exCluster3cds_subsetMatrix_exprs) = row.names(CD8exCluster3cds_subset)

cd8exClusterMatrixForGSVA<-cbind(CD8exCluster0cds_subsetMatrix_exprs,CD8exCluster1cds_subsetMatrix_exprs,CD8exCluster2cds_subsetMatrix_exprs,CD8exCluster3cds_subsetMatrix_exprs)
cd8exClusterMatrixForGSVA<-as.data.frame(cd8exClusterMatrixForGSVA)
write.csv(cd8exClusterMatrixForGSVA,file = "cd8exClusterMatrixForGSVA.csv")

#SCT data extraction
data <- as(as.matrix(respondersSeurat@assays$SCT@counts+0.001), 'sparseMatrix')
pd <-  respondersSeurat@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

#ITGAE
ITGAEviolin<-read.csv(file = "CD8exRespondersClusters(ligand-receptor clustered).csv",header = T)
ITGAEviolin1<-colnames(expr)
ITGAEviolin1<-as.data.frame(ITGAEviolin1)
ITGAEviolin<-ITGAEviolin[,-1]
fix(ITGAEviolin)
fix(ITGAEviolin1)
ITGAEviolin1<-dplyr::left_join(ITGAEviolin1,ITGAEviolin,by="ID")
expr@meta.data[["exCluster"]]<-ITGAEviolin1$CD8exRespondersClusters
VlnPlot(expr,features = "ITGAE",group.by = "exCluster",pt.size = 0)
Idents(expr)<-"exCluster"
ITGAEseurat<-SubsetData(expr,ident.use = c("0","1","2","3"))
VlnPlot(ITGAEseurat,features = "ITGAE",group.by = "exCluster",pt.size = 0)
exClusterDEG<-FindAllMarkers(ITGAEseurat,min.pct = 0.1,logfc.threshold = 0.25)

VlnPlot(ITGAEseurat,features = c("ENTPD1","ITGAE","CD96","KLRC1","CSF1","CD52","HAVCR2","CD74"),
        pt.size = 0,group.by = "exCluster",stack = TRUE,fill.by = "ident")

VlnPlot(ITGAEseurat,features = c("GZMB","NKG7","GNLY","CCL5","TNFRSF9","FASLG","CLEC2B"),
        pt.size = 0,group.by = "exCluster",stack = TRUE,fill.by = "ident",flip =TRUE )

FeaturePlot(respondersSeurat,features = c("ENTPD1"),pt.size = 0.1)
FeaturePlot(respondersSeurat,features = c("ITGAE"),pt.size = 0.1)
FeaturePlot(nonrespondersSeurat,features = c("ITGAE"),pt.size = 0.1)
FeaturePlot(nonrespondersSeurat,features = c("ENTPD1"),pt.size = 0.001)

#ENTPD1 ITGAE co-expression
Idents(respondersSeurat)<-"Cell.Type"
CD8respondersSeurat<-SubsetData(respondersSeurat,ident.use = c("CD8_act_T_cells","CD8_ex_T_cells","CD8_mem_T_cells"))
coEntpd1ItgaeRes<-FeatureScatter(CD8respondersSeurat, feature1 = "ENTPD1", feature2 = "ITGAE",group.by = "Cell.Type",pt.size = 2)+scale_color_manual(values = c("#9370DB","#98FB98","#F08080"))
Idents(nonrespondersSeurat)<-"Cell.Type"
CD8nonrespondersSeurat<-SubsetData(nonrespondersSeurat,ident.use = c("CD8_act_T_cells","CD8_ex_T_cells","CD8_mem_T_cells"))
coEntpd1ItgaeNonRes<-FeatureScatter(CD8nonrespondersSeurat, feature1 = "ENTPD1", feature2 = "ITGAE",group.by = "Cell.Type",pt.size = 2)+scale_color_manual(values = c("#9370DB","#98FB98","#F08080"))
coEntpd1ItgaeRes+coEntpd1ItgaeNonRes
Idents(expr)<-"Cell.Type"
CD8Seurat<-SubsetData(expr,ident.use = c("CD8_act_T_cells","CD8_ex_T_cells","CD8_mem_T_cells"))
DotPlot(CD8nonrespondersSeurat, features = c("ENTPD1","ITGAE"), group.by = "Cell.Type",dot.scale = 20) + RotatedAxis()
DotPlot(CD8Seurat, features = c("ENTPD1","ITGAE"), group.by = "Cell.Type",dot.scale = 20) + RotatedAxis()

#Ligand-receptor markers Profiles Fig 1L demo
marker_test_res <- top_markers(cds_subset, group_cells_by="Cell.Type", 
                               reference_cells=2000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(4, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds_subset,
                    top_specific_marker_ids,
                    group_cells_by="Cell.Type",
                    ordering_type="maximal_on_diag",
                    flip_percentage_mean=TRUE,
                    max.size=6)

plot_genes_by_group <- function(cds_subset,
                                top_specific_marker_ids,
                                group_cells_by="Cell.Type",
                                reduction_method = "UMAP",
                                norm_method = c("log", "size_only"),
                                lower_threshold = 0,
                                max.size = 10,
                                # maybe be also do the maximum color on the
                                # diagonal; the axis change be switched too
                                ordering_type = c('cluster_row_col',
                                                  'maximal_on_diag',
                                                  'none'),
                                axis_order = c('group_marker', 'marker_group'),
                                flip_percentage_mean = FALSE,
                                pseudocount = 1,
                                scale_max = 3,
                                scale_min = -3) {
  
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  
  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", "partition") |
                              group_cells_by %in% names(colData(cds)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table."))
  }
  
  norm_method = match.arg(norm_method)
  
  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>%
    dplyr::pull(rowname)
  if(length(gene_ids) < 1)
    stop(paste('Please make sure markers are included in the gene_short_name",
               "column of the fData!'))
  
  if(flip_percentage_mean == FALSE){
    major_axis <- 1
    minor_axis <- 2
  } else if (flip_percentage_mean == TRUE){
    major_axis <- 2
    minor_axis <- 1
  }
  
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c('Cell', 'Gene', 'Expression')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  
  
  if (group_cells_by == "cluster"){
    cell_group <- tryCatch({clusters(cds,
                                     reduction_method = reduction_method)},
                           error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group <- tryCatch({partitions(cds,
                                       reduction_method = reduction_method)},
                           error = function(e) {NULL})
  } else{
    cell_group <- colData(cds)[,group_cells_by]
  }
  
  if (length(unique(cell_group)) < 2) {
    stop(paste("Only one type in group_cells_by. To use plot_genes_by_group,",
               "please specify a group with more than one type. "))
  }
  
  names(cell_group) = colnames(cds)
  
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
    dplyr::summarize(mean = mean(log(Expression + pseudocount)),
                     percentage = sum(Expression > lower_threshold) /
                       length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)
  
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']
  
  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene,
                         value.var = colnames(ExpVal)[2 + major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id
  
  if(ordering_type == 'cluster_row_col') {
    row_dist <- stats::as.dist((1 - stats::cor(t(res)))/2)
    row_dist[is.na(row_dist)] <- 1
    
    col_dist <- stats::as.dist((1 - stats::cor(res))/2)
    col_dist[is.na(col_dist)] <- 1
    
    ph <- pheatmap::pheatmap(res,
                             useRaster = T,
                             cluster_cols=TRUE,
                             cluster_rows=TRUE,
                             show_rownames=F,
                             show_colnames=F,
                             clustering_distance_cols=col_dist,
                             clustering_distance_rows=row_dist,
                             clustering_method = 'ward.D2',
                             silent=TRUE,
                             filename=NA)
    
    ExpVal$Gene <- factor(ExpVal$Gene,
                          levels = colnames(res)[ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group,
                           levels = row.names(res)[ph$tree_row$order])
    
  } else if(ordering_type == 'maximal_on_diag'){
    
    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for(i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]
    
    if(major_axis == 1){
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene ,
                            levels = dimnames(res)[[2]][max_ind_vec])
    }
    else{
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)),
                                            max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group,
                             levels = dimnames(res)[[1]][max_ind_vec])
    }
  } else if(ordering_type == 'none'){
    ExpVal$Gene <- factor(ExpVal$Gene, levels = markers)
  }
  
  if(flip_percentage_mean){
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) +
      geom_point(aes(colour = mean,  size = mean)) +
      viridis::scale_color_viridis(name = 'mean') +
      scale_size(name = 'log(mean + 0.1)', range = c(0, max.size))
  } else {
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) +
      geom_point(aes(colour = mean,  size = mean)) +
      viridis::scale_color_viridis(name = 'log(mean + 0.1)') +
      scale_size(name = 'mean', range = c(0, max.size))
  }
  
  if (group_cells_by == "cluster"){
    g <- g + xlab("Cluster")
  } else if (group_cells_by == "partition") {
    g <- g + xlab("Partition")
  } else{
    g <- g + xlab(group_cells_by)
  }
  
  g <- g + ylab("Gene") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  if(axis_order == 'marker_group') {
    g <- g + coord_flip()
  }
  
  g
}



#celltype frequency cd8
celltypefrequencyforplot<-read.csv("celltypefrequencyCD8.csv",header = T)

CairoPDF(file="celltypefrequencyplot.pdf",width=13,height=10) 
showtext_begin()
ggplot(data=celltypefrequencyforplot,aes(status,counts,fill=cluster))+
  geom_bar(stat="identity",position="fill",color="white",width = 1,size=0.25)+
  scale_fill_manual(values = c("#9370DB","#98FB98","#F08080")) + theme(axis.line = element_line(linetype = "solid"), 
                                                axis.ticks = element_line(size = 1), 
                                                panel.grid.major = element_line(colour = "black", 
                                                                                linetype = "blank"), panel.grid.minor = element_line(colour = "black", 
                                                                                                                                     linetype = "blank"), axis.title = element_text(family = "Times"), 
                                                axis.text = element_text(family = "Times", 
                                                                         size = 11, colour = "black"), axis.text.x = element_text(family = "Times", 
                                                                                                                                  size = 11, colour = "black"), axis.text.y = element_text(family = "Times", 
                                                                                                                                                                                           size = 11, colour = "black"), plot.title = element_text(family = "Times"), 
                                                legend.text = element_text(family = "Times"), 
                                                legend.title = element_text(family = "Times"), 
                                                panel.background = element_rect(fill = NA), 
                                                legend.background = element_rect(fill = NA)) + theme(axis.title = element_text(size = 13)) +labs(x = "Treatment and Response Status", 
                                                                                                                                                 y = "Proportion", fill = "Cell Type")
showtext_end()
dev.off() 


#co-expression network
cds_exprs<-SingleCellExperiment::counts(CD8cds_subset)
markers_exprs = matrix(cds_exprs, nrow=nrow(CD8cds_subset))
colnames(markers_exprs) = colnames(SingleCellExperiment::counts(CD8cds_subset))
row.names(markers_exprs) = row.names(CD8cds_subset)
TF_matrix<-markers_exprs
TF_matrix <- CreateSeuratObject(counts = markers_exprs,project = "TF_matrix")
TF_matrix <- PercentageFeatureSet(TF_matrix, pattern = "^MT-", col.name = "percent.mt")
TF_matrix <- SCTransform(TF_matrix, vars.to.regress = "percent.mt", verbose = FALSE)
TF_matrix <- RunPCA(TF_matrix, verbose = FALSE)
TF_matrix <- RunUMAP(TF_matrix, dims = 1:30, verbose = FALSE)

TF_matrix <- FindNeighbors(TF_matrix, dims = 1:30, verbose = FALSE)
TF_matrix <- FindClusters(TF_matrix, verbose = FALSE,resolution = 0.3)
DimPlot(TF_matrix, label = TRUE) + NoLegend()

res <- FindAllMarkers(TF_matrix, only.pos = T, test.use = "wilcox", min.diff.pct = 0, print.bar = F, do.print = T)
CD8exRespondersClusters<-TF_matrix@meta.data[["seurat_clusters"]]
CD8exRespondersClusters<-as.character(CD8exRespondersClusters)
CD8exRespondersClustersId<-colnames(TF_matrix)
CD8exRespondersClusters<-cbind(CD8exRespondersClustersId,CD8exRespondersClusters)
CD8exRespondersClusters<-as.data.frame(CD8exRespondersClusters)
write.csv(CD8exRespondersClusters,file = "CD8exRespondersClusters(ligand-receptor clustered).csv")
write.csv(res, "TF_markers.csv")

tf.names<-ciliated_genes

aux <- dcast(res, gene ~ cluster, value.var = "p_val")
rownames(aux) <- aux[,1]
aux <- aux[,-1]
DF  <-  as.matrix(as.data.frame(lapply(aux, as.numeric)))
DF <- -log(DF,10)
DF[is.na(DF)] <- 0
dd <- dist(DF, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
heatmap.2(DF)
aux <- t(aux)
colnames(aux) <- aux[1,]
aux <- aux[-1,]

mat<-TF_matrix@assays[["RNA"]]@counts
groups<-TF_matrix@meta.data[["seurat_clusters"]]
group_averages <- function(mat, groups){
  group_names = unique(groups)
  means = matrix(0, dim(mat)[1], length(group_names))
  colnames(means) = group_names
  rownames(means) = rownames(mat)
  for(group in group_names){
    means[,group] = Matrix::rowMeans(mat[,groups == group,drop=FALSE])
  }
  means
}

output.tf.cor <- function(data, pval_cutoff, avg_logFC_cutoff, nmarkers.per.celltype, markerfile="TF_markers.csv",height=10){
  
  # read in markers
  markers <- read.csv(markerfile, row.names = 1)
  markers$gene <- rownames(markers)
  markers$avg_logFC <- as.numeric(as.character(markers$avg_logFC))
  markers$p_val <- as.numeric(as.character(markers$p_val))
  
  print(head(markers))
  
  # Get top genes by cell type
  genes.use <- as.character(unique(markers %>% dplyr::group_by(cluster) %>% 
                                     filter(gene %in% tf.names) %>% 
                                     dplyr::top_n(n=nmarkers.per.celltype, wt=avg_logFC) %>% 
                                     pull(gene)))
  genes.use<-ciliated_genes
  print(paste("genes plotted: ",length(genes.use)))
  # Calculate correlations
  # mat.use <- t(as.matrix(tiss@data[genes.use, ]))
  
  # use cell type averages instead (changed for revision)
  annotation.2.means <- group_averages(TF_matrix@assays[["RNA"]]@counts, TF_matrix@meta.data[["seurat_clusters"]])
  mat.use <- t(annotation.2.means[genes.use, ])
  tf.cor    <- cor(mat.use)
  
  return(tf.cor)
}
# cluster genes based on expression correlation and plot correlogram
# returns output of barb.cormap, which includes the list of genes in the order they appear in the correlogram
plot.correlogram <- function(data, pval_cutoff, avg_logFC_cutoff, nmarkers.per.celltype, markerfile="TF_markers.csv",height=10){
  
  # read in markers
  markers <- read.csv(markerfile, row.names = 1)
  markers$gene <- rownames(markers)
  markers$avg_logFC <- as.numeric(as.character(markers$avg_logFC))
  markers$p_val <- as.numeric(as.character(markers$p_val))
  
  print(head(markers))
  
  # Get top genes by cell type
  genes.use <- as.character(unique(markers %>% dplyr::group_by(cluster) %>% 
                                     filter(gene %in% tf.names) %>% 
                                     dplyr::top_n(n=nmarkers.per.celltype, wt=avg_logFC) %>% 
                                     pull(gene)))
  print(paste("genes plotted: ",length(genes.use)))
  # Calculate correlations
  # mat.use <- t(as.matrix(tiss@data[genes.use, ]))
  
  # use cell type averages instead (changed for revision)
  annotation.2.means <- group_averages(data@assays[["RNA"]]@counts, data@meta.data[["seurat_clusters"]])
  mat.use <- t(annotation.2.means[genes.use, ])
  tf.cor    <- cor(mat.use)
  print(colnames(tf.cor))
  enrich.score <- dcast(markers, gene ~ cluster, value.var = 'avg_logFC')
  rownames(enrich.score) <- enrich.score[,"gene"]
  enrich.score <- enrich.score[, 2:ncol(enrich.score)]
  # change rownames to the cell type and tissue that each gene is enriched in (ordering by avg_logFC above)
  topIDenriched   <- sapply(colnames(tf.cor), function(x) {
    names(sort(t(enrich.score)[, x],  decreasing = T))[1]}) 
  colnames(tf.cor) <- topIDenriched
  
  correlo.out <- barb.cormap(tf.cor,  'TF_cormap1.pdf', height =height, width = height)
  
  # Generate row colors of heatmap corresponding to highest-expressing cell type
  topID.plotorder <- correlo.out[[2]]
  meta.summary <- data@meta.data %>% distinct(seurat_clusters)
  print(head(meta.summary))
  ntypes=length(unique(meta.summary$cell_ontology_class))
  tmp <- colorRampPalette(brewer.pal(min(ntypes, 11), 'Paired'))(ntypes)
  annot_colors <- data.frame(annot.colors = tmp, cell_ontology_class = unique(meta.summary$cell_ontology_class))
  write.csv(annot_colors, file=paste0('annotColors.csv'))
  
  print(head(annot_colors))
  meta.summary <- merge(meta.summary, annot_colors, by = 'cell_ontology_class')
  print(head(meta.summary))
  plot.colors <- data.frame(annotation.2 = topID.plotorder)
  print(head(plot.colors))
  plot.colors <- merge(plot.colors, meta.summary %>% select(annotation.2, annot.colors, tiss.color), by = 'annotation.2')
  print((plot.colors))
  plot.colors <- plot.colors[match(topID.plotorder, plot.colors$annotation.2), ]
  print((plot.colors))
  plot.colors$ymin <- 0.1*(0:(nrow(plot.colors)-1))
  plot.colors$ymax <- 0.1*(1:(nrow(plot.colors)))
  print((plot.colors))
  plot.colors$tiss.color <- toupper(plot.colors$tiss.color)
  plot.colors$annot.colors <- toupper(plot.colors$annot.colors)
  require(grDevices)
  pdf( 'TF_cormap_rowcolors_AnnotRight_TissueLeft.pdf', height = 10, width = 4)
  plot(c(0, 2), c(0, max(plot.colors$ymax) + 1), type = "n", xlab = "", ylab = "",
       main = "plot colors")
  rect(0,plot.colors$ymin, 1 , plot.colors$ymax, col = plot.colors$annot.colors, border = NA)
  rect(1,plot.colors$ymin, 2 , plot.colors$ymax, col = plot.colors$tiss.color, border = NA)
  dev.off()
  return(correlo.out)
}
barb.cormap <- function(mat.cor,fname,width=12,height=12,method="complete",cex=0.5,mincor=-1,maxcor=1){
  require(lattice)
  require(cba)
  rowdist <- dist(mat.cor)
  coldist <- dist(mat.cor, by_rows = F)
  hc.cor <- hclust(coldist, method=method)
  hr.cor <- hclust(rowdist, method=method)
  optimal.row <- order.optimal(rowdist,hr.cor$merge)
  optimal.col <- order.optimal(coldist,hc.cor$merge)
  
  ord.row <- optimal.row$order
  ord.col <- optimal.col$order
  
  plt = levelplot(mat.cor[ord.row,ord.col],xlab=NULL,ylab=NULL,
                  at=do.breaks(c(mincor-0.01,maxcor+0.01),19),scales=list(x=list(rot=90),cex=cex),
                  colorkey=list(space="top"),
                  col.regions=colorRampPalette(c("dodgerblue4", "dodgerblue", "white", "lightcoral", "firebrick4"), space="Lab"))
  pdf(fname,width=width,height=height)
  print(plt)
  dev.off()
  
  return(list(rownames(mat.cor[ord.row, ]), colnames(mat.cor[, ord.col]), plt, hc.cor, hr.cor))
}



get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


tf.cor.all <- output.tf.cor(data = TF_matrix, pval_cutoff = 10^-3.5, avg_logFC_cutoff = 0.15, nmarkers.per.celltype = 10, height = 30)
reorder_tf.cor.all <- reorder_cormat(tf.cor.all)
upper_tf.cor.all <- get_upper_tri(reorder_tf.cor.all)

melted_tf.cor.all <- melt(upper_tf.cor.all)
# ggplot(data = melted_tf.cor.all, aes(x=Var1, y=Var2, fill=value)) + 
# geom_tile()
ggplot(data = melted_tf.cor.all, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "blue", mid = "white", na.value="white",
                       midpoint = 0,  limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 0, hjust = 1))+
  coord_fixed()+ 
  
  theme(text=element_text(size=4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = rel(1), angle = 0),
        legend.position="none")
ggsave("tf_cor_all.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 20, height = 20, units = "cm",
       dpi = 150, limitsize = FALSE)


correlo.all <- plot.correlogram(data = TF_matrix, pval_cutoff = 10^-3.5, avg_logFC_cutoff = 0.15, nmarkers.per.celltype = 20, height = 30)



#exCluster1 pre vs post
Idents(expr)<-"Cell.Type"
exCluster1<-SubsetData(expr,ident.use = "CD8_ex_T_cells")
Idents(exCluster1)<-"exCluster"
exCluster1<-SubsetData(exCluster1,ident.use = "1")
Idents(exCluster1)<-"Treatment.Status"
exCluster1PrePostDEG<-FindAllMarkers(exCluster1,only.pos = TRUE)
exCluster1PrePostDELigandReceptors<-FindAllMarkers(exCluster1,only.pos = TRUE,features = ciliated_genes)
write.csv(exCluster1PrePostDEG,file = "exCluster1PrePostDEG.csv")
#top10
top10exclusterPrePostDEG <- exCluster1PrePostDEG %>%
  filter(avg_logFC >= 0.84)
top10exclusterPrePostDEG<-top10exclusterPrePostDEG[top10exclusterPrePostDEG$cluster=="post",]
top10exclusterPrePostDEGname<-rownames(top10exclusterPrePostDEG)
top10exclusterPrePostDEGname<-c(top10exclusterPrePostDEGname,"BST2","IL16","IFNG","KLRC1","ADGRE5","CSF1")
Idents(exCluster1)<-"Treatment.Status"
VlnPlot(exCluster1,features = top10exclusterPrePostDEGname,pt.size = 0,group.by = "Treatment.Status",stack = TRUE,fill.by = "ident")




#Fig. S1 demo
PretreatResponderVsNonrespondersCombineAll<-read.csv(file = "PretreatResponderVsNonrespondersCombineAll.csv",header = TRUE)
ggplot(PretreatResponderVsNonrespondersCombineAll,aes(x=interacting_cell,y=interacting_pair,fill=LogRatioRespondersVsNonresponders))+
  geom_point(colour="black",shape=22,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(PretreatResponderVsNonrespondersCombineAll$LogRatioRespondersVsNonresponders),0,max(PretreatResponderVsNonrespondersCombineAll$LogRatioRespondersVsNonresponders)))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                        fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))

Fig. S2 demo
#significant_meansReshapePreEffectiveVsNonEffective
significant_meansReshapePreEffectiveVsNonEffective<-read.csv("significant_meansReshapePreEffectiveVsNonEffective.csv",header = T)
significant_meansReshapePreEffectiveVsNonEffective$interacting_cell<-gsub("\\.","-",significant_meansReshapePreEffectiveVsNonEffective$interacting_cell)
#all
CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffective.pdf",height=60,width=65) 

ggplot(significant_meansReshapePreEffectiveVsNonEffective,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                       axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                       axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                       axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                  colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                       legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                            fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))

dev.off()

#Obvious
significant_meansReshapePreEffectiveVsNonEffective<-significant_meansReshapePreEffectiveVsNonEffective[significant_meansReshapePreEffectiveVsNonEffective$FoldChange==">2 or <0.5",]
CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffectiveObvious.pdf",height=35,width=35) 

ggplot(significant_meansReshapePreEffectiveVsNonEffective,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values = 23)+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffective$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 30)) +labs(title = "Comparison of Cell-Cell Crosstalk between Responsers and Non-Responsers (Pre-Treatment Status)")+
  theme(axis.title = element_text(size = 30), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20), 
    legend.key = element_rect(size = 5), 
    legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
    vjust = 2))
dev.off()

#significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq
significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq<-significant_meansReshapePreEffectiveVsNonEffective[(significant_meansReshapePreEffectiveVsNonEffective$gene_aName %in% significant_meansReshapePreEffectiveVsNonEffective$anti.PD1RNASEQDEG.preEffectiveVsNonEffective.GSE78220)|(significant_meansReshapePreEffectiveVsNonEffective$gene_bName %in% significant_meansReshapePreEffectiveVsNonEffective$anti.PD1RNASEQDEG.preEffectiveVsNonEffective.GSE78220),]

CairoPDF(file="significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq.pdf",height=20,width=25) 

ggplot(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq,aes(x=interacting_cell,y=interacting_pair,fill=log(Ratio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$Ratio)),0,max(log(significant_meansReshapePreEffectiveVsNonEffectiveIntersectAntipd1Rnaseq$Ratio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 25)) +labs(title = "Intersecting Ligand-Receptor Pairs with Hugo et al. Database")+
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()


#Fig. S3 demo
PretreatVsPosttreatResponderVsNonrespondersCombineAll<-read.csv(file = "significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.csv",header = TRUE)
ggplot(PretreatVsPosttreatResponderVsNonrespondersCombineAll,aes(x=interacting_cell,y=interacting_pair,fill=log2(EffectiveVsNonEffectiveRatio)))+
  geom_point(colour="black",shape=22,size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),na.value = NA,values = rescale(c(min(log2(PretreatVsPosttreatResponderVsNonrespondersCombineAll$EffectiveVsNonEffectiveRatio)),0,max(log2(PretreatVsPosttreatResponderVsNonrespondersCombineAll$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                 axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                 axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                 axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                            colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                 legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                      fill = "Log2(Ratio of Score)") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 15,family = "Times"), 
        axis.text = element_text(size = 13,family = "Times"), legend.text = element_text(size = 10,family = "Times"), panel.background = element_rect(size = 0), plot.background = element_rect(size = 0.2))+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 1, vjust = 0)) +labs(title = "Change of Cell-Cell Crosstalk in Responders (Post Vs. Pre)")+ theme(plot.title = element_text(hjust = +0.5, 
                                                                                                                                                                      vjust = 2))



                                                                                                          y = "Proportion", fill = "Cell Type")




Fig S4 demo
#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll
#Ligand Scene(all)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll<-read.csv("significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.csv",header = T)
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell<-gsub("\\.","-",significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$interacting_cell)
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll.pdf",height=60,width=60) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio),shape=FoldChange))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  coord_equal()+
  scale_shape_manual(values=c(23,22))+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 20,family = "Times"),legend.title = element_text(size = 20,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))

dev.off()

#significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq
significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq<-significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll[(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_aName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061|significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$gene_bName %in% significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreAll$IntersectingAnti.PD1RNASEQDEG.EffectiveVsNonEffectiveofPostVsPre.GSE91061),]
CairoPDF(file="significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq2.pdf",height=25,width=60) 

ggplot(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq,aes(x=interacting_cell,y=interacting_pair,fill=log(EffectiveVsNonEffectiveRatio),shape=FoldChange.EffectiveCombineNonEffective.))+
  geom_point(colour="black",size=5)+
  scale_x_discrete()+
  scale_size(limits = 0.5)+
  scale_shape_manual(values=c(23,22))+
  coord_equal()+
  scale_fill_gradientn(colors=c("#2e68b7","#ffffff","#ce0000"),values=rescale(c(min(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio)),0,max(log(significant_meansReshapeEffectiveVsNonEffectiveOfPostVsPreALLIntersectingAntiPd1RnaSeq$EffectiveVsNonEffectiveRatio))))) + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                   axis.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                   axis.text = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                   axis.text.x = element_text(family = "Times", 
                                                                                                                                                                                                                                                                                                                                                                              colour = "black", angle = 90), plot.title = element_text(family = "Times"), 
                                                                                                                                                                                                                                                                                                                                                   legend.title = element_text(family = "Times")) +labs(x = "Ligand-Receptor Cell Set", y = "Ligand-Receptor", 
                                                                                                                                                                                                                                                                                                                                                                                                        fill = "Log2 Ratio of Change") + theme(axis.text = element_text(colour = "black")) + theme(axis.text = element_text(size = 4.5)) + theme(axis.text.x = element_text(vjust = -0.0005))+
  theme(panel.grid.major = element_line(size = 0.3), 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 11), legend.text = element_text(size = 17,family = "Times"),legend.title = element_text(size = 17,family = "Times"), panel.background = element_rect(size = 0), 
        plot.background = element_rect(size = 0.2))+ 
  theme(plot.title = element_text(size = 40)) +labs(title = "Change Comparison between Responsers and Non-Responsers Intersecting Riaz et al. Database")+
  theme(axis.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.key = element_rect(size = 5), 
        legend.background = element_rect(size = 5)) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                      vjust = 2))

dev.off()



#Featureplot Myofibroblasts and CAFs
CairoPDF(file="Featureplot Myofibroblasts and CAFs.pdf",height=8,width=24) 
showtext_begin()
FeaturePlot(expr,features=c("COL1A2","ACTA2","FAP"),ncol = 3 ,pt.size = 0.05)
showtext_end()
dev.off()

#Fig S5 GSVA CAF vs MYO demo 
CAFvsMYOgsvaFilter<-read.csv("CAFvsMYOgsvaFilter1(CAFsVsMyo).csv",header = T)
order<-sort(CAFvsMYOgsvaFilter$logFC,index.return=T,decreasing = T)
CAFvsMYOgsvaFilter$Geneset<-factor(CAFvsMYOgsvaFilter$Geneset,levels=CAFvsMYOgsvaFilter$Geneset[order$ix])
CairoPDF(file="CAFvsMYOgsvaFilter(CAFsVsMyo)2.pdf",width=20,height=35) 
ggplot(CAFvsMYOgsvaFilter, aes(x = Geneset, y = logFC,fill=Fill)) +
  geom_bar(stat="identity",width = 0.8,size=0.1,alpha=1) +
  coord_fixed() +
  coord_flip()+
  
  theme_light()  + theme(axis.title = element_text(family = "Times"), 
    axis.text = element_text(family = "Times", 
        size = 18, colour = "black"), axis.text.x = element_text(family = "Times",size = 18), 
    axis.text.y = element_text(family = "Times",size = 18), 
    plot.title = element_text(family = "Times", 
        size = 15), legend.text = element_text(size = 10, 
        family = "Times"), legend.title = element_text(size = 12, 
        family = "Times")) +labs(title = "GSVA between CAFs and Myofibroblasts") + theme(panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(size = 13), 
    axis.text = element_text(size = 18), 
    plot.title = element_text(size = 20, 
        hjust = -4)) + theme(axis.title = element_text(size = 18), 
    axis.text.x = element_text(colour = "black"), 
    axis.text.y = element_text(colour = "black"), 
    plot.title = element_text(size = 30), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20)) + theme(legend.position = c(0.8, 0.98)) +labs(fill = "Cell Type", size = 30) + theme(axis.title = element_text(size = 30)) +labs(x = "GO Terms")+ 
  theme(axis.text = element_text(size = 50, 
    face = "bold")) + 
  theme(axis.text = element_text(size = 20), 
    axis.text.y = element_text(size = 23)) + theme(legend.text = element_text(colour = NA), 
    legend.title = element_text(colour = NA), 
    legend.position = "none")
dev.off() 


library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(pheatmap)
DefaultAssay(object = expr) <- "RNA"
Idents(object = expr) <- 'Cell.Type'
myo<-SubsetData(expr,ident.use = "Myofibroblasts")
CAF<-SubsetData(expr,ident.use = "CAFs")
myomatrix<-myo[["RNA"]]@counts
myomatrix<-as.matrix(myomatrix)
CAFmatrix<-CAF[["RNA"]]@counts
CAFmatrix<-as.matrix(CAFmatrix)
myoCAFmatrix<-cbind(myomatrix,CAFmatrix)



type_1<-rep("Myofibroblasts",378)
type_2<-rep("CAFs",1066)
Type=c(type_1,type_2)

adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

design <- model.matrix(~ factor(Type))
colnames(design) <- c("Myofibroblasts", "CAFs vs. Myofibroblasts")
design[,2]<-ordercount

#DEGSET
geneSets <- getGmt("KEGG gene sets.gmt")

res_es <- gsva(myoCAFmatrix, geneSets,kcdf="Poisson", min.sz=1, max.sz=500, verbose=TRUE, parallel.sz=1)

fit <- lmFit(res_es, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="CAFs vs. Myofibroblasts", number=Inf)
DEgeneSets <- topTable(fit, coef="CAFs vs. Myofibroblasts", number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)

DEgeneSetsG<-DEgeneSets


#DEG
myoCAFmatrix<-voom(myoCAFmatrix,design,plot=T)
fit <- lmFit(myoCAFmatrix, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="CAFs vs. Myofibroblasts", number=Inf)
DEgenes <- topTable(fit, coef="CAFs vs. Myofibroblasts", number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)
write.csv(DEgenes,file = "CAFsMyoDEG.csv")



myoCAFmatrix <- apply(myoCAFmatrix, 2, function(x) (x/sum(x))*10000)
myoCAFmatrix <- t(myoCAFmatrix)
violinMYOCAFtype<-read.csv("violinMYOCAFtype.csv",header = T)
myoCAFmatrix<-cbind(violinMYOCAFtype,myoCAFmatrix)
myoCAFmatrix<-myoCAFmatrix[,c("CellType","CXCL12","CXCL14","IL6","S100A4","SFRP2","MMP1","MMP2","MMP3","MMP9","MMP14","WNT5A","WNT2")]
myoCAFmatrix<-cbind(rownames(myoCAFmatrix),myoCAFmatrix)
colnames(myoCAFmatrix)[1]<-"id"
myoCAFmatrix<-melt(myoCAFmatrix,id=(c("id","CellType")))
myoCAFmatrix<-myoCAFmatrix[,-1]
myoCAFmatrix$value1<-log2(myoCAFmatrix$value+1)


CairoPDF(file="violinMYOCAFDEG.pdf",height=4,width=10)
ggplot(myoCAFmatrix,aes(x=variable,y=value1))+
         geom_boxplot(outlier.size = 0,aes(fill=factor(CellType)),position = position_dodge(0.7),size=0.02,width=0.6)+
         guides(fill=guide_legend(title="Cell Type"))+
         theme_light() + theme(axis.line = element_line(linetype = "solid"), 
    axis.ticks = element_line(colour = "black", 
        size = 1), panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(family = "Times", 
        size = 12, face = "bold"), axis.text = element_text(family = "Times", 
        size = 10, face = "bold", colour = "black"), 
    plot.title = element_text(family = "Times", 
        size = 13, face = "bold", hjust = 0.5), 
    legend.text = element_text(size = 10, 
        face = "bold", family = "Times"), 
    legend.title = element_text(size = 10, 
        face = "bold", family = "Times"), 
    panel.background = element_rect(fill = NA), 
    plot.background = element_rect(size = 0)) +labs(title = "Expression Level of Genes between CAFs and Myofibroblasts", 
    x = "Gene", y = "Log2(Expression Level)", 
    fill = "Cell Type")
dev.off()


CairoPDF(file="violinMYOCAFDEG.pdf",height=2.5,width=9.5)
ggplot(myoCAFmatrix,aes(variable,value1,fill=CellType))+
  geom_bar(stat="identity",position=position_dodge(),width = 0.7,size=0.25) + theme(axis.line = element_line(linetype = "solid"), 
    axis.ticks = element_line(colour = "black"), 
    panel.grid.major = element_line(colour = "black", 
        linetype = "blank"), panel.grid.minor = element_line(colour = "black", 
        linetype = "blank"), axis.title = element_text(family = "Times", 
        size = 10, face = "bold"), axis.text = element_text(family = "Times", 
        size = 9, face = "bold", colour = "black"), 
    plot.title = element_text(family = "Times", 
        size = 14, face = "bold", hjust = 0.5), 
    legend.text = element_text(size = 9, 
        face = "bold", family = "Times"), 
    legend.title = element_text(size = 9, 
        face = "bold", family = "Times"), 
    panel.background = element_rect(fill = NA), 
    legend.key = element_rect(fill = NA), 
    legend.background = element_rect(fill = NA)) +labs(title = "Expression Level of Genes between CAFs and Myofibroblasts", 
    x = "Genes", y = "Log2(Expression Level)", 
    fill = "Cell Type") + theme(legend.position = c(0.9, 0.85))
dev.off()


#Chi-square test demo
install.packages("reshape")
library(reshape)
celltypefrequencyTcellsNK <- read_csv("celltypefrequencyTcellsNK.csv")
celltypefrequencyTcellsNK<-celltypefrequencyTcellsNK[,-1]

ab<-celltypefrequencyTcellsNK[c(1,2),]
ac<-celltypefrequencyTcellsNK[c(1,3),]
cd<-celltypefrequencyTcellsNK[c(3,4),]
bd<-celltypefrequencyTcellsNK[c(2,4),]

abt<-chisq.test(ab)
act<-chisq.test(ac)
cdt<-chisq.test(cd)
bdt<-chisq.test(bd)
