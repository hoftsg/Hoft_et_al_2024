setwd("/Users/macbook/Documents/Stella/Human_AIG_vs_Hp")
options(future.globals.maxSize = 100000 * 1024^2)
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(metap)
library(gridExtra)

#Identifying Clusters by Gene Expression
combined<-readRDS('human_hp_aig_v5.rds')
All_markers<-FindAllMarkers(combined, only.pos=TRUE, min.pct=0.25, logfc.threshold = 1)
write.csv(All_markers, file="./CSV/findallmarkers_v5.csv")
write.csv(table(combined@meta.data$seurat_clusters,combined@meta.data$condition), file="./CSV/number_by_path_v5.csv")
combined<-RenameIdents(combined,'0'="B",'1'="B",'2'="B",'3'="Neck",'4'="T",
                       '5'="B",'6'="Pit",'7'="T",'8'="SPEM",'9'="Pit",'10'="Cancer",
                       '11'="Prolif Epi",'12'="Endo",'13'="B",'14'="Fibro",'15'="Neuro",
                       '16'="Neck",'17'="Fibro",'18'="T",'19'="DC",'20'="Chief",
                       '21'="SM",'22'="Neuro",'23'="B",'24'="Mast",'25'="Prolif",
                       '26'="T",'27'="SM",'28'="Parietal",'29'="B",'30'="IM",
                       '31'="Neuro",'32'="Neuro",'33'="Neuro",'34'="B",'35'="B")
combined$epi_clusters<-Idents(combined)
DimPlot(combined, reduction="umap", label=FALSE)


#Applying Metaplasia signature
meta_sig <- list(c("TFF2","AQP5","CD44","MUC6","TFF3","DMBT1","CLU","CFTR","OLFM4","WFDC2","CDX1","CDX2","MUC4","MUC2","MAL2"))
combined <- AddModuleScore(object = combined, features = meta_sig, name = "meta_sig_score")
FeaturePlot(object = combined, features = "meta_sig_score1", pt.size=0.4, label=TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PRGn")))
#Seperate Metaplasia sig by disease
Idents(combined)<-"condition"
AIG <- subset(combined, ident = "AIG")
Idents(AIG)<-"epi_clusters"
Hp <- subset(combined, ident = "Hp")
Idents(Hp)<-"epi_clusters"
FeaturePlot(object = AIG, features = "meta_sig_score1", pt.size=0.4, label=TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PRGn")))
FeaturePlot(object = Hp, features = "meta_sig_score1", pt.size=0.4, label=TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PRGn")))

#Metaplasia subclustering
quantile(combined@meta.data[["meta_sig_score1"]], probs = c(0.9, 0.95))
# 90% = 0.2394622 
Meta<-subset(combined, subset = meta_sig_score1 > 0.2394622 & PTPRC <=0 & CD79A<=0 & CD3E<=0 & EPCAM>0)
DimPlot(object = Meta, reduction = "umap",label=TRUE, pt.size = 1, group.by = "epi_clusters", shuffle=FALSE) + NoLegend()
write.csv(table(Meta@meta.data$epi_clusters,Meta@meta.data$condition), file="./CSV/number_by_path_meta_subset.csv")


#proportion Metaplasia of total
Idents(Meta)<-"condition"
head(table(Idents(Meta))) 
Idents(combined)<-"condition"
head(table(Idents(combined))) 
#generating Metaplasia subset object
Meta[["RNA"]] <- split(Meta[["RNA"]], f = Meta$condition)
Meta
Meta<- NormalizeData(Meta)
Meta<- FindVariableFeatures(Meta)
Meta<- ScaleData(Meta)
Meta<- RunPCA(Meta)
integrated.meta <- IntegrateLayers(
  object = Meta, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
integrated.meta[["RNA"]] <- JoinLayers(integrated.meta[["RNA"]])
Meta<-FindNeighbors(integrated.meta, reduction = "integrated.rpca", dims=1:30)
Meta<-FindClusters(Meta, resolution = 0.4)
Meta<-RunUMAP(Meta, dims=1:30, reduction = "integrated.rpca")
DimPlot(object = Meta, reduction = "umap",label=TRUE, pt.size = 1)
DimPlot(object = Meta, reduction = "umap",label=TRUE, pt.size = 1, split.by="condition") + NoLegend()
saveRDS(Meta, file="Meta_aig+hp_v5_int.rds")


#DEG and proportions of metaplasia subtypes
Meta<-readRDS('Meta_aig+hp_v5_int.rds')
Meta$celltype<- Idents(Meta)
write.csv(table(Meta@meta.data$condition,Meta@meta.data$seurat_clusters), file="./CSV/number_meta_by_path_int.csv")
All_markers<-FindAllMarkers(Meta, only.pos=TRUE, min.pct=0.5, logfc.threshold = 1)
write.csv(All_markers, file="./CSV/findallmarkers_meta_int.csv")


#Heatmap
DefaultAssay(Meta)<-"RNA"
features<-c("REG3A","PRR4","MUC6","BPIFB1","MSMB",
            "LCN2","CXCL2","CXCL3","SOCS3","EFNA1",
            "LIPF","PGA3","PGA4","PGA5","PDIA2",
            "MUC5AC","GKN1","GKN2","TFF1","MUCL3",
            "CHRM3","SHISA6","RORA","CFTR","OLFM4",
            "ANPEP","DMBT1","PRAP1","FABP1","ALDOB",
            "CPE","PCSK1N","TTR","CHGA","SCGN",
            "PCLAF","UBE2C","CENPK","STMN1","TUBA1B",
            "TFPI","ODAM","CLDN2","ZG16B","SOD3",
            "TFF3","MUC2","ITLN1","SPINK4","REG4",
            "SFTPB","SFTA2","NAPSA","HOPX","WFDC2",
            "CHIA","CBLIF","APOA1","MT1F","MT1G")
DoHeatmap(subset(Meta, downsample=50), features = features, angle = 90, size=4, lines.width = 1) + 
  scale_fill_gradientn(colours = c("dark blue","gray","yellow")) + theme(axis.text.y = element_text(size = 10))


#Expected metaplasia transcript violin plot
Meta$celltype<- Idents(Meta)
meta<-VlnPlot(Meta, features=c("TFF2","MUC6","TFF3","TACSTD2"), group.by = "celltype", pt.size = 0,combine=FALSE)
CombinePlots(plots=meta, ncol=4, legend='none')


#Meta subtype signature
library(msigdbr)
msigdbr_show_species()
homo_df<- msigdbr(species = "Homo sapiens")
homo_sets<- homo_df %>% split(x = .$gene_symbol, f = .$gs_name)

Chief<-homo_sets$BUSSLINGER_GASTRIC_CHIEF_CELLS
chief_list <- list(Chief)
Goblet<-homo_sets$BUSSLINGER_DUODENAL_GOBLET_CELLS
gob_list <- list(Goblet)
Imm_entero<-homo_sets$BUSSLINGER_DUODENAL_LATE_IMMATURE_ENTEROCYTES
imm_entero_list <- list(Imm_entero)
ecl<-homo_sets$BUSSLINGER_GASTRIC_OXYNTIC_ENTEROCHROMAFFIN_LIKE_CELLS
ecl_list <- list(ecl)

combined<-readRDS('human_hp_aig_v5.rds')
combined<-RenameIdents(combined,'0'="B",'1'="B",'2'="B",'3'="Neck",'4'="T",
                       '5'="B",'6'="Pit",'7'="T",'8'="SPEM",'9'="Pit",'10'="Cancer",
                       '11'="Prolif",'12'="Endo",'13'="B",'14'="Fibro",'15'="Neuro",
                       '16'="Neck",'17'="Fibro",'18'="T",'19'="DC",'20'="Chief",
                       '21'="SM",'22'="Neuro",'23'="B",'24'="Mast",'25'="Prolif",
                       '26'="T",'27'="SM",'28'="Parietal",'29'="B",'30'="IM",
                       '31'="Neuro",'32'="Neuro",'33'="Neuro",'34'="B",'35'="B")
combined$epi_clusters<-Idents(combined)
Prolif_markers<-FindMarkers(combined, ident.1="Prolif", only.pos=TRUE, min.pct=0.5, logfc.threshold = 1)
prolif_list<-list(rownames(Prolif_markers))

Meta<-readRDS('Meta_aig+hp_v5_int.rds')
Meta$celltype<- Idents(Meta)

Meta <- AddModuleScore(object = Meta, features = chief_list, name = "chief_list_score")
chief_plot<-VlnPlot(Meta, features="chief_list_score1", 
                    group.by = "celltype", pt.size = 0, combine=FALSE, y.max = 1.5)

Meta <- AddModuleScore(object = Meta, features = gob_list, name = "gob_list_score")
gob_plot<-VlnPlot(Meta, features="gob_list_score1", 
                  group.by = "celltype",pt.size = 0, combine=FALSE, y.max = 3)

Meta <- AddModuleScore(object = Meta, features = imm_entero_list, name = "imm_entero_list_score")
imm_plot<-VlnPlot(Meta, features="imm_entero_list_score1", 
                  group.by = "celltype",pt.size = 0, combine=FALSE, y.max = 1.5)

Meta <- AddModuleScore(object = Meta, features =prolif_list, name = "prolif_list_score")
prolif_plot<-VlnPlot(Meta, features="prolif_list_score1", 
                     group.by = "celltype",pt.size = 0, combine=FALSE, y.max=1)

Meta <- AddModuleScore(object = Meta, features =ecl_list, name = "ecl_list_score")
ecl_plot<-VlnPlot(Meta, features="ecl_list_score1", 
        group.by = "celltype",pt.size = 0, combine=FALSE, y.max=1)

CombinePlots(plots=c(prolif_plot,chief_plot,imm_plot,ecl_plot,gob_plot), ncol=5, legend='none')


#cancer sig on meta
meta<-readRDS('Meta_aig+hp_v5_int.rds')
DefaultAssay(meta)<-"RNA"
cancer_sig <- list(c("EPCAM","TACSTD2","MMP7","CEACAM5","REG4","CLDN7","S100A6","APOA1","MUC13","CDH17",
                     "CLDN3","CDX1","OLFM4","ANPEP","ONECUT2","EFNA2","ANXA13","GUCY2C","DMBT1","HKDC1"))
meta <- AddModuleScore(object = meta, features = cancer_sig, name = "cancer_sig_score")
FeaturePlot(object = meta, features = "cancer_sig_score1", pt.size=0.4, label=TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdGy")))

quantile(meta@meta.data[["cancer_sig_score1"]], probs = c(0.9, 0.95))
hist(meta@meta.data[["cancer_sig_score1"]], 10)

Idents(meta)<-"condition"
a_meta<- subset(meta, idents="AIG")
Idents(a_meta)<-"seurat_clusters"
h_meta<- subset(meta, idents="Hp")
Idents(h_meta)<-"seurat_clusters"

a_meta <- AddModuleScore(object = a_meta, features = cancer_sig, name = "cancer_sig_score")
p1<-FeaturePlot(object = a_meta, features = "cancer_sig_score1", pt.size=0.4, label=TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdGy")))
h_meta <- AddModuleScore(object = h_meta, features = cancer_sig, name = "cancer_sig_score")
p2<-FeaturePlot(object = h_meta, features = "cancer_sig_score1", pt.size=0.4, label=TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdGy")))
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PRGn")))
grid.arrange(p2,p1, ncol = 2)


#Cell type projection mapping
meta<-readRDS('Meta_aig+hp_v5_int.rds')
DefaultAssay(meta)<-"RNA"

combined<-readRDS('/Users/macbook/Documents/Stella/Human_Hp_Meta/pos_pat/human_hp_adj_canc.rds')
DefaultAssay(combined)<-"RNA"
Idents(combined)<-"seurat_clusters"
combined<-RenameIdents(combined,'0'="B",'1'="B",'2'="B",'3'="T",'4'="B",
                       '5'="Fibro",'6'="Meta/Neck",'7'="T",'8'="Endo",'9'="Cancer",'10'="T",
                       '11'="B",'12'="Pit",'13'="Meta",'14'="Mast",'15'="SM",'16'="Mast",
                       '17'="Pit",'18'="Cancer",'19'="SM",'20'="Prolif",'21'="B",'22'="Neuro",
                       '23'="Pit",'24'="B",'25'="Mac",'26'="Mast",'27'="Stromal",'28'="B",'29'="B",
                       '30'="B",'31'="T",'32'="Chief/Parietal")
combined$epi_clusters<-Idents(combined)
epi<-subset(combined, idents = c("Pit","Meta","Meta/Neck","Chief/Parietal","Cancer","Neuro","Fibro","Endo","Prolif","SM"))

VlnPlot(combined, features="ANPEP", group.by = "epi_clusters", pt.size = 0,combine=FALSE)
VlnPlot(meta, features="ANPEP", group.by = "seurat_clusters", split.by="condition", cols=c("red","blue"), pt.size = 0.1,combine=FALSE)

anchors <- FindTransferAnchors(reference = epi, query = meta,
                               dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = epi$epi_clusters,
                            dims = 1:30)
meta.pred.combined <- AddMetaData(meta, metadata = predictions) #Adds predictions object
#visualize predicted mouse cluster IDs based on human clusters
DimPlot(meta.pred.combined, reduction = "umap", group.by = "predicted.id", label =TRUE, repel = TRUE, raster=FALSE,cols=c("#B0172B","#71797E","#C0C0C0","#D3D3D3","#818589","#A9A9A9","#36454F"))
DimPlot(meta.pred.combined, reduction = "umap", group.by = "predicted.id", split.by="condition",label = FALSE, repel = TRUE, raster=FALSE, cols=c("#B0172B","#71797E","#C0C0C0","#D3D3D3","#818589","#A9A9A9","#36454F"))
