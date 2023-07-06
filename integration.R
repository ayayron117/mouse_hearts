library(Seurat)

input_path <- "/Users/asm/Creighton/core_lab/scRNAseq_Ren"

# matrices_path <- file.path(input_path, "matrices")
# dir.create(matrices_path)
seurats_path <- file.path(input_path, "seurats")
dir.create(seurats_path)

################################################################################

wild <- Read10X(file.path(input_path, "Wildtype", "filtered_feature_bc_matrix"),  
                gene.column = 2,  
                cell.column = 1,  
                unique.features = TRUE,  
                strip.suffix = FALSE)

knock <- Read10X(file.path(input_path, "Knockout", "filtered_feature_bc_matrix"),  
                 gene.column = 2,  
                 cell.column = 1,  
                 unique.features = TRUE,  
                 strip.suffix = FALSE)

################################################################################

wild_s <- CreateSeuratObject(wild, project = "wildtype")
knock_s <- CreateSeuratObject(knock, project = "knockout")

################################################################################

wild_s[["sample_id"]] <- "wildtype"
knock_s[["sample_id"]] <- "knockout"

################################################################################

seurat_obj <- merge(x = wild_s, y = knock_s,
                add.cell.ids = c("wildtype","knockout"),
                project = "mouse_hearts")

saveRDS(seurat_obj, file.path(seurats_path, "merged_seurat.rds"))

################################################################################

ifnb.list <- SplitObject(seurat_obj, split.by = "sample_id")

################################################################################

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

################################################################################

features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = TRUE)
})

saveRDS(ifnb.list, file.path(seurats_path, "ifnb_list.rds"))

################################################################################

anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
saveRDS(anchors, file.path(seurats_path, "anchors.rds"))

################################################################################

# anchors <- readRDS(file.path(seurats_path, "anchors.rds"))

combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"

saveRDS(combined, file.path(seurats_path, "integrated_mouse_seurat.rds"))















