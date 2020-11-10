##' Perform standard Seurat processing
##' @rdname std_seu_processing
##' @name std_seu_processing
##' @export std_seu_processing
##' @import Seurat
##' @import Matrix
##' @import M3Drop
##' @import ParetoTI
std_seu_processing = function(seu, scale.factor = 10000,
                              nfeatures = 5000, # N highly variable genes
                              M3Drop_threshold = 0.1, # FDR-corrected p-value cutoff for selecting genes with M3Drop
                              vst.mean_filter = 0, # filtering genes
                              percent.mt_filter = 15, # filtering cells
                              npcs = 75, # number of PCs to compute
                              nn.eps = 0.5, resolution = 2, # SNN and clustering parameters
                              min.dist = 0.75, # UMAP parameter
                              species = 9606
) {


  # count UMI and cells per gene
  seu@assays$RNA@meta.features$nCount_RNA = Matrix::rowSums(seu@assays$RNA@counts)
  seu@assays$RNA@meta.features$nCells_RNA = Matrix::rowSums(seu@assays$RNA@counts > 0)

  # remove practically undetected genes
  seu = seu[seu@assays$RNA@meta.features$nCells_RNA > 3,]

  # normalise ============================================================
  seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = scale.factor)

  # compute highly variable genes ============================================================
  # select genes with surprisingly high dropout rate for their mean expression level
  m3d = M3Drop::M3DropConvertData(seu@assays$RNA@counts, is.counts=T)
  m3d = as(m3d, "dgCMatrix")
  #rownames(m3d) = rownames(seu)
  #colnames(m3d) = colnames(seu)
  res = M3Drop::M3DropFeatureSelection(m3d, mt_method="fdr",
                                       mt_threshold=M3Drop_threshold, suppress.plot=FALSE,
                                       xlim=NA)
  #hist(log10(res$effect.size))
  m3d_ind = rownames(res)[log10(res$effect.size) >
                            mean(log10(res$effect.size)) + sd(log10(res$effect.size))]

  # select gene with surprisingly high difference in expression between cells
  seu = FindVariableFeatures(seu, nfeatures = nfeatures) #using vst
  sc_ind = seu@assays$RNA@var.features

  # plot overlap in two selection strategies
  grid::grid.newpage()
  print(grid::grid.draw(VennDiagram::draw.pairwise.venn(area1 = length(sc_ind), area2 = length(m3d_ind),
                                                        cross.area = length(intersect(sc_ind, m3d_ind)),
                                                        category = c("Seurat VST", "M3Drop"))))
  grid::grid.newpage()

  # combine M3Drop and deviance genes
  seu@assays$RNA@var.features = union(sc_ind, m3d_ind)
  # find genes with high poisson deviance
  #seu = compute_gene_deviance(seu)

  # filter by expression ============================================================
  seu = seu[seu@assays$RNA@meta.features$vst.mean > vst.mean_filter,]

  # mitochondrial reads ============================================================
  # get mitochondria_located_genes (mitochondria- & nucleus-encoded):
  if(species != 9606) {
    col = c("GOALL", "ENSEMBL", "SYMBOL")
  } else {
    col = c("GOALL", "ENSEMBL", "SYMBOL", "MAP")
  }
  go_annot = ParetoTI::map_go_annot(taxonomy_id = species, keys = "GO:0005739",
                                    columns = col, keytype = "GOALL",
                                    ontology_type = c("CC"))
  annot_dt = go_annot$annot_dt[!is.na(ENSEMBL)]
  annot_dt = annot_dt[ENSEMBL %in% rownames(seu),]

  # compute the fraction of mitochondria-encoded reads (have no chromosome position "MAP")
  if(species != 9606) {
    mt_genes = seu@assays$RNA@meta.features$ENSEMBL[grepl("^mt-", seu@assays$RNA@meta.features$SYMBOL)]
    mitochondria_located_genes = unique(annot_dt[, ENSEMBL])
    mitochondria_located_genes = mitochondria_located_genes[!mitochondria_located_genes %in% mt_genes]
  } else {
    mt_genes = annot_dt[is.na(MAP), ENSEMBL]
    mitochondria_located_genes = unique(annot_dt[!is.na(MAP), ENSEMBL])
  }

  if(length(mt_genes) > 0) { # if MT-genes found
    seu$percent.mt = PercentageFeatureSet(seu, features = mt_genes) / 100
  } else seu$percent.mt = 0

  # compute the fraction of nucleus-encoded reads
  mitochondria_located_genes = mitochondria_located_genes[mitochondria_located_genes %in% rownames(seu)]
  seu$nucl_mito_genes = Matrix::colSums(seu@assays$RNA@counts[mitochondria_located_genes, ]) / seu$nCount_RNA

  #cut_off_line = geom_path(aes(x, y), color = "orange", data = data.table(x = seq(0, 1, 0.01), y = sqrt((seq(0, 1, 0.01) - 0.0) * 0.05)))
  cut_off_line = geom_vline(xintercept = percent.mt_filter/100, color = "orange")
  print(qplot(seu$percent.mt, seu$nucl_mito_genes, geom = "bin2d", bins = c(100, 25)) +
          scale_fill_viridis_c(trans = "log10") + scale_color_viridis_c(trans = "log10") +
          xlim(0, 1) + ylim(0, 0.25) +
          xlab("Mitochondria-encoded") + ylab("Nucleus-encoded") +
          coord_fixed() + theme_bw() +
          ggtitle("all cells") + cut_off_line)
  print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nucl_mito_genes"), ncol = 3))

  seu = seu[, seu$percent.mt * 100 < percent.mt_filter]

  # scaling data  ============================================================
  seu = ScaleData(seu, vars.to.regress = NULL)

  # PCA ============================================================
  seu = RunPCA(seu, npcs = npcs, ndims.print = 1:5, nfeatures.print = 5)
  print(ElbowPlot(seu, ndims = npcs))

  # Graph-based clustering ============================================================
  seu = FindNeighbors(seu, reduction = "pca", dims = 1:npcs, nn.eps = nn.eps)
  seu = FindClusters(seu, resolution = resolution, n.start = 10)

  # compute UMAP ============================================================
  seu = RunUMAP(seu, dims = 1:npcs, min.dist = min.dist)

  seu
}
