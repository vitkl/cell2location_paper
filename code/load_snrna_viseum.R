##' @rdname load_snrna_viseum
##' @name load_snrna_viseum
##' @description \code{load_snrna_viseum}: Load snRNA-seq of mouse cortex that matches trialed Viseum chips
##' @export load_snrna_viseum
##' @import data.table
##' @import Seurat
##' @import Matrix
##' @import DropletUtils
load_snrna_viseum = function(nfeatures = 5000, M3Drop_threshold = 0.1, emptyDrops_fdr_p_val = 0.01, scrublet_score = 0.2,
                             data_dir = "./data/mouse_viseum_snrna",
                             dir10x = c("raw_feature_bc_matrix", "filtered_feature_bc_matrix")){

  folders = list.dirs(paste0(data_dir, "/rawdata/", dir10x[1]), recursive = F) #paste0("./data/mouse_viseum_snrna/rawdata/", manifest$Sample, "/")

  cache_file = paste0(data_dir, "/rawdata/all.mtx")

  if(file.exists(cache_file)) { # read processed matrix (cells selected with emptydrops)

    # read data
    data = Matrix::readMM(cache_file)
    # read cell names
    cell_names = fread(paste0(data_dir, "/rawdata/all_cells.txt"),
                       header = T, stringsAsFactors = F)
    colnames(data) = cell_names$cell_id
    # read gene names
    gene_names = fread(paste0(data_dir, "/rawdata/all_genes.txt"),
                       header = T, stringsAsFactors = F)
    rownames(data) = gene_names$gene_id

  } else { # select cells with emptydrops

    data = foreach(i = seq_along(folders), .combine = cbind) %do% {

      files = paste0(folders[i], "/", list.files(folders[i]))

      # read data
      data = Matrix::readMM(files[3])

      # read gene names
      row_data = fread(files[2], stringsAsFactors = F, header = F)[, 1:2]
      colnames(row_data) = c("ENSEMBL", "SYMBOL")

      # read cell barcodes
      col_data = fread(files[1], stringsAsFactors = F, header = F)

      # assign row and col names
      rownames(data) = row_data$ENSEMBL
      sample_name = substr(folders[i], nchar(folders[i]) - nchar("5705STDY8058280") + 1, nchar(folders[i]))
      colnames(data) = paste0(sample_name, "_", col_data$V1)

      if(dir10x != "filtered_feature_bc_matrix") {

        # remove empty droplets
        dr = DropletUtils::emptyDrops(data, lower=100, retain=NULL)

        is.cell = dr$FDR <= emptyDrops_fdr_p_val
        #sum(is.cell, na.rm=TRUE)

        # Check if p-values are lower-bounded by 'niters'
        # (increase 'niters' if any Limited==TRUE and Sig==FALSE)
        message(sample_name)
        table(Sig=is.cell, Limited=dr$Limited)

        print(plot_ngenes_vs_cum_cells(data, is.cell, sample_name))

        # return cells
        return(data[, is.cell & !is.na(is.cell)])

      } else {

        return(data)
      }

    }

    Matrix::writeMM(data, cache_file)
    fwrite(data.table(cell_id = colnames(data)),
           paste0(data_dir, "/rawdata/all_cells.txt"))
    fwrite(data.table(gene_id = rownames(data)),
           paste0(data_dir, "/rawdata/all_genes.txt"))

  }

  files = paste0(folders[1], "/", list.files(folders[1]))
  row_data = fread(files[2], stringsAsFactors = F, header = F)[, 1:2]
  colnames(row_data) = c("ENSEMBL", "SYMBOL")
  row_data = as.data.frame(row_data)
  rownames(row_data) = row_data$ENSEMBL

  #remove non-expressed genes (less than 3 UMI total)
  data = data[Matrix::rowSums(data) > 3,]
  row_data = row_data[rownames(data),]

  # create metadata
  meta_data = data.table(cell_id = colnames(data))
  meta_data$sample = gsub("_[[:alpha:]]{16}-1", "", meta_data$cell_id)
  #meta_data = merge(meta_data, manifest, by.x = "sample", by.y = "Sample", all.x = T, all.y = F)
  meta_data = as.data.frame(meta_data)
  rownames(meta_data) = meta_data$cell_id

  # create Seurat object
  seu = CreateSeuratObject(counts = data,
                           meta.data = meta_data,
                           project = "snrna_viseum")

  # add gene metadata
  seu@assays$RNA@meta.features = row_data

  # scrublet doublet detection pipeline
  seu = detect_scrublets(seu, sample_col = "sample",
                         doublet_cutoff=scrublet_score, expected_doublet_rate = 0.06,
                         sim_doublet_ratio = 2, min_counts=3, min_cells=3,
                         min_gene_variability_pctl=75, n_prin_comps=50
  )
  # seu = recall_doublets(seu, doublet_cutoff=0.19)

  # standard Seurat processing pipeline
  seu = std_seu_processing(seu, nfeatures = nfeatures, M3Drop_threshold = M3Drop_threshold,
                           percent.mt_filter = 16, resolution = 3,
                           species = 10090)

  # # plot diagnostic plots
  # p1 = FeaturePlot(seu, reduction = "umap", pt.size = 0.1, features = "percent.mt") +
  #   ggtitle(label = "percent mito genes")
  # p2 = FeaturePlot(seu, reduction = "umap", pt.size = 0.1, features = "nCount_RNA") +
  #   ggtitle(label = "nUMI")
  # p3 = DimPlot(seu, reduction = "umap", pt.size = 0.1, group.by = "sample") +
  #   ggtitle(label = "Sample")
  # p4 = DimPlot(seu, reduction = "umap", pt.size = 0.1, group.by = "Sorting") +
  #   ggtitle(label = "Sorting")
  # p5 = DimPlot(seu, reduction = "umap", pt.size = 0.1, group.by = "Region") +
  #   ggtitle(label = "Region")
  # p6 = FeaturePlot(seu, reduction = "umap", pt.size = 0.1, features = "doublet_scores") +
  #   ggtitle(label = "scrublet scores")
  # print(cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, align = "hv"))

  seu

}
