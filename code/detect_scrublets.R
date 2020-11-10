detect_scrublets = function(seu, # Seurat object
                            sample_col, # column name of 10X (or other droplet protocol) run ids telling which cells can potentially be doublets
                            doublet_cutoff=0.3, #replace the default scrublet cutoff
                            expected_doublet_rate = 0.06, # Details: https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb
                            sim_doublet_ratio = 2, # the number of doublets to simulate, relative to the number of observed transcriptomes
                            min_counts=3, 
                            min_cells=3, 
                            min_gene_variability_pctl=75, 
                            n_prin_comps=50
) {
  
  # import scrublet python module
  scr = reticulate::import("scrublet")
  
  sample_ids = unique(seu@meta.data[,sample_col])
  
  doublets = foreach(i = sample_ids, .combine = rbind) %dopar% {
    
    seu_i = seu[, seu@meta.data[,sample_col] == i]
    mat = Matrix::t(seu_i@assays$RNA@counts)
    
    scrub = scr$Scrublet(mat, expected_doublet_rate=expected_doublet_rate, 
                         sim_doublet_ratio = sim_doublet_ratio)
    
    doublet_scores = scrub$scrub_doublets(min_counts=as.integer(min_counts), 
                                          min_cells=as.integer(min_cells), 
                                          min_gene_variability_pctl=min_gene_variability_pctl, 
                                          n_prin_comps=as.integer(n_prin_comps))
    names(doublet_scores) = c("doublet_scores", "predicted_doublets_scr")
    
    doublet_scores = as.data.frame(doublet_scores)
    doublet_scores$predicted_doublets = doublet_scores$doublet_scores > doublet_cutoff
    doublet_scores$cell_id = rownames(mat)
    rownames(doublet_scores) = doublet_scores$cell_id
    
    doublet_scores
  }
  
  seu$doublet_scores = doublets[colnames(seu), "doublet_scores"]
  seu$predicted_doublets = doublets[colnames(seu), "predicted_doublets"]
  seu$predicted_doublets_scr = doublets[colnames(seu), "predicted_doublets_scr"]
  
  plot_list = lapply(sample_ids, function(sample_id) {
    qplot(seu[, seu@meta.data[,sample_col] == sample_id]$doublet_scores,
          geom = "histogram", bins = 50) +
      scale_y_log10() +
      geom_vline(xintercept = doublet_cutoff) +
      xlab("doublet_scores") + ggtitle(sample_id) +
      xlim(-0.01, 1.01)
  })
  print(cowplot::plot_grid(plotlist = plot_list, ncol = ceiling(length(sample_ids) / 2), align = "hv"))
  
  seu
  
}

recall_doublets = function(seu, doublet_cutoff=0.2){ # reassign doublets based on the new cutoff
  
  seu$predicted_doublets = seu$doublet_scores > doublet_cutoff
  
  plot_list = lapply(sample_ids, function(sample_id) {
    qplot(seu[, seu@meta.data[,sample_col] == sample_id]$doublet_scores,
          geom = "histogram", bins = 50) +
      scale_y_log10() +
      geom_vline(xintercept = doublet_cutoff) +
      xlab("doublet_scores") + ggtitle(sample_id) +
      xlim(-0.01, 1.01)
  })
  print(cowplot::plot_grid(plotlist = plot_list, ncol = ceiling(length(sample_ids) / 2), align = "hv"))
  
  seu
  
}