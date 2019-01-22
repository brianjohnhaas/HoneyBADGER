#!/usr/bin/env Rscript

options(error = function() traceback(2))

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="hb.example.matrix",
                                    annotations_file="hb.example.cell_annots",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("normal") )
                                   



saveRDS(infercnv_obj, file='infercnv.hb_example.obj')

