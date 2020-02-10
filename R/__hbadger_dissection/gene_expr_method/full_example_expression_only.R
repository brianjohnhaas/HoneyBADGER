#!/usr/bin/env Rscript

library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)

mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
hb <- new('HoneyBADGER', name='MGH31')

hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)

hb$setMvFit(verbose=TRUE) ## model variance
hb$setGexpDev(verbose=TRUE) ## model necessary expression deviation to identify CNVs
hb$calcGexpCnvBoundaries(init=TRUE, verbose=FALSE) ## HMM

hb$retestIdentifiedCnvs(retestBoundGenes = TRUE, retestBoundSnps = FALSE, verbose=FALSE)

results <- hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)

write.table(results, "hb.gene_expr_only.results.tsv", quote=F, sep="\t")

