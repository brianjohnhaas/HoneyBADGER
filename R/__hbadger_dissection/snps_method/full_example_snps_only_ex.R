#!/usr/bin/env Rscript

library(HoneyBADGER)
data(gexp) ## tumor cells
data(ref) ## reference
require(biomaRt) ## for gene coordinates
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
data(r) ## alternate allele
data(cov.sc) ## total coverage
library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## in order to map SNPs to genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

hb <- new('HoneyBADGER', name='MGH31')
hb$setAlleleMats(r.init=r, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=detectCores())
hb$setGeneFactors(txdb) ## map SNPs to genes
#hb$plotAlleleProfile()

hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
bsf <- get('bound.snps.final', slot(hb, '.xData'))
snps <- get('snps', slot(hb, '.xData'))
regions.snp <- range(snps[unlist(bsf)])
print(regions.snp)



hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=TRUE)
results <- hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
