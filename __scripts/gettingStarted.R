#!/usr/bin/env Rscript

## from:  https://jef.works/HoneyBADGER/Getting_Started.html

library(HoneyBADGER)

data(gexp) ## tumor cells
data(ref) ## reference

require(biomaRt) ## for gene coordinates
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")

print(gexp[1:5,1:5])


print(ref[1:5])

print(mart.obj)



hb <- new('HoneyBADGER', name='MGH31')
hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)

pdf("hbadger.plots")

hb$plotGexpProfile() ## initial visualization


hb$setMvFit(verbose=TRUE)
hb$setGexpDev(verbose=TRUE)
hb$calcGexpCnvBoundaries(init=TRUE, verbose=FALSE)


## double check what CNVs were identified
bgf <- hb$bound.genes.final
genes <- hb$genes
regions.genes <- range(genes[unlist(bgf)])
print(regions.genes)

## Indeed, our initial HMM has identified a number of candidate CNVs to test. We can now retest all identified CNVs on all cells to derive the final posterior probability of each CNV in each cell. We can cluster cells on these posterior probabilities and visualize them as a heatmap.

hb$retestIdentifiedCnvs(retestBoundGenes = TRUE, retestBoundSnps = FALSE, verbose=FALSE)

## look at final results
results <- hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
print(head(results[,1:5]))

## visualize as heatmap 
trees <- hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(25,15))

## order cells
hc <- trees$hc
order <- hc$labels[hc$order]
## plot all chromosomes
hb$plotGexpProfile(cellOrder=order)


## plot just identified cnvs
hb$plotGexpProfile(cellOrder=order, region=hb$cnvs[['gene-based']][['amp']])

hb$plotGexpProfile(cellOrder=order, region=hb$cnvs[['gene-based']][['del']])

dev.off()

