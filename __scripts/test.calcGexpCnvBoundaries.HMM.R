#!/usr/bin/env Rscript

library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
gexp.norm = gexp.mats$gexp.norm
genes = gexp.mats$genes
m=0.15
chrs=paste0('chr', c(1:22))
min.traverse=3
t=1e-6
min.num.genes=3
trim=0.1
verbose=TRUE


genes <- genes[rownames(gexp.norm)]
## order
gos <- as.data.frame(genes)
rownames(gos) <- names(genes)
gos <- gos[rownames(gexp.norm),]
tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
    na.omit(gexp.norm[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],])
})
tl <- tl[chrs]  # gexp broken out by chromosome, restricted to main chrs 


gexp.norm <- do.call(rbind, lapply(tl, function(x) x))

## smooth
k = 101
mat.smooth <- apply(gexp.norm, 2, caTools::runmean, k)
d <- dist(t(mat.smooth))
d[is.na(d)] <- 0
d[is.nan(d)] <- 0
d[is.infinite(d)] <- 0
hc <- hclust(d, method="ward.D2")

## iterative HMM
heights <- seq_len(min(min.traverse, ncol(gexp.norm)))
## cut tree at various heights to establish groups

h=3 # just to explore one example
ct <- cutree(hc, k = h)
cuts <- unique(ct)

## just grab cut 1 for example purposes.
group=1
mat.smooth <- apply(gexp.norm[, ct==1], 1, mean) # just getting the mean for each group.
## change point
delta <- c(0, 1, 0)
t <- t
pd <- -m
pn <- 0
pa <- m
sd <- sd(mat.smooth)
z <- HiddenMarkov::dthmm(mat.smooth,
                         matrix(c(1-2*t, t, t, t, 1-2*t, t, t, t, 1-2*t),
                                byrow=TRUE, nrow=3),
                         delta,
                         "norm",
                         list(mean=c(pd, pn, pa),
                              sd=c(sd,sd,sd)))

results <- HiddenMarkov::Viterbi(z)

print(results)


