#!/usr/bin/env Rscript

library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
hb <- new('HoneyBADGER', name='MGH31')
hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)


gexp.norm <- hb$gexp.norm

#####


#setGexpDev=function(alpha=0.05, n=100, seed=0, plot=FALSE, verbose=FALSE) {

alpha=0.05
n=100
seed=0
plot=FALSE
verbose=FALSE


k = 101
set.seed(seed)
gexp.sd <- sd(gexp.norm)
devs <- seq_len(10)/10

pvs <- unlist(lapply(devs, function(dev) {
    mean(unlist(lapply(seq_len(n), function(i) {
        pv <- ks.test(rnorm(k, 0, gexp.sd), rnorm(k, dev, gexp.sd))
        pv$p.value
    })))
}))


if(plot) {
    plot(pvs, devs, xlab="p-value", ylab="deviation", xlim=c(0,1))
}
fit <- lm(devs ~ pvs)
optim.dev <- predict(fit, newdata=data.frame(pvs=alpha))
if(verbose) {
    cat('Optimal deviance: ')
        cat(optim.dev)
}
#dev <<- optim.dev
message("Dev: ", optim.dev)


