#!/usr/bin/env Rscript

library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
hb <- new('HoneyBADGER', name='MGH31')
hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)


gexp.norm <- hb$gexp.norm




#' Model expected gene expression variance as a function of number of genes
#'
#' @name HoneyBADGER_setMvFit
#' @param num.genes Number of random genes sampled (default: seq(5, 100, by=5))
#' @param rep Number of repeats/resampling (default: 50)
#' @param plot Whether to plot (default: FALSE)
#' @param verbose Verbosity (default: TRUE)
#'
#HoneyBADGER$methods(
#                setMvFit=function(

num.genes = seq(5, 100, by=5)
rep = 50
plot=FALSE
verbose=TRUE

if(verbose) {
    cat('Modeling expected variance ... ')
}
## Models variance of the sample mean

mean.var.comp <- lapply(num.genes, function(ng) {
    set.seed(0)
    m <- do.call(rbind, lapply(1:rep, function(i) {
        nrmchr.sub <- gexp.norm[sample(1:nrow(gexp.norm), ng),] # sample ng genes from cancer cell expression matrix
        nm <- apply(nrmchr.sub, 2, mean) # get mean of gene expression per cell for a sample of genes.
        nm
    }))
    return(m)
})
names(mean.var.comp) <- num.genes

# mena.var.comp: list of matrices, each matrix is (sample mean) x (cell), and one matrix for each num.genes sample size


fits <- lapply(1:ncol(gexp.norm), function(k) {
    ## for one cell
    mean.comp <- do.call(cbind, lapply(mean.var.comp, function(x) x[,k]))

    if(plot) {
        par(mfrow=c(1,3), mar=rep(2,4))
        perf.test <- function(mat) {
            require(ggplot2)
            require(reshape2)
            m <- melt(mat)
            p <- ggplot(m) + geom_boxplot(aes(x = factor(Var2), y = value))
            return(p)
        }
        perf.test(mean.comp)
    }

    if(plot) {
        plot(log10(num.genes),log10(apply(mean.comp, 2, var)), type="l")
    }
    df <- data.frame('x'=num.genes, 'y'=apply(mean.comp, 2, var))
    fit <- lm(log10(y)~log10(x), data=df)

    if(plot) {
        x2 <- log10(num.genes)
        y2 <- predict(fit, x=x2, interval="predict")[, 'fit']
        plot(x2, y2)
        plot(log10(df$x), log10(df$y), type="l")
        points(x2, y2, type="l", col="red")
        plot(df$x, df$y, type="l")
        points(10^x2, 10^y2, type="l", col="red")
    }

    return(fit)
})
names(fits) <- colnames(gexp.norm)

## each cell gets fit separately.



