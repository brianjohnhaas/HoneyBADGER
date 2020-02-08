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

gexp.norm <- hb$gexp.norm
genes <- hb$genes
dev <- hb$dev
mvFit <- hb$mvFit


#' Calculate posterior probability of CNVs using normalized expression data
#'
#' @name HoneyBADGER_calcGexpCnvProb
#' @param gexp.norm.sub Optional normalized gene expression matrix. If not provided, internal normalized gene expression matrix is used.
#' @param m Expected mean deviation due to copy number change (default: 0.15)
#' @param region Optional GenomicRanges region of interest such as expected CNV boundaries. (default: NULL)
#' @param quiet Boolean for whether to suppress progress display (default: TRUE)
#' @param verbose Verbosity (default: FALSE)
#'
#HoneyBADGER$methods(

calcGexpCnvProb=function(gexp.norm.sub=NULL, m=0.15, region=NULL, quiet=TRUE, verbose=FALSE) {
    if(!is.null(gexp.norm.sub)) {
        gexp.norm <- gexp.norm.sub
        genes <- genes[rownames(gexp.norm.sub)]
        mvFit <- mvFit[colnames(gexp.norm.sub)]
    }

    gexp <- gexp.norm
    gos <- genes
    fits <- mvFit

    if(!is.null(region)) {
        ## removing this code for now... doesn't seem to be used w/ current call here.
        stop("region is specified...  I was wrong.")
    }

    mu0 <- apply(gexp, 2, mean)  # for each cell, get the mean across genes.
    ng <- nrow(gexp) # num genes in this region.

    # get the cell-specific variance in expression according to number of genes in the block.
    sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))


    ## Model
    if(verbose) {
        cat('Aggregating data to list ... \n')
    }
    data <- list(
        'K' = length(mu0), # number of cells
        'JJ' = nrow(gexp), # number of genes in block
        'gexp' = gexp,
        'sigma0' = sigma0,
        'mag0' = m  # dev
    )
    modelFile <-  system.file("bug", "expressionModel.bug", package = "HoneyBADGER")

    if(verbose) {
        cat('Initializing model ... \n')
    }

    inits <- list(
        list(S = rep(0, ncol(gexp)), dd = 0),
        list(S = rep(1, ncol(gexp)), dd = 0),
        list(S = rep(0, ncol(gexp)), dd = 1),
        list(S = rep(1, ncol(gexp)), dd = 1)
    )

    require(rjags)
    model <- jags.model(modelFile, data=data, inits=inits, n.chains=4, n.adapt=100, quiet=quiet)
    update(model, 100, progress.bar=ifelse(quiet,"none","text"))

    parameters <- c('S', 'dd', 'mu')
    samples <- coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
    samples <- do.call(rbind, samples) # combine chains

    if(verbose) {
        cat('...Done!')
    }

    snpLike <- samples
    v <- colnames(snpLike)
    S <- snpLike[,grepl('S', v)]
    dd <- snpLike[,grepl('dd', v)]
    mu <- snpLike[,grepl('mu', v)]
    ##plot(mu0, colMeans(mu))
    delcall <- apply(S*(1-dd), 2, mean)
    ampcall <- apply(S*dd, 2, mean)
    ##plot(mu0, delcall)
    ##plot(mu0, ampcall)
    names(ampcall) <- names(delcall) <- colnames(gexp)

    return(list('posterior probability of amplification'=ampcall,
                'posterior probability of deletion'=delcall,
                'estimated mean normalized expression deviation'=mu0))
}


### test:

## test amp:
bound.genes.cont = read.table("test.amp.3_regions.dat")
genenames = rownames(bound.genes.cont)
bound.genes.cont = bound.genes.cont[,1]
names(bound.genes.cont) = genenames
region_name = "1732"
bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == region_name]

prob <- calcGexpCnvProb(gexp.norm.sub=gexp.norm[bound.genes.new, ], m=dev, verbose=TRUE)

cat("AMPLIFICATION PROBABILITY: ")
cat(prob[[1]])

cat("DELETION PROBABILITY: ")
cat(prob[[2]])


## test del:
bound.genes.cont = read.table("test.del.3_regions.dat")
genenames = rownames(bound.genes.cont)
bound.genes.cont = bound.genes.cont[,1]
names(bound.genes.cont) = genenames
region_name = "2667"
bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == region_name]

prob <- calcGexpCnvProb(gexp.norm.sub=gexp.norm[bound.genes.new, ], m=dev, verbose=TRUE)

cat("AMPLIFICATION PROBABILITY: ")
cat(prob[[1]])

cat("DELETION PROBABILITY: ")
cat(prob[[2]])




