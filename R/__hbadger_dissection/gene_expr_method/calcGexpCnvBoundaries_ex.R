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


pred.genes.r <- matrix(0, nrow(gexp.norm), ncol(gexp.norm))
rownames(pred.genes.r) <- rownames(gexp.norm)
colnames(pred.genes.r) <- colnames(gexp.norm)
bound.genes.old <- c()
bound.genes.final <- list()

#' Recursive HMM to identify CNV boundaries using normalized gene expression data
#'
#' @name HoneyBADGER_calcGexpCnvBoundaries
#' @param gexp.norm.sub Optional normalized gene expression matrix. Useful for manual testing of restricted regions. If NULL, gexp.norm within the HoneyBADGER object will be used.
#' @param chrs List of chromosome names. Genes not mapping to these chromosomes will be excluded. Default autosomes only: paste0('chr', c(1:22))
#' @param min.traverse Depth traversal to look for subclonal CNVs. Higher depth, potentially smaller subclones detectable. (default: 3)
#' @param min.num.genes Minimum number of genes within a CNV. (default: 5)
#' @param t HMM transition parameter. Higher number, more transitions. (default: 1e-6)
#' @param init Initialize recursion (default: FALSE)
#' @param verbose Verbosity (default: FALSE)
#' @param ... Additional parameters for \code{\link{calcGexpCnvProb}}
#'
#' Operation summary:
#'   start w/ init, using all cells and all genes.
#'   smooth expr vals along chrs, and cluster cells
#'   iterative HMM on segmented clusters
#'      get lists of amp and del genes for each cluster segment
#'   for each list of (amp or del) regions:
#'      collapse del and amp genes from above iterations into genomic regions of min length
#'   for each collapsed del or amp region:
#'      get cell probabilities for each region according to del or amp.
#'   make combined table of (row=region, col=cell, val=pvalue)
#'   binarize table
#'   pick region with highest clonality for amp or deletion
#'   extract cells into two groups: cell_G1: those with high confidence amp/del, and cell_G1: those with no amp/del
#'   for cell_G1:
#'      store the region genes
#'      record these genes were assigned to these cells.
#'   recurse on cell_G1, minus assigned genes.
#'      further partition cell_G1 based on other genes.
#'   recurse on cell_G2, minus assigned genes.
#'      further partition cell_G2 based on other genes.
#'
#'
#'


#HoneyBADGER$methods(

calcGexpCnvBoundaries=function(gexp.norm.sub=NULL, chrs=paste0('chr', c(1:22)), min.traverse=5, min.num.genes=10, t=1e-6, init=FALSE, verbose=FALSE, recursion_count = 0, G1orG2 = "BOTH",...) {
    if(!is.null(gexp.norm.sub)) {
        gexp.norm <- gexp.norm.sub
        genes <- genes[rownames(gexp.norm)]
    }
    if(init) {
        pred.genes.r <<- matrix(0, nrow(gexp.norm), ncol(gexp.norm))
        rownames(pred.genes.r) <<- rownames(gexp.norm)
        colnames(pred.genes.r) <<- colnames(gexp.norm)
        bound.genes.old <<- c()
        bound.genes.final <<- list()
    }
    if(is.null(bound.genes.final)) {
        cat('ERROR! USE init=TRUE! ')
        return()
    }

    message("Dim(gexp): ", paste(dim(gexp.norm), sep=",", collapse=","), ", recursion count: " , recursion_count, ", G1orG2: ", G1orG2 )

    ## remove old bound genes
    vi <- !(rownames(gexp.norm) %in% bound.genes.old)
    gexp.norm <- gexp.norm[vi,]

    ## order
    gos <- as.data.frame(genes)
    rownames(gos) <- names(genes)
    gos <- gos[rownames(gexp.norm),]
    tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
        na.omit(gexp.norm[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],])
    })
    tl <- tl[chrs]
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
    heights <- seq_len(min(min.traverse, ncol(gexp.norm))) # default min.traverse=5 ... so 1 to 5 clusters max.
    ## cut tree at various heights to establish groups

    boundgenes.pred <- lapply(heights, function(h) {
        ct <- cutree(hc, k = h)
        cuts <- unique(ct)  # assigns each cell to their cut grouping number
        ## look at each group, if deletion present
        boundgenes.pred <- lapply(cuts, function(group) {
            if(sum(ct==group)>1) { # have more than one cell in that group

                # get the mean expr per gene for this group of cells.
                mat.smooth <- apply(gexp.norm[, ct==group], 1, mean)

                # mat.smooth is now just a numeric vector of mean gene expr vals.

                ## run HMM for segmentationm
                ## change point
                delta <- c(0, 1, 0)
                t <- t
                pd <- -dev
                pn <- 0
                pa <- dev
                sd <- sd(mat.smooth)
                z <- HiddenMarkov::dthmm(mat.smooth,
                                         matrix(c(1-2*t, t, t, t, 1-2*t, t, t, t, 1-2*t),
                                                byrow=TRUE, nrow=3),
                                         delta,
                                         "norm",
                                         list(mean=c(pd, pn, pa), # using +/- dev to define means for amp and deltions.
                                              sd=c(sd,sd,sd) ) # shared sd, set based on sd for this set of genes.
                                         )
                results <- HiddenMarkov::Viterbi(z)

                ampgenes <- names(mat.smooth)[which(results==3)]
                delgenes <- names(mat.smooth)[which(results==1)]
                boundgenes <- list('amp'=ampgenes, 'del'=delgenes)

                return(boundgenes)
            }
        })

    })


    boundgenes.pred <- unlist(boundgenes.pred, recursive=FALSE) # flatten to just get all the various amp and del gene sets.


    ## get per-cell amp and del probabilites

    getTbv <- function(boundgenes.pred, amp_or_del) { # restricted to just the amp or dels, as invoked separately below.

        # define the list of genes
        foo <- rep(0, nrow(gexp.norm));
        names(foo) <- rownames(gexp.norm)
        foo[unique(unlist(boundgenes.pred))] <- 1  # vector of genes, zero=no CNV, 1=yes CNV

        ## vote
        vote <- rep(0, nrow(gexp.norm))
        names(vote) <- rownames(gexp.norm)
        lapply(boundgenes.pred, function(b) {
            vote[b] <<- vote[b] + 1
        })
        vote[bound.genes.old] <- 0 ## do not want to rediscover old bounds

        if(verbose) {
            cat(paste0('max vote:', max(vote), '\n'))
        }
        if(max(vote)==0) {
            if(verbose) {
                cat('Exiting; no new bound genes found.\n')
            }
            return() ## exit iteration, no more bound genes found
        }

        vote[vote > 0] <- 1
        mv <- 1 ## at least 1 vote
        cs <- 1
        bound.genes.cont <- rep(0, length(vote))
        names(bound.genes.cont) <- names(vote)
        for(i in 2:length(vote)) {
            ##if(vote[i] == mv & vote[i] == vote[i-1]) {
            if(vote[i] >= mv & vote[i] == vote[i-1]) {
                bound.genes.cont[i] <- cs
            } else {
                # increment chain number
                cs <- cs + 1
            }
        }
        tb <- table(bound.genes.cont)
        tbv <- as.vector(tb);
        names(tbv) <- names(tb)
        tbv <- tbv[-1] # get rid of 0 (the no cnv preds)
        ## tbv contains the genomic region names and their counts of genes.
        ## all detected deletions have fewer than 5 genes...reached the end
        tbv[tbv < min.num.genes] <- NA
        tbv <- na.omit(tbv)
            if(length(tbv)==0) {
                if(verbose) {
                    cat(paste0('Exiting; fewer than ', min.num.genes, ' new bound genes found.\n'))
                }
                return()
            }

        ## Have tbv now: just a vector of gene counts labeled by region name.

        cache_amp_del_genes_filename = paste0(amp_or_del, '.', length(tbv), "_regions", ".", length(bound.genes.cont), "_genes.dat")
        message("-writing table: ", cache_amp_del_genes_filename)
        write.table(bound.genes.cont, file=cache_amp_del_genes_filename, quote=F, sep="\t")


        ## test each of these highly confident deletions
        prob.info <- lapply(names(tbv), function(ti) {
            ## get the list of genes in that cluster
            bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == ti]
            if(verbose) {
                cat('GENES POTENTIALLY AFFECTED BY CNV: ')
                cat(bound.genes.new)
            }
            ## now that we have boundaries, run on all cells
            prob <- calcGexpCnvProb(gexp.norm.sub=gexp.norm[bound.genes.new, ], m=dev, verbose=verbose, ...)

            if(verbose) {
                cat("AMPLIFICATION PROBABILITY: ")
                cat(prob[[1]])
                cat("\n")
                cat("DELETION PROBABILITY: ")
                cat(prob[[2]])
                cat("\n")
            }

            return(list('ap'=prob[[1]], 'dp'=prob[[2]], 'bs'=bound.genes.new))
        })
        ##cat(prob.info)
        return(prob.info)
    }

    ## amp and del probabilities for each cell
    amp.prob.info <- getTbv(lapply(boundgenes.pred, function(x) x[['amp']]), amp_or_del='amp')
    del.prob.info <- getTbv(lapply(boundgenes.pred, function(x) x[['del']]), amp_or_del='del')
    ## above info has:
    ##   $ap: amplification cell probs
    ##   $dp: deletion cell probs
    ##   $bs:  the list of genes that define each genomic region.

    ## get tables of amplification region and deletion region probabilities for cells.  column = cell, row=region, val=pvalue
    del.prob <- do.call(rbind, lapply(seq_along(del.prob.info), function(i) del.prob.info[[i]]$dp))
    amp.prob <- do.call(rbind, lapply(seq_along(amp.prob.info), function(i) amp.prob.info[[i]]$ap))

    ## head(del.prob)
    ##  MGH31_A02 MGH31_A03 MGH31_A04 MGH31_A05 MGH31_A06 MGH31_A07 MGH31_A08
    ## [1,]   0.00325   0.00000    0.0295   0.00000     0.000   0.00000         0
    ## [2,]   1.00000   0.08275    0.0000   0.93625     0.905   0.98225         0
    ## [3,]   0.11725   0.00000    0.0000   1.00000     1.000   1.00000         0
    ##  MGH31_A10 MGH31_A11 MGH31_B01 MGH31_B02 MGH31_B04 MGH31_B05 MGH31_B06
    ## [1,]         0    0.0000    0.0005   0.00000         0     0.002   0.00000
    ## [2,]         1    0.7505    0.9990   0.02775         1     0.000   0.98775
    ## [3,]         1    0.0000    0.0000   0.99125         1     0.000   0.00050
    ## ...
    ## rows are the different amp or del regions.

    ## capture the lists of genes per region
    del.bound.genes.list <- lapply(seq_along(del.prob.info), function(i) del.prob.info[[i]]$bs)
    amp.bound.genes.list <- lapply(seq_along(amp.prob.info), function(i) amp.prob.info[[i]]$bs)

    ## merge amp and del prob and gene list info
    prob.bin <- prob <- rbind(del.prob, amp.prob)  #table of cell p-values;  prob.bin gets binarized below as yes/no/unsure.
    bound.genes.list <- c(del.bound.genes.list, amp.bound.genes.list) # lists of gene lists

    ## make the probs.bin binary reflecting having or not having amps, or NA for unsure regions.
    if(sum(prob < 0.25) > 0) {
        ## cells unlikely to have amp or del
        prob.bin[prob < 0.25] <- 0
    }
    if(sum(prob > 0.75) > 0) {
        ## cells likely to have amps or dels
        prob.bin[prob > 0.75] <- 1
    }
    if(sum(prob <= 0.75 & prob >= 0.25) > 0) {
        ## unsure regions for cells
        prob.bin[prob <= 0.75 & prob >= 0.25] <- NA
    }

    if(is.null(prob.bin)) {
        return()
    }

    if(ncol(prob.bin)>1) {
        ## have multiple cells w/ amps or dels
        dps <- rowSums(prob.bin, na.rm=TRUE) # number of cells w/ amp or del per region

        ## ####################################################################
        ## ** select the region that has the greatest number of cells w/ CNV **
        ## ####################################################################

        dpsi <- which(dps == max(dps))[1] ## pick one of the most clonal
        ## gets the region that has the largest number of cells containing the amp or del

        ## store gene list and cell probabilities.
        prob.fin <- t(as.matrix(prob[dpsi,]))
        bound.genes.new <- bound.genes.list[[dpsi]]
    } else {
        ## a single cell
        prob.fin <- prob
        bound.genes.new <- bound.genes.list
    }

    if(verbose) {
        print("CNV SNPS:")
        print(bound.genes.new)
        print(prob.fin)
    }

    ## need better threshold
    g1 <- colnames(prob.fin)[prob.fin > 0.75]
    if(verbose) {
        print("GROUP 1 (have CNV):")
        print(g1)
    }
    g2 <- colnames(prob.fin)[prob.fin <= 0.25]
    if(verbose) {
        print("GROUP 2 (dont have CNV):")
        print(g2)
    }
    ##clafProfile(r[, g1], n.sc[, g1], l, n.bulk)
    ##clafProfile(r[, g2], n.sc[, g2], l, n.bulk)

    ## record
    if(length(g1) > 0) {

        ## ##################################################
        ## Store the genes and cells w/ predicted CNV region.

        pred.genes.r[bound.genes.new, g1] <<- 1
        bound.genes.final[[length(bound.genes.final)+1]] <<- list(bound.genes.new)
    }

    ## all discovered
    bound.genes.old <<- unique(c(bound.genes.new, bound.genes.old))

    ## Recursion
    ## print('Recursion for Group1')
    if(length(g1)>=3) {
        tryCatch({
            calcGexpCnvBoundaries(gexp.norm.sub=gexp.norm[, g1], recursion_count = recursion_count + 1, G1orG2 = "G1")
        }, error = function(e) { cat(paste0("ERROR: ", e)) })
    }

    ##print('Recursion for Group2')
    if(length(g2)>=3) {
        tryCatch({
            calcGexpCnvBoundaries(gexp.norm.sub=gexp.norm[, g2], recursion_count = recursion_count + 1, G1orG2 = "G2")
        }, error = function(e) { cat(paste0("ERROR: ", e)) })
    }

}


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

calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)



## report the CNVs that were identified

print("CNV regions found:\n")
regions.genes <- range(genes[unlist(bound.genes.final)])
print(regions.genes)

