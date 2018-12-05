#!/usr/bin/env Rscript

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' calcGexpCnvProb
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate posterior probability of CNVs using normalized expression data
#' calcGexpCnvProb(gexp.norm.sub=NULL, : Optional normalized gene expression matrix. If not provided, internal normalized gene expression matr
ix is used.
#'                   m=dev,            : Expected mean deviation due to copy number change (default: 0.15) 
#'                   region=NULL,      : Optional GenomicRanges region of interest such as expected CNV boundaries. (default: NULL)
#'                   quiet=TRUE,       : Boolean for whether to suppress progress display (default: TRUE)
#'                   verbose=FALSE)    : Verbosity (default: FALSE)
#'
#' dev = setGexpDev == Set needed absolute gene expression deviance to be able to distinguish neutral from amplified or deletion regions
#' mvFit = Linear regression model of variance in mean of N sampled genes (modeled between 5 to 100 random genes)
#'         Each cell gets a separate linear model.

#'
#' Seek to find the posterior distribution of [S^k, dd^k | gexp^k]
#' k = cell 
#' S^k = 1 if cell has cnv
#' S^k = 0 cell dose not has cnv
#' dd^k = 1 copy number gain 
#' dd^k = 0 copy number loss 
#' gexp^k = observed average normalized gene expression 
#' 
#'
#' gexp^k = theta_k + epsilon, 
#'     epsilon ~ Normal(0, sigma^2_k)
#'
#' sigma^2_k = cell specific variance to account for expected noise 
#'
#' S^k ~ Bernoulli(alpha_k) 
#' dd^k ~ Bernoulli(beta)
#'           aplpha, beta ~ uniform(0,1) uniform priors 



library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
mvFit <- setMvFit(gexp.mats$gexp.norm) # cell-specific prediction of var(gexp) given n_genes for a given cell.
gexp.norm = gexp.mats$gexp.norm # ref-subtracted gexp
genes = gexp.mats$genes
m=0.15
region=GenomicRanges::GRanges('chr10', IRanges::IRanges(0,1e9))
verbose=TRUE


gexp <- gexp.norm
gos <- genes[rownames(gexp.norm)]
fits <- mvFit[colnames(gexp.norm)]
quiet <- !verbose

# restrict to genes in region
overlap <- IRanges::findOverlaps(region, gos)
hit <- rep(FALSE, length(gos))
names(hit) <- names(gos)
hit[S4Vectors::subjectHits(overlap)] <- TRUE
vi <- hit
gexp <- gexp[vi,] #  gene expression restricted to genes in cnv block


mu0 <- apply(gexp, 2, mean)  # mean expr across all cells
ng <- nrow(gexp) # number of genes in cnv block
sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit']))) # predicted var(gexp) for each cell given number of genes.

data <- list(
    'K' = length(mu0), # number of cells
    'JJ' = nrow(gexp), # number of genes in cnv block
    'gexp' = gexp,     # normalized gene expression
    'sigma0' = sigma0, # standard deviations per cell
    'mag0' = m         # expression deviation due to copy number change
)
modelFile <-  system.file("bug", "expressionModel.bug", package = "HoneyBADGER")

## model {
##    ## Likelihood
##    ## Single cell
##    for( k in 1:K ) {
##
##        for ( j in 1:JJ ) {
##            ## Expression level
##            gexp[j, k] ~ dnorm(mu[k], sigma0[k])
##        }
##
##        # 0 if neutral
##        # + mag if amplification
##        # - mag if deletion
##        mu[k] <- 0 * ( 1 - S[k] ) +
##                -mag0 * ( S[k] * (1 - dd)) +
##                mag0 * ( S[k] * dd)
##
##        ## Cell level
##        S[k] ~ dbern(alpha[k]) # cnv or not
##        alpha[k] ~ dunif(0,1) # cell specific hyper-parameter prior to allow for better mixing
##    }
##    dd ~ dbern(beta) # direction of cnv
##    beta ~ dunif(0,1) # cell specific hyper-parameter
## }


inits <- list(
    list(S = rep(0, ncol(gexp)), dd = 0),
    list(S = rep(1, ncol(gexp)), dd = 0),
    list(S = rep(0, ncol(gexp)), dd = 1),
    list(S = rep(1, ncol(gexp)), dd = 1)
)
model <- rjags::jags.model(modelFile, data=data, inits=inits, n.chains=4, n.adapt=100, quiet=quiet)
update(model, 100, progress.bar=ifelse(quiet,"none","text"))

parameters <- c('S', 'dd')
library(rjags)
samples <- coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
samples <- do.call(rbind, samples) # combine chains

snpLike <- samples
v <- colnames(snpLike)
S <- snpLike[,grepl('S', v)]
dd <- snpLike[,grepl('dd', v)]
delcall <- apply(S*(1-dd), 2, mean)
ampcall <- apply(S*dd, 2, mean)
names(ampcall) <- names(delcall) <- colnames(gexp)

print(list('posterior probability of amplification'=ampcall,
           'posterior probability of deletion'=delcall))




