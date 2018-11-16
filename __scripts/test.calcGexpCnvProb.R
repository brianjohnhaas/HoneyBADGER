#!/usr/bin/env Rscript

library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
mvFit <- setMvFit(gexp.mats$gexp.norm)
gexp.norm = gexp.mats$gexp.norm
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
gexp <- gexp[vi,]


## smooth data
mu0 <- apply(gexp, 2, mean)  # mean expr in each cell
ng <- nrow(gexp)
sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))

data <- list(
    'K' = length(mu0), # number of cells
    'JJ' = nrow(gexp), # number of genes
    'gexp' = gexp,     # normalized gene expression
    'sigma0' = sigma0, # standard deviations per cell
    'mag0' = m         # expression deviation due to copy number change
)
modelFile <-  system.file("bug", "expressionModel.bug", package = "HoneyBADGER")
        
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




