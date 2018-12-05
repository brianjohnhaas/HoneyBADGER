
## hacked together from other bits of hbadger code

plot_HMM_states <- function (states_matrix,
                             chrs = paste0("chr", c(1:22)),
                             region = NULL, 
                             window.size = 101,
                             zlim = c(-2, 2),
                             cellOrder = NULL,
                             widths = NULL,
                             recluster_cells = TRUE) {

    ## cluster cells
    if (recluster_cells) {
        d = dist(t(states_matrix))
        h = hclust(d, method='ward.D')
        states_matrix = states_matrix[,h$order]
    }
    
    require(biomaRt)
    mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",     dataset = 'hsapiens_gene_ensembl',     host = "jul2015.archive.ensembl.org")
    
    id="hgnc_symbol"
    gos <- getBM(values=rownames(states_matrix),attributes=c(id, "chromosome_name","start_position","end_position"),filters=c(id),mart=mart.obj)
    rownames(gos) <- make.unique(gos[[id]])
    gos <- gos[rownames(states_matrix),]
    
    if(nrow(gos)==0) {
        cat('ERROR! WRONG BIOMART GENE IDENTIFIER. Use ensembl_gene_id instead of hgnc_symbol? \n')
    }
    
    if(length(grep('chr', gos$chromosome_name))==0) {
        gos$chromosome_name <- paste0('chr', gos$chromosome_name)
    }
    gos <- na.omit(gos)
    gs <- with(gos, GenomicRanges::GRanges(as.character(chromosome_name),
                                           IRanges::IRanges(as.numeric(as.character(start_position)),
                                                            as.numeric(as.character(end_position)))))
    names(gs) <- rownames(gos)
    gvi <- intersect(rownames(states_matrix), names(gs))
    genes <- gs[gvi]
    
    ## remove genes with no position information
    states_matrix <- states_matrix[gvi,]
    
    ## cluster cells:
    
    
    
    genes <- genes[rownames(states_matrix)]

    ## restrict to region?
    if (!is.null(region)) {
        overlap <- IRanges::findOverlaps(region, genes)
        hit <- rep(FALSE, length(genes))
        hit[S4Vectors::subjectHits(overlap)] <- TRUE
        if (sum(hit) < 10) {
            cat(paste0("WARNING! ONLY ", sum(hit), " GENES IN REGION! \n"))
        }
        vi <- hit
        states_matrix <- states_matrix[vi, ]
        genes <- genes[rownames(states_matrix)]
    }

    
    
    gos <- as.data.frame(genes)
    rownames(gos) <- names(genes)
    mat <- states_matrix
    tl <- tapply(1:nrow(gos), as.factor(gos$seqnames), function(ii) {
        na.omit(mat[rownames(gos)[ii[order((gos[ii, ]$start + 
            gos[ii, ]$end)/2, decreasing = F)]], , drop = FALSE])
    })
    tl <- tl[chrs]
    if (!is.null(region)) {
        vi <- unlist(lapply(tl, function(x) {
            if (is.null(x)) {
                return(FALSE)
            }
            else {
                return(nrow(x) > 1)
            }
        }))
        tl <- tl[vi]
    }
    if (is.null(widths)) {
        widths <- rep(1, length(tl))
    } else if (widths[1] == "set") {
        widths <- sapply(tl, nrow)
        widths <- widths/max(widths) * 100
    }
    l <- layout(matrix(seq(1, length(tl)), 1, length(tl), byrow = TRUE), 
        widths = widths)
    if (is.null(cellOrder)) {
        cellOrder <- colnames(states_matrix)
    } else if (cellOrder[1] == "set") {
        avgd <- do.call(rbind, lapply(names(tl), function(nam) {
            d <- tl[[nam]]
            d <- apply(d, 2, caTools::runmean, k = window.size, 
                align = "center")
            d
        }))
        hc <- hclust(dist(t(avgd)))
        cellOrder <- hc$order
    }
    tlsub <- tl
    pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 
        "RdBu")))(256)
    tlsmooth <- lapply(names(tlsub), function(nam) {
        d <- tlsub[[nam]]
        d <- apply(d, 2, caTools::runmean, k = window.size, align = "center")
        d[d < zlim[1]] <- zlim[1]
        d[d > zlim[2]] <- zlim[2]
        d <- d[, cellOrder]
        par(mar = c(0.5, 0.2, 3, 0.2), mgp = c(2, 0.65, 0), cex = 0.8)
        image(seq_len(nrow(d)), seq_len(ncol(d)), d, col = pcol, 
            zlim = zlim, xlab = "", ylab = "", axes = FALSE, 
            main = nam)
        box()
        return(d)
    })
    #return(list(tl = tl, tlsmooth = tlsmooth, cellOrder = cellOrder))
}
