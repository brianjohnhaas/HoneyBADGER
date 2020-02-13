#!/usr/bin/env Rscript

library(HoneyBADGER)

library(HoneyBADGER)
data(gexp)
data(ref)
require(biomaRt)

mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
hb <- new('HoneyBADGER', name='MGH31')

hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)

gexp.norm = hb$gexp.norm

genes = hb$genes
gos <- as.data.frame(genes)
rownames(gos) <- names(genes)
mat <- gexp.norm
## organize into chromosomes
##   gos$seqnames = chromosome names


gos$genename = rownames(gos)
gexp.melt = melt(gexp.norm)
colnames(gexp.melt) = c('genename', 'cell', 'exp.norm')

data = left_join(gexp.melt, gos, key='genename')
data$chr = str_replace(string=data$seqnames, pattern="chr", replacement="")

data = data %>% filter(chr %in% 1:22)


data = data %>% mutate(chr = ordered(chr, levels=1:22))

data$pos = as.numeric( (data$start + data$end)/2)


## get chr bounds for plotting later.
chr_maxpos = data %>% group_by(chr) %>% summarize(maxpos = max(pos))
chr_maxpos$minpos = 1

## set up base plot
## define chr boundaries based on max coordinates for now.


p = ggplot(data=data) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
          ) +

    geom_vline(data=chr_maxpos, aes(xintercept=minpos), color=NA) +
    geom_vline(data=chr_maxpos, aes(xintercept=maxpos), color=NA) +

    geom_point(aes(x=pos, y=cell, color=exp.norm), alpha=0.6) +

    scale_colour_gradient2(low = "blue", mid = "white",
                           high = "red", midpoint = 0, space = "Lab",
                           na.value = "grey50", guide = "colourbar", aesthetics = "colour")

plot(p)





