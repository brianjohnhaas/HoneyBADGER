#!/usr/bin/env Rscript

library(HoneyBADGER)
data(r)
data(cov.sc)


r_melt = melt(r)
colnames(r_melt) = c('chrpos', 'cell', 'alleleCov')

sc_cov_melt = melt(cov.sc)
colnames(sc_cov_melt) = c('chrpos', 'cell', 'totCov')

data = full_join(r_melt, sc_cov_melt, by=c('chrpos','cell'))

## include chr and postion separately.
data = data %>% separate(chrpos, ":", into=c('chr', 'pos'), remove=FALSE)


data = data %>% mutate(chr = ordered(chr, levels=1:22))
data$pos = as.numeric(data$pos)


## get chr bounds for plotting later.
chr_maxpos = data %>% group_by(chr) %>% summarize(maxpos = max(pos))
chr_maxpos$minpos = 1

## select snps where the fraction of cells containing an allelic site is between 0.1 and 0.9 of cells with read coverage.
het_snps = data %>% group_by(chrpos) %>%  mutate(l=sum(alleleCov>0), n.bulk=sum(totCov>0), E=l/n.bulk) %>% filter(E>0.1 & E<0.9) %>% pull('chrpos')
data = data %>% filter(chrpos %in% het_snps)



## filter snps, require at least 3 cells have coverage
min.cells = 3
snps_min_cells = data %>%  filter(totCov > 0) %>% group_by(chrpos) %>% tally() %>% filter(n>=min.cells) %>% pull(chrpos)
data = data %>% filter(chrpos %in% snps_min_cells)


## compute allele freq and mBAF
data = data %>% mutate(AF = alleleCov/totCov)

data = data %>% mutate(mBAF = pmax(AF, 1-AF))

## further restrict to het snps based on allele frequency data
het_snps = data %>% filter(AF > 0.25 & AF < 0.75) %>% group_by(chrpos) %>% tally() %>% filter(n>=3) %>% pull(chrpos)

data = data %>% filter(chrpos %in% het_snps)  ## note, want het_snps from normal!!!  but dont have it here.



data = data %>% filter(totCov > 0)  ## just here


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

    geom_point(aes(x=pos, y=cell, color=mBAF, size=totCov), alpha=0.6) + scale_radius()

plot(p)





