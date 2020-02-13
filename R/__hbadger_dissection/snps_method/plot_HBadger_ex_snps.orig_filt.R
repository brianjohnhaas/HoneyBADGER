#!/usr/bin/env Rscript

library(HoneyBADGER)
data(r)
data(cov.sc)



hb <- new('HoneyBADGER', name='MGH31')
hb$setAlleleMats(r.init=r, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=detectCores())

pdf("orig_snp_plot.pdf", width=8, height=6)
hb$plotAlleleProfile()
dev.off()


pdf("new_snp_plot.pdf", width=8, height=6)
r.maf = hb$r.maf
n.sc = hb$n.sc

r_melt = melt(r.maf)
colnames(r_melt) = c('chrpos', 'cell', 'alleleCov')

sc_cov_melt = melt(n.sc)
colnames(sc_cov_melt) = c('chrpos', 'cell', 'totCov')

data = full_join(r_melt, sc_cov_melt, by=c('chrpos','cell'))

## include chr and postion separately.
data = data %>% separate(chrpos, ":", into=c('chr', 'pos', 'pos2'), remove=FALSE)
data$chr = sub(x=data$chr, pattern="chr", replacement="")


data = data %>% mutate(chr = ordered(chr, levels=1:22))
data$pos = as.numeric(data$pos)


## get chr bounds for plotting later.
chr_maxpos = data %>% group_by(chr) %>% summarize(maxpos = max(pos))
chr_maxpos$minpos = 1


## compute allele freq and mBAF
data = data %>% mutate(AF = alleleCov/totCov)

data = data %>% mutate(mBAF = pmax(AF, 1-AF))

## restrict to het snps
#het_snps = data %>% filter(AF > 0.25 & AF < 0.75) %>% group_by(chrpos) %>% tally() %>% filter(n>=3) %>% pull(chrpos)

#data = data %>% filter(chrpos %in% het_snps)  ## note, want het_snps from normal!!!  but dont have it here.



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

    geom_point(aes(x=pos, y=cell, color=mBAF, size=totCov))

plot(p)


dev.off()





