#######################################################################
# Note that this note can directly be run in R.
# Version: GeneCycle 1.0.3 (February 2008)
#######################################################################

#
# EXAMPLE SESSION FOR CELL CYCLE ANALYSIS
# 

# for details see:
#
# Wichert, S., K. Fokianos, and K. Strimmer. 2004.
# Identifying periodically expressed transcripts in microarray
# time series data. Bioinformatics 20:4-20 



# load GeneCycle library
library("GeneCycle")

#######################################################################

# THE DATA:

# the normalized data need to be ready in time series format, i.e. in
# a matrix where each *column* corresponds to a gene, and where the
# *rows* correspond to the individual measurements (time points).

# Example: the Caulobacter data set
data(caulobacter)

# how many samples (11) and how many genes (1444)?
dim(caulobacter)
summary(caulobacter)
get.time.repeats(caulobacter)

# plot first nine time series
plot(caulobacter, 1:9)



#######################################################################

# IDENTIFYING PERIODICALLY EXPRESSED GENES:

# A statistical test developed by Fisher is used to detect 
# periodically expressed genes, and the average periodogram
# is used to visualize the dominant frequencies

# compute and plot average periodogram
avgp.caulobacter <- avgp(caulobacter, "Caulobacter")
avgp.caulobacter

# p-values from Fisher's g test
pval.caulobacter <- fisher.g.test(caulobacter)
pval.caulobacter


#######################################################################
# multiple testing 

# proportion of null p-values for different methods
pval.estimate.eta0(pval.caulobacter, method="conservative")
pval.estimate.eta0(pval.caulobacter, method="adaptive")
pval.estimate.eta0(pval.caulobacter, method="bootstrap")
pval.estimate.eta0(pval.caulobacter, method="smoother")


# (local) false discovery rates (using bootstrap)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
sum(fdr.out$qval < 0.05) # 52
sum(fdr.out$lfdr < 0.2) # 94


