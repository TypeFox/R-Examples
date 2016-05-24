## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
library(alakazam)

# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Partitions the data based on the SAMPLE column
clones <- countClones(df, groups="SAMPLE")
head(clones, 5)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Partitions the data based on both the SAMPLE and ISOTYPE columns
# Weights the clone sizes by the DUPCOUNT column
clones <- countClones(df, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
head(clones, 5)

## ---- eval=TRUE, results='hide', warning=FALSE, fig.width=6, fig.height=4----
# Partitions the data on the SAMPLE column
# Calculates a 95% confidence interval via 200 bootstrap realizations
clones <- estimateAbundance(df, group="SAMPLE", ci=0.95, nboot=200)

## ---- eval=TRUE, warning=FALSE, fig.width=6, fig.height=4----------------
head(clones, 5)

# Plots a rank abundance curve of the relative clonal abundances
p1 <- plotAbundance(clones, legend_title="Sample")

## ---- eval=TRUE, results='hide'------------------------------------------
# Compare diversity curve across values in the "SAMPLE" column
# q ranges from 0 (min_q=0) to 32 (max_q=32) in 0.05 incriments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 2000 resampling realizations are performed (nboot=200)
sample_div <- rarefyDiversity(df, "SAMPLE", min_q=0, max_q=32, step_q=0.05, 
                                 ci=0.95, nboot=200)

# Compare diversity curve across values in the "ISOTYPE" column
# Analyse is restricted to ISOTYPE values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_div <- rarefyDiversity(df, "ISOTYPE", min_n=30, min_q=0, max_q=32, 
                                  step_q=0.05, ci=0.95, nboot=200)

## ---- eval=TRUE, fig.width=6, fig.height=4-------------------------------
# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
p2 <- plotDiversityCurve(sample_div, main_title=sample_main, 
                         legend_title="Sample", log_q=TRUE, log_d=TRUE)

# Plot isotype diversity using default set of Ig isotype colors
isotype_main <- paste0("Isotype diversity (n=", isotype_div@n, ")")
p3 <- plotDiversityCurve(isotype_div, colors=IG_COLORS, main_title=isotype_main, 
                         legend_title="Isotype", log_q=TRUE, log_d=TRUE)

## ---- eval=TRUE, results='hide'------------------------------------------
# Test diversity at q=0 (species richness) across values in the "SAMPLE" column
# 2000 bootstrap realizations are performed (nboot=200)
sample_test <- testDiversity(df, 0, "SAMPLE", nboot=200)

## ---- eval=TRUE----------------------------------------------------------
sample_test

## ---- eval=TRUE, results='hide'------------------------------------------
# Test diversity across values in the "ISOTYPE" column
# Analyse is restricted to ISOTYPE values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_test <- testDiversity(df, 2, "ISOTYPE", min_n=30, nboot=200)

## ---- eval=TRUE----------------------------------------------------------
isotype_test

