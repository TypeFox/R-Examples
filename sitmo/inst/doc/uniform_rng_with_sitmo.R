## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----correlation---------------------------------------------------------
# Load packages
library(ggplot2)
library(sitmo)

# Number of Observations to Generate
n = 1e6

# Number of seeds to try (1 ... S)
nseeds = 30

# Storage for seed number and the correlation of the realizations between generators. 
cppdf = data.frame(s1=numeric(nseeds), s2=numeric(nseeds),
                   cor=numeric(nseeds), stringsAsFactors = F)

# Generate observations under the seeds
count = 0
for(i in 1:nseeds){
  for(j in i:nseeds){
    u1 = runif_sitmo(n, 0.0, 1.0, i)
    u2 = runif_sitmo(n, 0.0, 1.0, j)
    count = count + 1
    cppdf[count,] = c(i, j, cor(u1,u2))
  }
}

## ----corr_plot, fig.width = 7, fig.height = 4----------------------------
# Create Correlation Plot
ggplot(cppdf) + geom_tile(aes(x = s1, y = s2, fill=cor)) +
  xlab("Seed 1") + ylab("Seed 2") + 
  ggtitle("Correlations between seed realizations using `sitmo`")

