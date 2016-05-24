## ----short_example-------------------------------------------------------
library(qualvar)
set.seed(1)

# create a vector of frequencies for four categories
x <- rmultinom(1, 100, rep_len(0.25, 4))
x <- as.vector(t(x))
# compute the DM index
DM(x)

# Now let's compute DM indices for each row of a data frame where each column represents a category
df <- rmultinom(10, 100, rep_len(0.25, 4))
df <- as.data.frame(t(df))
names(df) <- c("a", "b", "c", "d")
apply(df, 1, DM)

## ----replication, results='asis'-----------------------------------------
library(DT)
data(wilcox1973)

wilcox1973$MDA <- apply(wilcox1973[,2:4], 1, MDA)
wilcox1973$DM <- apply(wilcox1973[,2:4], 1, DM)
wilcox1973$ADA <- apply(wilcox1973[,2:4], 1, ADA)
wilcox1973$VA <- apply(wilcox1973[,2:4], 1, VA)
wilcox1973$HREL <- apply(wilcox1973[,2:4], 1, HREL)
wilcox1973$B <- apply(wilcox1973[,2:4], 1, B)

wilcox1973[,5:10] <- apply(wilcox1973[,5:10], 2, function(x) round(x, digits = 3))

datatable(wilcox1973, options = list(pageLength = 60))


## ----correlation, fig.width=8, fig.height=8, fig.cap="Scatterplots, kernel density and correlation between all six indices."----
library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)

wilcox1973 %>%
  ggpairs(5:10) +
  theme_bw()


