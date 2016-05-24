# hilbe.NBR2.F9.4.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Figure 9.4     user may amend as required for own model
# Table 9.32
#
load("c://source/azpro.RData")
attach(azpro)
myTable <- function(x) {
myDF <- data.frame( table(x) )
myDF$Prop <- prop.table( myDF$Freq )
myDF$CumProp <- cumsum( myDF$Prop )
myDF
}
by(los, procedure, myTable)
windows(record=TRUE)      # set so graphs both show
by(los, procedure, hist)  # produces histograms
# OR
library("ggplot2")
qplot(los,geom="histogram",facets=procedure~.)


