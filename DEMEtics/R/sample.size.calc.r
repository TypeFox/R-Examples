# Function called within allelefreq.r to calculate locus-wise sample sizes
sample.size.calc <- function(y2){

y2.without.zeros <- y2[y2$fragment.length!=0,]
x <- (table(y2.without.zeros[,2]))/2

tab <- cbind(as.data.frame.table(x, row.names = NULL, optional = FALSE),as.character((y2$locus)[1]))
# The sample sizes for the actual locus and the corresponding
# populations are combined in a data frame.
}
