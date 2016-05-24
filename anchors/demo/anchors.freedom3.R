#######################################################################
##
## Demo    : anchors.freedom3
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created : 2007-02-07
##  
## DESCRIPTION: Replication of Wand et al (2007)
##              Figure 2 (4 variations on C using 3 vignettes)
##
## MODIFIED:
##   2008-04-04 : JW
##   - updated to anchors 3.0 syntax
##   - estimates cpolr models separately by country
##
################# making a bar chart...
cat("Repl Wand et al (2007) Fig 2 histo with 3 vign\n")
data(freedom)

op <- par( no.readonly = TRUE)
par(mfrow=c(2,2))
ylim <- c(0,.5)

fo <- list(self = self ~ 1,
           vign = cbind(vign1,vign3,vign6) ~ 1,
           cpolr= ~ sex+age+educ )

## rank
ra1e <- anchors(fo, freedom, subset=country=="Eastasia", method = "C")
ra1o <- anchors(fo, freedom, subset=country=="Oceania"  , method = "C")

barplot(ra1e, ra1o, ties="omit"      , ylim=ylim, main = "Omit Tied Cases")
barplot(ra1e, ra1o, ties="uniform"   , ylim=ylim, main = "Uniform Allocation")
barplot(ra1e, ra1o, ties="cpolr"     , ylim=ylim, main = "Censored Ordered Probit Allocation")
barplot(ra1e, ra1o, ties="minentropy", ylim=ylim, main = "Minimum Entropy Allocation")

summary(ra1e, ties = c("omit","uniform","minentropy","cpolr"))
summary(ra1o, ties = c("omit","uniform","minentropy","cpolr"))

par(op)
