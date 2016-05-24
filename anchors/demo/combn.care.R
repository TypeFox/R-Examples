#######################################################################
##
## Demo    : combn.care.R
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created : 2006-11-03
##  
##
## MODIFIED:
##   2008-04-04 : JW
##   - updated to anchors 3.0 syntax
##
################# making a bar chart...
cat("Repl from King and Wand (2007) Fig2\n")

data(selfcare)
dim(selfcare)

selfcare$single  <- (selfcare$q1008==1)*1
selfcare$married <- (selfcare$q1008==2)*1
selfcare$height  <- selfcare$q1006/100

fo <- list(self = q2020 ~ 1,
           vign = cbind( q2101d,q2119d,q2113d,q2105d,q2117d ) ~ 1,
           cpolr= ~ q1001+q1002+q1004+q1010+married)

z0 <- anchors( fo,  data=selfcare, method="C", combn=TRUE)
summary(z0)


#fname <- "fig_entropy_self.ps"
#postscript( fname , height=8, width=8)

plot(z0, cex=.8 )
plot(z0, xy = c("minimum","estimated"))

# dev.off()
# system(paste("epstopdf",fname))

