#######################################################################
##
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created : 2006-11-01
##  
##
## MODIFIED:
##   2008-04-04 : JW
##   - updated to anchors 3.0 syntax
##
################# making a bar chart...
cat("Repl from King and Wand (2007) Fig3\n")

data(mexchn)

fo <- list(self = xsayself ~ 1,
           vign = cbind(xsay5,xsay4,xsay3,xsay2,xsay1) ~ 1,
           cpolr= ~ china + age + male + educyrs)

z0 <- anchors( fo, mexchn, method="C", options=anchors.options(verbose=TRUE), combn=TRUE)
summary(z0)

# fname <- "fig_entropy_pol.ps"
# postscript( fname , height=8, width=8)

plot(z0, cex=.8 )
plot(z0, xy = c("minimum","estimated"))


# dev.off()
# system(paste("epstopdf",fname))

