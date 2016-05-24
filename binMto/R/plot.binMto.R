"plot.binMto" <-
function(x, ltyH0=3, H0line=0, ltyCI=2, main=NULL, xlab=NULL,...)
{
esti <- rev(x$conf.int[,1])
lower <- rev(x$conf.int[,2])
upper <- rev(x$conf.int[,3])
compnames <- rownames(x$conf.int)

conf.level <- x$conf.level
alt <- x$alternative
k<- length(esti)
method=x$method

if(alt=="two.sided")
 {
  plotlim<-range(lower, upper, H0line)
 }

if(alt=="less")
 {
  plotlim<-range(esti, upper, H0line)
 }

if(alt=="greater")
 {
  plotlim<-range(esti, lower, H0line)
 }

if(x$adj=="Dunnett"){adjI<-"Dunnett-adjusted "}
if(x$adj=="Dunnettappr"){adjI<-"Approximative Dunnett-adjusted "}
if(x$adj=="Bonf"){adjI<-"Bonferroni-adjusted "}
if(x$adj=="Unadj"){adjI<-"Unadjusted "}

if (is.null(main))
 {main<-paste(adjI, method,"-intervals")}

if(is.null(xlab))
 {xlab<-paste(round(conf.level, digits=5)*100, "-% confidence intervals")}


plot.new()
args <- list(...) 

# adjust the left plot margin depending on the size of the names and the given par(cex.axis)

cex.axis <- args$cex.axis
if (!is.null(cex.axis)) 
par(cex.axis = cex.axis)
# save the old graphics settings
oldmai <- mymai <- par("mai")
ywidth <- max(strwidth(compnames, units = "inches", cex = par("cex.axis"))) * 1.2
if (mymai[2] < ywidth) 
 {mymai[2] <- ywidth}
par(mai = mymai, new = TRUE)

pr = rbind(c(plotlim[1], 1), c(plotlim[2], k))
pargs = c(list(x = pr[, 1]), list(y = pr[, 2]), type = "n", axes = FALSE, xlab = xlab, ylab = "", main = main, args)


do.call("plot", pargs)
axis(1, ...)

axis(2, 1:k, compnames[k:1], las = 1, ...)
box(...)

for (i in 1:k)
 {

segments(x0=lower[i], y0=i, x1=upper[i], y1=i, lty=ltyCI,...)

points(lower[i], i, pch = "(", ...)
points(upper[i], i, pch = ")", ...)
points(esti[i], i, pch = 19, ...)
       
 }
    

abline(v = H0line, lty = ltyH0, ...)

# restore the old plot settings
par(mai = oldmai)

}

