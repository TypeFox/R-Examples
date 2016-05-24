granova.ds <- function(data, revc = FALSE, sw = 0.4, ne=.5, ptpch=c(19,3), ptcex=c(.8,1.4), labcex = 1, ident = FALSE, 
            colors = c(1,2,1,4,2,'green3'), pt.lab = NULL,
            xlab = NULL, ylab = NULL, main = NULL, sub = NULL, par.orig = TRUE){

# Given dependent sample data (X,Y) data points are plotted (blue, small closed circles).
# The X and Y scores are used to generate differences X - Y (default) or as Y - X (if revc = T; see below).
# The main 45 degree diagonal line (black) corresponds to X = Y, a reference (or identity) diagonal.
# Mean of differences corresponds to heavy (red) dashed line, parallel to identity line.
# The marginal distribution of the differences is plotted after making (parallel) projections to line segment
# at the lower left corner of the figure. Marginal distributions of Y & X are given as rug
# plots, top and right side; the means of these marginal distributions are shown as vertical and horizontal lines.
# Note that the mean D corresponds to the intersection of means for X and Y.

# 'xdata' should be an n X 2 dataframe.
# 'revc = TRUE' reverses X,Y axes;
# 'sw' extends axes on lower ends, effectively moving circles to lower left, or southwest.
# 'ne' extends ases on upper ends, effectively moving circles to upper right, or northeast.
# Making both sw and ne smaller moves points farther apart, both larger moves them closer together.
# 'ptpch' controls the X,Y point and marginal dotplot pch.
# 'ptcex' controls the X,Y point and marginal dotplot cex.
# 'labcex' controls the size of the axes labels.
# 'ident' allows user to identify specific points on the plot
# 'colors' is a vector of colors found in the plot:  round points, dashed mean lines, light diagonal dashed lines, 
#   marginal plot points, dashed diagonal difference mean line, confidence interval. 
# 'pt.lab' allows user to provide labels for points, else the rownames of xdata are used (if defined), or if not labels are 1:n.
# 'xlab', 'ylab', 'main', 'sub' are optional, with axes labels taken from column names otherwise.
# 'par.orig' returns par to original settings; if trellis plots are desired it may be best to set to 'FALSE'.
# Please REPORT experiences/suggestions for editing: 'rmpruzek@yahoo.com' or 'james.helmreich@marist.edu' ..w/ thanks!

#Setting margins, forcing square plot to aid interpretations (saving current settings; restore settings at end).
op <- par(no.readonly = TRUE)
if(par.orig){on.exit(par(op))}
par(mai = c(1, 1.7, 1, 1.7), pty="s")

xdata <- data

col.dim <- dim(as.matrix(xdata))[2]
if(!col.dim == 2){stop("Input data must be a n X 2 dataframe or matrix.")}

n <- dim(as.matrix(xdata))[1]

#Defining variables (reversed if called);
#difference values for each x,y combination; number of strata
x <- xdata[,1]
y <- xdata[,2]
if(revc){x <- xdata[,2]
         y <- xdata[,1]}
d <- x - y

#Upper and lower bound of plot;
#Constants 'sw' and 'ne' may need modification to assure a pleasing appearance
xr <- range(x)
yr <- range(y)
min.xy <- min(xr[1], yr[1])
max.xy <- max(xr[2], yr[2])
lwb <- min.xy - 1.2 * ne * (max.xy - min.xy)
upb <- max.xy + .5 * sw * (max.xy - min.xy)

#Weighted mean difference; standard deviation of strata estimates
mn.diff <- mean(d)
sd.diff <- sd(d)

#Main Plot
if(is.null(main)){main<-"Dependent Sample Assessment Plot"}
plot(x, y, xlim = c(lwb, upb), ylim = c(lwb,upb), pch = ptpch[1], col = colors[1],
        cex = ptcex[1], xlab = "", ylab = "", main = main, sub = sub, lwd=1.5)

#Y=X line
abline(0, 1, lwd = 2)

#Vertical and horizontal lines through X Y means
mnx <- mean(x)
mny <- mean(y)
mnxy <- (mnx + mny)/2
abline(h = mny, v = mnx, lty = 3, col = colors[2])

#Dashed lines from points to perpendicular cross line in lower left corner, and '+' signs
kp <- (3*lwb+upb)/4
ex <- .008*(upb-lwb)
for(i in 1:n){segments(kp + d[i]/2 + ex,kp - d[i]/2 + ex,x[i],y[i],lty = 3, lwd = .8, col = colors[3])
             points(kp + d[i]/2 + ex, kp - d[i]/2 + ex, pch= ptpch[2], lwd=1.4, col = colors[4], cex=ptcex[2])}

#Perpendicular cross segment: location fixed to lower quarter of plot
ext <- .025*(upb-lwb)
segments((lwb+upb)/2+ext,lwb-ext,lwb-ext,(lwb+upb)/2+ext,lwd=2)

#Rug plots of data
rug(x, side = 3, lwd = 0.6)
rug(y, side = 4, lwd = 0.6)

#Default labels from axes taken from column names, or if none and none given explicitly are X, Y.
if(is.null(colnames(xdata))){colnames(xdata)<-c("X","Y")}
if(is.null(xlab)){xlab<-colnames(xdata)[1]}
if(is.null(ylab)){ylab<-colnames(xdata)[2]}

#Labels for axes
if(revc){xla <- ylab
         yla <- xlab
         xlab <- xla
         ylab <- yla
        }
title(xlab = xlab)
title(ylab = ylab)

X.mn <- mnx
Y.mn <- mny

#Statistics for output Summary
es.d<-(mn.diff)/sd.diff  #standardized effect size (ES(D))
ttest<-t.test(x,y,paired=TRUE)
pval<-ttest$p.value
ci<-ttest$conf.int
tstat <- ttest$statistic
df<-ttest$parameter
se.diff<-sd.diff/sqrt(n-1)
llm.ci<-ci[1]
ulm.ci<-ci[2]
df.t<-df
r.x.y<-cor(x,y)
r.xplsy.d<-cor(x+y,d)
legend(x="topleft",legend = list(paste("Mean diff. =", round(mn.diff,2)),"95% CI",paste("t =", round(tstat,2))),
        pch =c("_","_",""),col=c(colors[c(5,6)],0), bty = "n")
Summary<-(matrix(c(n, X.mn, Y.mn, mn.diff ,sd.diff ,es.d ,r.x.y ,r.xplsy.d ,llm.ci ,ulm.ci,tstat ,df.t ,pval), 13,1))
dimnames(Summary) <- list(c("n", "mean(x)","mean(y)","mean(D=x-y)","SD(D)","ES(D)","r(x,y)", "r(x+y,d)",
"LL 95%CI", "UL 95%CI",  "t(D-bar)","df.t","pval.t"), "Summary Stats")
Summary<-round(Summary,3)

#Putting the CI below the cross plot
xl<-.75*lwb+.25*upb-ext+ci[1]/2
yl<-.75*lwb+.25*upb-ext-ci[1]/2
xu<-.75*lwb+.25*upb-ext+ci[2]/2
yu<-.75*lwb+.25*upb-ext-ci[2]/2
segments(xl,yl,xu,yu,lwd=5,col=colors[6])
segments(xl-.05*ext+.7*ext/2,yl+.05*ext+.7*ext/2,xl-.05*ext-.7*ext/2,yl+.05*ext-.7*ext/2,lwd=2,col=colors[6])
segments(xu+.05*ext+.7*ext/2,yu-.05*ext+.7*ext/2,xu+.05*ext-.7*ext/2,yu-.05*ext-.7*ext/2,lwd=2,col=colors[6])

#Red dotted line representing weighted mean difference;
#Redraw points to be sure they are on top
abline(-mn.diff, 1, lwd = 3, col = colors[5], lty = 3)
points(x, y, pch = ptpch[1], col = colors[1], cex = ptcex[1])

if(ident){
         if(is.null(pt.lab) & !is.null(rownames(xdata))){pt.lab<-rownames(xdata)}
         if(is.null(pt.lab) & is.null(rownames(xdata))){pt.lab<-1:n}
         identify(x,y,labels = pt.lab)
         }

return(Summary)
}
