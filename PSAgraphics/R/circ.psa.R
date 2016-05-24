circ.psa <- function(response, treatment = NULL, strata = NULL, 
    summary = FALSE, statistic = "mean", trim = 0,
    revc = FALSE, confint = TRUE, sw = 0.4, ne = .5, inc = 0.25, 
    pw = 0.4, lab = TRUE, labcex = 1, xlab = NULL, ylab = NULL, 
    main = NULL){

# Given 'treatment' (0=Control, 1=Treatment), function calculates 'statistic' 
# for response within strata and treatment levels. Within strata, 'statistic' is 
# plotted as circles centered at (X,Y) = (C-stat,T-stat). The radii of the circles correspond to
# the sizes of the strata. The X and Y scores are used to generate 
# differences X - Y (default) or as Y - X (if revc = T; see below).
# The main 45 degree diagonal line (black) corresponds to X = Y, a reference diagonal.
# Mean of weighted differences corresponds to heavy (blue) dashed line, also with
# slope = 1. The marginal distribution of the weighted differences is
# plotted after making (parallel) projections to line segment at the lower
# part left of the figure.  Marginal distributions of Y & X are given as rug
# plots, top and right side; the weighted means of these marginal distributions
# are shown as vertical and horizontal lines. Note that the mean D corresponds
# to the intersection of means for X and Y.

# 'statistic' is any univariate function of the response data that returns a single 
#       value, eg "mean", "median" or 
#       "tr.mean" where user has predefined the latter as (say): 
#       tr.mean<-function(x){mean(x,trim=.1)}.  Default is 'mean'
# 'revc = TRUE' reverses X,Y axes;
# 'sw' extends axes on lower ends, effectively moving circles to lower left, or southwest.
# 'ne' extends ases on upper ends, effectively moving circles to upper right, or northeast.  
#       Making both sw and ne smaller moves circles farther apart, both larger moves them closer together.
# 'inc' controls circle sizes, but RELATIVE circle sizes are controlled via 'pw'. In general one wants
#       circle areas to appear SUBJECTIVELY to be sized in accordance w/ strata sizes.
# 'pw'  controls relative circle sizes.
# 'lab' labels circles with strata number.
# 'labcex' controls the size of the axes labels.

#Setting margins, forcing square plot to aid interpretations (saving current settings; restore settings at end).
op <- par(no.readonly = TRUE)
on.exit(par(op))     
par(mai = c(1, 1.7, 1, 1.7), pty="s")

#If "response" has three columns, treat as r, t, s.
if(dim(as.data.frame(response))[2]==3){treatment   <- response[,2]
                                       strata      <- response[,3]
                                       response    <- response[,1]}

tr.mean<-function(x){mean(x,trim=trim)}
if(!trim==0){statistic<-tr.mean}

if(!summary){                                       
n <- length(response)
nstrat <- dim(table(strata))
ct.means <- tapply(response,list(strata,treatment),statistic)
ncontrol <- as.data.frame(table(strata,treatment))[1:nstrat,3]
ntreat <- as.data.frame(table(strata,treatment))[(nstrat+1):(2*nstrat),3]
summary.strata <- cbind(ncontrol,ntreat,ct.means)}else{
summary.strata <- response
ct.means <- response[,3:4]
n <- sum(response[,1:2])
nstrat <- length(response[,1])}
wts<-rowSums(summary.strata[,1:2])

#Defining variables (reversed if called); 
#difference values for each x,y combination; number of strata
x <- ct.means[,1]
y <- ct.means[,2]
if(revc){x <- ct.means[,2]
         y <- ct.means[,1]}
d <- y - x  

#Upper and lower bound of plot; 
#Constants 'sw' and 'ne' may need modification to assure a pleasing appearance
xr <- range(x)
yr <- range(y)
min.xy <- min(xr[1], yr[1])
max.xy <- max(xr[2], yr[2])
lwb <- min.xy - 1.2 * ne * (max.xy - min.xy)
upb <- max.xy + .5 * sw * (max.xy - min.xy)

#Weighted mean difference; standard deviation of strata estimates 
diff.wtd <- sum(d * wts)/n

#Calculating the weighted st.dev using Conniffe's definition
#Note: 7/30/8: really calculating the standard error here, though I've called it st.dev.
if(!summary){
o<-order(treatment)
ord.strata<-strata[o]
nc<-table(treatment)[1]
nt<-table(treatment)[2]
ord.response<-response[o]

var.0<-tapply(ord.response[1:nc],ord.strata[1:nc],var)
ni.0<-table(ord.strata[1:nc])
frac.0<-var.0/ncontrol

ncp1<-nc+1
ncpnt<-nc+nt
var.1<-tapply(ord.response[ncp1:ncpnt],ord.strata[ncp1:ncpnt],var)
ni.1<-table(ord.strata[ncp1:ncpnt])
frac.1<-var.1/ntreat

se.wtd<-((sum(frac.0) + sum(frac.1))^.5)/nstrat
             
strat.labels<-sort(unique(strata))             
             }
                          
#If data are summarized, can't calculate se.wtd.
if(summary){
se.wtd<-NULL
strat.labels<-rownames(response)
}


#Finding radii and plotting circles
wtt <- wts/max(wts)
dfrng <- diff(c(lwb, upb))
radii <- (abs(dfrng/20)) * wtt
symbols(x, y, circles = radii^pw, inches = inc, xlim = c(lwb, upb), ylim = 
            c(lwb,upb), cex = 0.86, xlab = "", ylab = "", main = main, lwd=1.5)
if(lab){for(i in 1:nstrat){
        legend(x[i], y[i], strat.labels[i], bty="n", xjust=.71, yjust=.48, cex=labcex)}}

#Y=X line and parallel line representing weighted mean difference
abline(0, 1, lwd = 2)
abline(diff.wtd, 1, lwd = 3, col = 4, lty = 3)

#Vertical and horizontal lines through X Y means
wtss <- wts/n
mnx <- sum(x * wtss)
mny <- sum(y * wtss)
mnxy <- (mnx + mny)/2
abline(h = mny, v = mnx, lty = 3, col = 2)

#Dashed lines from circles to perpendicular cross line in lower left corner
kp <- (3*lwb+upb)/4
ex <- .008*(upb-lwb)
for(i in 1:nstrat){segments(kp - d[i]/2 + ex,kp + d[i]/2 + ex, x[i], y[i], lty = 2, lwd = .8)
                   points(kp - d[i]/2 + ex, kp + d[i]/2 + ex, pch = 3,lwd = 1.4, col = 4, cex = 1.7)
                  }
              
#Perpendicular cross segment: location fixed to lower quarter of plot
ext <- .025*(upb-lwb)
segments((lwb + upb)/2 + ext, lwb - ext,lwb - ext, (lwb + upb)/2 + ext, lwd = 2)

#Rug plots of centers of circles, ie strata means
rug(x, side = 3)
rug(y, side = 4) 

#Labels for axes
   dimnamez <- NULL
   dimnamezz <- unlist(dimnames(ct.means)[2])
   dimnamez <- unlist(dimnames(ct.means)[2])
   if(revc) dimnamez <- dimnamez[c(2, 1)]
   if(is.null(xlab)){title(xlab = unlist(dimnamez)[1])}else{title(xlab = xlab)}
   if(is.null(ylab)){title(ylab = unlist(dimnamez)[2])}else{title(ylab = ylab)}

Cname<-unlist(dimnamez)[1]
Tname<-unlist(dimnamez)[2]

C.wtd <- mnx; T.wtd <- mny

#if(!revc){diff.wtd<-(-1*diff.wtd)}

approx.t<-diff.wtd/se.wtd
df<-n-2*nstrat

if(!summary){colnames(summary.strata)<-c(paste("n.",dimnamezz[1],sep=""),paste("n.",dimnamezz[2],sep=""),
              paste("means.",dimnamezz[1],sep=""),paste("means.",dimnamezz[2],sep=""))}

#Note: changed out name of diff.wtd to ATE, did not change internal name diff.wtd.

out <- list(summary.strata, C.wtd, T.wtd, diff.wtd, se.wtd, approx.t = approx.t, df = df)
names(out)<-c("summary.strata",paste("wtd.Mn.",Cname,sep=""),paste("wtd.Mn.",Tname,sep=""),"ATE","se.wtd","approx.t","df") 
           
#Putting the CI below the cross plot
if(confint){
ci.diff<- -diff.wtd
#if(revc==FALSE)ci.diff<- -ci.diff
ci<-c(ci.diff-qt(.975,df)*se.wtd,ci.diff+qt(.975,df)*se.wtd)
ci.out<--ci[c(2,1)]
xl<-.75*lwb+.25*upb-ext+ci[1]/2
yl<-.75*lwb+.25*upb-ext-ci[1]/2
xu<-.75*lwb+.25*upb-ext+ci[2]/2
yu<-.75*lwb+.25*upb-ext-ci[2]/2
segments(xl,yl,xu,yu,lwd=5,col="green3")
segments(xl-.05*ext+.7*ext/2,yl+.05*ext+.7*ext/2,xl-.05*ext-.7*ext/2,yl+.05*ext-.7*ext/2,lwd=2,col="green3")
segments(xu+.05*ext+.7*ext/2,yu-.05*ext+.7*ext/2,xu+.05*ext-.7*ext/2,yu-.05*ext-.7*ext/2,lwd=2,col="green3")
out <- list(summary.strata, C.wtd, T.wtd, diff.wtd, se.wtd, approx.t, df, ci.out)
names(out)<-c("summary.strata",paste("wtd.Mn.",Cname,sep="",collapse=""),paste("wtd.Mn.",Tname,sep=""),"ATE","se.wtd","approx.t","df","CI.95")

           }           
                      

return(out)
}

