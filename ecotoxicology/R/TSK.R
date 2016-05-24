
# Trimmed Spearman-Karber method, as per Hamilton and EPA.
# Written by Brenton R. Stone, last revised May 18 2010
# Translated to R by Jose Gama 2013

#' Inverse error function
#' @description Returns the inverse error function
#' @param x numeric vector
#' @return the inverse error function
#' @references Abramowitz and Stegun 29.2.29
#' \url{http://stat.ethz.ch/R-manual/R-devel/library/stats/html/Normal.html}
#' @author Jose Gama
#' @examples
#' erfinv(1:10)
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)# inverse error function (see Abramowitz and Stegun 29.2.29) http://stat.ethz.ch/R-manual/R-devel/library/stats/html/Normal.html

#Functions to calculate the variance, as per the appendix of Hamilton.
V1<-function(data,A,L,U){
((data[L+1,'x'] - data[L,'x']) * (data[L+1,'p'] - A)^2 / (data[L+1,'p'] - data[L,'p'])^2)^2 * data[L,'p'] * (1 - data[L,'p']) / data[L,'n']
}
V2<-function(data,A,L,U){
((data[L,'x'] - data[L+2,'x']) + (data[L+1,'x'] - data[L,'x']) * (A - data[L,'p'])^2 / (data[L+1,'p'] - data[L,'p'])^2)^2 * data[L+1,'p'] *
(1-data[L+1,'p']) / data[L+1,'n']
}
V3<-function(data,A,L,U){
v3 = (data[(L+1):(U-3),'x'] - data[(L+3):(U-1),'x'])^2 * data[(L+2):(U-2),'p'] * (1 - data[(L+2):(U-2),'p']) / data[(L+2):(U-2),'n']
sum(v3)
}
V4<-function(data,A,L,U){
((data[U-2,'x'] - data[U,'x']) + (data[U,'x'] - data[U-1,'x']) * (data[U,'p'] - 1 + A)^2 / (data[U,'p'] - data[U-1,'p'])^2)^2 * data[U-1,'p'] * 
(1-data[U-1,'p']) / data[U-1,'n']
}
V5<-function(data,A,L,U){
((data[U,'x'] - data[U-1,'x']) * (1 - A - data[U-1,'p'])^2 / (data[U,'p'] - data[U-1,'p'])^2)^2 * data[U,'p'] * (1-data[U,'p']) / data[U,'n']
}
V6<-function(data,A,L,U){
((data[U,'x'] - data[L+1,'x']) * (1 - A - data[U,'p'])^2 / (data[U,'p'] - data[L+1,'p'])^2 - (data[L+1,'x'] - data[L,'x']) * (A - data[L,'p'])^2 /
(data[L+1,'p'] - data[L,'p'])^2 + (data[L,'x'] - data[U,'x']))^2 * data[L+1,'p'] * (1 - data[L+1,'p']) / data[L+1,'n']
}


#' Trimmed Spearman-Karber method, as per Hamilton and EPA
#' @description Returns the Trimmed Spearman-Karber (TSK) method, as per Hamilton and EPA
#' @param x numeric vector
#' @param r numeric vector
#' @param n numeric vector
#' @param A numeric vector
#' @param conf numeric vector
#' @return mu=mu,gsd=gsd,left=left,right=right
#' @references 
#' Hamilton,M.A.,Russo,R.L.,Thurston,R.V.,1977.
#' Trimmed Spearman–Karber method for estimating median lethal concentrations.
#' Environ. Sci. Tech. 11,714–719.
#' @author Jose Gama
#' @examples
#' x<-c(15.54,20.47,27.92,35.98,55.52)
#' n1<-c(20,20,20,19,20)
#' r<-c(0,0,0,5.26,100)/100*n1
#' n<-c(20,20,20,19,20)
#' TSK(x,r,n)
TSK<-function(x,r,n,A=0,conf=0.95)
 {
 
 #x=c(0,6.25,12.5,25,50,100);r=c(0,0,1,0,0,16);n=20;A=0.2;conf=erf(2/sqrt(2))
 #x=x;r=p;n=n;A=0.05;conf=conf
 #x=c(15.54,20.47,27.92,35.98,55.52);r=c(0,0,0,5.26,100)/100*n1;n=c(20,20,20,19,20);A=0;conf=0.95
 
 #---CHECK THAT THE INPUT IS VALID------------------------------------------
 if (any(is.na(x))) stop('x must not contain NA')
 if (any(is.na(r))) stop('r must not contain NA')
 if (any(is.na(n))) stop('n must not contain NA')
 if (is.na(A)) stop('A must not be NA')
 if (is.na(conf)) stop('conf must not be NA')
 if (!is.numeric(x)) stop('x must be numeric')
 if (!is.numeric(r)) stop('r must be numeric')
 if (!is.numeric(n)) stop('n must be numeric')
 if (!is.numeric(A)) stop('A must be numeric')
 if (!is.numeric(conf)) stop('conf must be numeric')
# Check if xr, r, and nr match
N<-length(x)
if (N != length(r)) stop('Different numbers of doses and responses were given')
if (length(n)== 1) n = rep(n,N)
if (N != length(n)) stop('Different numbers of doses and subject groups were given')
# Check if x is in the right range and order
if (any(x < 0)) stop('All doses must be positive or zero. Please input concentrations, not logs of concentrations.')
if (!all(x==x[order(x)])) stop('Data is not in order of increasing dose.')
# Check if r is in the right range
if (any(r < 0)) stop('Responses must be positive or zero')
# Check if n is nonnegative
if (any(n < 0)) stop('Population sizes should be positive.')
# Check if A is in the right range
if ((A<0) | (A>=0.5)) stop ('The trim A is not a valid value.')
# Check if conf is in the right range
if ((conf>=1) | (conf<=0.5)) stop('The confidence is not a valid value.')
#---MASSAGE THE DATA-------------------------------------------------------
# If first concentration is 0, scale everything else to account for that. 
# Also, transform to log scale.
titl = 'TSK'
if (x[1] == 0){
#mydata <- cbind('x'=log10(x[-1]), 'p'=0, 'n'=0)
if (r[1]==0) {p <- r[-1]/n[-1];n <- n[-1];titl <- sprintf('Responses, trim of %d%%', A*100)
} else {
warning('Control dose has a non-zero response. Scaling everything else to correct for that.')
p <- (r[-1]-r[1])/(n[-1]-n[1])
n <- n[-1]-r[1]
titl = sprintf('Responses scaled for control, trim of %d%%', A*100)
}
mydata <- cbind('x'=log10(x[-1]), 'p'=p, 'n'=n)
N = length(n)
} else mydata <- cbind('x'=log10(x), 'p'=r/n, 'n'=n)

# Smooth by fixing nondecreasingness, as per the first step in Hamilton.
# There's probably a more efficient way to do this but this is what's in
# Hamilton. It does not seem to add much to the execution time.

datasmooth <- mydata
if (any(datasmooth[-N,'p'] > datasmooth[-1,'p'])) warning('Responses are not monotonically increasing with dose. Smoothing the data by averaging adjacent nonincreasing responses.')
while (any(datasmooth[-N,'p'] > datasmooth[-1,'p']))
for (i in 1:(N-1))
if (datasmooth[i,'p'] > datasmooth[i+1,'p'])
{
pave <- (datasmooth[i,'p'] * datasmooth[i,'n'] + datasmooth[i+1,'p'] * datasmooth[i+1,'n']) / (datasmooth[i,'n'] + datasmooth[i+1,'n'])
datasmooth[i,'p'] <- pave
datasmooth[i+1,'p'] <- pave
}

#---SCALE AND TRIM---------------------------------------------------------
# As per steps 2 and 3 in Hamilton

datascale <- datasmooth
if (A != 0) datascale[,'p'] <- (datasmooth[,'p'] - A) / (1 - 2 * A)

if (datascale[1,'p'] > 0 | datascale[N,'p'] < 1){
    #Cannot do it with this trim. Determine smallest trim that analysis can
    #be done with.
    Amaybe <- max(datasmooth[1,'p'], 1 - datasmooth[N,'p'])
    stop(paste('Responses do not cover the space from A to 1-A. Consider using a value of A that is ',Amaybe,' or larger.',sep=''))
}
# Linearly interpolate points where the lines meet the trim.
# Check if any points lie inside the trim
#   x: [1.4459 1.5561 1.7444]
#    p: [0 0.0526 1]
keepers = ((datascale[,'p'] > 0) & (datascale[,'p'] < 1))
if (any(keepers)){
wkeep<-which(keepers)
    Uscale = wkeep[length(wkeep)]
    Lscale = wkeep[1]
    trimx = datascale[keepers,'x']
    trimp = datascale[keepers,'p']
    if (datascale[Lscale-1,'p'] != 0){
        i=Lscale
        xhead = approx(c(datascale[i-1,'p'], datascale[i,'p']), c(datascale[i-1,'x'], datascale[i,'x']),0) #interp1q([datascale.p(i-1); datascale.p(i)], [datascale.x(i-1); datascale.x(i)],0)
    } else xhead = datascale[Lscale-1,'x']
    if (datascale[Uscale+1,'p'] != 1){
        i = Uscale
        xtail = approx(c(datascale[i,'p'], datascale[i+1,'p']), c(datascale[i,'x'], datascale[i+1,'x']),1)
        #interp1q([datascale.p(i); datascale.p(i+1)], [datascale.x(i); datascale.x(i+1)],1);

    } else xtail = datascale[Uscale+1,'x']
    datatrim = list('x'=c(ifelse(is.list(xhead),xhead$y,xhead), trimx, ifelse(is.list(xtail),xtail$y,xtail)),'p'=c(0, trimp, 1))
} else # When all datapoints lie outside the trim
{
    i = which(datascale[,'p'] > 0)[1] -1
    xhead = approx(c(datascale[i,'p'], datascale[i+1,'p']), c(datascale[i,'x'], datascale[i+1,'x']),0)
    #interp1q([datascale.p(i); datascale.p(i+1)], [datascale.x(i); datascale.x(i+1)],0);
    xtail = approx(c(datascale[i,'p'], datascale[i+1,'p']), c(datascale[i,'x'], datascale[i+1,'x']),1)
    #interp1q([datascale.p(i); datascale.p(i+1)], [datascale.x(i); datascale.x(i+1)],1);
    datatrim = list('x'=c(ifelse(is.list(xhead),xhead$y,xhead), ifelse(is.list(xtail),xtail$y,xtail)),'p'=c(0, 1))
}

#---FIND THE MEAN----------------------------------------------------------
# Step 4 in Hamilton
trimN = length(datatrim$x)
mids  = (datatrim$x[1:trimN-1] + datatrim$x[2:trimN]) / 2
delp  = datatrim$p[2:trimN] - datatrim$p[1:trimN-1]
logmu = sum(mids * delp)
mu    = 10^logmu

#---FIND THE VARIANCE------------------------------------------------------
# Appendix in Hamilton
w<-which(datasmooth[,'p'] <= A)
L = w[length(w)]#find(datasmooth.p <= A, 1, 'last');
U = which(datasmooth[,'p'] >= 1-A)[1]#find(datasmooth.p >= 1-A, 1, 'first');
s = U - L
if (s==0){
        #If this case occurs it should have failed before now but just in
        #case...
        Var = NaN
} else if (s==1){
        Var = (datasmooth[U,'x'] - datasmooth[L,'x'])^2 * 
            ((0.5-datasmooth[U,'p'])^2 / (datasmooth[U,'p']-datasmooth[L,'p'])^4 * 
            datasmooth[L,'p'] * (1-datasmooth[L,'p']) / datasmooth[L,'n'] 
            + (0.5-datasmooth[L,'p'])^2 / (datasmooth[U,'p']-datasmooth[L,'p'])^4 * 
            datasmooth[U,'p'] * (1-datasmooth[U,'p']) / datasmooth[U,'n'])
} else if (s==2){
        Var = (V1(datasmooth,A,L,U) + V5(datasmooth,A,L,U) + V6(datasmooth,A,L,U)) / (2 - 4*A)^2
} else if (s==3){
        Var = (V1(datasmooth,A,L,U) + V2(datasmooth,A,L,U) + V4(datasmooth,A,L,U) + V5(datasmooth,A,L,U)) / (2 - 4*A)^2
} else {
        Var = (V1(datasmooth,A,L,U) + V2(datasmooth,A,L,U) + V3(datasmooth,A,L,U) + V4(datasmooth,A,L,U) + V5(datasmooth,A,L,U)) / (2 - 4*A)^2
}

gsd=10^(sqrt(Var))
#Calculate confidence interval
v     = sqrt(2) * erfinv(conf)
left  = mu * gsd^(-v)
right = mu * gsd^v

#---PLOT DATA--------------------------------------------------------------, log = "y"
plot((10^datatrim$x), (100*(datatrim$p * (1 - 2*A) + A)), type = "p",xlab = "Dose", ylab = "Response, %",
main =titl, ylim=c(0,100), xlim=c(0,100),pch='x')
par(new=T)
plot((10^datatrim$x), (100*(datatrim$p * (1 - 2*A) + A)), type = "l",xlab = "", ylab = "",
main ='', ylim=c(0,100), xlim=c(0,100),lty=2)
par(new=T)
plot(mu, 50, ylim=c(0,100), xlim=c(0,100),xlab = '',ylab = '',main = '',pch='o')
par(new=T)

#Error bars for the confidence interval on the mean. HERRORBAR doesn't do
#the ends right since we're only plotting one point, so we have to do this
#by hand. This is based on the code from herrorbar.m.
tee = 0.02 # make tee .02 x-distance for error bars
ptop = 0.5 + tee
pbot = 0.5 - tee

# build up nan-separated vector for bars
xb = matrix(0,9)
xb[1] = left
xb[2] = left
xb[3] = NaN
xb[4] = left
xb[5] = right
xb[6] = NaN
xb[7] = right
xb[8] = right
xb[9] = NaN

yb = matrix(0,9)
yb[1] = ptop
yb[2] = pbot
yb[3] = NaN
yb[4] = 0.5
yb[5] = 0.5
yb[6] = NaN
yb[7] = ptop
yb[8] = pbot
yb[9] = NaN

plot(xb,100 * yb,pch='|', ylim=c(0,100), xlim=c(0,100),xlab = '',ylab = '',main = '')
par(new=T)
plot(xb,100 * yb, type = "l",lty=1,ylim=c(0,100), xlim=c(0,100),xlab = '',ylab = '',main = '')
list(mu=mu,gsd=gsd,left=left,right=right)
 }
