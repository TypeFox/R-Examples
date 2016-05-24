
.packageName <- 'EMMIXskew'

EmSkew<-function(dat, g, distr = "mvn",ncov= 3,   
clust   = NULL,init= NULL,  
itmax   = 1000,epsilon = 1e-6,nkmeans = 0, nrandom = 10, 
nhclust = FALSE,debug = TRUE,initloop=20)
{
dat<- as.matrix(dat)

if(!is.null(init)  | !missing(init))
obj<-EmSkewfit2(dat,g, init, distr,ncov,itmax,epsilon)
else
{
if(is.null(clust) | missing(clust)) {

init<- try(init.mix(dat,g,distr,ncov,nkmeans,nrandom,nhclust,initloop))

if(!is.null(init)) 
obj<-EmSkewfit2(dat,g, init, distr,ncov,itmax,epsilon)
else
{
warning("not find initial values")
obj<-list()
obj$error=20
}
}
else {
obj<-EmSkewfit1(dat,g,clust, distr,ncov,itmax,epsilon,initloop)
}
}

error <- obj$error;

ret<-NULL

msg <- switch(tolower(error),
'1' = paste("stopped at (did not converge within) ", itmax, " iterations"),
'2' = paste("density fails at initial steps! "),
'3' = paste("allocation fails at initial steps"),
'12' = paste("density fails at estps! "),
'13' = paste("allocation fails at esteps"),
'20' =paste("not find initials"))


if(error>1)
{
cat('\n-----------------------\n')
warning("error code=",error,'\n',msg,"\n")
cat('\n-----------------------\n\n')
}

# summarize the results

if(error<=1) {

ICL <- getICL(dat,nrow(dat),ncol(dat),g,distr,ncov,obj$pro,obj$mu,obj$sigma,obj$dof,obj$delta,obj$clust)

ret<-obj

ret$ICL<-ICL$ICL

# mode point for each component

ret$modpts<- EmSkewMOD(ncol(dat),g,distr,ret$mu,ret$sigma,ret$dof,ret$delta)


if(debug) {

msg<-switch(tolower(distr), 
'mvn'=paste(g,"- Component Multivariate Normal Mixture Model"),
'mvt'=paste(g,"- Component Multivariate t      Mixture Model"),
'msn'=paste(g,"- Component Multivariate Skew Normal Mixture Model"),
'mst'=paste(g,"- Component Multivariate Skew-t Mixture Model"))


cat('\n-----------------------\n\n')
cat(msg,"\n")
cat('\n-----------------------\n\n')
 
switch(tolower(distr), 'mvn' = print(obj[1:8]),'mvt' = print(obj[1:9]),
'msn' = print(obj[c(1:8,10)]),'mst' = print(obj[1:10]))
print(ICL)
cat('\n-----------------------\n')
}

}
ret
}

getICL<-function(x,n,p,g,distr,ncov,pro,mu,sigma,dof,delta,clust)
{

x<-as.matrix(x);loglik<-0;nc=0

#if(length(table(clust<- unclass(as.ordered(clust))))!=g)
#stop(paste("labels should be of ", g, "levels"))

nn <- sum(clust>0) #outliers are numbered zero.

lnden<-(as.matrix((ddmix(x,n,p,g, distr, mu, sigma, dof, delta))))

for(h in 1:g) loglik=loglik+sum(ifelse(clust==h,log(pro[h])+lnden[,h],0) )

nc<-switch(tolower(ncov),
    '1' = (g-1)+g*p+p*(1+p)/2,  #common covariance
    '2' = (g-1)+g*p+p,          #common diagonal covariance
    '3' = (g-1)+g*(p+p*(1+p)/2),#general covariance
    '4' = (g-1)+g*(p+p),        #general diagonal covariance
    '5' = (g-1)+g*(p+1)  )      #sigma(h)*I_p

nc <- switch(tolower(distr),"mvn" = nc,"mvt" = nc + g,"msn" = nc + g*p,"mst" = nc + g*p + g)

ICL = loglik - nc/2*log(nn)

list(ICL=ICL)
}




EmSkewMOD<- function(p,g,distr,mu,sigma,dof,delta,nrand=10000){
distr<-tolower(distr)
mvrand<-function(n,p,distr,mean,cov,nu,del)
{
switch(distr,
'msn' = rdmsn(n,p,mean=mean,cov=cov,      del=del),
'mst' = rdmst(n,p,mean=mean,cov=cov,nu=nu,del=del),NULL)
}

if(!(distr %in% c("mvn","mvt","msn","mst")))
stop("model specified not available yet")

mu    <- array(mu,c(p,g))
sigma <- array(sigma,c(p,p,g))
dof   <- c(dof)
delta <- array(delta,c(p,g))

modpts <- mu

if(distr %in% c("mvn","mvt"))
return(modpts)


for(h in 1:g) {
dat <-  cbind(mvrand(nrand,p,distr,mu[,h],sigma[,,h],dof[h],delta[,h]))
den <-  ddmix(dat,nrand,p,1, distr,mu[,h],sigma[,,h],dof[h],delta[,h])
id  <-  which.max(den)
modpts[,h] <- dat[id,]
}
return(modpts)
}

EmSkewfit1<-function(dat,g,clust,distr,ncov,itmax,epsilon,initloop=20)
{
ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)
if(ndist>4|ncov<1|ncov>5) 
stop("the model specified is not available yet")

dat<-as.matrix(dat)
n <- nrow(dat)
p <- ncol(dat)

if(n <= 20*g)
stop("sample size is too small")

if(missing(clust) | is.null(clust))
stop("initial clust must be given")

clust <- unclass(as.ordered(clust))

if(max(clust)!=g) stop(paste("The levels of cluster should be g=",g))


obj<-.C('emskewfit1',PACKAGE="EMMIXskew",
as.double(dat),as.integer(n),as.integer(p),
as.integer(g),as.integer(ncov),as.integer(ndist),#6
pro   = double(g),mu  = double(p*g),sigma  = double(p*p*g),
dof   = double(g),delta=double(p*g), #11
tau   = double(n*g),double(n*g),double(n*g),double(n*g),double(n*g),
sumtau=double(g), sumvt=double(g),
sumzt=double(g), sumlnv=double(g), 
ewy=double(p*g),ewz=double(p*g),ewyy = double(p*p*g),#23
loglik= double(1),lk= double(itmax),aic= double(1),bic= double(1),
clust = as.integer(clust),#28
error = integer(1),as.integer(itmax),
as.double(epsilon),as.integer(initloop))[c(7:12,24:29)]

lk<-obj$lk;lk<-lk[lk!=0]

list(distr=distr,error=obj$error,loglik=obj$loglik,bic=obj$bic,aic=obj$aic,
pro=obj$pro,mu= array(obj$mu,c(p,g)),
sigma=array(obj$sigma,c(p,p,g)),dof=obj$dof,delta=array(obj$delta,c(p,g)),
clust=obj$clust,tau=array(obj$tau,c(n,g)),lk=lk)
}


EmSkewfit2<-function(dat,g, init, distr,ncov,itmax,epsilon)
{
ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4|ncov<1|ncov>5) 
stop("the model specified is not available yet")


dat<-as.matrix(dat)
#
n <- nrow(dat)
p <- ncol(dat)

if(n <= 20*g)
stop("sample size is too small")



if(is.null(init)|missing(init))
stop("init should be provided")

pro   <- init$pro
mu    <- init$mu
sigma <- init$sigma
dof   <- init$dof
delta <- init$delta

obj<-.C('emskewfit2',PACKAGE="EMMIXskew",
as.double(dat),as.integer(n),as.integer(p),
as.integer(g),as.integer(ncov),as.integer(ndist),#6
pro   = as.double(pro),mu  = as.double(mu),sigma  = as.double(sigma),
dof   = as.double(dof),delta= as.double(delta), #11
tau   = double(n*g),double(n*g),double(n*g),double(n*g),double(n*g),
sumtau=double(g), sumvt=double(g),
sumzt=double(g), sumlnv=double(g), #20
loglik= double(1),lk= double(itmax),aic= double(1),bic= double(1),
clust = integer(n),#25
error = integer(1),as.integer(itmax),as.double(epsilon))[c(7:12,21:26)]

lk<-obj$lk;lk<-lk[lk!=0]

list(distr=distr,error=obj$error,loglik=obj$loglik,bic=obj$bic,aic=obj$aic,
pro=obj$pro,mu= array(obj$mu,c(p,g)),
sigma=array(obj$sigma,c(p,p,g)),dof=obj$dof,delta=array(obj$delta,c(p,g)),
clust=obj$clust,tau=array(obj$tau,c(n,g)),lk=lk)
}





#*****************************************
#initial values
#*****************************************


initEmmix<-function(dat,g,clust,distr,ncov,maxloop=20)
{
ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4|ncov<1|ncov>5) 
stop("the model specified is not available yet")

dat<-as.matrix(dat)
n <- nrow(dat)
p <- ncol(dat)

clust<- unclass(as.ordered(clust))


obj<-.C('initfit_',PACKAGE="EMMIXskew",
as.double(dat),as.integer(n),as.integer(p),
as.integer(g),as.integer(ncov),as.integer(ndist),
pro   = double(g),mu  = double(p*g),sigma  = double(p*p*g),
dof   = double(g),delta=double(p*g), #11
tau   = double(n*g),double(n*g),double(n*g),double(n*g),double(n*g),
sumtau=double(g), sumvt=double(g),
sumzt=double(g), sumlnv=double(g), 
ewy=double(p*g),ewz=double(p*g),ewyy = double(p*p*g),#23
loglik= double(1),as.integer(clust),
error = integer(1),as.integer(maxloop)) [c(7:11,24,26)]


error <- obj$error
ret   <-NULL

if(error==0) {
ret<-list(distr=distr,error=error,loglik=obj$loglik,
pro=obj$pro,mu= array(obj$mu,c(p,g)),
sigma=array(obj$sigma,c(p,p,g)),
dof=obj$dof,delta=array(obj$delta,c(p,g)))
} else warning("error:",error)

ret

}

init.mix<-function(dat,g,distr,ncov,nkmeans,nrandom,nhclust,maxloop=20)
{
found<-list()
found$loglik<--Inf

n<-nrow(dat)
clust<-rep(1,n)
mclust<-NULL


if(g>1){
if(nkmeans>0) {


for(i in 1:nkmeans)
{
clust<-kmeans(dat,g,nstart=5)$cluster

if(min(table(clust))<10) next

obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))

if(length(obj)!=8 | obj$error) next

if(obj$loglik>found$loglik)
{
found<-obj
mclust<-clust
}

}

if(is.null(mclust))
nrandom <- 10
}



if(nrandom>0)
for(i in 1:nrandom) {
clust<-sample(1:g,n,replace=TRUE)

obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))
if(length(obj)!=8 | obj$error!=0) next
if(obj$loglik>found$loglik)
{
found<-obj 
mclust<-clust
}

} 

#methods<-c( "ward", "single", "complete", "average", "mcquitty", "median","cen")
methods<-c("complete")
#
if(nhclust)
{
dd <- as.dist((1 - cor(t(dat)))/2)  
#Use correlations between variables ``as distance''

for(j in methods){	
clust<- cutree(hclust(dd, j), k = g)

obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))
if(length(obj)!=8 | obj$error!=0) next

if(obj$loglik>found$loglik)
{
found<-obj;
mclust<-clust
}

} #end of for j
} #end if
#------------------------------------------------
}
else
{
obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))

if(length(obj)==8)  
{found<-obj;mclust<-clust}

}

if(is.null(mclust)) {
found <-NULL
warning("failed to find a initial values!")
}

found
}


# distance functions


# calculate the Scale free Weighted (mahalonobis distance) Ratio (SWR)and UWR
getSWR<-function(dat,g,sigma,clust, tau)
{

ret <- NULL

if(g>1)
{

intra <- intradist(dat,g,sigma, clust, tau) 
if(intra$error) stop("error")

inter <- interdist(dat,g,sigma, clust, tau) 
if(inter$error) stop("error")

ret <- list(SWR=sqrt(intra$OverallIntraDist1/inter$OverallInterDist1),
UWR=sqrt(intra$OverallIntraDist2/inter$OverallInterDist2))
}

ret

}



mypanel2<- function(x,y,...) {
par(new=TRUE);
smoothScatter(x,y,..., nrpoints=0)
}

mypanel3<- function(x,y,...) {
par(new=TRUE);
Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
smoothScatter(x,y, colramp = Lab.palette)

}

mypanel4<- function(x,y,...) {
xy <- cbind(x,y)
par(new=TRUE);
plot(xy, col = densCols(xy), pch=20)
}

panel.density <- function(x, col=1,...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )

oo = density(x)
y = oo$y
lines(oo$x,y/max(y))

}

conplot <-function(x,y, pro, mu, sigma, dof, delta,distr,
grid=300, nrand=6000,levels=seq(5,95,by=20),col ='white')
{

ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)
if(ndist>4) 
stop("the model specified is not available yet")  

ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, dof, delta)
{
ret<-ddmix(dat,n,p,g, distr, mu,sigma,dof,delta)
c(exp(ret)%*%pro )
} #joint density

    xlim = range(x)+c(-1,0)
    ylim = range(y)+c(-1,0)

g<-length(pro);p<-2

x1 <- seq(xlim[1], xlim[2], length=grid) 
y1 <- seq(ylim[1], ylim[2], length=grid) 

nx <- length(x1)
ny <- length(y1)
xoy <- cbind(rep(x1,ny), as.vector(matrix(y1,nx,ny,byrow=TRUE)))
X <- matrix(xoy, nx*ny, 2, byrow=FALSE)

dens     <- ddemmix(X, nx*ny,2,g, distr, pro, mu, sigma, dof, delta)
dens.mat <- matrix(dens,nx,ny)

n <- table(sample(1:g, nrand, replace = TRUE, prob = pro))
nn <- n
if(length(n) <g) {
nn <- rep(0,g)
for(i in as.numeric(names(n)))
nn[i] <- n[paste(i)]
}
rand     <-  rdemmix(nn,p,g,distr,mu,sigma,dof,delta)
rand.den <-  ddemmix(rand, nrand,2,g, distr, pro, mu, sigma, dof, delta)
cont     <-  quantile(rand.den, prob=levels/100)
contour(x1, y1, dens.mat, levels=cont, add=TRUE,drawlabels=FALSE,lty=1,col =col)
}

conplot2 <- function (x, y, pro, mu, sigma, dof, delta, distr, grid = 300, 
    nrand = 6000, levels = seq(5, 95, by = 20)) 
{
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3, mst = 4, 
        5)
    if (ndist > 4) 
        stop("the model specified is not available yet")
    ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, 
        dof, delta) {
        ret <- ddmix(dat, n, p, g, distr, mu, sigma, dof, delta)
        c(exp(ret) %*% pro)
    }
    g <- length(pro)
    p <- 2

#make mesh

    xlim = range(x)+c(-1,0)
    ylim = range(y)+c(-1,0)

    x1 <- seq(xlim[1], xlim[2], length= grid)
    y1 <- seq(ylim[1], ylim[2], length= grid)

    nx <- length(x1)
    ny <- length(y1)

    xoy <- cbind(rep(x1, ny), as.vector(matrix(y1, nx, ny, byrow = TRUE)))
    X <- matrix(xoy, nx * ny, 2, byrow = FALSE)
    dens <- ddemmix(X, nx * ny, p, g, distr, pro, mu, sigma, 
        dof, delta)

    dens.mat <- matrix(dens, nx, ny)

#
    n <- table(sample(1:g, nrand, replace = TRUE, prob = pro))
    nn <- n
    if (length(n) < g) {
        nn <- rep(0, g)
        for (i in as.numeric(names(n))) nn[i] <- n[paste(i)]
    }

    rand <- rdemmix(nn, p, g, distr, mu, sigma, dof, delta)
    rand.den <- ddemmix(rand, nrand, 2, g, distr, pro, mu, sigma,dof, delta)
    cont <- quantile(rand.den, prob = 1-levels/100)

    samp <- cbind(x,y)
    samp.den <-ddemmix(samp, length(x), 2, g, distr, pro, mu, sigma,dof, delta)

select <-  which(samp.den>cont)

clust  <-  ifelse(samp.den>cont,2,1)

list(select,clust,x1=x1, y1=y1, dens.mat=dens.mat, cont=cont)
}


conplot3 <- function (x, y, pro, mu, sigma, dof, delta, modpts,distr, grid =300, 
    nrand = 10000, levels = seq(5, 95, by = 20)) 
{
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3, mst = 4, 
        5)
    if (ndist > 4) 
        stop("the model specified is not available yet")
    ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, 
        dof, delta) {
        ret <- ddmix(dat, n, p, g, distr, mu, sigma, dof, delta)
        c(exp(ret) %*% pro)
    }

    g <- length(pro)
    p <- 2

#----------------------------------------------------
#mesh
#----------------------------------------------------

#make mesh

    xlim = range(x)+c(-1,0)
    ylim = range(y)+c(-1,0)

    x1 <- seq(xlim[1], xlim[2], length=grid)
    y1 <- seq(ylim[1], ylim[2], length=grid)

    nx <- length(x1)
    ny <- length(y1)

    xoy <- cbind(rep(x1, ny), as.vector(matrix(y1, nx, ny, byrow = TRUE)))

    X <- matrix(xoy, nx * ny, 2, byrow = FALSE)

#--------------------------------------------------



for(h in 1:g) { #do each component one by one

    dens <- ddemmix(X, nx * ny, 2, 1, distr, c(1), mu[,h], sigma[,,h],dof[h], delta[,h])

    dens.mat <- matrix(dens, nx, ny)

#  randon sample

    rand <- rdemmix(c(nrand), 2, 1, distr, mu[,h], sigma[,,h], dof[h], delta[,h])

    rand.den <- ddemmix(rand, nrand, 2, 1, distr, c(1), mu[,h], sigma[,,h], dof[h], delta[,h])

#-----------------------------------------

    cont <- quantile(rand.den, prob = 1-levels/100)

    contour(x1, y1, dens.mat, levels = cont, add = TRUE, drawlabels = FALSE,lty = 1, col = h)
    
        if (!is.null(modpts)) 
        points(t(modpts[, h]), col = h,pch=3)
} #end of h loop


}

EmSkew.filter <- function (S, g=1,distr="mst", diag.panel = TRUE, upper.panel = "type2", 
    lower.panel = "type3", levels = 90, attop = FALSE,title="",path="",plot=TRUE) 
{
    
S <- as.matrix(S)

dat <- S[,1:2]

obj <- EmSkew(dat,g,distr,ncov=3,itmax=200,debug=0)

ppp <- conplot2(c(dat[,1]), c(dat[,2]), obj$pro, obj$mu, obj$sigma, 
obj$dof, obj$delta, obj$distr, nrand=10000,levels = levels)

clust <- ppp[[2]]
select<- ppp[[1]]

#---------------------------------

mypanel <- function(x, y, ...) {

        par(new = TRUE)
        points(x, y, ..., col = clust)

        st <- pmatch(c(x[1], y[1]), S[1, ])
#st <- c(current.row(),current.column())

if(st[1]==1&st[2]==2) {

a=contourLines(ppp$x1, ppp$y1, ppp$dens.mat, levels=ppp$cont)[[1]]

ax <- a$x
#ax[ax<ppp$x1[1]]=ppp$x1[1]+1

ay <- a$y

lines(ax,ay,lty = 2, col = 'blue')

}

}

    if (diag.panel) {
        diag.panel <- panel.density
    }
    else diag.panel <- NULL


    upper.panel <- switch(upper.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)

    lower.panel <- switch(lower.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)

#-------------------------------

if(plot) {

dev.new()
pairs(S, pch = ".", panel=mypanel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("Before Filtering (EmSkew):",toupper(distr)))


dev.new()
pairs(S[select,], upper.panel=upper.panel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("After Filtering (EmSkew):",toupper(distr)))
}

#-------------------------------

if(path!='') {

png(paste(path,'/',title,"-before.png",sep=''),width=512,height=512)
 
pairs(S, pch = ".", panel=mypanel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("Before Filtering (EmSkew):",toupper(distr)))

dev.off()

png(paste(path,'/',title,"-after.png",sep=''),width=512,height=512)
 
pairs(S[select,], upper.panel=upper.panel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("After Filtering (EmSkew):",toupper(distr)))
dev.off()

}

list(subset=select,clust=clust,filter=obj)
}





EmSkew.contours <- function (S, obj = NULL, clust = NULL,distr="",diag.panel = TRUE, upper.panel = "type2", 
    lower.panel = "type3", levels = seq(5, 95, by = 20), plot=TRUE, title="",path='',attop = FALSE) 
{
    mypanel <- function(x, y, ...) {
        par(new = TRUE)
        smoothScatter(x, y, ..., nrpoints = 0)
        g <- length(obj$pro)
        st <- pmatch(c(x[1], y[1]), S[1, ])
#st <- c(current.row(),current.column())

	
	conplot3(x, y, obj$pro, obj$mu[st, ], obj$sigma[st, st, 
            ], obj$dof, obj$delta[st, ],obj$modpts[st,], obj$distr, levels = levels)
    }
    if (diag.panel) {
        diag.panel <- panel.density
    }
    else diag.panel <- NULL
    upper.panel <- switch(upper.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)
    lower.panel <- switch(lower.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)


if(plot) {

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of Components using EmSkew:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EmSkew:",toupper(distr), "Distribution") )

}


if(path!='') {

png(paste(path,'/',title,".png",sep=''),width=512,height=512)

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of Components using EmSkew:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EmSkew:",toupper(distr), "Distribution") )

dev.off()

}

}



EmSkew.flow <-function(S,obj=NULL,distr="",diag.panel=TRUE,
upper.panel="type2",lower.panel="type3",
levels=seq(5,95,by=20),attop=FALSE,clust=NULL,title="",path="",plot=TRUE) {

mypanel<- function(x,y,...) {
par(new=TRUE);
smoothScatter(x,y,..., nrpoints=0)
g <- length(obj$pro)

st <- pmatch(c(x[1],y[1]),S[1,])
#st <- c(current.row(),current.column())


conplot(x,y, obj$pro,obj$mu[st,],obj$sigma[st,st,],
obj$dof,obj$delta[st,],obj$distr,levels=levels,nrand=10000)

if(!is.null(obj$modpts))
points(t(obj$modpts[st,]),col=1:g,pch=3)
}



if(diag.panel) {
diag.panel<- panel.density}
else diag.panel<- NULL

upper.panel<-switch(upper.panel,
                    "type2"=mypanel2,
		    "type3"=mypanel3,
		    "type4"=mypanel4,
		    NULL)


lower.panel<-switch(lower.panel,
                    "type2"=mypanel2,
		    "type3"=mypanel3,
		    "type4"=mypanel4,
		    NULL)

if(plot){

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of Mixture using EmSkew:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EmSkew:",toupper(distr), "Distribution") )



}


if(path!='') {

png(paste(path,'/',title,".png",sep=''),width=512,height=512)

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of EmSkew:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EmSkew:",toupper(distr), "Distribution") )

dev.off()
}

}


ddmvn<-function(dat,n,p,mean = rep(0,p),cov = diag(p)                   )
{
exp(ddmix(dat,n,p,1,"mvn", mean,cov,0,rep(0,p)))

}

ddmvt<-function(dat,n,p,mean = rep(0,p),cov=diag(p),nu=4                )
{
exp(ddmix(dat,n,p,1, "mvt", mean,cov,nu,rep(0,p)))
}

ddmsn<-function(dat,n,p,mean=rep(0,p),cov=diag(p),       del = rep(0,p))
{
exp(ddmix(dat,n,p,1, "msn", mean,cov,0,del))
}

ddmst<-function(dat,n,p,mean = rep(0,p),cov=diag(p),nu=4,del = rep(0,p))
{
exp(ddmix(dat,n,p,1, "mst", mean,cov,nu,del))
}



ddmix <- function(dat,n,p,g,distr, mu, sigma, dof=NULL, delta=NULL)
{

if(is.null(dof))
dof <- rep(4,g)

if(is.null(delta))
delta <- array(0,c(p,g))

ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4) 
stop("the model specified is not available yet")

dat<-as.matrix(dat);

if(n == 1 & (ncol(dat) ==1))
dat<-t(dat)

if(nrow(dat)!=n | ncol(dat)!=p )
stop("dat does not match n and p.")

#is mu,sigma,dof,delta specified correctly?

if(length(c(mu)) != (p*g))
stop(paste("mu should be a ",p, 'by', g, "matrix!"))

if(length(c(sigma)) != (p*p*g))
stop(paste("sigma should be a ",p, 'by', p,'by', g, " array!"))

if(length(c(dof)) != g)
stop(paste("dof should be a ",g, " vector!"))

if(length(c(delta)) != (p*g))
stop(paste("delta should be a ",p, 'by', g, " array!"))

obj<-.C('ddmix',PACKAGE="EMMIXskew",
as.double(dat),as.integer(n),
as.integer(p),as.integer(g),as.integer(ndist),
as.double(mu),as.double(sigma),
as.double(dof),as.double(delta),
den = double(n*g),error = integer(1))[10:11] 

if(obj$error) stop("error")

(matrix(obj$den,ncol=g))

}




# rdmvn is a wrapper of the function rmvnorm from R package "mvtnorm"


rdmvn<-function (n, p,mean = rep(0,p), cov = diag(p)) 
{
cov<-as.matrix(cov)
    
if (nrow(cov) != ncol(cov)) {
        stop("cov must be a square matrix")
    }
    if (length(mean) != nrow(cov)) {
        stop("mean and cov have non-conforming size")
    }

rmvnorm(n, mean = mean, sigma = cov,method="chol")

}



rdmvt<-function(n,p,mean = rep(0,p),cov=diag(p),nu=3)
{
cov<-as.matrix(cov)
u<-rgamma(n,nu/2,nu/2)
t(t(rdmvn(n,p,cov=cov)/sqrt(u))+mean)
}


rdmsn<-function(n,p,mean=rep(0,p),cov=diag(p),del=rep(0,p))
{

x<-rdmvn(n,p,mean,cov)
z<-abs(rnorm(n))
as.matrix(z%*%t(del)+x)

}



rdmst<-function(n,p,mean = rep(0,p), cov=diag(p),nu=10,del = rep(0,p))
{
u<-rgamma(n,nu/2,nu/2)
x<-t(t(rdmvn(n,p,cov=cov)/sqrt(u))+mean)
z<-abs(rnorm(n)/sqrt(u))
as.matrix(z%*%t(del)+x)
}


rdemmix2<-function(n,p,g,distr,pro,mu,sigma,dof=NULL,delta=NULL)
{

n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
nn <- n0
if(length(nn) <g) {
nn <- rep(0,g)
for(i in as.numeric(names(n0)))
nn[i] <- n0[paste(i)]
}

names(nn) <- NULL

rdemmix(nn,p,g,distr,mu,sigma,dof,delta)

}

rdemmix3<-function(n,p,g,distr,pro,mu,sigma,dof=NULL,delta=NULL)
{

if(length(pro) != g)
stop(paste("pro should be a ",g, " vector!"))

n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
nn <- n0
if(length(nn) <g) {
nn <- rep(0,g)
for(i in as.numeric(names(n0)))
nn[i] <- n0[paste(i)]
}

names(nn) <- NULL

dat <- rdemmix(nn,p,g,distr,mu,sigma,dof,delta)

list(data = dat, cluster = rep(1:g,nn) )

}


rdemmix<-function(nvect,p,g,distr,mu,sigma,dof=NULL,delta=NULL)
{

if(length(c(nvect))!=g) stop("nvect should be a vector")

ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4) 
stop("the model specified is not available yet")

if(is.null(dof))
dof <- rep(4,g)

if(is.null(delta))
delta <- array(0,c(p,g))

if(length(c(mu)) != (p*g))
stop(paste("mu should be a ",p, 'by', g, "matrix!"))

if(length(c(sigma)) != (p*p*g) )
stop(paste("sigma should be a ",p, 'by', p,'by', g, " array!"))

if(length(c(dof)) != g)
stop(paste("dof should be a ",g, " vector!"))

if(length(c(delta)) != (p*g) )
stop(paste("delta should be a ",p, "by", g, " array!"))

# to fix the "g=1" bug,

mu    = array(mu, c(p,g))
sigma = array(sigma, c(p,p,g))
delta = array(delta, c(p,g))

dat<-array(0,c(10,p))

mvrand<-function(n,p,ndist,mean,cov,nu,del)
{

switch(ndist,
'1' = rdmvn(n,p,mean=mean,cov=cov              ),
'2' = rdmvt(n,p,mean=mean,cov=cov,nu=nu        ),
'3' = rdmsn(n,p,mean=mean,cov=cov,      del=del),
'4' = rdmst(n,p,mean=mean,cov=cov,nu=nu,del=del))
}


if(g>=1)
for(h in 1:g)
{
if(nvect[h]>0)
dat<-rbind(dat,mvrand(nvect[h],p,ndist,mu[,h],sigma[,,h],dof[h],delta[,h]))

}

dat[-(1:10),]
}


# BOOTSTRAP functions

bootstrap <- function(x,n,p,g,distr,ncov,popPAR,B=99, replace=TRUE,itmax=1000,epsilon=1e-5)
{
x<-as.matrix(x);


if(missing(popPAR))
stop("please run the function EmSkew() first")

counter <- 0
nnn <- g*(1 + p + p*p + 1 + p) 

ret <- array(0, c(B,nnn)) 

dimnames(ret) <- list(1:B, c(
paste("pi",1:g,sep=''),
paste("mu",rep(1:p,g),rep(paste(1:g,sep=''),rep(p,g)),sep=''),
paste("sigma",rep(paste(rep(1:p,rep(p,p)),rep(1:p,p),sep=''),g),rep(paste(',',1:g,sep=''),rep(p*p,g)),sep=''),
paste("dof",1:g,sep=''),
paste("delta",rep(1:p,g),rep(paste(1:g,sep=''),rep(p,g)),sep='')))



for(i in 1:(2*B) )
{

if(replace)
dat <- x[sample(1:n,n,replace=TRUE),]
else
dat <- rdemmix3(n,p,g,distr,popPAR$pro,popPAR$mu,popPAR$sigma,popPAR$dof,popPAR$delta)


obj <- EmSkewfit2(dat,g, popPAR, distr,ncov,itmax,epsilon)

if(obj$error > 1) next

counter <- counter +1 

ret[counter,] <- c(obj$pro,obj$mu,obj$sigma,obj$dof,obj$delta)

if(counter >= B) break 

}

std <- sqrt(apply(ret[1:counter,],MARGIN=2,FUN= "var"))

names(std) <- dimnames(ret)[[2]]
std
}

bootstrap.noc <- function(x,n,p,g1,g2,distr,ncov,B=99, replace=TRUE,itmax=1000,epsilon=1e-5)
{

x<-as.matrix(x);

if(g1 >= g2)
stop("g1 should be less than g2")

if(g1 < 1)
stop("g1 should be greater than 0")

counter <- 0

vlk <- rep(0,g2-g1+1)

ret <- array(0,c(B,g2-g1))

dimnames(ret) <- list(1:B,paste(1+g1:(g2-1),"vs",(g1:(g2-1)),sep=' '))

lk0 <- -Inf
# start

clust <- rep(1,n)

for( g in g1:g2) {

counter <- 0

lk1 <- -Inf

while(counter < 10) {

if(g>1)
clust <- kmeans(x,g,nstart=5)$cluster

emobj <- EmSkewfit1(x,g,clust, distr,ncov,itmax,epsilon)

if(emobj$error>1) next

if(emobj$loglik > lk1) 
lk1 <- emobj$loglik 

counter = counter +1

}

#save the results for g

dput(emobj,paste("ReturnOf_g_",g,".ret",sep=''))

#----------------------

counter <- 0

lk0 <- lk1

vlk[g-g1+1] <- lk0

if(g < g2) {


for(i in 1:(2*B) )
{

if(replace)
dat <- x[sample(1:n,n,replace=TRUE),]
else
dat <- rdemmix2(n,p,g,distr,emobj$pro,emobj$mu,emobj$sigma,emobj$dof,emobj$delta)

if(is.null(dat)) stop("I can not generate the data!")

obj <- EmSkewfit2(dat,g, emobj, distr,ncov,itmax,epsilon)

if(obj$error > 1) next

ii <- 0

lk2 <- -Inf

while(ii<10) {

clust <- kmeans(dat,g+1,nstart=5)$cluster

obj2 <- EmSkewfit1(dat,g+1,clust, distr,ncov,itmax,epsilon)

ii <- ii + 1

if(obj2$error>1) next

if(obj2$loglik > lk2)
lk2 <- obj2$loglik 

} #end ii loop


counter <- counter +1 

ret[counter,g-g1+1] <- -2*(obj$loglik-lk2)

if(counter >= B) break 

} #end i loop
} # end g loop

}# end if

pvalue <- rep(0,g2-g1)

for(i in 1:(g2-g1))
{
pvalue[i] <- sum(ret[,i] < 2*(vlk[i+1]-vlk[i]))/B
}

list(ret=ret,vlk=vlk,pvalue=pvalue)
}
# mahalonobis distance

mahalonobis<-function(p, g, mu, sigma) 
{


obj<-.C('mahalonobis_',PACKAGE="EMMIXskew",
as.integer(p),as.integer(g),as.double(mu),as.double(sigma), 
dist = double(g*g), error = integer(1)) 

if(obj$error) stop("") 

matrix(obj$dist, ncol=g)


}


intradist<-function(dat,g,sigma, clust, tau) 
{
dat<-as.matrix(dat)

intraobj<-.C('intradist_',PACKAGE="EMMIXskew",
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(g),as.integer(clust),as.double(sigma),as.double(tau),
dist1=double(g+1),dist2 = double(g+1), error = integer(1)) 

list(error=intraobj$error,dist1 = intraobj$dist1[1:g],dist2 = intraobj$dist2[1:g],
OverallIntraDist1=intraobj$dist1[g+1],OverallIntraDist2=intraobj$dist2[g+1])
}


interdist<-function(dat,g,sigma, clust, tau) 
{

dat<-as.matrix(dat)


interobj<-.C('interdist_',PACKAGE="EMMIXskew",
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(g),as.integer(clust),as.double(sigma),as.double(tau),
dist1=double(g*g+1),dist2 = double(g*g+1), error = integer(1)) 



list(error=interobj$error,dist1 = matrix(interobj$dist1[1:(g*g)],ncol=g),
dist2 = matrix(interobj$dist2[1:(g*g)],ncol=g), 
OverallInterDist1 = interobj$dist1[(g*g)+1],
OverallInterDist2 = interobj$dist2[(g*g)+1]   )

}


#----------------------------------
# U\V |  V_1  V_2 ... V_C  | sums
# ---------------------------------
# U_1 |  n_11 n_12... n_1c | a_1
# U_2 |  n_21 n_22... n_2c | a_2
# .
# .
# .
# U_R |  n_r1 n_r2... n_rc | a_r
#---------------------------------
#sums  |  b_1  b_2 ... b_c  | N
#---------------------------------

rand.index<- function(LabelA,LabelB) {


u <- unclass(as.ordered(LabelA))
v <- unclass(as.ordered(LabelB))

if((N <- length(u)) != length(v))
stop("Label A and B does not match!")

#Adjusted Rand Index (ARI)

row <- max(u)
col <- max(v)

nvect <- array(0,c(row,col))

for(i in 1:row) {
 for(j in 1:col) {
   nvect[i,j]<-sum(u==i&v==j)

}}

SumsA <- rowSums(nvect)
SumsB <- colSums(nvect)

a = 0
for(i in 1:row)
a=a+choose(SumsA[i],2)

b = 0
for(j in 1:col)
b=b+choose(SumsB[j],2)

c <- a*b/choose(N,2)

d = 0
for(i in 1:row) {
 for(j in 1:col) {
   d=d+choose(nvect[i,j],2)
}}

#Adjusted Rand INdex
arj <- (d-c)/((a+b)/2-c)

#Rand Index (RI)

a=d

b=0

for(l in 1:row) {

for(i in 1:(col-1)) {
 for(j in (i+1):col) {
   b=b+ nvect[l,i]*nvect[l,j]

}}

}

c=0

for(l in 1:col) {

for(i in 1:(row-1)) {
 for(j in (i+1):row) {
   c=c+ nvect[i,l]*nvect[j,l]

}}

}

#d= choose(N,2)-a-b-c

#rad= (a+d)/choose(N,2)

rad= (choose(N,2)-b-c)/choose(N,2)

ind <- c(rad,arj)

names(ind) <- c("Rand Index (RI)","Adjusted Rand Index (ARI)")

ind
}


inverse <-function(sigma,p)
{
if(length(c(sigma))!=p*p | ncol(sigma)!=p)
stop("sigma should be p by p matrix")

obj <- .Fortran('inverse3',PACKAGE="EMMIXskew",
as.double(sigma),inv=double(p*p), 
det=double(1),as.integer(p), error = integer(1),
count = integer(1),index = integer(p))

if(obj$error) stop("")
a <- array(obj$inv,c(p,p))
as.matrix(t(a)%*%a)
}

tau2clust<-function(tao)
{
apply(tao,FUN=which.max,MARGIN=1)
}

getcov <-function(msigma,sumtau,n,p,g,ncov)
{
sigma<-array(0,c(p,p))
if( (ncov==1)|(ncov==2))
{
for(h in 1:g)
sigma<-sigma+sumtau[h]*msigma[,,h]
sigma<-as.matrix(sigma/n)

if(ncov==2)
sigma<-diag(c(diag(sigma)),p)
for(h in 1:g)
msigma[,,h]=sigma
}

if(p>1)
{
if(ncov==4)
for(h in 1:g)
msigma[,,h]<-diag(c(diag(msigma[,,h])),p)

if(ncov==5)
for(h in 1:g)
msigma[,,h]<-diag(sum(diag(msigma[,,h]))/p,p)
}

msigma
}


mvt.dof <-
function(sumtau,sumlnv,lx=2+1e-4,ux=200)
{

if(sumtau <=2)
return(4L)

f<-function(v,sumlnv,sumtau) 
{
sumtau*( log(v/2)-digamma(v/2)+1)+ sumlnv
}

if(f(lx,sumlnv,sumtau)*f(ux,sumlnv,sumtau)>0)
return(ux)
else
(uniroot(f,c(lx,ux),sumlnv=sumlnv,sumtau=sumtau)$root)
}



error.rate<-function(clust1,clust2)
{

clust1 <- unclass(as.ordered(clust1))
clust2 <- unclass(as.ordered(clust2))

if((n=length(clust1))!=length(clust2))
{warning("error: length not equal");return}

if( (g=length(table(clust1)))!=length(table(clust2)))
{warning("the number of clusters are not equal");return}

permute<-function(a)
{
n<-length(a)
if(n==1)
f<-a
else
{
nm<-gamma(n)
f<-array(0,c(n,n*nm))
j<-1

for(i in a)
{
 f[1, (j-1)*nm+1:nm]<-i
 f[-1,(j-1)*nm+1:nm]<-permute(setdiff(a,i))
 j<-j+1
}
}

f
}


#
id<-1:n

cmb<-permute(1:g)

nperm<-ncol(cmb)

rate<-rep(0,nperm)

#
for(i in 1:nperm)
{

tmp<-rep(0,g)

tc<-rep(0,n)

for(j in 1:g)
tc[clust2==j]=cmb[j,i]

for(j in 1:g)
{  
tmp1<-0 

for(k in (1:g)[-j])
        tmp1<-tmp1+length(intersect(id[clust1==j],id[tc==k]))

tmp[j]<-tmp1
}

rate[i]<-sum(tmp)/n
}

min(rate)
}


#end


