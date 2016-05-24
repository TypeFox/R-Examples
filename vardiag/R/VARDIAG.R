#assign("p11", NULL, envir = .GlobalEnv)
#assign("p12", NULL, envir = .GlobalEnv)
#assign("p21", NULL, envir = .GlobalEnv)
#assign("p22", NULL, envir = .GlobalEnv)

#p11 <<- par() p22 <<- par() p12 <<- par() p21 <<- par()
#p11 <- par() p22 <- par() p12 <- par() p21 <- par()
.varDiagOptions <- new.env(FALSE, globalenv())
assign("p11", NULL, envir = .varDiagOptions)
assign("p12", NULL, envir = .varDiagOptions)
assign("p21", NULL, envir = .varDiagOptions)
assign("p22", NULL, envir = .varDiagOptions)

gamsph<-function(h,th=rbind(1,1,1)){(0<h)*(h<=th[3])*(th[1]+th[2]*(3/2*(h/th[3])-1/2*(h/th[3])^3))+(h>th[3])*(th[1]+th[2])}

fth<-function(th,y,h1,w1=1){(y-gamsph(h1,th))/w1}

ftc<-function(th,y,h1,w1){(y-gamsph(h1,th))/gamsph(h1,th)}

ftg<-function(th,y,h1,cv1){cv1%*%(y-gamsph(h1,th))}

fts <- function(th, y, h1, cv1) {
	cv1 %*% (y-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(h1,th)^0.25)
}

ftsOpt <- function(th, y, h1, cv1) {
	ret = cv1 %*% (y-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(h1,th)^0.25)
	mean(ret^2)
}

gamsph1<-function(h,th=rbind(1,1,1)){1}

gamsph2<-function(h,th=rbind(1,1,1)){(0<h)*(h<=th[3])*(3/2*(h/th[3])-1/2*(h/th[3])^3)+(h>th[3])}

gamsph3<-function(h,th=rbind(1,1,1)){(0<h)*(h<=th[3])*3/2*th[2]/th[3]*((h/th[3])^3-h/th[3])}

hyperg<-function(r){
f<-1+1.125*r+1.1484375*r^2+1.158007813*r^3+1.16317749*r^4+1.166408539*r^5;
#a<-0.75;
#b<-0.75;
#c<-0.5;
#k<-1;
#f<-1;
#n<-ceiling(10+exp((max(r)-1)*30)*500);
#n<-10;
#for (i in 2:50){
#	k<-k*(a+i-2)*(b+i-2)/(c+i-2)*r/(i-1);
#	f<-f+k
#	}
f}

ficorr<-function(r){gamma(0.75)^2/(sqrt(pi)-gamma(0.75)^2)*((1-r^2)*hyperg(r^2)-1)}

estvar <- function(h0, y, iter=50, tolerance=0.0002, trace=1, th0=rbind(0,1,1)) 
	{

	#EJP added:
	#stop("this function requires nlregb (an S-Plus proprietary function) to work")
	
	n<-ceiling(sqrt(2*length(h0)))
	
	#Vorbereitung fuer covgamma
	n1<-n*(n-1)/2
	#1. index der gamma[i,j] matrix
	i1<-matrix(1:n,n,n)
	#1. teil des zeilenindex der covgamma gamma matrix
	k1<-matrix(i1[row(i1)<col(i1)],n1,n1)
	#2. teil des zeilenindex der covgamma gamma matrix
	k2<-matrix(t(i1)[row(i1)<col(i1)],n1,n1)
	#1. teil des spaltenindex der covgamma gamma matrix
	k3<-t(k1)
	#2. teil des spaltenindex der covgamma gamma matrix
	k4<-t(k2)
	
	if(!missing(th0)) {
		#EJP outcommented:
		#opt<-nlregb(n*(n-1)/2,cbind(0,max(y/2),max(h0)),fts,y=y^0.25,h1=h0,cv1=diag(n1),lower=cbind(0,0,0))
		opt<-optim(par = c(0,max(y/2),max(h0)), ftsOpt,
			lower=cbind(0,0,0), method = "L-BFGS-B",
			y=y^0.25, h1=h0, cv1=diag(n1))
		th1 <- opt$par
	}
	else
		th1<-th0
	th1<-cbind(0,max(y/2),max(h0))
	#th0<-th1_c(3.72635248595876, 15.5844183738953, 1.22109233789852)
	#th1<-c(0.0000000,7.6516077,0.7808538)
	for (i in 1:iter) {
		if(trace>0) 
			print(i)
		gg<-sqrt(2*gamsph(h0,th1))
		#Spalte 1, Spalte 2, ...
		#gamma vektor wird als matrix dargestellt
		tt<-matrix(gg[(t(i1)-2)*(t(i1)-1)/2+i1],n,n)
		#symmetrisierung
		tt1<-tt
		tt1[row(tt1)>col(tt1)]<-t(tt)[row(tt1)>col(tt1)]
		#diagonale loeschen
		tt1[row(tt1)==col(tt1)]<-0
		#covgamma wird berechnet
		cg<-matrix(tt1[(k4-1)*n+k1]+tt1[(k2-1)*n+k3]-tt1[(k3-1)*n+k1]-tt1[(k4-1)*n+k2],n1,n1)
		cgcg<-outer(gg,gg,"*")
		corg<-sqrt(cgcg)*ficorr((cg*lower.tri(cg))/cgcg)
		corg<-sqrt(2)*(sqrt(pi)-gamma(0.75)^2)/pi*(corg+t(corg)+diag(gg))
		infm<-solve(corg);
		cv<-chol((infm+t(infm))/2);
		#sc<-cbind(1/th1[2],1/th1[2],1/th1[3])
		#EJP outcommented:
		#opt<-nlregb(n*(n-1)/2,th1,fts,y=y^0.25,h1=h0,cv1=cv,lower=cbind(0,0,0))
		opt <- optim(par = th1, ftsOpt,
			lower=cbind(0,0,0), method = "L-BFGS-B",
			y=y^0.25,h1=h0,cv1=cv)
		if(trace>0) print(opt$par)
		if(sum(abs((th1-opt$par)/(th1+0.00001)))<=tolerance)
			break
		th1<-opt$par
	}
	print("Fertig")
	v<-list(pars=opt$par)
	v$cg<-corg
	v$res<-y^0.25-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(h0,v$pars)^0.25
	v$lof<-t(v$res)%*%solve(corg,v$res)
	v
}

varobj<-function(m,iter=50,tolerance=0.0002,trace=1,loo=FALSE){
n<-dim(m)[1]

#a1<-t(m[,3]-t(matrix(m[,3],n,n)))
#b1<-t(m[,1]-t(matrix(m[,1],n,n)))
#c1<-t(m[,2]-t(matrix(m[,2],n,n)))
a1<-outer(m[,3],m[,3],FUN="-")
b1<-outer(m[,1],m[,1],FUN="-")
c1<-outer(m[,2],m[,2],FUN="-")
#d1<-cbind(sqrt(b1[row(b1)<col(b1)]^2+c1[row(c1)<col(c1)]^2),a1[row(a1)<col(a1)]^2)
y<-a1[row(a1)<col(a1)]^2/2
h0<-sqrt(b1[row(b1)<col(b1)]^2+c1[row(c1)<col(c1)]^2)

v<-estvar(h0,y,iter,tolerance,trace)

XM<-cbind(gamsph1(h0,v$pars),gamsph2(h0,v$pars),gamsph3(h0,v$pars))*(gamsph(h0,v$pars))^(-0.75)/4
v$info<-solve(t(XM)%*%solve(v$cg,XM))

loores<-matrix(0,n,n)
tha<-matrix(0,n,3)
lofa<-matrix(0,n,1)
cda<-matrix(0,n,1)

v$h<-h0
v$y<-y

if(loo==TRUE){
  for (i in 1:n){
	print(i)
	m1<-m[-i,]

	a1<-t(m1[,3]-t(matrix(m1[,3],n-1,n-1)))
	b1<-t(m1[,1]-t(matrix(m1[,1],n-1,n-1)))
	c1<-t(m1[,2]-t(matrix(m1[,2],n-1,n-1)))
	y<-a1[row(a1)<col(a1)]^2/2
	h0<-sqrt(b1[row(b1)<col(b1)]^2+c1[row(c1)<col(c1)]^2)

	z<-estvar(h0,y,iter,tolerance,trace,th0=v$pars)
	lofa[i,1]<-v$lof-z$lof
	tha[i,]<-z$pars
	cda[i,1]<-t(v$pars-z$pars)%*%v$info%*%(v$pars-z$pars)
	mm2<-m[i,]
	mm3<-t(t(m)-mm2)^2/2
	h<-sqrt(mm3[,1]+mm3[,2])
	loores[i,]<-mm3[,3]^0.25-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(h,z$pars)^0.25
	}
}
v$loores<-loores
v$tha<-tha
v$lofa<-lofa
v$cda<-cda
v$data<-m

class(v)<-"varobj"
v
}

print.varobj<-function(x,...){print(x$pars); print(x$lof);invisible(x)}

#r[row(r)<col(r)]<-v$res
#r<-r+t(r)

PlotDiag.varobj<-function(v, region = NULL, xyi = 0, zmv = 0) {

  palette(c("black","cyan","magenta","green3","yellow","blue","white","red"))

  n<-length(v$h)
  infm<-solve(v$cg);
  cv<-chol((infm+t(infm))/2);
  XM<-cbind(gamsph1(v$h,v$pars),gamsph2(v$h,v$pars),gamsph3(v$h,v$pars))*(gamsph(v$h,v$pars))^(-0.75)/4
  Vare<-v$cg-XM%*%solve(t(XM)%*%solve(v$cg,XM),t(XM))
  #sig<-mean(sqrt(diag(Vare)))
  e<-v$res
  sig<-sqrt(sum(e^2)/(n-3))
  gdd<-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(v$h,v$pars)^0.25+e*sig/sqrt(diag(Vare))
  
  r1<-v$loores[row(v$loores)<col(v$loores)]
  tloores<-t(v$loores)
  r2<-tloores[row(tloores)<col(tloores)]
  resi<-v$loores-v$loores
  resi[row(resi)<col(resi)]<-v$res
  resi<-resi+t(resi)

  n0<-length(v$lofa)
  xn<-v$data[,c(2,1)]
  xy<-matrix(0,n,4)
  te<-crossprod(matrix(1,1,n0),t(xn[,1]))
  xy[,1]<-te[row(te)<col(te)]
  te<-crossprod(t(xn[,1]),matrix(1,1,n0))
  xy[,2]<-te[row(te)<col(te)]
  te<-crossprod(matrix(1,1,n0),t(xn[,2]))
  xy[,3]<-te[row(te)<col(te)]
  te<-crossprod(t(xn[,2]),matrix(1,1,n0))
  xy[,4]<-te[row(te)<col(te)]
  
  if(!missing(xyi)){
   ix<-ceiling(sqrt((2*(xyi)+0.25))-0.5)+1
   iy<-(xyi)-ceiling(sqrt((2*(xyi)+0.25))-0.5)/2*(ceiling(sqrt((2*(xyi)+0.25))-0.5)-1)
   nl<-n
   #*(n-1)/2
   ind1<-ceiling(sqrt((2*(1:nl)+0.25))-0.5)+1
   ind2<-(1:nl)-ceiling(sqrt((2*(1:nl)+0.25))-0.5)/2*(ceiling(sqrt((2*(1:nl)+0.25))-0.5)-1)
  }

  paro<-par(no.readonly=TRUE)
  par(mfrow=c(2,2), mar = c(3,2,2,1)+.1, bg="white")

#EJP: moved (1,1) to be first plotted:
  #graphsheet(win.width=0.5,win.height=0.5,win.left=0,win.top=0,Name="1) MAP")
  # plot map view as left plot in first row
  dg<-1:length(v$lofa)
  for(i in 1:length(v$lofa)){
    dg[i]<-sum((0.822*(gamsph(v$h,c(v$tha[i,1],v$tha[i,2],v$tha[i,3]))^0.25-gamsph(v$h,v$pars)^0.25))^2)
  }
  #if(! exists("zmv")) zmv<-0
  #if(zmv==0) plot(xn[,2],xn[,1],xlim=c(-0.6,1.1),ylim=c(-0.28, 0.28))
  #if(zmv==0) plot(xn[,2],xn[,1],xlim=c(-1.1,1.1),ylim=c(-0.28, 0.28),lwd=3)
  if (zmv==0 && !is.null(region))
    plot(xn[,2],xn[,1],xlim=c(min(region[,1]), max(region[,1])),
		ylim=c(min(region[,2]), max(region[,2])),lwd=3, 
		asp = 1, xlab = "", ylab = "")
  else
    plot(xn[,2],xn[,1],xlim=c(min(xn[,2]), max(xn[,2])),
		ylim=c(min(xn[,1]), max(xn[,1])),lwd=3, 
		asp = 1, xlab = "", ylab = "")
  #lolo<-lof[1,1]-lofa
  z<-xn[,1]
  if(zmv>0){
    z<-switch(zmv,v$data[,3],v$cda,v$lofa,dg)
    inc<-0.25
    rmin<-0.03
    epsi<-(max(z)-min(z))/(inc/rmin-1)
  #  symbols(xn[,2],xn[,1],circles=z-min(z)+epsi,inches=inc,xlim=c(-0.63,1.14),ylim=c(-0.3, 0.28))
    if(zmv>0 && !is.null(region))
      symbols(xn[,2],xn[,1],circles=z-min(z)+epsi,inches=inc,
	  	xlim=c(min(region[,1]), max(region[,1])),
		ylim=c(min(region[,2]), max(region[,2])),lwd=3)
    else
      symbols(xn[,2],xn[,1],circles=z-min(z)+epsi,inches=inc,
	  	xlim=c(min(xn[,2]), max(xn[,2])),
		ylim=c(min(xn[,1]), max(xn[,1])),lwd=3)
    #  symbols(xn[,2],xn[,1],circles=z-min(z)+epsi,inches=inc,lwd=3)
  }
  if(!is.null(region)) polygon(region[,1],region[,2],density=0,col=2)
    title(paste("Map View",switch(zmv+1,'','(y)',"(Cook's Distance)",'(Mahalanobis Distance)',"(Cook's Distance)")))
  #gsdmpv<-dev.cur()

  if(!missing(xyi)){
   segments(xy[xyi,3],xy[xyi,1],xy[xyi,4],xy[xyi,2],pch=16,col=3,lwd=3)
   points(xy[xyi,3],xy[xyi,1],pch=16,col=6)
   points(xy[xyi,4],xy[xyi,2],pch=16,col=8)
#   identify(xn[,2],xn[,1],plot=T,pts=cbind(xn[ix,2],xn[ix,1]))
   text(xn[ix,2],xn[ix,1]-(max(z)-min(z))/10,paste(ix))
#   identify(xn[,2],xn[,1],plot=T,pts=cbind(xn[iy,2],xn[iy,1]))
   text(xn[iy,2],xn[iy,1]-(max(z)-min(z))/10,paste(iy))
  }
  assign("p11", par(no.readonly=TRUE), envir = .varDiagOptions)
  # EJP-end

  #graphsheet(win.width=0.8,win.height=1,win.left=0,win.top=0,Name="Interactive Variogram Plot")
  #windows(width = 8, height = 5.5,rescale="R")
  #windows()
  par(mfg=c(2,1))
  #graphsheet(win.width=0.5,win.height=0.5,win.left=0,win.top=0.5,Name="3) LOO")
  plot(matrix(cbind(v$res,v$res),n*2,1), matrix(cbind(r1,r2),n*2,1),
  	pch=1, xlab="", ylab="", lwd=1)
  lines(c(min(v$res),max(v$res)),c(min(v$res),max(v$res)),col=8)
  segments(v$res,r1,v$res,r2)
  title("Leave One Out Residuals")
  if(!missing(xyi)){
    print("xyi")
    print(xyi)
   points(v$res[xyi],r1[xyi],pch=18,col=3)
   points(v$res[xyi],r2[xyi],pch=18,col=5)
   points(t(resi[ix,-ix]),t(v$loores[ix,-ix]),pch=16,col=6)
   points(t(resi[iy,-iy]),t(v$loores[iy,-iy]),pch=16,col=8)
   segments(v$res[xyi],r1[xyi],v$res[xyi],r2[xyi],col=3,lwd=5)
  }
  assign("p21", par(no.readonly=TRUE), envir = .varDiagOptions)

  cv1<-cv
  i<-1:n
  di<-dim(v$cg)[1]
  if(!missing(xyi)){
#   di<-dim(v$cg)[1]
#   pm<-diag(di)
#   pm[xyi,]<-diag(di)[di,]
#   pm[di,]<-diag(di)[xyi,]
#   cg1<-pm%*%v$cg%*%pm
#   i[n]<-xyi
#   i[xyi]<-n
#   print(max(abs(cv1-cv)))
   i<-c(sample(seq(di)[-xyi]),xyi)
   cg1<-v$cg[i,i]
   infm<-solve(cg1);
   cv1<-chol((infm+t(infm))/2);
  }
  par(mfg=c(2,2))
  #graphsheet(win.width=0.5,win.height=0.5,win.left=0.5,win.top=0.5,Name="4) DCR")
  x<-((2^0.25*gamma(0.75)/sqrt(pi))*gamsph(v$h,v$pars)^0.25)[i]
  y<-v$res[i]
  cv1<-cv1/cv1[di,di]
  plot(cv1%*%x,cv1%*%y,xlab="",ylab="",lwd=1)
  if(!missing(xyi))
    points(x[n],y[n],pch=16,col=3)
  #sm<-lowess(cv1%*%x,cv1%*%y)
  #lines(sm$x,sm$y,lwd=3)
  glu<-min(cv1%*%x)
  glo<-max(cv1%*%x)
  lines(c(glu,glo),c(0,0))
  title("Decorrelated Residuals")
  assign("p22", par(no.readonly=TRUE), envir = .varDiagOptions)

  xv<-seq(0.0001,max(v$h),0.01)
  par(mfg=c(1,2))
  #graphsheet(win.width=0.5,win.height=0.5,win.left=0.5,win.top=0,Name="2) SVC")
  plot(v$h,gdd,xlab="",ylab="",lwd=1)
  lines(xv,(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(xv,v$pars)^0.25,lwd=3)
  title("Studentized Square Root Cloud")
  if(!missing(xyi)){
   points(v$h[ind1==ix | ind2 == ix],gdd[ind1==ix | ind2 == ix],pch=16,col=6)
   points(v$h[ind1==iy | ind2 == iy],gdd[ind1==iy | ind2 == iy],pch=16,col=8)
   points(v$h[xyi],gdd[xyi],pch=16,col=3)
  }

  assign("p12", par(no.readonly=TRUE), envir = .varDiagOptions)

  par(paro)
  n
}

CookRLF.varobj<-function(v){
n<-length(v$lofa)
lofa<-matrix(0,n,1)
i1<-matrix(1:n,n,n)
for (k in 1:n){
    ii<-(i1[row(i1)<col(i1)]==k)|(t(i1)[row(t(i1))<col(t(i1))]==k)
    cgt<-v$cg[!ii,!ii]
    rt<-v$y[!ii]^0.25-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(v$h[!ii],v$tha[k,])^0.25
    lofa[k]<-rt%*%solve(cgt,rt)
}
dg<-1:length(v$lofa)
for(i in 1:length(v$lofa)){
    dg[i]<-sum((0.822*(gamsph(v$h,c(v$tha[i,1],v$tha[i,2],v$tha[i,3]))^0.25-gamsph(v$h,v$pars)^0.25))^2)
}
plot((v$lof[1]-lofa)/v$lof[1]*187/19,dg/(3*v$lof[1])*187,ylab="Cook's Distance",xlab="Reduktion im Lack of Fit")
identify((v$lof[1]-lofa)/v$lof[1]*187/19,dg/(3*v$lof[1])*187)
}


QQVarcloud.varobj<-function(v){

n<-length(v$h)
infm<-solve(v$cg);
cv<-chol((infm+t(infm))/2);
plot(qnorm(seq(from=1/(2*n),length=n,by=1/n)),sort(cv%*%v$y),
	xlab="quantile of standard normal distribution",
	ylab="orderd decorrelated residual")
lines(c(-3,3),c(-3,3),col=8,lwd=3)
apply(t(apply(apply(matrix(rnorm(n*100),ncol=100),2,sort),1,quantile,probs=c(0.05/n,1-0.05/n))),2,lines,x=qnorm(seq(from=1/(2*n),length=n,by=1/n)))
}


QQDecorr.varobj<-function(v){

n<-length(v$h)
infm<-solve(v$cg);
cv<-chol((infm+t(infm))/2);
plot(qchisq(seq(from=1/(2*n),length=n,by=1/n),1),sort(v$y/gamsph(v$h,v$pars)),
	xlab="quantile of Chi-square distribution",
	ylab="ordered value of [Z(s)-Z(s')]^2/(2g(s-s'))")
apply(t(apply(apply(matrix(rchisq(n*100,1),ncol=100),2,sort),1,quantile,probs=c(0.05/n,1-0.05/n))),2,lines,x=qchisq(seq(from=1/(2*n),length=n,by=1/n),1))
lines(c(0,8),c(0,8),col=8,lwd=3)
}

interact.varobj<-function(v,region=NULL,g="s",pchi=0.05,zmv=0){
#Identifikation in studentisierter VC

palette(c("black","cyan","magenta","green3","yellow","blue","white","red"))

n<-length(v$h)
infm<-solve(v$cg);
cv<-chol((infm+t(infm))/2);
XM<-cbind(gamsph1(v$h,v$pars),gamsph2(v$h,v$pars),gamsph3(v$h,v$pars))*(gamsph(v$h,v$pars))^(-0.75)/4
Vare<-v$cg-XM%*%solve(t(XM)%*%solve(v$cg,XM),t(XM))
#sig<-mean(sqrt(diag(Vare)))
e<-v$res
sig<-sqrt(sum(e^2)/(n-3))
gdd<-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(v$h,v$pars)^0.25+e*sig/sqrt(diag(Vare))
xn<-v$data[,c(2,1)]
r1<-v$loores[row(v$loores)<col(v$loores)]
tloores<-t(v$loores)
r2<-tloores[row(tloores)<col(tloores)]
resi<-v$loores-v$loores
resi[row(resi)<col(resi)]<-v$res
resi<-resi+t(resi)

n0<-length(v$lofa)
xn<-v$data[,c(2,1)]
xy<-matrix(0,n,4)
te<-crossprod(matrix(1,1,n0),t(xn[,1]))
xy[,1]<-te[row(te)<col(te)]
te<-crossprod(t(xn[,1]),matrix(1,1,n0))
xy[,2]<-te[row(te)<col(te)]
te<-crossprod(matrix(1,1,n0),t(xn[,2]))
xy[,3]<-te[row(te)<col(te)]
te<-crossprod(t(xn[,2]),matrix(1,1,n0))
xy[,4]<-te[row(te)<col(te)]

if(g=="l"){
   par(mfrow=c(2,2), mfg=c(2,1))
   par(get("p21", envir = .varDiagOptions))
   par(fig=c(0,0.5,0,0.5))
   xyi<-identify(matrix(cbind(v$res,v$res),n*2,1),matrix(cbind(r1,r2),n*2,1),plot=FALSE,n=1)
   if(xyi>n) xyi<-xyi-n
}

if(g=="m"){
   par(mfrow=c(2,2), mfg=c(1,1))
   par(get("p11", envir = .varDiagOptions))
   par(fig=c(0,0.5,0.5,1))
   ix0<-identify(xn[,2],xn[,1],plot=TRUE,n=1)
   points(xn[ix0,2],xn[ix0,1],pch=16,col=6)

      ind1<-ceiling(sqrt((2*(1:n)+0.25))-0.5)+1
      ind2<-(1:n)-ceiling(sqrt((2*(1:n)+0.25))-0.5)/2*(ceiling(sqrt((2*(1:n)+0.25))-0.5)-1)
      par(get("p12", envir = .varDiagOptions))
      par(mfrow=c(2,2),mfg=c(1,2),fig=c(0.5,1,0.5,1))
#      par(mfg=c(1,2,2,2))
#      par(fig=c(0.5,1,0.5,1))
      points(v$h[ind1==ix0 | ind2 == ix0],gdd[ind1==ix0 | ind2 == ix0],pch=16,col=6)
      par(mfrow=c(2,2), mfg=c(2,1))
      #par(p21)
      par(get("p21", envir = .varDiagOptions))
      par(fig=c(0,0.5,0,0.5))
      points(t(resi[ix0,-ix0]),t(v$loores[ix0,-ix0]),pch=16,col=6)

      par(mfrow=c(2,2), mfg=c(1,1))
   par(get("p11", envir = .varDiagOptions))
   par(fig=c(0,0.5,0.5,1))
   iy0<-identify(xn[,2],xn[,1],plot=FALSE,n=1)
   if(length(iy0)>0){
      ix<-max(ix0,iy0)
      iy<-min(ix0,iy0)
      xyi<-(ix-1)*(ix-2)/2+iy}
   else{
      xyi<-0
#      dev.off()
#      PlotDiag.varobj(v,region,zmv=zmv)
#      par(mfrow=c(2,2))
#      par(mfg=c(1,1,2,2))
#      par(p11)
#      par(fig=c(0,0.5,0.5,1))
#      points(xn[ix0,2],xn[ix0,1],pch=16,col=6)
#      identify(xn[,2],xn[,1],plot=T,pts=cbind(xn[ix0,2],xn[ix0,1]))
      ind1<-ceiling(sqrt((2*(1:n)+0.25))-0.5)+1
      ind2<-(1:n)-ceiling(sqrt((2*(1:n)+0.25))-0.5)/2*(ceiling(sqrt((2*(1:n)+0.25))-0.5)-1)
      par(mfrow=c(2,2), mfg=c(1,2))
      #par(p12)
      par(get("p12", envir = .varDiagOptions))
      par(fig=c(0.5,1,0.5,1))
      points(v$h[ind1==ix0 | ind2 == ix0],gdd[ind1==ix0 | ind2 == ix0],pch=16,col=6)
      par(mfrow=c(2,2), mfg=c(2,1))
      #par(p21)
      par(get("p21", envir = .varDiagOptions))
      par(fig=c(0,0.5,0,0.5))
      points(t(resi[ix0,-ix0]),t(v$loores[ix0,-ix0]),pch=16,col=6)
      }
}

if(g=="s"){
   par(mfrow=c(2,2), mfg=c(1,2))
   #par(p12)
   par(get("p12", envir = .varDiagOptions))
   par(fig=c(0.5,1,0.5,1))
   xyi<-identify(v$h,gdd,plot=FALSE,n=1)
}

if(g=="t"){
   par(mfrow=c(2,2), mfg=c(1,2))
   #par(p12)
   par(get("p12", envir = .varDiagOptions))
   par(fig=c(0.5,1,0.5,1))
   p<-locator(n=500,type="l",col=4)
   m<-length(p$x)
   lines(p$x[c(m,1)],p$y[c(m,1)],col=2)
   i<-t(outer(v$h,p$x,FUN="-"))/(p$x[c(2:m,1)]-p$x)
   gt<-apply(t((i*(p$y[c(2:m,1)]-p$y)+p$y))>=gdd&t((i>=0)&(i<=1)),1,"sum")
   s<-apply(t((i*(p$y[c(2:m,1)]-p$y)+p$y))<=gdd&t((i>=0)&(i<=1)),1,"sum")
   i0<-(s%%2)|(gt%%2)
   dev.off()
   PlotDiag.varobj(v,region)
   par(mfg=c(1,2))
   #par(p12)
   par(get("p12", envir = .varDiagOptions))
   par(fig=c(0.5,1,0.5,1))
   points(v$h[i0],gdd[i0],pch=16,col=3)
   polygon(p,density=0,col=4)

   par(mfg=c(1,1))
   par(get("p11", envir = .varDiagOptions))
   par(fig=c(0,0.5,0.5,1))
   segments(xy[i0,3],xy[i0,1],xy[i0,4],xy[i0,2],pch=16,col=3,lwd=3)

   xyi<-0
}

if(g=="x"){
   par(mfrow=c(2,2), mfg=c(1,2))
   #par(p12)
   par(get("p12", envir = .varDiagOptions))
   par(fig=c(0.5,1,0.5,1))

   i0<-(gdd-(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(v$h,v$pars)^0.25)/sig>qnorm(1-pchi/2)
   i0a<-(-gdd+(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(v$h,v$pars)^0.25)/sig>qnorm(1-pchi/2)
   dev.off()
   PlotDiag.varobj(v,region,zmv=zmv)
   par(mfg=c(1,2))
   #par(p12)
   par(get("p12", envir = .varDiagOptions))
   par(fig=c(0.5,1,0.5,1))
   points(v$h[i0],gdd[i0],pch=16,col=3)
   points(v$h[i0a],gdd[i0a],pch=16,col=4)

   xv<-seq(0.0001,max(v$h),0.01)

#   lines(xv,gamsph(xv,v$pars)*qchisq(1-pchi,1),lty=4,lwd=2)
   lines(xv,(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(xv,v$pars)^0.25+sig*qnorm(1-pchi/2),lty=4,lwd=2)
   lines(xv,(2^0.25*gamma(0.75)/sqrt(pi))*gamsph(xv,v$pars)^0.25-sig*qnorm(1-pchi/2),lty=4,lwd=2)
   par(get("p11", envir = .varDiagOptions))

   par(mfg=c(1,1))
   par(get("p11", envir = .varDiagOptions))
print(xy[i0,])
   par(fig=c(0,0.5,0.5,1))
   segments(xy[i0,3],xy[i0,1],xy[i0,4],xy[i0,2],pch=16,col=3,lwd=2)
   segments(xy[i0a,3],xy[i0a,1],xy[i0a,4],xy[i0a,2],pch=16,col=4,lwd=2)

   xyi<-0
}

if(g=="n"){
   par(mfrow=c(2,2), mfg=c(1,1))
   par(get("p11", envir = .varDiagOptions))
   par(fig=c(0,0.5,0.5,1))
   p<-locator(n=500,type="l",pch=16,col=4)
   m<-length(p$x)
   lines(p$x[c(m,1)],p$y[c(m,1)],col=2)
   i<-t(outer(xn[,2],p$x,FUN="-"))/(p$x[c(2:m,1)]-p$x)
   gt<-apply(t((i*(p$y[c(2:m,1)]-p$y)+p$y))>=xn[,1]&t((i>=0)&(i<=1)),1,"sum")
   s<-apply(t((i*(p$y[c(2:m,1)]-p$y)+p$y))<=xn[,1]&t((i>=0)&(i<=1)),1,"sum")
   i0<-(s%%2)|(gt%%2)
   nl<-length(v$h)
   ind1<-ceiling(sqrt((2*(1:nl)+0.25))-0.5)+1
   ind2<-(1:nl)-ceiling(sqrt((2*(1:nl)+0.25))-0.5)/2*(ceiling(sqrt((2*(1:nl)+0.25))-0.5)-1)
   i00<-match(ind1,(1:n0)[i0],nomatch=FALSE)&match(ind2,(1:n0)[i0],nomatch=FALSE)
   dev.off()
   PlotDiag.varobj(v,region)
   par(mfg=c(1,2))
   #par(p12)
   par(get("p12", envir = .varDiagOptions))
   par(fig=c(0.5,1,0.5,1))
   points(v$h[i00],gdd[i00],pch=16,col=3)

   par(mfg=c(1,1))
   par(get("p11", envir = .varDiagOptions))
   par(fig=c(0,0.5,0.5,1))
   polygon(p,density=0,col=4)
   segments(xy[i00,3],xy[i00,1],xy[i00,4],xy[i00,2],pch=16,col=3,lwd=3)
   xyi = 0
}

#print(xyi)

if(g!="t"&g!="x"&g!="n"& xyi>0){
   dev.off()
   PlotDiag.varobj(v,region,xyi=xyi,zmv=zmv)
}

xyi}

