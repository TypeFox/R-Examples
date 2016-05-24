#####################################
##  OTOLITH SHAPE ANALYSIS 
##      Morphological functions
##      by Lisa Anne Libungan
#####################################


.shapeR.check.outline.list = function(object)
{
  if( length(object@outline.list) == 0 || (length(object@outline.list) ==1  && length(object@outline.list[[1]]) ==0 ))
    stop("No outlines. Need to run detect.outline?")
}

.shapeR.check.shape.coefs = function(object)
{
  if( length(object@shape.coef.raw) == 0 )
    stop("No shape coefficients. Need to run getShapeCoefficients?")
}

.shapeR.check.shape.coefs.std = function(object)
{
  if( length(object@wavelet.coef.std) == 0 || length(object@fourier.coef.std) == 0)
    stop("Coefficients have not been standardized. Need to run stdCoefs?")
}

.shapeR.rbind.fill.matrix = function (...) 
{
  matrices <- list(...)
  if (length(matrices) == 0) 
    return()
  if (is.list(matrices[[1]]) && !is.matrix(matrices[[1]])) {
    matrices <- matrices[[1]]
  }
  tmp <- unlist(lapply(matrices, is.factor))
  if (any(tmp)) {
    stop("Input ", paste(which(tmp), collapse = ", "), " is a factor and ", 
         "needs to be converted first to either numeric or character.")
  }
  matrices[] <- lapply(matrices, as.matrix)
  lcols <- lapply(matrices, function(x) dimnames(x)[[2]])
  cols <- unique(unlist(lcols))
  rows <- unlist(lapply(matrices, nrow))
  nrows <- sum(rows)
  output <- matrix(NA, nrow = nrows, ncol = length(cols))
  colnames(output) <- cols
  pos <- matrix(c(cumsum(rows) - rows + 1, rows), ncol = 2)
  for (i in seq_along(rows)) {
    rng <- seq(pos[i, 1], length = pos[i, 2])
    output[rng, lcols[[i]]] <- matrices[[i]]
  }
  output
}


#######################################################################
#
# Based on Claude, J. 2008. Morphometrics with R. Springer, New York.
# Error in the function in the book, wrote a new function
# See: https://sites.google.com/site/raduiovita/morphometrics-software
# f5.10. this function has been slightly modified 
#
######################################################################

.shapeR.NEF<-function(M, n=dim(M)[1]/2,start=F)
{
  ef<-.shapeR.efourier(M,n)
  A1<-ef$an[1]; B1<-ef$bn[1]
  C1<-ef$cn[1]; D1<-ef$dn[1]
  theta<-(0.5*atan(2*(A1*B1+C1*D1)/(A1^2+C1^2-B1^2-D1^2)))%%pi
 
  phaseshift<-matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
  M2<-matrix(c(A1,C1,B1,D1),2,2)%*%phaseshift
  v<-apply(M2^2,2, sum)
  if (v[1]<v[2]){theta<-theta+pi/2} 
  theta<-(theta+pi/2)%%pi-pi/2
 
  Aa<-A1*cos(theta)+B1*sin(theta)
  Cc<-C1*cos(theta)+D1*sin(theta)
  scale<-sqrt(Aa^2+Cc^2)
  psi<-atan(Cc/Aa)
  size<-(1/scale)
  rotation<-matrix(c(cos(psi),-sin(psi),sin(psi),cos(psi)),2,2)
  A<-B<-C<-D<-numeric(n)
  if (start){theta<-0}
  for (i in 1:n){
    mat<-size*rotation%*%matrix(c(ef$an[i],ef$cn[i],ef$bn[i],ef$dn[i]),2,2)%*%matrix(c(cos(i*theta),sin(i*theta),-sin(i*theta),cos(i*theta)),2,2) #this resizes the first harmonic so its major axis is equal to 1; this is the classical way to standardize (normalize) EF coefficients, but may not be appropriate 
    #mat<-rotation%*%matrix(c(ef$an[i],ef$cn[i],ef$bn[i],ef$dn[i]),2,2)%*%matrix(c(cos(i*theta),sin(i*theta),-sin(i*theta),cos(i*theta)),2,2) #without the size correction, just the rotation
    A[i]<-mat[1,1]
    B[i]<-mat[1,2]
    C[i]<-mat[2,1]
    D[i]<-mat[2,2]}
  list(A=A,B=B,C=C,D=D,size=scale,theta=theta,psi=psi,ao=ef$ao,co=ef$co)
}

########### Elliptic Fourier ######################

#Based on Claude, J. 2008. Morphometrics with R. Springer, New York.
.shapeR.efourier<-function(M, n=dim(M)[1]/2)
{
  p<-dim(M)[1]
  Dx<-M[,1]-M[c(p,(1:p-1)),1]
  Dy<-M[,2]-M[c(p,(1:p-1)),2]
  Dt<-sqrt(Dx^2+Dy^2)
  # Error in this function. Same coordinate, side by side, divided by zero
  # Fixed with this
  Dt[which(Dt==0)] = 1
  t1<-cumsum(Dt)
  t1m1<-c(0, t1[-p])
  T<-sum(Dt)
  an<-bn<-cn<-dn<-numeric(n)
  for (i in 1:n){
    an[i]<- (T/(2*pi^2*i^2))*sum((Dx/Dt)*
                                   (cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
    bn[i]<- (T/(2*pi^2*i^2))*sum((Dx/Dt)*
                                   (sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))
    cn[i]<- (T/(2*pi^2*i^2))*sum((Dy/Dt)*
                                   (cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
   dn[i]<- (T/(2*pi^2*i^2))*sum((Dy/Dt)*
                                  (sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))
  }
  ao<-2*sum(M[,1]*Dt/T)
  co<-2*sum(M[,2]*Dt/T)
  list(ao=ao,co=co,an=an,bn=bn,cn=cn,dn=dn)
}

#As in Claude, J. 2008. Morphometrics with R. Springer, New York.
.shapeR.iefourier<-function(an,bn,cn,dn,k,n,ao=0,co=0)
{
  theta<-seq(0,2*pi, length=n+1)[-(n+1)]
  harmx <- matrix (NA, k, n)
  harmy <- matrix (NA, k, n)
  for (i in 1:k){
    harmx[i,]<-an[i]*cos(i*theta)+bn[i]*sin(i*theta)
    harmy[i,]<-cn[i]*cos(i*theta)+dn[i]*sin(i*theta)
  }
  x<-(ao/2) + apply(harmx, 2, sum)
  y<-(co/2) + apply(harmy, 2, sum)
  list(x=x, y=y)
}


.shapeR.inverse.wavelet = function(coefs,m.o,plotDetail=T,return.polar=F,num.points=512*2,n.levels=5)
{
  wwd = wd(rep(0,num.points))
  
  for(level in 0:(log2(num.points)-1))
    wwd = putD(wwd,level=level,v=rep(0,2^level))
  
  wwd = putD(wwd,level=0,v=coefs[1])
  wwd = putD(wwd,level=1,v=coefs[2:3])
  
  for(level in 2:n.levels)
    wwd = putD(wwd,level=level,v=coefs[(sum(2^(c(1:(level-1)))):sum(2^(c(1:level))))[-1]+1])
  
  rad=wr(wwd)+m.o
  angle = seq(-pi,pi,by=2*pi/num.points)[-1]
  
  if(return.polar==F)
    list(X=rad*cos(angle),Y=rad*sin(angle))
  else
    list(angle=angle,radii=rad)
}


.shapeR.resizeImage = function(im, size.ratio) {
  # function to resize an image 
  # im = input image, w.out = target width, h.out = target height
  # Bonus: this works with non-square image scaling.
  
  # initial width/height
  w.in = nrow(im)
  h.in = ncol(im)
  
  w.out = floor(w.in*size.ratio)
  h.out = floor(h.in*size.ratio)
  # Create empty matrix
  im.out = matrix(rep(0,w.out*h.out), nrow =w.out, ncol=h.out )
  
  
  # Do resizing -- select appropriate indices
  im.out <- im[ floor(size.ratio* 1:w.out), floor(size.ratio* 1:h.out)]

  return(im.out)
}

.shapeR.find.outline=function(skra,threshold=.15,mouse.click=F,main=NA,display.images=F)
{
  #M<- read.pnm(skra) # skra innan ""
  M<- readJPEG(skra) # skra innan ""
  if(length(dim(M))>2){
    M<- suppressWarnings(pixmapRGB(M[,,1:3]))
    M<- as(M, "pixmapGrey")
  }
  else{
    M<- suppressWarnings(pixmapGrey(M))
  }
  M@grey[which(M@grey<=threshold)]<-0
  M@grey[which(M@grey>threshold)]<-1
  
  if(mouse.click== T || display.images==T)
    plot(M,main=main)
  
  start=list(x=NA,y=NA)
  if(mouse.click==T)
  {
    # uses coordinates to start, need to select one spot 
    start<-locator(1)
  }
  else{
    start$x=M@size[2]/2
    start$y=M@size[1]/2
  }

  Rc<-.shapeR.Conte(c(round(start$x),round(start$y)),M@grey)
  if(mouse.click == T || display.images==T)
    lines(Rc$X, Rc$Y, lwd=2,col="red")
  return(Rc)
}

# From Claude, J. 2008. Morphometrics with R. Springer, New York.
.shapeR.Conte=function(x, imagematrix)
{
  I<-imagematrix
  x<-rev(x)
  x[1]<-dim(I)[1]-x[1]
  while (abs(I[x[1],x[2]]-I[x[1],(x[2]-1)])<0.1)
  {x[2]<-x[2]-1;
  #if(x[2]<1)
  #  simpleError("Problem in detecting outline when collecting coordinates.")  
  }
  a<-1
  M<-matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),2,8,byrow=T)
  M<-cbind(M[,8],M,M[,1])
  X<-0; Y<-0;
  x1<-x[1]; x2<-x[2]
  SS<-NA; S<-6

  while ((any(c(X[a],Y[a])!=c(x1,x2) ) | length(X)<3))
  {
    
    if (abs(I[x[1]+M[1,S+1],x[2]+M[2,S+1]]-I[x[1],x[2]])<0.1)
    {
      a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,S+1]
      SS[a]<-S+1; S<-(S+7)%%8
    }
    else if (abs(I[x[1]+M[1,S+2],x[2]+M[2,S+2]]
               -I[x[1],x[2]])<0.1)
    {
      a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,S+2]
      SS[a]<-S+2; S<-(S+7)%%8
    }
    else if (abs(I[x[1]+M[1,(S+3)],x[2]+M[2,(S+3)]]
               -I[x[1],x[2]])<0.1)
    {
      a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,(S+3)]
      SS[a]<-S+3; S<-(S+7)%%8
    }
    else S<-(S+1)%%8
    
    if(a>(dim(I)[1]+dim(I)[2])*1e2)
    {
      X[a]=x1
      Y[a]=x2
    }
  }
  list(X=(Y[-1]), Y=((dim(I)[1]-X))[-1])
}

#############################################################
#
#  Function to get the maximum length, height, area etc. 
#  Use gpc.poly to get perimeter
#
#############################################################
.shapeR.rotate.svd=function(outline)
{
  M=cbind(outline$X,outline$Y)
  Ma=M%*%svd(var(M))$u
  Ma[,1]=-Ma[,1]
  return(Ma)
}

.shapeR.otolith.image.parameters=function(x)
{
  xbil=max(diff(range(x[,1])))
  ybil=max(diff(range(x[,2]))) 
  area = .shapeR.area(x)
  perimeter=.shapeR.fourier(x,16)$perimeter
  return(list(width=xbil,height=ybil,area=area,perimeter=perimeter))
}

.shapeR.area = function(m){
 n = dim(m)[1]
 
 x1 = m[,1]
 y1 = m[,2]

 x2 = c(m[-1,1],m[1,1])
 y2 = c(m[-1,2],m[1,2])

 tarea = abs(sum((x2-x1)*(y2+y1))/2)

 return(tarea)
}


.shapeR.get.otolith.image.parameters = function(outline){
  x = .shapeR.rotate.svd(outline)
  .shapeR.otolith.image.parameters(x)
}


.shapeR.fourier=function(M,n)
{
  p<-dim(M)[1]
  an<-numeric(n)
  bn<-numeric(n)
  tangvect<-M-rbind(M[p,],M[-p,])
  perim<-sum(sqrt(apply((tangvect)^2, 1, sum)))
  v0<-(M[1,]-M[p,])
  tet1<-Arg(complex(real=tangvect[,1],
                    imaginary = tangvect [,2]))
  tet0<-tet1[1]
  t1<-(seq(0, 2*pi, length=(p+1)))[1:p]
  phi<-(tet1-tet0-t1)%%(2*pi)
  ao<- 2* sum(phi)/p
  for (i in 1:n){
    an[i]<- (2/p) * sum( phi * cos (i*t1))
    bn[i]<- (2/p) * sum( phi * sin (i*t1))
  }
  list(ao=ao, an=an, bn=bn, phi=phi, t=t1,perimeter=perim, thetao=tet0)
}

#Claude, J. 2008. Morphometrics with R. Springer, New York., page 53
.shapeR.regularradius<-function(Rx, Ry, n)
{
  le<-length(Rx)
  M<-matrix(c(Rx, Ry), le,2)
  M1<-matrix(c(Rx-mean(Rx), Ry-mean(Ry)), le,2)
  V1<-complex(real=M1[,1], imaginary=M1[,2])
  M2<-matrix(c(Arg(V1), Mod(V1)), le,2)
  V2<-NA
  #The following code finds the indices of the nearest pixel on the outline using the angul
  #increment.
  for (i in 0:(n-1))
  {
    V2[i+1]<-which.max((cos(M2[,1]-2*i*pi/n)))
  }
  V2<-sort(V2)
 
  ord = order(M2[V2,1])
  list("pixindices"=V2[ord],"radii"=M2[V2[ord],2],"coord"=M1[V2[ord],],"angle"=M2[V2[ord],1])
}


####### Center images, find centroid
.shapeR.centroid=function(outline)
{
  M=outline
  M$X=outline$X-mean(outline$X)
  M$Y=outline$Y-mean(outline$Y)
  #plot(utl2$X,utl2$Y,type="l")
  return(M)
}


###########
#LLEONART (2000)
.shapeR.standard.coef=function(y.in,in.class,x,std.type="mean",is.na.classes,p.crit)
{
  if(is.na.classes)
    in.class = rep("1",length(x))
  #y=y.in+abs(y.in)+.1
  y=y.in+ifelse(min(y.in)<0,-min(y.in),0)+.1
  unique(in.class)->class
  mean.x=tapply(x,factor(in.class),mean)
  nc = length(class)
  if(std.type=="mean")
    mean.x = rep(mean(mean.x),nc)
  if(std.type=="min")
    mean.x = rep(min(mean.x),nc)
  if(std.type=="max")
    mean.x = rep(max(mean.x),nc)
  
  
  lm.y=list()
  for(i in 1:nc) lm.y[[i]]=lm(log(y[in.class==class[i]])~log(x[in.class==class[i]]))
  #Default regression set as 0, i.e. no standardization.
  b=rep(0,nc)
  for(i in 1:nc)
  {
    p.value.b = summary(lm.y[[i]])$coefficients["log(x[in.class == class[i]])","Pr(>|t|)"]
    #If p value from regression is less than 0.05, we keep the regression coefficient
    if(is.finite(p.value.b) &&  p.value.b<p.crit){
      #print(paste(class[i],"p.value=",p.value.b,"b=",lm.y[[i]]$coefficients[2]))    
      b[i] = lm.y[[i]]$coefficients[2]
    }
  }
  
  yy=y
  for(i in 1:nc)
  {
    ind = in.class==class[i]
    yy[ind]=y[ind]*(mean.x[i]/x[ind])^b[i]
  }
  return(yy)
  
  
}



.shapeR.coef.standardize.f = function(coef,classes,std.by,std.type="mean",is.na.classes=F,p.crit=0.05,bonferroni=TRUE,...)
{
  
  coef.standard = matrix(NA,nrow=dim(coef)[1],ncol=dim(coef)[2])

  k=1
  if(bonferroni)
    k = dim(coef)[2]
  
  r.y = c()
  for(i in 1:ncol(coef)) 
  {
    #ANCOVA
    if(is.na.classes){
      t.ancova = aov(coef[,i]~std.by)
      if(summary(t.ancova)[[1]]["std.by","Pr(>F)"]<p.crit/k)
      {
        r.y = c(r.y,i)
      }
    }
    else{
      t.ancova = aov(coef[,i]~classes*std.by)
      if(summary(t.ancova)[[1]]["classes:std.by","Pr(>F)"]<p.crit/k)
      {
        r.y = c(r.y,i)
      }
    }
    
    coef.standard[,i]= .shapeR.standard.coef(coef[,i],classes,std.by,std.type,is.na.classes,p.crit)
  }
  
  if(is.null(r.y))
  {
    message("No coefficients removed")
    return(list(coef.standard=coef.standard,r.y))
  }

  message("Removed coefficients: ",appendLF=F)
  message(paste(r.y,collapse=","))
  
  return(list(coef.standard=coef.standard[,-r.y],r.y=r.y))
  
}

.shapeR.variation.among = function(in.classes,coef)
{
  classes = factor(in.classes)

  mean.var=function(x) mean(tapply(x,classes,var))
  mean.var.among=function(x) var(tapply(x,classes,mean))
  
  a = length(levels(classes))
  na = tapply(coef[,1],classes,length)
  
  n0 = 1/(a-1)*(dim(coef)[1] - sum(na^2)/sum(na))

  var.within=apply(coef,2,mean.var)
  var.among=n0*apply(coef,2,mean.var.among)
  
  add.var.component=(var.among-var.within)/n0
  
  return(add.var.component/(add.var.component+var.within))
}

#####################################################################
#
# Based on Claude, J. 2008. Morphometrics with R. Springer, New York.
#
######################################################################
.shapeR.smoothout<-function(M.in, n)
{
  M = cbind(M.in$X,M.in$Y)
  p<-dim(M)[1]
  a<-0
  while (a<=n)
  {
    a<-a+1
    Ms<-rbind(M[p,],M[-p,])
    Mi<-rbind(M[-1,],M[1,])
    M<-M/2+Ms/4+Mi/4
  }
  return(data.frame(X=M[,1],Y=M[,2]))
}
