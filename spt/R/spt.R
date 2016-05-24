spt <- function(A,B){
  if(missing(A)||missing(B))
    stop("Sizes of at least two angles need to be specified.")
  if(A<=0 || B<=0 || A>=180 || B >=180 || A+B>=180)
    stop("Invaid angle size(s) specified: 0<A<180, 0<B<180, and 0<A+B<180.")
  A1 = A; B1 = B; C1 = 180.0-(A1+B1)
  triname = paste("SPT(",as.character(format(A1,20)),",",
    as.character(format(B1,20)),",",
    as.character(format(180-A1-B1)),")",sep='')
  angles = sort(c(A1,B1,C1),decreasing=TRUE)

  minA = angles[3]; maxA = angles[1];
  r = pi/180;  h=100;
  C = angles[3]*r; B = angles[1]*r; A = angles[2]*r;
  xa = 0; ya=0;
  xb = h/tan(A)+h/tan(B); yb=0;
  xc = h/tan(A); yc = h;
  A0 = 0; B0 = pi-B; C0 = pi+A;

  xmin = min(xa,xb,xc); xmax=max(xa,xb,xc); 
  ymin = 0; ymax=h;
  if(maxA > 90){
    delta = abs(h * sin(maxA))*.25
    xmin = xmin - delta;
    delta = (xmax-xmin-ymax+ymin)*.8;
    ymin = ymin - delta;
    ymax = ymax + delta;
  }
  viewport =c(xmin,ymin,xmax,ymax);
  if(maxA > 90) sptdim = NA
  else if(maxA==90) sptdim = 2.0
  else sptdim = .Fortran(.F_SptDimFP, as.double(angles),as.double(0))[[2]]

  out = structure(list(angles=angles, viewport=viewport,
    ABC = cbind(c(xa,xb,xc), c(ya,yb,yc)),
    tri = list(A=A,B=B,C=C,xa=xa,xb=xb,xc=xc,ya=ya,yb=yb,yc=yc,A0=A0,B0=B0,C0=C0),
    Dim = sptdim, data.name = triname), class = "spt")
}


print.spt <- function(x,...){
  tmp = x$ABC
  Angles = c(x$angles[3],x$angles[1],x$angles[2])
  tmp = data.frame(Angles=Angles,x=tmp[,1],y=tmp[,2])
  row.names(tmp) = c("A","B","C")
  cat("\n\t", x$data.name,"\n\n")
  cat("Dimension = ", x$Dim,
      "\tViewport=(",x$viewport[1],",",
      x$viewport[2],",",
      x$viewport[3],",",
      x$viewport[4],")\n\n")
  print(tmp)
  cat("\n")
}


plot.spt <- function(x,iter,tol=0.0001,main=NULL,...){
  xmin = x$viewport[1];
  ymin = x$viewport[2];
  xmax = x$viewport[3];
  ymax = x$viewport[4];
  maxA = x$angles[1]
  if(is.null(main)) main = x$data.name
  plot(0,0, type='n', bty='n', xaxt='n',yaxt='n',
       xlab='', ylab='', asp=1, main = main,
       xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  polygon(x$ABC[,1], x$ABC[,2]);
  if(missing(iter)) iter = 11
  if(iter < 1) stop("Iteration number can not be zero/negative.")
  Tri0 = x$ABC;
  tol = tol * (xmax-xmin) * (ymax-ymin)
  trilist = NULL # to save all resulted SPT (the triangle to be pasted)
  if(maxA == 90)
    .sptpedal2(iter,Tri0,tol)
  else if(maxA < 90)
    .sptpedal3(iter,Tri0,tol)
  else trilist = .sptpedal4(min(12,iter),Tri0,tol)
  invisible(trilist)
}
