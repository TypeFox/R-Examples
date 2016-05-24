circum <- function(x,y){
  # interface to:
  #        SUBROUTINE CIRCUM (X1,Y1,X2,Y2,X3,Y3,RATIO, XC,YC,CR,
  #     .                   SA,AR)

  if(length(y)!=3 || length(x)!=3)
    stop("need exactly three points!")
  
  ret <- .Fortran("circum",
                  as.double(x[1]),
                  as.double(y[1]),
                  as.double(x[2]),
                  as.double(y[2]),
                  as.double(x[3]),
                  as.double(y[3]),
                  ratio=logical(1),
                  xc=double(1),
                  yc=double(1),
                  cr=double(1),
                  sa=double(1),
                  ar=double(1),
                  PACKAGE="tripack")
  list(x=ret$xc,y=ret$yc,radius=ret$cr,signed.area=ret$sa,aspect.ratio=ret$ar)
}       
