"mintvmon" <-
function(y,sigma=-1,DYADIC=TRUE,thresh=-1,method=2,MONCONST=TRUE,CONVCONST=FALSE)
{
n <- length(y)
if(thresh<0)
  {
  if(sigma < 0)
    sigma <- mad(y[-1]-y[-n])/sqrt(2)
  print(sigma)
  thresh <- sqrt(2*log(n))*sigma
  }
tmp <- .C("mintvmon", as.integer(n),
f=double(n),derivsign=integer(n),secsign=integer(n),as.integer(method),as.integer(DYADIC),as.double(thresh),as.double(y),as.integer(MONCONST),as.integer(CONVCONST),jact=integer(n),kact=integer(n),signact=integer(n),nact=integer(1),piecesleft=integer(n),piecesright=integer(n),PACKAGE="ftnonpar")

list(y=tmp$f,derivsign=tmp$derivsign,secsign=tmp$secsign,jact=tmp$jact[1:tmp$nact],kact=tmp$kact[1:tmp$nact],signact=tmp$signact[1:tmp$nact],pl=tmp$piecesleft[1:tmp$nact],pr=tmp$piecesright[1:tmp$nact])

}
