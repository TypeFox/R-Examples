`RectDense` <-
function( INx, INy, icut=1, u = par("usr"), ndivs=10)
{
#######  given a set of points (x,y)
  ##  find rectangles that have greater than icut
###  in those rectangles
  if(missing(icut)) { icut=1 }
  if(missing(ndivs)) { ndivs=10 }
  if(missing(u)) { u = par("usr") }
  
  ex1 = seq(from=u[1], to=u[2], length=ndivs)
  why1 = seq(from=u[3], to=u[4], length=ndivs)

  intX = findInterval( INx, ex1)
  intY = findInterval( INy, why1)

  flag = INx>ex1[1] & INx<ex1[length(ex1)] & INy>why1[1] &  INy<why1[length(why1)]

  ucorns = unique( cbind(intX[flag], intY[flag]))

  ccorns = cbind(ex1[ucorns[,1]], why1[ucorns[,2]], ex1[ucorns[,1]+1], why1[ucorns[,2]+1])

  LENS = apply(ccorns, 1, function(x) length(INx[INx>x[1]&INy>x[2] & INx<x[3]& INy<x[4]]))

  ipass = which(LENS>=icut)

  return(list( icorns=ccorns[ipass, ], ilens=LENS[ipass] , ipass=ipass  ,corners=ccorns, lens=LENS))
}

