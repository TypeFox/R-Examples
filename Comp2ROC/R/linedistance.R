linedistance <-
function(curve.fpr,curve.tpr,curve.segslope,curve.slope,line.slope,ref.point) {
  line.xint1=c()
  line.yint1=c()
  line.dist1=c()
  n.points=dim(curve.fpr)
  K=length(line.slope)
  for (i in 1:K)
  {
    for (j in 1:(n.points-1))
    {
      if (line.slope[i]<=curve.slope[n.points-j] && line.slope[i]>curve.slope[n.points-j+1])
      {
        if (curve.segslope[n.points-j]==-Inf)
          line.xint1[i]=curve.fpr[n.points-j]
        else
          line.xint1[i]=ref.point[1]+(curve.tpr[n.points-j]-ref.point[2]-curve.segslope[n.points-j]*curve.fpr[n.points-j]+curve.segslope[n.points-j]*ref.point[1])/(line.slope[i]-curve.segslope[n.points-j])
      }
    }
    line.yint1[i]=line.slope[i]*(line.xint1[i]-ref.point[1])-ref.point[2]
    line.dist1[i]=sqrt((line.xint1[i]-ref.point[1])^2+(line.yint1[i]-ref.point[2])^2)
  }
  answer=list(dist=line.dist1,x=line.xint1,y=line.yint1)
  return(answer)
  }
