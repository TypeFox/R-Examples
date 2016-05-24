wave.trans <- function(x, method="modwt",wf="la8",n.levels=4,boundary="periodic")
{

  if(method == "dwt"){
    x.wt <- dwt(x,wf,n.levels,boundary)
  }else{
    x.wt <- modwt(x,wf,n.levels,boundary)
  }

  x.bw <- brick.wall(x.wt,wf,method)
  return(x.bw)
}