
as.timedData <- function(timeseries)
{
  mat=as.matrix(timeseries)
  rownames(mat)=time(timeseries)
  attr(mat,"class")<-NULL
  attr(mat,"tsp")<-NULL
  mat
}

# eof
