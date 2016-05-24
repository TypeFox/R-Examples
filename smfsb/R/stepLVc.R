
stepLVc <- function(x0,t0,deltat,th=c(1,0.005,0.6))
{
  new=.C("stepLV",as.integer(x0),as.double(t0),as.double(deltat),as.double(th))[[1]]
  names(new)=names(x0)
  new
}

# eof

