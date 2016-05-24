desscop <-
function(copest,xname,yname,normalize=FALSE,axes=TRUE)
{
  dims=dim(copest$cop)
  dim=length(dims)
  subsampsize=dims[1]
  nabcisse=match(xname,copest$varnames)
  nordonnee=match(yname,copest$varnames)
  proj=as.numeric(apply(copest$cop,c(nabcisse,nordonnee),sum))
  if (normalize)
  {
    proj=(proj-min(proj))/(max(proj)-min(proj))
    ech=0.2
  } else ech=10
  if (axes)
  {
    labx=copest$varnames[nabcisse]
    laby=copest$varnames[nordonnee]
  } else
  {
    labx=''
    laby=''
  }
  symbols(rep(1:subsampsize,subsampsize),rep(1:subsampsize,each=subsampsize),
    circles=(proj)*ech,
    xlab=labx,ylab=laby,inches=FALSE)
}

