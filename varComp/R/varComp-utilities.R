is.formula=function(x)inherits(x, 'formula')

get.seed <-
function(){
#Retrieve the current random number seed the RNGkind() for reproducibility. 
  if(!exists('.Random.seed', envir=globalenv(), mode='numeric', inherits=FALSE)) runif(1L)
  seed=get('.Random.seed', envir=globalenv(), mode='numeric', inherits=FALSE)
  attr(seed, 'RNGkind')=RNGkind()
  seed
}

safeseq=function(from=1L, to=1L, by=1L,...)
{
  disc=by*(from-to)
  if(disc>0){
    vector(class(disc), 0L)
  }else seq(from=from, to=to, by=by, ...)
}

if(FALSE){# mget replaces this
char2list=function(Names)
{
	pfm=parent.frame(); pfm
	ans=lapply(Names, get, envir=pfm)
	names(ans) = Names
	ans
}

}
