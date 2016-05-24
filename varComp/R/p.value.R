p.value <-
function(x) UseMethod('p.value')

p.value.htest = function(x) x$p.value
p.value.default = p.value.htest

p.value.varComp.test <-
function(x)
{
#Extracting p-values from objects from varComp.test. 
#i.	x: an object of class varComp.test
  f=function(zz){
    if(inherits(zz, 'htest')) return(zz$p.value)
    if(is.list(zz)) return(unlist(sapply(zz, f)))
    NA
  }
  p=drop(unlist(sapply(x, f)))
  
  cn=function(zz){  
    if(inherits(zz, 'htest')) return(NULL)
    nn=names(zz)
    ans=character(0L)
    for(i in seq_along(nn)) {
      tmp=cn(zz[[i]])
      ans=c(ans, if(is.null(tmp)) nn[i] else paste(nn[i], tmp, sep='.'))
    }
    ans
  }
  
  names(p)=cn(x)
  p
}

p.value.varCompFixEf=function(x)
{
	Aov=attr(x, 'anova')
	ans=Aov[, 'Pr(>|t|)']
	attr(ans, 'Overall') = attr(Aov, 'Overall')[,'Pr(>F)']
	ans
}