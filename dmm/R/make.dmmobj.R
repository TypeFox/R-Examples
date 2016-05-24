make.dmmobj <-
function(p=NULL,components=c("VarG(Ia)"),...){
#  make.dmmobj() - construct an object of class 'dmm' from user-supplied
#                  variance/covariance matrices

# make component tables
  ctable <- make.ctable()
# check components in ctable
  if(any(is.na(match(components,ctable$all)))){
    print(components)
    stop("Component(s) not recognized:\n")
  }
  
  if(is.null(p)) {
    stop("make.dmmobj(): phenotypic covariance matrix must be specified\n")
  }
# p <- nearPD(p)$mat
  p <- as.matrix(nearPD(p,ensureSymmetry=T)$mat)
  
  traits <- dimnames(p)[[1]]
  traitpairs <- permpaste(traits)
  l <- length(traits)
  comp <- list(...)
  vcomp <- length(comp)
  v <- length(components)

  siga <- matrix(0,v,l*l)
  for(i in 1:v){
    siga[i, ] <- as.vector(comp[[i]])
  }
  dimnames(siga) <- list(components,traitpairs)
# kludge an am object for siga.posdef
  am <- list(v=v,l=l)
# check siga positive definite
  siga <- siga.posdef(siga,am,ctable)

  varcom <- matrix(0,v+1,l*l)
  varcom[1:v, ] <- siga[ , ]
  varcom[v+1, ] <- as.vector(p)
  dimnames(varcom) <- list(c(components,"VarP(I)"),traitpairs)

# kludge a b matrix
  b <- matrix(0,1,l)
  dimnames(b) <- list("dummy",traits)

  retobj <- list(aov=NULL,mdf=NULL,fixform=NULL,b=b,seb=NULL,vara=NULL,totn=NULL,
                 degf=NULL,dme.wmat=NULL,dme.correl=NULL,dme.fit=NULL,dmeopt=NULL,
                 siga=siga, sesiga=NULL,vard=NULL,degfd=NULL,component=components,
                 correlation=NULL,correlation.variance=NULL,correlation.se=NULL,
                 fraction=NULL,fraction.variance=NULL,fraction.se=NULL,
                 variance.components=varcom,variance.components.se=NULL,
                 phenotypic.variance=p,phenotypic.variance.se=NULL,
                 observed.variance=NULL,call=NULL)
  class(retobj) <- "dmm"
  return(retobj)
}
