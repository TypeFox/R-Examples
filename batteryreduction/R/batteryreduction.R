#' A function for data reduction
#'
#' Battery reduction is a method used in data reduction. It uses Gram-Schmidt orthogonal rotations to find out a subset of variables best representing the original set of variables.
#' @param vars a vector of variable names
#' @param numfact	numeric, number of variables to be selected
#' @param data a dataframe
#' @return If numfact=1, reda matrix and selected variable are returned. If numfact>=2, reda matrix, rota matrix, and selected variables are returned.
#' @references D'Agostino, R. B., & Zhang, Z. (1992). Data/variable reduction by principal components, battery reduction and variable clustering. MATRIX, 7(60), 06.
#' @references Hans Werner Borchers (2015). pracma: Practical Numerical Math Functions. R package version 1.8.6. http://CRAN.R-project.org/package=pracma
#' @examples ## Generate an example dataset
#' set.seed(1234)
#' data<-data.frame(x1=rnorm(n=100, mean=14, sd=7),
#' x2=rnorm(n=100, mean=3, sd=1),
#' x3=rpois(n=100, lambda=1),
#' x4=rpois(n=100, lambda=10),
#' x5=rgamma(n=100, shape=1),
#' x6=rgamma(n=100, shape=10))

#' ## Demonstrate batteryreduction
#' vars<-c('x1','x2','x3','x4','x5','x6')
#' numfact<-3
#' batteryreduction(vars, numfact, data)

batteryreduction<-function(vars, numfact, data){
  #1. Compute the correlation matrix
  corr.matrix<-cor(data[, vars], use = "everything",method = "pearson")
  
  #2. Create the loading matrix
  m<-eigen(corr.matrix)$values
  d<-diag(m)
  e<-eigen(corr.matrix)$vectors #Seems e estimate is different from SAS
  a<-e%*%sqrt(d)
  nfact<-numfact
  nvars<-nrow(a)
  reda<-as.matrix(a[, c(1:nfact)])
  rownames(reda)<-vars
  af<-matrix(nrow=nvars, ncol=1)
  af[,1]<-0
  
  #3. Performs the gram-schmidt rotation
  if (nfact>=2){
    for (i in 1:(nfact-1)){
      if (i==1){
        a1<-a[, c(i:nfact)]
      }else{
        colno<-ncol(a)
        a1<-a[, c(2:colno)]
      }
      
      sum.of.squares<-function(x) sum(x*x)
      sq<-apply(a1, 1, sum.of.squares)
      max.index<-which.max(sq)
      norm1<-a1[max.index,]
      sq1<-sq[[max.index]]
      sq1<-sqrt(sq1)
      normn1<-norm1/sq1
      nfact1<-nfact-i
      ident=diag(nfact1)
      values<-rep(0,times=nfact1)
      zero<-matrix(values,nrow=1,ncol=nfact1)
      tempy<-rbind(ident,zero)
      #normn1<-t(normn1)
      ty<-cbind(normn1,tempy)
      
      ayp1<-a1 %*% gramSchmidt(ty)$Q
      a<-ayp1
      colno<-ncol(a)
      
      if (i-(nfact-1)<0) {
        af<-cbind(af, a[, c(1)])
      }else {
        af<-cbind(af, a[, c(1:2)])
      }
      
    }
    
    rota<-af[,c(2:(nfact+1))]
    rownames(rota)<-vars
    per_var<-apply(rota, 2, sum.of.squares)
    per_var<-per_var/nvars
    
    selected<-rep(NA,nfact)
    for (i in 1:nfact){
      selected[i]<-rownames(rota)[which.max(abs(rota[,i]))] 
    }
    
    return(list(reda = reda, rota = rota, per_var=per_var, selected=selected))
  }else{
    selected<-vars[which.max(abs(reda))] 
    return(list(reda = reda, selected=selected))
  }
}