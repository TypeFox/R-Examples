

## Function to calculate the MLE of the hybrid design for multiple counties
## Input: MM (matrix), NN (matrix), cc (list), aprximation
## Output: MLE
## Notes: cc should be "NULL" if no case-control data observed
hybdesEco = function(MM.Z, MM.W, NN, cc, aprx='binom', start.mle=NA, group.int=FALSE){

  betafct = function(x){ return(c(x[1],x[2],x[3],x[2]+x[3],x[4]))}
  analytics = TRUE
  analytic.hessian = FALSE
  nlm.print.level = 0
  nlm.iterlim = 100
  nlm.hessian=TRUE

  if(!is.matrix(MM.Z)){stop("MM.Z must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.matrix(MM.W)){stop("MM.W must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.matrix(NN)){stop("NN must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.list(cc)){stop("cc must be a list, with elements equal to case-control data")}

  K = dim(NN)[1]
  J = dim(MM.Z)[2]*dim(MM.W)[2]-1

  if( (dim(MM.Z)[2]!=2) | (dim(MM.W)[2]!=2) ){stop("MM.Z and MM.W must be K x 2 matrices.")}
  
  nn = t(sapply(cc,function(x){if(!is.null(x)){ return(apply(x,2,sum))}; if(is.null(x)){ return(c(0,0))} } ) )
  mm.z = t(sapply(cc,function(x){if(!is.null(x)){ mm = apply(x,1,sum); mm.z = c(mm[1]+mm[3],mm[2]+mm[4]); return(mm.z)}; if(is.null(x)){ return(rep(0,dim(MM.Z)[2]))} } ) )
  mm.w = t(sapply(cc,function(x){if(!is.null(x)){ mm = apply(x,1,sum); mm.w = c(mm[1]+mm[2],mm[3]+mm[4]); return(mm.w)};  if(is.null(x)){ return(rep(0,dim(MM.W)[2])) }} ) )

  analytics.txt = paste("      attr(lik, 'gradient') = -f$gradient
      attr(lik, 'hessian') = -f$hessian",sep='')

  if(analytic.hessian == FALSE){
    analytics.txt = paste("      attr(lik, 'gradient') = -f$gradient",sep='')
  }

  if(analytics == FALSE){analytics.txt = ""}

   dummy = function(x){return(x)} # Fake function to appease CRAN warnings/notes; the real 'dummy' funciton is defined below.

  nsums = 0
  for(k in 1:K){
    nsums = nsums + enumerate.count(MM.Z[k,]-mm.z[k,],MM.W[k,]-mm.w[k,])
  }

  eval.text = paste("dummy = function(beta){ 
    if( group.int==FALSE ){
      beta.new = betafct(beta)
      beta.matrix = matrix(rep(beta.new, K), byrow=TRUE, nrow=K)
    }
    if( group.int==TRUE ){
      beta.int = beta[1:K]
      beta.or = beta[(K+1):length(beta)]
      beta.matrix = betafct(c(beta.int[1], beta.or))
      if(K > 1){
        for(k in 2:K){
          beta.matrix = rbind(beta.matrix, betafct(c(beta.int[k], beta.or)) )
      }}
    }
    f = hyblik.eco(beta.matrix, MM.Z, MM.W, NN, cc, Rk = NA, aprx ='",aprx,"', group.int=group.int)
    lik = -f$hyblik
    ", analytics.txt, "
    return(lik) }",sep='')
 
  eval(parse(text = eval.text))

  # Determine start value for MLE calculation

  mle = nlm(dummy, start.mle, hessian=nlm.hessian, print.level = nlm.print.level, iterlim = nlm.iterlim)
#  mle = nlm(dummy, start.mle, hessian=FALSE, print.level = nlm.print.level, iterlim = nlm.iterlim)
  
  return(list(mle=mle, start.mle = start.mle, nsums = nsums))

}

