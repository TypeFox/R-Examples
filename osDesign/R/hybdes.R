## Function to calculate the MLE of the hybrid design for multiple counties
## Input: MM (matrix), NN (matrix), cc (list), aprximation
## Output: MLE
## Notes: cc should be "NULL" if no case-control data observed

hybdes = function(MM, NN, cc, ntrue = 0, aprx='binom', start.mle=NA, group.int=FALSE, betafct = function(x){return(x[1] + c(0,x[-1]) )}, print.level = 0, iterlim = 100){

  if(!is.matrix(MM)){stop("MM must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.matrix(NN)){stop("NN must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.list(cc)){stop("cc must be a list, with elements equal to case-control data")}

  K = dim(MM)[1]
  J = dim(MM)[2] - 1

  nn = t(sapply(cc,function(x){if(!is.null(x)){ return(apply(x,2,sum))}; if(is.null(x)){ return(c(0,0))} } ) )
  mm = t(sapply(cc,function(x){if(!is.null(x)){ return(apply(x,1,sum))}; if(is.null(x)){ return(rep(0,J+1))} } ) )

#  analytics.txt = paste("      attr(lik, 'gradient') = -f$gradient
#      attr(lik, 'hessian') = -f$hessian",sep='')

#  if(analytic.hessian == FALSE){
#    analytics.txt = paste("      attr(lik, 'gradient') = -f$gradient",sep='')
#  }

#  if(analytics == FALSE){analytics.txt = ""}
   analytics.txt = ""

   dummy = function(x){return(x)} # Fake function to appease CRAN warnings/notes; the real 'dummy' funciton is defined below.
    
  if(is.na(aprx)){
    ntrue = K
    nsums = 0
    for(k in 1:K){
      nsums = nsums + enumerate.count(MM[k,]-mm[k,],NN[k,]-nn[k,])
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

      f = hyblik(beta.matrix, MM, NN, cc, aprx = NA, ntrue = ntrue, group.int=group.int)
      lik = -f$hyblik
      ", analytics.txt, "
      return(lik) }",sep='')
  }
  else{
    nsums = NA 
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
      f = hyblik(beta.matrix, MM, NN, cc, aprx ='",aprx,"', ntrue = ntrue, group.int=group.int)
      lik = -f$hyblik
      ", analytics.txt, "
      return(lik) }",sep='')
  } 

  eval(parse(text = eval.text))

  if(is.na(start.mle)[1]){ 
    if(identical(betafct(1:(J+1)),1 + c(0,2:(J+1))) ){
	  nn = t(sapply(cc, function(x){
          if(length(x)==0){
           return(c(0,0))
           } 
           else{return(apply(x,2,sum))} 
           }))
	  nk = apply(nn, 1, sum)
	  sampled.groups = (nk > 0)
      K.sampled = sum(sampled.groups)
      cc.data = unlist(cc[sampled.groups])       
      sampled.data = data.frame(group = rep(1:sum(sampled.groups), nk[sampled.groups]), Y = rep(rep(rep(0:1, each = (J+1)),K.sampled), cc.data), X = rep(rep(0:J,2*K.sampled), cc.data))
	  if(group.int==TRUE){
         tps.fit = tps(Y ~ -1 + as.factor(group) + as.factor(X), group = sampled.data$group, nn0 = NN[sampled.groups,1], nn1 = NN[sampled.groups,2], cohort=TRUE, method="PL", data = sampled.data)
         intpos = grep('as.factor(group)',names(tps.fit$coef), fixed = TRUE)
             
         start.mle = rep(NA, J + K)     
         start.mle[(1:K)[sampled.groups]] = tps.fit$coef[intpos]
         start.mle[(1:K)[!sampled.groups]] = mean(tps.fit$coef[intpos])
         start.mle[(K+1):(K+J)] = tps.fit$coef[-intpos]
      }
	  else{
	     tps.fit = tps(Y ~ as.factor(X), group = sampled.data$group, nn0 = NN[sampled.groups,1], nn1 = NN[sampled.groups,2], cohort=TRUE, method="PL", data = sampled.data)
		 start.mle = tps.fit$coef
	  }
	}
	else{
  		len = length(betafct(rep(1,J+1)))
		if(group.int==TRUE){
		  start.mle = rep(0, len + K-1)
		}
		else{
		  start.mle = rep(0, len)
		}
	}
  }
  
  mle = nlm(dummy, start.mle, hessian=TRUE, print.level = print.level, iterlim = iterlim)
  
  return(list(mle=mle, start.mle = start.mle, nsums = nsums))

}



#class(hybdes) <- "hybdes"



















