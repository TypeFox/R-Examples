 
## Function to calculate the hybrid likelihood for multiple counties
## Input: beta, MM (matrix), NN (matrix), cc (list), Rk, aprximation
## Output: Likelihood
## Notes: cc should be "NULL" if no case-control data observed

hyblik = function(beta.matrix, MM, NN, cc, aprx = 'binom', ntrue = 0, group.int=FALSE){
  if(!is.matrix(MM)){stop("MM must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.matrix(NN)){stop("NN must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.list(cc)){stop("cc must be a list, with elements equal to case-control data")}
  if(!(aprx %in% c('binom','norm','pois',NA))){stop("aprx must be either 'binom', 'norm', 'pois', or NA")}

  K = dim(MM)[1]

  if(!is.na(aprx)){aprx.text = paste("'",aprx,"'",sep='') }
  else{
    aprx.text = NA
	ntrue = K
  }
  aprx.text = rep(aprx.text, K)

  ntrue = min(ntrue, K)

  mm = t(sapply(cc, function(x){
       if(length(x)==0){
         return(c(0,0))
       } 
         else{return(apply(x,1,sum))} 
       }))
  nn = t(sapply(cc, function(x){
       if(length(x)==0){
         return(c(0,0))
       } 
         else{return(apply(x,2,sum))} 
       }))

  MM.star = MM - mm
  NN.star = NN - nn

  NN.ord = order(NN.star[,2], apply(NN,1,sum))

  if(ntrue > 0){ 
    aprx.text[NN.ord[1:ntrue]] = NA
  }

#  print(aprx.text)

  hyblik = 0

  if(length(dim(beta.matrix))==0){
    if(length(beta.matrix)!=dim(MM)[2]){stop("number of parameters in beta vector not correct")}
    beta.matrix = matrix(rep(beta.matrix),K, byrow=TRUE, nrow=K)
  }
  if(sum(dim(beta.matrix)==dim(MM))!=2){stop("dimensions of beta.matrix not correct")}

  for(k in 1:K){
    beta = beta.matrix[k,]
    if(is.null(cc[[k]])){
      eval.text = paste('log.ecportion.beta(beta,MM[k,],NN[k,],aprx=',aprx.text[k],')',sep='')
    }
    else{
      eval.text = paste('log.altdecomp.beta(beta,MM[k,],NN[k,],cc[[k]],aprx=',aprx.text[k],')',sep='')
    }
    contributions = eval(parse(text = eval.text))
    hyblik = hyblik + contributions$hyblik
  }
  return(list("hyblik" = hyblik))
}


 



## Function to calculate the hybrid likelihood for multiple counties
## Input: beta, MM (matrix), NN (matrix), cc (list), Rk, aprximation
## Output: Likelihood
## Notes: cc should be "NULL" if no case-control data observed

hyblik.eco = function(beta.matrix, MM.Z, MM.W, NN, cc, Rk = NA, aprx = NA, group.int=FALSE){

#  if(!is.matrix(MM)){stop("MM must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.matrix(NN)){stop("NN must be a matrix (rows = groups, columns = margin totals)")}
  if(!is.list(cc)){stop("cc must be a list, with elements equal to case-control data")}
  if(!is.list(Rk)){if(!is.na(Rk)){stop("Rk must be a list")}}

  if(is.list(Rk)){Rk.text = "Rk[[k]]"}
  else{Rk.text = NA}

  if(!is.na(aprx)){aprx.text = paste("'",aprx,"'",sep='') }
  else{aprx.text = NA}

  K = dim(NN)[1]

  hyblik = rep(0, K)
  if(group.int==FALSE){
    lenbeta = dim(beta.matrix)[2] - 1
  }
  else{
#   No-intercept model:
    lenbeta = dim(beta.matrix)[2] + dim(beta.matrix)[1] - 2
  }
  hybgrad = rep(0, lenbeta )
  hybhess = matrix(0, nrow = lenbeta, ncol = lenbeta )   


  K = dim(NN)[1]
  J = dim(cc[[1]])[1] - 1

  nn = t(sapply(cc,function(x){if(!is.null(x)){ return(apply(x,2,sum))}; if(is.null(x)){ return(c(0,0))} } ) )
  mm.z = t(sapply(cc,function(x){if(!is.null(x)){ mm = apply(x,1,sum); mm.z = c(mm[1]+mm[3],mm[2]+mm[4]); return(mm.z)}; if(is.null(x)){ return(rep(0,dim(MM.Z)[2]))} } ) )
  mm.w = t(sapply(cc,function(x){if(!is.null(x)){ mm = apply(x,1,sum); mm.w = c(mm[1]+mm[2],mm[3]+mm[4]); return(mm.w)};  if(is.null(x)){ return(rep(0,dim(MM.W)[2])) }} ) )

  MM.Z.star = MM.Z-mm.z
  MM.W.star = MM.W-mm.w

#  if(length(dim(beta.matrix))==0){
#    if(length(beta.matrix)!=dim(MM)[2]){stop("number of parameters in beta vector not correct")}
#    beta.matrix = matrix(rep(beta.matrix),K, byrow=TRUE, nrow=K)
#  }
#  if(sum(dim(beta.matrix)==dim(MM))!=2){stop("dimensions of beta.matrix not correct")}

  length.beta = dim(beta.matrix)[2]
  for(k in 1:K){
    or.est = exp(beta.matrix[k,length.beta])
    beta = beta.matrix[k,1:(length.beta-1)]

    margins <- list(M0=MM.Z.star[k,1], M1=MM.Z.star[k,2], N0=MM.W.star[k,1], N1=MM.W.star[k,2])
    possible.semieco <- enumerate(MM.Z.star[k,], MM.W.star[k,])

    if(is.null(cc[[k]])){
      eval.text = paste('log.ecportionECO.beta(beta,MM,NN[k,],nn.var=',Rk.text,',aprx=',aprx.text,')',sep='')
    }
    else{
      eval.text = paste('log.altdecompECO.beta(beta,MM,NN[k,],cc[[k]],nn.var=',Rk.text,',aprx=',aprx.text,')',sep='')
    }

    nadd = dim(possible.semieco)[2]
    lcontrib = 
    pcontrib = rep(NA, nadd)

#    gcontrib = rep(NA, lenbeta)
#    hcontrib = matrix(0, nrow=lenbeta, ncol=lenbeta)

    gcontrib = matrix(NA, nrow=nadd, ncol=lenbeta-1)
    hcontrib = vector("list", nadd)

    if(is.null(cc[[k]])){
      MM = rbind(MM.Z.star[k,] - possible.semieco, possible.semieco)
    }
    else{
      MM = rbind(MM.Z.star[k,] - possible.semieco, possible.semieco) + apply(cc[[k]],1,sum)
    }

    contributions = eval(parse(text = eval.text))
    lcontrib = contributions$hyblik
    gcontrib = t(contributions$hybgrad)
    hcontrib = t(contributions$hybhess)

    for(semieco in 1:nadd){
      pcontrib[semieco] = dXhyper( or.est, margins, N1x = possible.semieco[,semieco] )
    }

    gnu = modeXhyper(or.est, margins)
    dgnu = dXhyper( or.est, margins, N1x = c(MM.W.star[k,2]-gnu,gnu))
    ddgnu = ddXhyper( or.est, margins, N1x = c(MM.W.star[k,2]-gnu,gnu))
    dddgnu = dddXhyper( or.est, margins, N1x = c(MM.W.star[k,2]-gnu,gnu))

    gnuw = dgnu/ddgnu*gnu
    if(gnu == 0){gnuw = 0}

    gnu1 = dgnu/ddgnu*gnu/or.est
    if(gnu == 0){gnu1 = 0}
    gnu2 = dgnu/dddgnu*gnu*(gnu-1)/or.est^2
    if(gnu %in% c(0,1)){gnu2 = 0}


    dpsivals = pcontrib*(possible.semieco[2,] - gnuw)
#    dpsivals = pcontrib*(possible.semieco[2,] - gnu1*or.est)

    ddpsivals = pcontrib*or.est*(possible.semieco[2,]*(possible.semieco[2,]-1)/or.est - 
                2*possible.semieco[2,]*gnu1 +
                possible.semieco[2,]/or.est - 
                or.est*gnu2 -
                gnu1 +
                2*or.est*gnu1^2)


#    ddpsivals = dpsivals + dpsivals^2/pcontrib + 
#                 or.est*pcontrib*(possible.semieco[2,]/or.est - gnu2 + gnu1^2)



    xMM = log(pcontrib) + lcontrib
    hyblik[k] = xMM[which.max(xMM)[1]] + log(sum(exp(xMM - xMM[which.max(xMM)[1]])))

    gvals = apply(gcontrib * (exp(lcontrib) * pcontrib), 2, sum)/exp(hyblik[k])
    dpsi = sum(dpsivals*exp(lcontrib))/exp(hyblik[k])

    ddpsi = sum(ddpsivals*exp(lcontrib))/exp(hyblik[k])

#    llike = lcontrib[which.max(lcontrib)[1]] + log(sum(exp(lcontrib - lcontrib[which.max(lcontrib)[1]])))
#    gvals = apply(exp(lcontrib)*gcontrib,2,sum)/exp(llike)
    if(group.int==FALSE){
      hybgrad = hybgrad + c(gvals,dpsi)
#      hybhess = hybhess + contributions$hybhess
    }
    else{
      hybgrad[c(k,(K+1):(lenbeta))] = hybgrad[c(k,(K+1):(lenbeta))] + c(gvals,dpsi)

      hybhess[c(k,(K+1):(lenbeta-1)),c(k,(K+1):(lenbeta-1))] = 
           hybhess[c(k,(K+1):(lenbeta-1)),c(k,(K+1):(lenbeta-1))] - 
           matrix(gvals,ncol=1)%*%matrix(gvals,nrow=1) +
           matrix(apply(pcontrib*exp(hyblik[k])*( t(apply(gcontrib, 1, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)})) + hcontrib),2,sum), nrow=3)/exp(hyblik[k])

      hybhess[c(k,(K+1):(lenbeta-1)),lenbeta] =
      hybhess[lenbeta, c(k,(K+1):(lenbeta-1))] = 
        hybhess[c(k,(K+1):(lenbeta-1)),lenbeta] - gvals*dpsi + 
        apply(gcontrib * (exp(lcontrib) * dpsivals), 2, sum)/exp(hyblik[k])

      hybhess[lenbeta,lenbeta] = hybhess[lenbeta,lenbeta] - dpsi^2 + ddpsi
    }

  }

  HYBLIK = sum(hyblik)

#  hybgrad = hybgrad[-lenbeta]
#  hybhess = hybhess[-lenbeta,][,-lenbeta]

  return(list("hyblik" = HYBLIK, "gradient" = hybgrad, "hessian" = hybhess))
#  return(list("hyblik" = HYBLIK, "gradient" = hybgrad))
}




