
# In this file we develop the procedure described in Brockwell, Davis and Yang (2012) 
# for estimation of the underlying noise once the parameters of carma(p,q) are previously 
# obtained
 

yuima.eigen2arparam<-function(lambda){
  MatrixLamb<-diag(lambda)
  
  matrixR<-matrix(NA,length(lambda),length(lambda))
  for(i in c(1:length(lambda))){
    matrixR[,i]<-lambda[i]^c(0:(length(lambda)-1))
  }
  
  AMatrix<-matrixR%*%MatrixLamb%*%solve(matrixR)
  
  acoeff<-AMatrix[length(lambda),]
}



yuima.carma.eigen<-function(A){
  diagA<-eigen(A)
  diagA$values<-diagA$values[order(diagA$values, na.last = TRUE, decreasing = TRUE)]
  n_eigenval<-length(diagA$values)
  diagA$vectors<-matrix(diagA$values[1]^(c(1:n_eigenval)-1),n_eigenval,1)
  if(n_eigenval>=2){
    for (i in 2:n_eigenval){
      diagA$vectors<-cbind(diagA$vectors, 
                           matrix(diagA$values[i]^(c(1:n_eigenval)-1),n_eigenval,1))
    }
  }
  return(diagA)
}


StateVarX<-function(y,tt,X_0,B,discr.eul){
  # The code obtains the first q-1 state variable using eq 5.1 in Brockwell, Davis and Yang 2011
  Time<-length(tt) 
  q<-length(X_0)
  e_q<-rep(0,q)
  e_q[q]<-1
  X_0<-matrix(X_0,q,1)
  e_q<-matrix(e_q,q,1)
  X<-matrix(0,q,Time)
  int<-matrix(0,q,1)
  for (i in c(3:Time)){
    if(discr.eul==FALSE){
      int<-int+(expm(B*(tt[i]-tt[(i-1)]))%*%(e_q*y[i-1]))*(tt[i]-tt[(i-1)])
      X[,i]<-as.matrix(expm(B*tt[i])%*%X_0+int)
    }else{
      X[,i]<-X[,i-1]+(B%*%X[,i-1])*(tt[i]-tt[(i-1)])+(e_q*y[i-1])*(tt[i]-tt[(i-1)])
    }
  }
  return(X)
}


# StateVarXp<-function(y,X_q,tt,B,q,p){
#   # The code computes the  state variable X using the derivatives of eq 5.2 
#   # see Brockwell, Davis and Yang 2011
#   
#   diagMatB<-yuima.carma.eigen(B)
#   if(length(diagMatB$values)>1){
#     MatrD<-diag(diagMatB$values)
#   }else{
#     MatrD<-as.matrix(diagMatB$values)
#   }
#   MatrR<-diagMatB$vectors
#   idx.r<-c(q:(p-1))
#   elem.X <-length(idx.r)
#   YMatr<-matrix(0,q,length(y))
#   YMatr[q,]<-y
#   OutherX<-matrix(NA,elem.X,length(y))
#   # OutherXalt<-matrix(NA,q,length(y))
#   for(i in 1:elem.X){
#     OutherX[i,]<-((MatrR%*%MatrD^(idx.r[i])%*%solve(MatrR))%*%X_q+
#       (MatrR%*%MatrD^(idx.r[i]-1)%*%solve(MatrR))%*%YMatr)[1,]  
#   }
#   X.StatVar<-rbind(X_q,OutherX)
#   return(X.StatVar)
# }


StateVarXp<-function(y,X_q,tt,B,q,p,nume.Der){
  # The code computes the  state variable X using the derivatives of eq 5.2 
  # see Brockwell, Davis and Yang 2011
  if(nume.Der==FALSE){
    yuima.warn("We need to develop this part in the future for gaining speed")
    return(NULL)
#     diagMatB<-yuima.carma.eigen(B)
#     if(length(diagMatB$values)>1){
#       MatrD<-diag(diagMatB$values)
#     }else{
#       MatrD<-as.matrix(diagMatB$values)
#     }
#     MatrR<-diagMatB$vectors
#     idx.r<-c(q:(p-1))
#     elem.X <-length(idx.r)
#     YMatr<-matrix(0,q,length(y))
#     YMatr[q,]<-y
#     OutherX<-matrix(NA,elem.X,length(y))
#     # OutherXalt<-matrix(NA,q,length(y))
#     for(i in 1:elem.X){
#       OutherX[i,]<-((MatrR%*%MatrD^(idx.r[i])%*%solve(MatrR))%*%X_q+
#                       (MatrR%*%MatrD^(idx.r[i]-1)%*%solve(MatrR))%*%YMatr)[1,]  
#     }
#     X.StatVar<-rbind(X_q,OutherX)
#     return(X.StatVar)
  }else{
    X_q0<-as.numeric(X_q[1,])
    qlen<-length(X_q0)
    X.StatVar<-matrix(0,p,qlen)
    X.StatVar[1,]<-X_q0[1:length(X_q0)]
    diffStatVar<-diff(X_q0)
    difftime<-diff(tt)[1]
    for(i in 2:p){
      in.count<-(p-i+1)
      fin.count<-length(diffStatVar)
      dummyDer<-(diffStatVar/difftime)
      X.StatVar[i,c(1:length(dummyDer))]<-(dummyDer)
      diffStatVar<-diff(dummyDer)
#       if(in.count-1>0){
#         difftime<-difftime[c((in.count-1):fin.count)]
#       }else{difftime<-difftime[c((in.count):fin.count)]}
    }
    X.StatVar[c(1:dim(X_q)[1]),]<-X_q
    return(X.StatVar[,c(1:length(dummyDer))])
  }
}


bEvalPoly<-function(b,lambdax){
  result<-sum(b*lambdax^(c(1:length(b))-1))
  return(result)
}

aEvalPoly<-function(a,lambdax){
  p<-length(a)
  a.new<-c(1,a[c(1:(p-1))])
  pa.new<-c(p:1)*a.new
  result<-sum(pa.new*lambdax^(c(p:1)-1))
  return(result)
}


CarmaNoise<-function(yuima, param, data=NULL,NoNeg.Noise=FALSE){
  if( missing(param) ) 
    yuima.stop("Parameter values are missing.")
  
  if(!is.list(param))
    yuima.stop("Argument 'param' must be of list type.")
  
  vect.param<-as.numeric(param)
  name.param<-names(param)
  names(vect.param)<-name.param
  
  if(is(yuima,"yuima")){
    model<-yuima@model
    if(is.null(data)){
      observ<-yuima@data
    }else{observ<-data}
}else{
  if(is(yuima,"yuima.carma")){
    model<-yuima
    if(is.null(data)){
      yuima.stop("Missing data")
    }
    observ<-data
  }
}

  if(!is(observ,"yuima.data")){
   yuima.stop("Data must be an object of class yuima.data-class")  
  }
  
  info<-model@info
  
  numb.ar<-info@p
  name.ar<-paste(info@ar.par,c(numb.ar:1),sep="")
  ar.par<-vect.param[name.ar]
  
  numb.ma<-info@q
  name.ma<-paste(info@ma.par,c(0:numb.ma),sep="")
  ma.par<-vect.param[name.ma]
  
  loc.par=NULL
  if (length(info@loc.par)!=0){
    loc.par<-vect.param[info@loc.par]
  }
  
  scale.par=NULL
  if (length(info@scale.par)!=0){
    scale.par<-vect.param[info@scale.par]
  }
  
  lin.par=NULL
  if (length(info@lin.par)!=0){
    lin.par<-vect.param[info@lin.par]
  }
  
  
  
  ttt<-observ@zoo.data[[1]]
  tt<-index(ttt)
  y<-coredata(ttt)
  
  levy<-yuima.CarmaNoise(y,tt,ar.par,ma.par, loc.par, scale.par, lin.par,NoNeg.Noise)
  inc.levy<-diff(as.numeric(levy))
  return(inc.levy[-1]) #We start to compute the increments from the second observation. 
}



yuima.CarmaNoise<-function(y,tt,ar.par,ma.par, 
                                loc.par=NULL, 
                                scale.par=NULL, 
                                lin.par=NULL, 
                                NoNeg.Noise=FALSE){
    
  if(!is.null(loc.par)){
      y<-y-loc.par
#       yuima.warn("the loc.par will be implemented as soon as possible")
#       return(NULL)
    } 
    if(NoNeg.Noise==TRUE){
      if(length(ar.par)==length(ma.par)){
        yuima.warn("The case with no negative jump needs to test deeply")
        #mean.y<-tail(ma.par,n=1)/ar.par[1]*scale.par
        mean.y<-mean(y)
        #y<-y-mean.y
        mean.L1<-mean.y/tail(ma.par,n=1)*ar.par[1]/scale.par
      }
    }  
    if(!is.null(scale.par)){
      ma.parAux<-ma.par*scale.par
      names(ma.parAux)<-names(ma.par)
      ma.par<-ma.parAux
#       yuima.warn("the scale.par will be implemented as soon as possible")
#       return(NULL)
    }
    if(!is.null(lin.par)){
      yuima.warn("the lin.par will be implemented as soon as possible")
      return(NULL)
    }
    
    p<-length(ar.par)
    q<-length(ma.par)
    
    A<-MatrixA(ar.par[c(p:1)])
    b_q<-tail(ma.par,n=1)
    ynew<-y/b_q
    
    if(q==1){
      yuima.warn("The Derivatives for recovering noise have been performed numerically.")
      nume.Der<-TRUE
      X_qtot<-ynew
      X.StVa<-matrix(0,p,length(ynew))  
      X_q<-X_qtot[c(1:length(X_qtot))]
      
      X.StVa[1,]<-X_q
      if(p>1){
        diffY<-diff(ynew)
        diffX<-diff(tt)[1]
        for (i in c(2:p)){          
          dummyDerY<-diffY/diffX
          X.StVa[i,c(1:length(dummyDerY))]<-(dummyDerY)
          diffY<-diff(dummyDerY)
        
        }
        X.StVa<-X.StVa[,c(1:length(dummyDerY))]
      }
      #yuima.warn("the car(p) process will be implemented as soon as possible")
      #return(NULL)
    }else{
    newma.par<-ma.par/b_q
    # We build the matrix B which is necessary for building eq 5.2
    B<-MatrixA(newma.par[c(1:(q-1))])
    diagB<-yuima.carma.eigen(B)
    
    e_q<-rep(0,(q-1))
    if((q-1)>0){
      e_q[(q-1)]<-1
      X_0<-rep(0,(q-1))
    }else{
      e_q<-1
      X_0<-0
    }
    
    
    discr.eul<-TRUE 
    # We use the Euler discretization of eq 5.1 in Brockwell, Davis and Yang
    X_q<-StateVarX(ynew,tt,X_0,B,discr.eul)
  
    
    nume.Der<-TRUE
    #plot(t(X_q))
    X.StVa<-StateVarXp(ynew,X_q,tt,B,q-1,p,nume.Der) #Checked once the numerical derivatives have been used    
    #plot(y)
    }
    
    diagA<-yuima.carma.eigen(A)

  
       
    
    
    BinLambda<-rep(NA,length(ma.par))
    for(i in c(1:length(diagA$values))){
      BinLambda[i]<-bEvalPoly(ma.par,diagA$values[i])
    }
    MatrBLam<-diag(BinLambda)
    
    # We get the Canonical Vector Space in eq 2.17
    Y_CVS<-MatrBLam%*%solve(diagA$vectors)%*%X.StVa #Canonical Vector Space
  
    # We verify the prop 2 in the paper "Estimation for Non-Negative Levy Driven CARMA process
    # yver1<-Y_CVS[1,]+Y_CVS[2,]+Y_CVS[3,]
    
    # plot(yver1)
    
    # plot(y)
    # Prop 2 Verified even in the case of q=0
    
    idx.r<-match(0,Im(diagA$values))
    lambda.r<-Re(diagA$values[idx.r])

    if(is.na(idx.r)){
      yuima.warn("all eigenvalues are immaginary numbers")
      idx.r<-1
      lambda.r<-diagA$values[idx.r]
    }
        int<-0
    
    derA<-aEvalPoly(ar.par[c(p:1)],lambda.r)
#     if(q==1){
#       tt<-tt[p:length(tt)]
#       # CHECK HERE 15/01 
#     }
    if(nume.Der==TRUE){
      tt<-tt[p:length(tt)]
      # CHECK HERE 15/01 
    }
    
    lev.und<-matrix(0,1,length(tt))
      
    Alternative<-FALSE
    if(Alternative==TRUE){
      incr.vect<-(X.StVa[,2:dim(X.StVa)[2]]-X.StVa[,1:(dim(X.StVa)[2]-1)])-(A%*%X.StVa[,1:(dim(X.StVa)[2]-1)]
                #                                                            +(matrix(c(0,mean.L1),2,(dim(X.StVa)[2]-1)))/diff(tt)[1]
                                                                            )*diff(tt)[1]
      #+(matrix(c(0,mean.L1),2,(dim(X.StVa)[2]-1)))
      lev.und<-as.matrix(cumsum(as.numeric(incr.vect[dim(X.StVa)[1],])),1,length(tt))
    }else{
      for(t in c(2:length(tt))){
        int<-int+Y_CVS[idx.r,t-1]*(tt[t]-tt[t-1])
        lev.und[,t]<-derA/BinLambda[idx.r]*(Y_CVS[idx.r,t]-Y_CVS[idx.r,1]-lambda.r*int)
      }
    }
#     if(NoNeg.Noise==TRUE){
#       if(length(ar.par)==length(ma.par)){
#         if(Alternative==TRUE){
#           mean.Ltt<-mean.L1*tt[c(2:length(tt))]
#         }else{mean.Ltt<-mean.L1*tt}
#         lev.und<-lev.und+mean.Ltt
#       
#       }
#     }
  return(Re(lev.und))
}


