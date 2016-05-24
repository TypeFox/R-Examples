
### Cook all the intertermediate results that will be used in various functions
ingredients<-function(object){

  

   
  Z<-object$x 
  Y<-object$y
  t.fail<-Y[,1]
  censor<-Y[,2]
  if(!is.matrix(Z))Z=as.matrix(Z)
  
  
### for  code testing purposes
   # failure.time<-test.data[,1]  
   #censor<-test.data[,2]
   #Z<-test.data[,c(3,4)]
      # Force the covariate input to be a matrix
    # Z<-as.matrix(Z)  
 
   ### number of parameters , number of observtaions  
   n.par<-object$npar
   n.obs<-object$nobs
   wts<-object$weights

   theta<-object$coef
   effects<-as.vector(Z%*%theta)
      ### Centralized the covariates
    # Z <- scale( Z,center= rep(1,n.par),scale=F)
   
   ### Devide the timeline to small intervals according to 
   ### the occurence of events or censoring  
 
  ordered<-sort(t.fail)
  t.unique<-unique(ordered)
  eta.ncol<-length(t.unique)
  t.diff<-t.unique-c(0,t.unique[-eta.ncol])

  # Time points where failures occur
  idx.c<-which(censor==1)
  t1.sorted<-sort(t.fail[idx.c])
  t1.unique<-unique(t1.sorted)

  # mark all the failure time on the timeline
  match.eta<-match(t.fail, t.unique)
  
  # mark all the event time on the timeline
  match.event<- match(t1.unique, t.unique)


  z.death <- eta.num <- eta.den <-matrix(NA, nrow=n.par, ncol=eta.ncol)
  n.death<-NULL
  eta.den.vec<-NULL

  for ( i in 1:eta.ncol){
    idx1<-(t.fail >= t.unique[i])
    idx2<-which(t.fail == t.unique[i])

    n.death[i]<-sum((censor*wts)[idx2])
    eta.den.vec[i]<-sum(wts[idx1])
    eta.den[,i]<-eta.den.vec[i]

    for (j in 1:n.par){
      eta.num[j,i] <-sum((Z*wts)[idx1,j])
      z.death[j,i]<-sum((Z[,j]*censor*wts)[idx2])
    }
  }

  eta<-eta.num/eta.den

  # because the eta for the first interval is 0, thus the eta.cum starts from 0 at t_1 and then eta*(t_2-t_1)... 
  eta.new <- cbind(rep(0, n.par), eta[, -eta.ncol]) 
  eta.cum <- as.matrix(apply(t(eta.new)*t.diff,2,cumsum)) # Note  eta.cum  is eta.ncol x n.par while eta is n.par x eta.ncol
  dLambda1 <- n.death/eta.den.vec
  lambda0 <- dLambda1-as.vector(t(eta.new)%*%theta*t.diff)
  # See equation (3.8) in Jie Hu's thesis. This is the results for the points at failure times in the original data
  # we do not need new data to caculated this quantity and its variance 
  Lambda0 <- cumsum(lambda0)






   return(list(Z=Z, t.fail=t.fail, censor=censor,  n.obs=n.obs, wts =wts, #observations
               t.unique=t.unique, t.diff=t.diff, t1.unique=t1.unique,  #time related
               n.par=n.par,theta=theta, effects=effects, #coef related
              eta.ncol=eta.ncol, match.eta=match.eta, eta=eta, eta.den.vec=eta.den.vec,  eta.cum=eta.cum,    #weighted predictors 
              idx.c=idx.c, match.event=match.event,  n.death=n.death, z.death=z.death, #events related results
              dLambda1=dLambda1, lambda0=lambda0, Lambda0=Lambda0 #baseline hazards related results
           ))
      invisible()
 }




### estimate the model-based variance for predicted values 
L.nvar<-function(t.pos= t.pos, pred = pred, newtime=newtime, n.par=n.par,n.death=n.death,z.death=z.death, t.unique =t.unique,
        eta.den.vec=eta.den.vec, eta=eta, dLambda1=dLambda1, eta.cum=eta.cum,  iA=iA, var = var, ...){  
        

  # See equation (3.42) on p.97 of Jie Hu's thesis  for the variance formula 
  # step 0:  find the location of s in the original timeline 

   t.mark <- t.unique[t.pos] # The closed time point in the original data timeline

  # calculate s - t.mark, the last intervla 
   last.interval <- newtime - t.mark  

  #step 1 
  #calculate \int_0^t.mark PdN(t)/(PY(t)^2)
  
  ## L.var1 = P\int_0^{s}\frac{1}{P^2Y(t)}Y(t)\left\{d\Lambda_0(t)+Z^T\theta_0dt\right\}.
  ## \{d\Lambda_0(t)+Z^T\theta_0dt is replaced by dN(t)
  L.var1<-cumsum((1/(eta.den.vec^2)*n.death)[t.pos])

  L.cov1 <- L.cov2 <- rep(0,n.par)
  L.var <- L.var2 <- L.cov <- NULL
      
  for ( i in 1:length(t.pos)){

    idx <- t.pos[i]

    #step 2
    #calculate \int_0^t.mark PZdN(t)/PY(t)
    #calculate \int_0^t.mark PdN(t)/PY(t) eta(t)

    #\{d\Lambda_0(t)+Z^T\theta_0dt is replaced by dN(t)
    # note a different approach is to remove the component involving $dLambda_0(t)$, since the two terms cancel out. 
    # z.death[j,i]<-sum((Z[,j]*censor*wts)[idx2])   
    # eta.den.vec[i]<-sum(wts[idx1])
    # the increment is = PZdN(t)/PY(t)
      
    L.cov1 <- L.cov1+z.death[,idx]/eta.den.vec[idx]

    # L.cov2 = P\int_0^{\tau}\frac{PY(t)Z}{PY(t)}\frac{Y(t)}{PY(t)}\{d\Lambda_0(t)+Z^T\theta_0dt\}
    # dLambda1<-n.death/eta.den.vec
    # eta = sum((Z*wts)[idx1,j])/sum(wts[idx1])
    # the increment is PdN(t)}/PY(t)}*\frac{PZY(t)}{PY(t)}
    L.cov2 <- L.cov2+ dLambda1[idx]*eta[,idx]
      
    #step 3 calculate the variance and covariance

    # calcuate D(s) where D(s) = D(t.mark[i]) + last.interval[s]*eta[,idx] 
    # Note  eta.cum  is eta.ncol x n.par while eta is n.par x eta.ncol
    Ds <- eta.cum[idx, ] + last.interval[i] * eta[, idx]

    #calcuate zs 
    # pred: n.obs x n.par
    zs <-pred * newtime[i]
    zs.scaled <- as.numeric(zs-Ds) 
    # L.var2(s) = \left\{zs-D(s)\right\}^TA^{-1}BA^{-1}\left\{zs-D(s)\right\}
    L.var2[i]<-sum((var %*% zs.scaled)* zs.scaled) 
    # L.cov[i](s) = \left\{zs-D(s)\right\}^TA^{-1}(L.cov1-L.cov2)
    L.cov[i]<-sum(iA%*%(L.cov1-L.cov2) *zs.scaled)

    L.var[i] <- L.var1[i] + L.var2[i] + 2*L.cov[i]
    
  }        

  return(L.var)
  invisible()

}

### This is not useful because it does't make sense to perform prediction if we don't assume the model
### However it will be used when estimating the two-phase estimators 
L.rvar<-function(newtime = newtime, pred = pred, t.pos=t.pos, wts=wts, t.unique=t.unique, eta.cum=eta.cum, eta=eta, resid=resid, L.resid=L.resid, ### resid is an output from fitted ah model 
         iA = iA, n.obs = n.obs, ...){

   t.mark <- t.unique[t.pos] # The closed time point in the original data timeline

  # calculate s - t.mark, the last intervla 
   last.interval <- newtime - t.mark  

    L.var<-NULL
    L.var1<-L.var2<-matrix(0, nrow=n.obs, ncol=length(t.pos))

    # See equantion before (3.42) for the formula for calculating the robust variance of predicted outcome   
    #  Let L.var1 = \int_0^{s}dM(t)/PY(t) and L.var2 = (zs-D(s)_A^{-1}\int_{0}^{\tau}[Z-\bar{Z}]dM(t)
    #  I calcuate the robust variance by (L.var1- L.var2)^2
    for ( i in 1:length(t.pos)){

      idx <- t.pos[i]
    
      L.var1[,i]<-L.resid[,i] #L.var1 is the L.resid at the closest t.fail on the left
    
      #calculate L.var2
      # calcuate D(s) where D(s) = D(t.mark[i]) + last.interval[s]*eta[,idx] 
      # Note  eta.cum  is eta.ncol x n.par while eta is n.par x eta.ncol
       Ds <- eta.cum[idx, ] + last.interval[i] * eta[, idx]

      #calcuate zs 
       # pred: n.obs x n.par
      zs <-pred * newtime[i]
      zs.scaled <- as.numeric(zs-Ds) 

      L.var2[,i]<- resid %*% iA%*%zs.scaled

      L.var[i]<-sum(((L.var1[,i]+L.var2[,i])*wts)^2)
    }
    return(L.var)
  }




  
  # Calculating \int_0^{s}dM(t)/PY(t)     for each individual  
  cook.L.resid<- function(time.pos=time.pos, match.eta=match.eta, n.obs=n.obs, 
    eta.den.vec=eta.den.vec, lambda0=lambda0,t1.unique =t1.unique, t.diff=t.diff,t.fail=t.fail, censor=censor,
    effects=effects,wts=wts,...){

        # I treat M(t) as a step function,  

        L.resid<-matrix(0, nrow=n.obs,ncol=length(time.pos))



        eta.den.lambda0<-cumsum(lambda0/eta.den.vec)
        eta.den.cum<-cumsum(t.diff/eta.den.vec)

        # match.eta<-match(t.fail, t.unique)
        eta.den.i<-eta.den.vec[match.eta]
        #eta.den.lambda0.i<-eta.den.lambda0[match.eta]
        
        #eta.den.cum.i<-eta.den.cum[match.eta]
        l<-length(time.pos)
        L1<-rep(0,l)
        L2<-eta.den.lambda0[time.pos]
        L3<-eta.den.cum[time.pos]

        # calculate L.resid for each individual at each observed time point 
        for ( k in 1: n.obs){
            L.resid1<-L1
            L.resid2<-L2
            L.resid3<-effects[k]*L3

            a<-which(t1.unique >=t.fail[k])[1]
            if(!is.na(a)){
               idx5<-a
               
               # L.resid1 <- int_0^{s}dN_i(t)/PY(t)
               L.resid1[idx5:l]<- censor[k]/eta.den.i[k] # when t< T_i, L.resid = 0; whent >=T_i L.resid is always 
                                                        #the same since N_i(t) is an indicator function 
               
               #L. resid2 <- int_0^{s} Y(t)d\Lambda(t)/PY(t) = int_0^{t.mark} Y(t)d\Lambda(t)/PY(t)
               L.resid2[idx5:l]<-L.resid2[idx5] # when s < T_i, L.resid2= int_0^{s} dLambda(t)/PY(t), 
                                                # based on our estimators for Lambda(t), dLambda(t) =0 for any t between two marks
                                                # when s > T_i  L.resid2= int_0^{T_i} dLambda(t)/PY(t), 
               #L. resid3<- int_0^{s} Y(t)Z^T\theta/PY(t)dt = nt_0^{s} Y(t)Z^T\theta/PY(t)dt 

               L.resid3[idx5:l]<-L.resid3[idx5] # think dt as Delta(t) which only have values at t.failures
               
            }
            
            
            L.resid[k,]<-(L.resid1-L.resid2-L.resid3)*wts[k]
        }
        return(L.resid)
        invisible()
  }




   #################################################################################################
   #################################  calculate the calibrated weight  #####################################
   ##################################################################################################
   ##################################################################################################

cook.wts.cal<-function(aux,aux.pha2,P,wts.pha2){
  
          

        
           if(!is.matrix(aux.pha2)) aux.pha2<-as.matrix(aux.pha2)
           aux.tot<-apply(aux,2,sum)
           aux.tot.pha2<-apply(aux.pha2*wts.pha2,2,sum)
           ### phase I total, 1 x q
          
           L0<-solve(P)%*%(aux.tot.pha2-aux.tot)

          
          
            
            
            model.calibration<-function(L){
              F<-NULL
              wts.fish<-as.vector(exp(-aux.pha2%*%L)*wts.pha2)
                 for (i in 1:dim(aux)[2] ){
                    F[i]<-sum(wts.fish*aux.pha2[,i])-aux.tot[i]
                  }
              F
                }

               eval<-rootSolve::multiroot(model.calibration,start=L0)
               L<-eval$root
              # est.acc<-eval$est.acc
               wts.cal<-as.vector(exp(-aux.pha2%*%L)*wts.pha2)
              return(wts.cal)
              invisible()
  }
       



#### Calculate the phase II calibrated variance for predicted values   
L.rvar.calibration<-function(pred = pred, newtime = newtime, t.pos=t.pos,  t.unique, eta.cum=eta.cum, eta = eta, resid=resid, L.resid=L.resid,
      iA=iA, n.obs =n.obs, new.wts=new.wts, aux.pha2=aux.pha2, P=P, wts.cal=wts.cal, ...){
           

           ## See page 109 of Jie Hu' thesis for the detailed formula 
           ## L.var2 = Q[\sqrt{\frac{1-\pi(V)}{\pi(V)} } * (\psi - Q {\psi \tilda{V}) * Q[\tilda{V}\tilda{V}^{-1}] * \tilda{V}})]^2
           ## I devide \psi to two parts \psi = L. resid + resid%*%iA%*% zs.scaled 
           ## Then sqrt(L.var2 )= new.weights * (L.var1 +L.var2) = new.weights [Q1 +Q2 *%*%iA%*% zs.scaled ]
           ## Q1 = L.resid - QL.resid \tilde{V} *   Q[\tilda{V}\tilda{V}^{-1}] * \tilda{V}})
           ## Q2 = resid - Qresid \tilde{V} *  Q[\tilda{V}\tilda{V}^{-1}] * \tilda{V}})
           
           ## Note expectation Qf= sum wts.cal*f

           # Q<- t(aux.pha2*sqrt(wts.cal))%*% (resid/sqrt(wts.cal))  
           # resid.adj<-resid-(aux.pha2%*%solve(P)%*% Q )*wts.cal     
            #temp<-resid.adj*sqrt((1-Pi.pha2)/(Pi.pha2*wts.pha2))

           t.mark <- t.unique[t.pos] # The closed time point in the original data timeline

           # calculate s - t.mark, the last intervla 
           last.interval <- newtime - t.mark  
 
           L.var<-NULL
           L.var1<-L.var2<-matrix(0, nrow=n.obs,ncol=length(t.pos))
           
           ## Q2: q x n  X n x1 = q x1 
           Q2<- t(aux.pha2*sqrt(wts.cal))%*% (resid/sqrt(wts.cal)) # resid are already mutiplied by wts.cal 
                                                                   # multiplied by sqrt(wts.cal) because Qf= sum wts.cal*f

           resid.adjust <- resid-(aux.pha2%*%solve(P)%*% Q2 )*wts.cal    ###  the 2nd term * wts.cal because resid are already weighted by wts.cal    

           for (i in 1:length(t.pos)){
                
                idx<-t.pos[i]
                # Calculate Q1
                 ## q x n   times nx1
                 Q1<- t(aux.pha2*sqrt(wts.cal))%*% (L.resid[,i]/sqrt(wts.cal)) # multiplied by sqrt(wts.cal) because Qf= sum wts.cal*f
                                                                               # L.resid are already weighted by wts.cal 

                ## L.resid[,i]: nx1, aux.pha2: n x q , P: q x q, Q1: qx1
                L.var1[,i] <- L.resid[,i]-(aux.pha2%*%solve(P)%*% Q1 ) * wts.cal  ###  the 2nd term * wts.cal because L.resid[,i] are already weighted by wts.cal 

                ### Calculate zs.scaled
                # calcuate D(s) where D(s) = D(t.mark[i]) + last.interval[s]*eta[,idx] 
                # Note  eta.cum  is eta.ncol x n.par while eta is n.par x eta.ncol
                Ds <- eta.cum[idx, ] + last.interval[i] * eta[, idx]

                #calcuate zs 
                # pred: n.obs x n.par
                zs <-pred * newtime[i]
                zs.scaled <- as.numeric(zs-Ds) 

                L.var2[,i]<- resid.adjust %*%iA%*%zs.scaled 

                    

                ## because Qf= sum wts.cal*f and L.var1 and L.var are a product something and wts.cal
                ## Thus new.wts=sqrt((1-Pi.pha2)/(wts.cal*Pi.pha2))

                L.var[i]<-sum(((L.var1[,i]+L.var2[,i])*new.wts)^2)
     



          }


        return(L.var)
        invisible()
}     



