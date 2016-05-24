

llk.C <- function(para,DT,IPRT,REint,AR,absTol) 
  {
    ## REint == 1 gamma
    ## REint == 2 log-normal
    if (AR){
      alpha <- exp(para[1])
      theta <- exp(para[2])
      delta <- ilgt(para[3])
      betas <- para[4:length(para)]
    }else{
      alpha <- exp(para[1])
      theta <- exp(para[2])
      betas <- para[3:length(para)]
      delta <- 0 
    }
    ## Add intercept term
    DT$x <- cbind(rep(1,length(DT$y)),DT$x)
    if(IPRT) cat("\n",para," ")
    nlk <- .Call("nllk_C",
                 as.numeric(DT$y), 
                 as.numeric(c(DT$x)),  
                 as.integer(DT$id - 1),
                 as.numeric(alpha), 
                 as.numeric(theta),
                 as.numeric(delta),
                 as.numeric(betas), 
                 as.integer(max(DT$ni)), 
                 as.integer(length(unique(DT$id))),
                 as.integer(REint),
                 as.numeric(DT$dif),
                 as.numeric(absTol),
                 package ="lmeNBBayes" ## Total number of patients
                 )[[1]]
    if (IPRT) cat("nlk ", nlk)
    return(nlk)
  }

fitSemiAR1 <- function(
                       formula,     ## an object of class "formula"
                       ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
                       data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
                       ## containing the variables in the model.
                       ID,          ## a vector of length n*ni containing patient IDs of observations in data
                       Vcode,
                       p.ini=NULL, # log(a), log(var(G)), logit(d), b0, b1,... if its NULL then initial values are set to 0.1,NULL,0.1,log(mean(y))+0.1,0.1,0.1,...
                       IPRT=TRUE,
                       ##FixN=-1, #see fitParaAR1
                       deps=1.e-3,      # stop iteration when max(bhat.new - bhat.old) < deps
                       maxit=100
                       )
{
  dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  DT <- getDT(dat$dat)
  
  ## DT = (id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames)

  ## DT$dif contains Vcode[j+1]-Vcode[j] except the initial scan of each patient (which is zero).
  DT$dif <- c(0,diff(Vcode))
  DT$dif[DT$ind] <- 0
  DT$dif <- DT$dif[1:DT$totN] #dif = scan lag 
  
  Dss <- getDjj(DT, Vcode)  # j-j' for all posible j > j'
  
  if (is.null(p.ini))
    {
      p.ini=rep(0, 4+DT$cn)
      p.ini[4]=mean(DT$y)   
    }
  p.new <- p.ini  #initial values for the first iteration
  Wh <- NULL       #weights for GLS

  dth <- 1  #dth = max(bhat.new-bhat.old)
  ahat <- 0 #initial value for alpha
  counter <- 1
  repeat
    {
      p.ini <- p.new #reset the initial values at each iteration
      b.ini <- p.ini[-(1:3)] ## b0,b1,...
      
      ## 
      ## Main Step 1 : estimate b using Generalized Least Square
      ## 
      if (DT$cn>0) #with covariates
        {
          temB <- optim(b.ini+0.1,## initial values of coefficients b0,b1,...
                        estb.ar1, ## weighted residuals function
                        hessian=TRUE,
                        dat=DT,## output of getDT
                        Wh=Wh## list of N elements, containing ni by ni weight matrix for each patient
                        )
          bhat <- temB$par 
        }else{  #without covariates
          tem <- optimize(estb.ar1, lower=-5, upper=5, dat=DT, Wh=Wh)
          bhat <- tem$min
        }

      ## Stopping criterion
      ## check if max(bhat.new - bhat.old) < deps
      dth <- max(abs(bhat-b.ini))
      ## if (IPRT) print(dth)
      if (dth < deps || counter > maxit) break

      ## 
      ## Main Step 2: estimate g_i i=1,...,N, g_i = w_i*(y_i+/mu_i+)+(1-w_i)
      ## 
      ## Substep 1: the estimation of mu_i+
      ##         hat.mu_i+ = exp( X_i^T bhat)
                                        # compute exp(b0+b1*x1+ ...)
      mu_ij <- rep(exp(bhat[1]), DT$totN)
      if (DT$cn>0) {
        tem <- exp(DT$x%*%bhat[-1]) 
        mu_ij <- tem*mu_ij ## vector of length Ntot containing hat.mu_ij
      } 
      ## mu_ip = u.i+ 
      mu_ip <- as.vector(tapply(mu_ij, DT$id, sum)) ## mu_i = sum_j mu_ij
      gh0 <- DT$ys/mu_ip ## vector of length n containing y_i+/mu_i+
      
      ## Substep 2: the estimation of w_i = sqrt(var(G_i)/var(Y_i+/mu_i+))
      ## weights for gi
      ## compute wi =sqrt(var(G)/Var(Y+/u+))
      if (ahat == 0) #wi =1 at the first iteration
        {
          wi <- rep(1, DT$np)
        }else{
          gsig <- by(data=cbind(mu_ij, Vcode), INDICES=DT$id, FUN=getWhs.ar1,
                     th=that, a=ahat, dt=dht, act="sum")/mu_ip^2 
          wi <- sqrt(that/gsig)  #that =var(G), gsig=var(Y+/u+) 
        }
      
      gh <- wi*gh0+1-wi

      ## normalized gh so that the average is 1
      gh <- gh/mean(gh)
      gh <- as.vector(gh)
      
      if (IPRT) {
        cat("\n iteration",counter)
        ## cat("\n The distribution of hat.g \n")
        ## print(summary(gh))
      }
      ## frequency table for gh
      tem <- sort(round(gh,6))
      gh1 <- unique(tem)
      ghf <- table(tem)
      gtb <- cbind(ghf, gh1)
      rownames(gtb) <- NULL
      gtb1 <- gtb
      gtb1[,1] <- gtb[,1]/DT$np #covert freq to prop
      ##estimate alpha and delta
      ## ~~~~~ Main Step 3: estimate a and d using profile likelihood ~~~~~ 
      ad.ini <- p.ini[c(1,3)]

      tt <- optim(ad.ini+0.1, ## c(log(alpha), logit(delta))
                  ar1.ad.lk,
                  hessian=TRUE,
                  dat=DT,
                  IPRT=FALSE,
                  gTB=gtb1,  ## frequency table for hat(gi): with proportions
                  uhat=mu_ij ## hat(u.ij)
                  )
      
      ahat <- exp(tt$par[1])
      dht <- ilgt(tt$par[2])

      ##  ~~~~~ Main Step 4: moment estimate of var(G) ~~~~~ 
      that <- estSg3.ar1(DT, mu_ip, mu_ij, ahat, dht, Dss$dis) 
      
      ##weights for GLS
      Wh <- by(data=cbind(mu_ij, Vcode), INDICES=DT$id,
              FUN=getWhs.ar1, th=that, a=ahat, dt=dht, act="inv") 

      ##update the estimate
      p.new <- c(tt$par[1], log(that), tt$par[2], bhat)

      if (IPRT){
        cat("\n log(alpha), log_var(G), lgt_d, log_u0", DT$xnames)
        cat("\n", p.new)
      }
      counter <- counter + 1
    }

  if (maxit==counter) warning("The maximum number of iterations occur!")
  vcm <- NULL
  if (DT$cn>0) vcm=solve(temB$hessian) #vcm=vcm.fun(he=temB$hessian, corr=F)
  p.new <- matrix(p.new,ncol=1)
  rownames(p.new)=c( "log_a", "log_var(G)", "lgt_d", "(Intercept)", DT$xnames)
  if (is.matrix(vcm)) colnames(vcm) <- rownames(vcm) <- c("(Intercept)", DT$xnames)
  re <- list(opt=tt, V=vcm, est=p.new, gi=as.vector(gh), RE="NoN",##idat = data.frame(dat),
             Vcode=Vcode,AR=TRUE,formula=formula)
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}


estb.ar1 <- function(b, # (b0, b1, ...)
                 dat, #data see fitParaIND
                 Wh=NULL ## inverse of Var(Yi)
                 ## a list of N elements.
                 ## each element is ni by ni matrix of weights
                 ## If Wh=NULL then the weight matrix are all regarded as identity
                 )
{
  ## GLS compute the weighted sum of the residuals
  mu.i <- exp(b[1]) #mu.i =exp(b0)
  ## If the number of covariate is nonzero:
  if (dat$cn>0) 
  {
    residual <- exp(dat$x%*%b[-1])
    mu.i <- residual*mu.i #mu.i = exp(b0+b1*x1 + ...) n*sn by 1 vector
  }
  residual <- dat$y-mu.i
  ## If weights are not provided then  W_i = diag(dt$vn) identity
  ## so the weighted sum of residuals is
  ## sum_i (y_i-X_i^T beta)^T (y_i-X_i^T beta) = sum(residual^2)
  if (is.null(Wh)){ 
    val <- sum(residual^2)
  }else{
    val <- 0
    for (ipat in 1 : dat$np)
    {
      ## pick is a vector of length ni containing the locations of
      ## repeated measures of patient ipat
      weightMat <- Wh[[ipat]]
      pick <- dat$ind[ipat]:(dat$ind[ipat+1]-1)
      val <- val+t(residual[pick])%*%weightMat%*%residual[pick]
    }
  }
  return(val)
}



getWhs.ar1 <- function(ud,  #cbind(mu_ij, Vcode)
                       th=exp(1.3), #Var(G)
                       a=exp(-0.5), dt=0.5, 
                       act="inv" #  transformation option: "none", "inverse", "sum"
                       )
{ mu_ij  <- ud[,1]
  Vcode <- ud[,2]

  w1 <- mu_ij%*%t(mu_ij)*th # u.i*u.j*Var(G)

  n <- length(mu_ij)
  Um <- diag(x=mu_ij,nrow=n)%*%matrix(1,n,n)
  mU <- pmin(Um, t(Um)) #pairwise min (u.j, u.j')

  Dm <- diag(x=Vcode,nrow=n)%*%matrix(1,n,n)
  mD <- abs(Dm-t(Dm))  

  Md <- dt^mD #Mdis

  UD <- mU*Md #min(uj, uj')*d^(distance)
  w2 <- (1+(th+1)*a)*UD   #[d^abs(j-j')]*u.l*(a*Var(G)+a+1); l=min(j,j')
  
  W <- w1+w2
  if (act=="inv") W=solve(W)
  if (act=="sum") W=sum(W)
  return(W)
}



##subroutines for mle.ar1.non1 and fitSemiAR1
                                        #likelihood function for estimating alpha and delta
ar1.ad.lk <- function(para=c(-0.5, -1), #c(log(alpha), logit(delta)) 
                      dat,  #see ar1.lk
                      IPRT=TRUE, 
                      gTB, #frequency table for hat(gi): with proportions
                      uhat #hat(u.ij)
                      ##,FixN=9  #see ar1.lk
                      )
{
  ## Li(a,d;bhat,ghat,ys)
  ## = 1/N* sum_{l=1}^N Pr(Y_i,1=y_i,1;ghat,bhat,a)*prod_{j=2}^ni Pr(Y_ij=y_ij;y_i(j-1),ghat,bhat,a) 
  if (IPRT) cat(para)
  ainv=exp(-para[1]) #ainv=1/a  
  dt=ilgt(para[2])   #delta
  
  th2=uhat*ainv  #r.ij
  Dl=dt^dat$dif  

  szm=c(0, th2[-dat$totN])
  U=Dl*szm
  V=szm-U
  sz2=th2-U
  tem=dat$ind[1:dat$np]
  sz2[tem]=th2[tem]
  if (any(sz2<=0)) return(1.e15)

  nllk=0
  lk0=0 #likelihood for a patient with all 0 counts
  
  Pr=ainv/(ainv+gTB[,2]) ## in this manuscript p=p
                                        #fq=gTB[,1]/dat$np      #proportions
  
  for (i in 1:dat$np)
    {

      ##if (lk0>0&dat$ni[i]==FixN&dat$ys[i]==0) lki=lk0
      ##else{
      ## Li(a,d;bhat,ghat,ys)
      ## = 1/N* sum_{l=1}^N Pr(Y_i,1=y_i,1;ghat,bhat,a)*prod_{j=2}^ni Pr(Y_ij=y_ij;y_i(j-1),ghat,bhat,a)
      ll=dat$ind[i]:(dat$ind[i+1]-1)
      tem=ar1.non(Pr=Pr, y=dat$y[ll], u=U[ll], v=V[ll], s2=sz2[ll])
      lki=sum(tem*gTB[,1])
      ##if (dat$ys[i]==0&dat$ni[i]==FixN) lk0=lki
      ##}
      nllk = nllk - log(lki)
    }
  if (IPRT) cat(" nllk=", nllk, "\n")
  return(nllk)
}

                                        #ar1.non=function(Pr=Pr, y=1:3, sz=rep(3,3), dt=0.3)
ar1.non=function(Pr, y=1:3, u=c(0,1,1), v=c(0, 2,2), s2=c(3,2,2))
{ #compute ar1.fun for a vector of Pr
                                        #also see ar1.intg
  tem=NULL
  for ( pr in  Pr)
    { #print(pr)
      tem=c(tem,ar1.fun(y, u, v, s2, pr)) 
    }
  return(tem)
}


                                        #moment estimate of var(G)
                                        #estSg3.ar1(DT, uhat, uh, ahat, dht, Dss$dis)
estSg3.ar1=function(dat, #see ar1.lk
  ui, uij,     #u.i+ and u.ij
  ahat, dhat,  #alpha hat and delta hat
  Dis          # j - j'
  ) 
{ 
  ssy=sum(dat$ys^2)-sum(dat$y^2)
  ssu=sum(ui^2)-sum(uij^2)
  
  DU=2*sum(dhat^Dis[,1]*uij[Dis[,2]]) 
  that=(ssy-DU)/(ssu+DU*ahat)
                                        #print(ssy/ssu)
  return(max(that-1, 0))
}

                                        # compute j-j' for all posible j > j'
getDjj <- function(DT,Vcode)
{
  kk <- Nj <- dj <- NULL
  for (i in 1:DT$np)
    {
      if (DT$ni[i] > 1)
        {
          ll <- DT$ind[i]:(DT$ind[i+1]-1)
          iVcode <- Vcode[ll]
          nj <- djj <- NULL
          id <- DT$id[DT$ind[i]]
          
          Ni <- DT$ni[i]-1
          for (j in 1:Ni )
            { 
              for (jj in (j+1):DT$ni[i])
                djj <- c(djj, iVcode[jj]-iVcode[j])
            }
          nj <- rep(1:Ni,Ni:1)
          kk <- c(kk,rep(id, length(nj)))
          Nj <- c(Nj, nj+DT$ind[i])
          dj <- c(dj, djj)
        }
    }
  return(list(id=kk, dis=cbind(dj, Nj-1)))
                                        # dis=([i]j'-[i]j, index for u.ij)
}

##route for weights
## act= inv => W = Var(Y)(^-1)
## act= sum => W= sum(Var(Y))
## Var(Y)=Cov(Yj, Yj')=u.j*u.j'*Var(G)+[d^abs(j-j')]*u.l*(a*Var(G)+a+1); 
## where l=min(i,j)

                                        #example ud.ex=cbind(rep(exp(.48),5), 1:5)
                                        #getWhs.ar1(ud=ud.ex, th=exp(1.46), a=exp(-0.825), dt=ilgt(-0.26), act="sum")
##==========  end SP =============



rbb <- function(n, N, u, v) {
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}

                                        #logit and inverse logit
lgt <- function(p)
{
  log(p/(1-p))
}

ilgt <- function(x)
{
  tem=exp(x)
  res=tem/(1+tem)
  return(res)
}


