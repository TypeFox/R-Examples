#' AME model fitting routine
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to relational data of various types
#' 
#' This command provides posterior inference for parameters in AME models of
#' relational data, assuming one of six possible data types/models:
#' 
#' "nrm": A normal AME model.
#' 
#' "bin": A binary probit AME model.
#' 
#' "ord": An ordinal probit AME model. An intercept is not identifiable in this
#' model.
#' 
#' "cbin": An AME model for censored binary data.  The value of 'odmax'
#' specifies the maximum number of links each row may have.
#' 
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#' 
#' "rrl": An AME model based on the row ranks. This is appropriate if the
#' relationships across rows are not directly comparable in terms of scale. An
#' intercept, row random effects and row regression effects are not estimable
#' for this model.
#' 
#' @usage ame(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(model=="rrl") ,
#' cvar = TRUE,  dcor = !symmetric, nvar=TRUE, R = 0, model="nrm",
#' intercept=!is.element(model,c("rrl","ord")),
#' symmetric=FALSE,
#' odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)), seed = 1, nscan =
#' 10000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE)
#' @param Y an n x n square relational matrix of relations. See model below for
#' various data types.
#' @param Xdyad an n x n x pd array of covariates
#' @param Xrow an n x pr matrix of nodal row covariates
#' @param Xcol an n x pc matrix of nodal column covariates
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)?  
#' @param dcor logical: fit a dyadic correlation (asymmetric case)?
#' @param nvar logical: fit nodal random effects (symmetric case)?
#' @param R integer: dimension of the multiplicative effects (can be zero)
#' @param model character: one of "nrm","bin","ord","cbin","frn","rrl" - see
#' the details below
#' @param intercept logical: fit model with an intercept? 
#' @param symmetric logical: Is the sociomatrix symmetric by design?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' models
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param plot logical: plot results while running?
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @return \item{BETA}{posterior samples of regression coefficients}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} \item{V}{posterior mean of multiplicative column effects v (asymmetric case)}
#' \item{UVPM}{posterior mean of UV (asymmetric case)} 
#' \item{ULUPM}{posterior mean of ULU (symmetric case)} 
#' \item{L}{posterior mean of L (symmetric case)} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' @author Peter Hoff
#' @examples
#' 
#' data(YX_frn) 
#' fit<-ame(YX_frn$Y,YX_frn$X,burn=5,nscan=5,odens=1,model="frn")
#' # you should run the Markov chain much longer than this
#' 
#'  
#' @export ame
ame<-function (Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
       rvar = !(model=="rrl") , cvar = TRUE, dcor = !symmetric, 
       nvar=TRUE, 
       R = 0,
       model="nrm",
       intercept=!is.element(model,c("rrl","ord")), 
       symmetric=FALSE,
       odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
       seed = 1, nscan = 10000, burn = 500, odens = 25,
       plot=TRUE, print = TRUE, gof=TRUE)
{ 
 
  # set random seed
  set.seed(seed)


  # set diag to NA
  diag(Y) <- NA 

  # force binary if binary model specified
  if(is.element(model,c("bin","cbin"))) { Y<-1*(Y>0) } 

  # observed and max outdegrees 
  if(is.element(model,c("cbin","frn","rrl")))
  {
    odobs<-apply(Y>0,1,sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y)) }
  }

  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }

  # construct design matrix
  n<-nrow(Y) 
  pr<-length(Xrow)/n
  pc<-length(Xcol)/n
  pd<-length(Xdyad)/n^2
  X<-design_array(Xrow,Xcol,Xdyad,intercept,nrow(Y)) 

  # design matrix warning for rrl 
  if( model=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0) 
                   & !any( apply(X,3,function(x){var(c(x))})==0) ) 
  {
    cat("WARNING: row effects are not estimable using this procedure ","\n")
  } 

  # design matrix warning for rrl and ord
  if( is.element(model,c("ord","rrl")) & 
      any( apply(X,3,function(x){var(c(x))})==0 ) )
  {
    cat("WARNING: an intercept is not estimable using this procedure ","\n")
  } 
  
  # construct matrix of ranked nominations for frn, rrl 
  if(is.element(model,c("frn","rrl")))
  {
    ymx<-max(apply(1*(Y>0),1,sum,na.rm=TRUE))
    YL<-NULL
    warn<-FALSE
    for(i in 1:nrow(Y))
    {
      yi<-Y[i,] ; rnkd<-which( !is.na(yi)&yi>0 )
      if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
      yi[rnkd]<-rank(yi[rnkd],ties.method="random")
      Y[i,]<-yi
      YL<-rbind(YL, match(1:ymx,yi))
    }
    if(warn){cat("WARNING: Random reordering used to break ties in ranks\n")}
  }

  # starting Z values
  if(model=="nrm") { Z<-Y }
  if(model=="ord") { Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) } 
  if(model=="rrl") { Z<-matrix(t(apply(Y,1,zscores)),nrow(Y),ncol(Y)) }  
  if(model=="bin")
  { 
    Z<-matrix(zscores(Y),nrow(Y),nrow(Y)) 
    z01<- .5* ( max(Z[Y==0],na.rm=TRUE) + min(Z[Y==1],na.rm=TRUE) ) 
    Z<-Z - z01
  } 

  if(is.element(model,c("cbin","frn")))
  {
    Z<-Y
    for(i in 1:nrow(Y))
    {
      yi<-Y[i,]
      zi<-zscores(yi)
      rnkd<-which( !is.na(yi) & yi>0 ) 
      if(length(rnkd)>0 && min(zi[rnkd])<0) 
      { 
        zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 
      }

      if(length(rnkd)<odmax[i]) 
      {
        urnkd<-which( !is.na(yi) & yi==0 ) 
        if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
      }

      Z[i,]<-zi
    } 
  }


  # starting values for missing entries 
  mu<-mean(Z,na.rm=TRUE) 
  a<-rowMeans(Z,na.rm=TRUE) ; b<-colMeans(Z,na.rm=TRUE)  
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  ZA<-mu + outer(a,b,"+") 
  Z[is.na(Z)]<-ZA[is.na(Z)] 
   

  # other starting values
  beta<-rep(0,dim(X)[3]) 
  s2<-1 
  rho<-0
  Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
  U<-V<-matrix(0, nrow(Y), R)  


  # output items
  BETA <- matrix(nrow = 0, ncol = dim(X)[3] - pr*symmetric)
  VC<-matrix(nrow=0,ncol=5-3*symmetric) 
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,nrow(Y))  
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y)
  GOF<-matrix(gofstats(Y),1,4)  
  rownames(GOF)<-"obs"
  colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
  names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y)
 
  # names of parameters, asymmetric case 
  if(!symmetric)
  {
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve") 
    colnames(BETA) <- dimnames(X)[[3]] 
  }

  # names of parameters, symmetric case
  if(symmetric)
  {
    colnames(VC) <- c("va", "ve")  
    rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
    bnames<-dimnames(X)[[3]]
    bni<-bnames[1*intercept] 
    bnn<-gsub("row",bnames[rb],replacement="node")       
    bnd<-bnames[-c(1*intercept,rb,cb)]
    colnames(BETA)<-c(bni,bnn,bnd) 
  }    

  # marginal means and regression sums of squares
  Xr<-apply(X,c(1,3),sum)            # row sum
  Xc<-apply(X,c(2,3),sum)            # col sum
  mX<- apply(X,3,c)                  # design matrix
  mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
  XX<-t(mX)%*%mX                     # regression sums of squares
  XXt<-t(mX)%*%mXt                   # crossproduct sums of squares

  # MCMC 
  have_coda<-suppressWarnings(
               try(requireNamespace("coda",quietly = TRUE),silent=TRUE))

  for (s in 1:(nscan + burn)) 
  { 

    # update Z 
    EZ<-Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
    if(model=="nrm"){ Z<-rZ_nrm_fc(Z,EZ,rho,s2,Y) } 
    if(model=="bin"){ Z<-rZ_bin_fc(Z,EZ,rho,Y) }
    if(model=="ord"){ Z<-rZ_ord_fc(Z,EZ,rho,Y) }
    if(model=="cbin"){Z<-rZ_cbin_fc(Z,EZ,rho,Y,odmax,odobs)}
    if(model=="frn"){ Z<-rZ_frn_fc(Z,EZ,rho,Y,YL,odmax,odobs)}
    if(model=="rrl"){ Z<-rZ_rrl_fc(Z,EZ,rho,Y,YL)} 

    # update s2
    if (model=="nrm") s2<-rs2_fc(Z-EZ,rho)  

    # update beta, a b
    tmp <- rbeta_ab_fc(Z-U%*%t(V), Sab, rho, X, mX, mXt,XX, XXt, Xr, Xc, s2)
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar 
    if(symmetric){ a<-b<-(a+b)/2 }

    # update Sab - full SRM
    if(rvar & cvar & !symmetric)
    {
      if(is.element(model,c("nrm","ord")))
      { 
        Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a, b))),3+nrow(Z))) 
      }

      if(model=="bin")
      {
        tmp<-raSab_bin_fc(Z,Y,a,b,Sab) ; Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }

      if(model=="cbin")
      {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
      if(model=="frn")
      { 
        tmp<-raSab_frn_fc(Z,Y,YL,a,b,Sab,odmax,odobs) 
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a 
      }
    }

    # update Sab - rvar only
    if (rvar & !cvar & !symmetric) 
    {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(a^2))/2)
    }

    # update Sab - cvar only
    if (!rvar & cvar & !symmetric) 
    {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(b^2))/2)
    }

    # update Sab - symmetric case
    if(symmetric & nvar)
    { 
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]   
    }

    # update rho
    if(dcor) 
    {
      rho<-rrho_mh(Z-(Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)), rho,s2)
    }

    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }

    # update U,V
    if (R > 0) 
    { 
      E<-Z-(Xbeta(X,beta)+outer(a,b,"+"))  ; if(symmetric){ E<-.5*(E+t(E)) }
      shrink<- (s>.5*burn)

      if(symmetric){ UV<-rUV_sym_fc(E, U, V, s2,shrink) }
      if(!symmetric){UV<-rUV_fc(E, U, V,rho, s2,shrink) }

      U<-UV$U ; V<-UV$V
    }    

    # burn-in countdown
    if(s%%odens==0&s<=burn){cat(round(100*s/burn,2)," pct burnin complete \n")}

    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    {  

       # save BETA and VC - symmetric case 
      if(symmetric)
      {
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2 
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA<-rbind(BETA,sbeta)

        VC<-rbind(VC,c(Sab[1,1],s2) )
      }

      # save BETA and VC - asymmetric case 
      if(!symmetric)
      {
        BETA<-rbind(BETA, beta)
        VC<-rbind(VC, c(Sab[upper.tri(Sab, diag = T)], rho,s2)) 
      }

      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b 

      # simulate from posterior predictive 
      EZ<-Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V) 
      if(symmetric){ EZ<-(EZ+t(EZ))/2 } 

      if(model=="bin") { Ys<-simY_bin(EZ,rho) }
      if(model=="cbin"){ Ys<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
      if(model=="frn") { Ys<-simY_frn(EZ,rho,odmax,YO=Y) }
      if(model=="rrl") { Ys<-simY_rrl(EZ,rho,odobs,YO=Y ) }
      if(model=="nrm") { Ys<-simY_nrm(EZ,rho,s2) }
      if(model=="ord") { Ys<-simY_ord(EZ,rho,Y) } 

      if(symmetric){ Ys[lower.tri(Ys)]<-0 ; Ys<-Ys+t(Ys)  }

      # update posterior sum
      YPS<-YPS+Ys

      # save posterior predictive GOF stats
      if(gof){ Ys[is.na(Y)]<-NA ; GOF<-rbind(GOF,gofstats(Ys)) }

      # print MC progress 
      if (print) 
      {
        cat(s,round(apply(BETA,2,mean),2),":",round(apply(VC,2,mean),2),"\n")
        if (have_coda & nrow(VC) > 3 & length(beta)>0) 
        {
          cat(round(coda::effectiveSize(BETA)), "\n")
        }
      }

      # plot MC progress
      if(plot) 
      { 
        # plot VC
        par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
        mVC <- apply(VC, 2, median)
        matplot(VC, type = "l", lty = 1)
        abline(h = mVC, col = 1:length(mVC)) 

        # plot BETA
        if(length(beta)>0) 
        {
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray") 
        } 
 
        # plot GOF 
        if(gof)
        {
          for(k in 1:4)
          {
            hist(GOF[-1,k],xlim=range(GOF[,k]),main="",prob=TRUE,
                 xlab=colnames(GOF)[k],col="lightblue",ylab="",yaxt="n")  
            abline(v=GOF[1,k],col="red") 
          }
        } 

      }


    }


  } # end MCMC   

  # output 

  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)
  UVPM<-UVPS/nrow(VC)
  YPM<-YPS/nrow(VC) 
  EZ<-Xbeta(X,apply(BETA,2,mean)) + outer(APM,BPM,"+")+UVPM  

  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
  dimnames(YPM)<-dimnames(EZ)<-dimnames(Y)
  rownames(BETA)<-NULL

  # asymmetric output 
  if(!symmetric) 
  {
    UDV<-svd(UVPM)
    U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    rownames(U)<-rownames(V)<-rownames(Y) 
    fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
                YPM=YPM,GOF=GOF)
  }

  # symmetric output
  if(symmetric) 
  {
    ULUPM<-UVPM
    eULU<-eigen(ULUPM)
    eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
    U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
    L<-eULU$val[eR]
    rownames(U)<-rownames(ULUPM)<-colnames(ULUPM)<-rownames(Y)
    EZ<-.5*(EZ+t(EZ)) ; YPM<-.5*(YPM+t(YPM)) 
    fit<-list(BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
              YPM=YPM,GOF=GOF)
  } 

  class(fit) <- "ame"
  fit
}



