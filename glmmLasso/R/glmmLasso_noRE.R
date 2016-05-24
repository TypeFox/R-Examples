est.glmmLasso.noRE<-function(fix,data,lambda,family,final.re,switch.NR,control)
{  
  if(grepl("\\*", fix[3]))
    stop("Usage of '*' not allowed in formula! Please specify the corresponding variables separately.")  
  
  fix.old<-fix
  ic.dummy<-attr(terms(fix),"intercept")

  y <- model.response(model.frame(fix, data))
  very.old.names<-attr(terms(fix),"term.labels")

  if(ic.dummy!=1 && sum(substr(very.old.names,1,9)=="as.factor")>0){
    fix.help <- update(fix,~ .+1)
    orig.names <- colnames(model.matrix(fix.help, data))[-1]
  }else{
    orig.names <- colnames(model.matrix(fix, data))
  }
  
  control<-do.call(glmmLassoControl, control)
  
#   if(!is.null(control$index))
#   {
#     order.vec<-order(control$index)
#     very.old.names<-very.old.names[order.vec]
#     control$index<-control$index[order.vec]
#   }else{
#     control$index<-1:length(very.old.names)
#   }
  

    if(is.null(control$index))
    control$index<-1:length(very.old.names)

  if(length(control$index)!=length(very.old.names))
    stop("Length of vector defining the grouping of the variables doesn't match with 
         the formula!")
  
  attr(control$index,"names")<-very.old.names
  
  fix<-formula(paste("y~-1+",paste(very.old.names,collapse="+"))) 
  
  if(ic.dummy==1)
  {
    fix<-update(fix,~ .+1) 
    control$index<-c(NA,control$index)
    names(control$index)[1]<-"(Intercept)"
  }
  
  index.new<-c()
  fac.variab<-logical()
  for(i in 1:length(control$index))
  {
    if(!grepl("as.factor",names(control$index)[i]))
    {
      index.new<-c(index.new,control$index[i]) 
      fac.variab<-c(fac.variab,F)
    }else{
      if(!grepl("\\:", names(control$index)[i]))
      {  
        fac.name<-strsplit(strsplit(names(control$index)[i],"\\(")[[1]][2],"\\)")[[1]][1]
      }else{
        fac.name<-unlist(strsplit(unlist(strsplit(unlist(strsplit(names(control$index)[i],"\\(")),"\\)")),"\\:"))
        fac.name<-paste(fac.name[2],":",fac.name[length(fac.name)],sep="")
      }
      if(!grepl("\\:", fac.name))
      {
        length.fac<-length(levels(as.factor(data[,fac.name])))-1
      }else{
        length.fac<-(length(levels(data[,strsplit(fac.name,":")[[1]][1]]))-1)*(length(levels(data[,strsplit(fac.name,":")[[1]][2]]))-1)
      }
      index.new<-c(index.new,rep(control$index[i],length.fac))
      fac.variab<-c(fac.variab,rep(T,length.fac))
    }
  }

  if(ic.dummy!=1 && sum(substr(very.old.names,1,9)=="as.factor")>0){
    fix.help <- update(fix,~ .+1)
    X <- model.matrix(fix.help, data)[,-1]
  }else{
    X <- model.matrix(fix, data)  
  }
 
  K <- NULL
  
  if(!is.null(family$multivariate)){
    y.fac <- as.factor(y)
    K <- length(levels(y.fac))-1
    if(family$family=="acat"){
      y <- c(t(model.matrix(~0+y.fac, contrasts = list(y.fac = "contr.treatment"))[,-length(levels(y.fac))]))
    }
    if(family$family=="cumulative"){
      get.resp <- function(x){as.numeric(as.numeric(x) <= 1:K)}
      y <- c(sapply(y.fac,get.resp))
    }
  }
  
  if(!is.null(family$multivariate)){
    if(all(X[,1]==1)){
      X <- X[,-1]
    }
    names.x <- colnames(X)
    theta <- matrix(rep(diag(1,K),nrow(X)),ncol=K,byrow=TRUE)
    X <- cbind(theta, matrix(rep(X,each=K),ncol=ncol(X)))
    colnames(X) <- c(paste0("theta",1:K),names.x)
    if(orig.names[1]=="(Intercept)")
    {  
      orig.names <- orig.names[-1]
      index.new <- c(rep(NA,K),index.new[-1]) 
    }else{
      index.new <- c(rep(NA,K),index.new) 
    }  
    orig.names <- c(paste0("theta",1:K),orig.names)
  }
  
  transf.names <- colnames(X)
  center <- control$center
  standardize <- control$standardize
  ####### Center & Standardization
  
  ## Which are the non-penalized parameters?
  any.notpen    <- any(is.na(index.new))
  inotpen.which <- which(is.na(index.new))
  nrnotpen      <- length(inotpen.which)
  
  intercept.which <- which(apply(X == 1, 2, all))
  has.intercept   <- length(intercept.which)
  
  ## Index vector of the penalized parameter groups
  if(any.notpen){
    ipen <- index.new[-inotpen.which]
    ipen.which <- split((1:ncol(X))[-inotpen.which], ipen)
  }else{
    if(has.intercept)
      warning("All groups are penalized, including the intercept.")
    ipen <- index.new
    ipen.which <- split((1:ncol(X)), ipen)
  }
  
  if(center){
    if(!has.intercept & is.null(family$multivariate)) ## could be removed; already handled above
      stop("Need intercept term when using center = TRUE")
    
    mu.x                 <- apply(X[,-intercept.which], 2, mean)
    X[,-intercept.which] <- sweep(X[,-intercept.which], 2, mu.x)
  }
  
  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
    ##warning("...Using standardized design matrix.\n")
    stand        <- blockstand(X, ipen.which, inotpen.which)
    X            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  
  if(ncol(X)==1)
  {
    if(colnames(X)=="(Intercept)")
      stop("No terms to select! Use glmer, glmmPQL or glmmML!")
  }

  old.names<-attr(X,"dimnames")[[2]]

  if(control$print.iter)
    message("Iteration 1")

  very.old.names<-very.old.names[!is.na(control$index)]
  
  block<-as.numeric(table(index.new[!is.na(index.new)]))
  
  BLOCK<-FALSE
  if(!all(block==1))
    BLOCK<-TRUE
  
  lin<-ncol(X)
  
  if(is.null(control$start))
    control$start<-c(rep(0,lin))
  
  if(family$family=="cumulative" & all(control$start==0))
    control$start[1:K] <- ((1:K)-mean(1:K))/K
    
  N<-length(y)
  
  beta_null<-control$start[1:lin]
  if(is.null(attr(beta_null,"names")))
    attr(beta_null,"names")<-orig.names
  beta_null<-beta_null[colnames(X)]
  
  Z_fastalles<-X
  
  if(!control$overdispersion && family$family=="gaussian")
    control$overdispersion<-T
 
  phi <- 1
  #######################################################################  
  ######################## allow switch to Newton Raphson ###############
  #######################################################################  
  if(switch.NR)
  {
    #######################################################################  
    ###########################  1. No Smooth #############################  
    #######################################################################  
    if(is.null(control$smooth))
    {  
      if(lin>1)
      {
        Eta_start<-X%*%beta_null
      }else{
        Eta_start<-rep(beta_null,N)
      }
      
      if(is.null(family$multivariate)){
        D<-family$mu.eta(Eta_start)
        Mu<-family$linkinv(Eta_start)
        SigmaInv <- 1/family$variance(Mu)
      }else{
        Eta_cat <- matrix(Eta_start, byrow = TRUE, ncol = K)
        Mu_cat <- family$linkinv(Eta_cat)
        D <- family$deriv.mat(Mu_cat)
        SigmaInv <- family$SigmaInv(Mu_cat)
        Mu <- c(t(Mu_cat))
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      q<-ncol(X)
      
      Z_alles<-cbind(X,U)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,lin)
      Delta[1,1:lin]<-beta_null[final.names]
      active_old<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      logLik.vec<-c()
      
      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      l=1

      if(is.null(family$multivariate)){
      score_vec<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)
      }else{
      score_vec<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))
      }
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1)
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:lin],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
          Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*score_vec
          if(t_opt>t_edge)
            Delta[1,crit.obj$whichmin+q]<-0  
          Eta<-Z_alles%*%Delta[1,]

          if(is.null(family$multivariate)){
            D<-family$mu.eta(Eta)
            Mu<-family$linkinv(Eta)
            SigmaInv <- 1/family$variance(Mu)
          }else{
            Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
            Mu_cat <- family$linkinv(Eta_cat)
            D <- family$deriv.mat(Mu_cat)
            SigmaInv <- family$SigmaInv(Mu_cat)
            Mu <- c(t(Mu_cat))
          }
          
          logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F,K=K)
          
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))

          betaact<-Delta[1,active]
          
      NRstep<-F; vorz <- T
        
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)
        
          if(is.null(family$multivariate)){
            score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)
          }else{
            score_old2<-score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))
          }
          
          
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda)
          }
          score_vec2<-c(score_vec2[1:q],grad.1)
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin])
          t_edge<-crit.obj$min.rate
          
          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l-1,1:lin],
                                             Grad=score_vec2,family=family,lower = 0, upper = Inf, K=K))
          
          t_opt<-optim.obj$par
          
          if(!NRstep && !(all(active_old==active)))
             NRstep <- T
            
          tryNR <- (t_opt<t_edge)  && NRstep  
          
          if(tryNR) 
          {
            lin_akt<-q+sum(!is.element(Delta[l-1,(q+1):lin],0))
            
            if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)
            }else{
            F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_aktuell)))
            }
            
            InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
            if(class(InvFisher)=="try-error")
              InvFisher<-try(solve(F_gross),silent=T)
            
            if(class(InvFisher)=="try-error")
            { 
              tryNR <- FALSE
            }else{
              Delta.test<-Delta[l-1,active]+nue*InvFisher%*%score_old2[active]
              if(lin_akt>q)
              {
                vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
              }else{
                vorz<-T  
              }
            }
          }
          
        score_old<-score_old2
        score_vec<-score_vec2
        
        if(!vorz)
          tryNR<-F
        
        if(tryNR)
        {
          Delta[l,active]<-Delta[l-1,active]+nue*InvFisher%*%score_old[active]
          NRstep<-T
          #       print("NR-step!!!")
        }else{
          Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*score_vec
          if(t_opt>t_edge)
            Delta[l,crit.obj$whichmin+q]<-0  
          NRstep<-F
        }
        Eta<-Z_alles%*%Delta[l,]
        
        if(is.null(family$multivariate)){
          D<-family$mu.eta(Eta)
          Mu<-family$linkinv(Eta)
          SigmaInv <- 1/family$variance(Mu)
        }else{
          Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          D <- family$deriv.mat(Mu_cat)
          SigmaInv <- family$SigmaInv(Mu_cat)
          Mu <- c(t(Mu_cat))
        }
        
        logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F,K=K)
        
              active_old<-active
              active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0))
              Z_aktuell<-Z_alles[,active]
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
             
          betaact<-Delta[l,active]
              
          Eta.ma[l+1,]<-Eta
              
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
        }}
      
      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }
      
      Delta_neu<-Delta[l,]
      Mu_opt<-Mu

      if(!final.re)
      {
        if(is.null(family$multivariate)){
          D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
          W_opt <- D*SigmaInv*D
          F_gross <- t(Z_aktuell)%*%(Z_aktuell*W_opt)
        }else{
          W_opt <- D%*%(SigmaInv%*%t(D))
          F_gross <- t(Z_aktuell)%*%(W_opt%*%Z_aktuell)
          W_inv_t <- chol(W_opt)
        }
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)!="try-error")
        {  
          if(is.null(family$multivariate)){
            FinalHat.df<-(Z_aktuell*sqrt(W_opt))%*%InvFisher2%*%t(Z_aktuell*sqrt(W_opt))
          }else{
            FinalHat.df<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
          }
          
        df<-sum(diag(FinalHat.df))
        if(control$overdispersion)
          phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-sum(diag(FinalHat.df)))
        
      }else{
        df <- NA
      }}
      
      aaa<-!is.element(Delta_neu[1:lin],0)
      
      if(final.re)
      {    
        ############ final re-estimation
        
          glmm_fin<-try(glmm_final_noRE(y,Z_fastalles[,aaa],K=K,
                                   Delta_start=Delta_neu[aaa],steps=control$maxIter,
                                   family=family,overdispersion=control$overdispersion,
                                   phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
          
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_noRE(y,Z_fastalles[,aaa],K=K,
                                      Delta_start=Delta_start[aaa],steps=control$maxIter,
                                      family=family,overdispersion=control$overdispersion,
                                      phi=control$phi,print.iter.final=control$print.iter.final,
                                      eps.final=control$eps.final),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }

        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {
          cat("Warning:\n")
          cat("Final Fisher scoring reestimation did not converge!\n")
        }
        #######
      }else{
        glmm_fin<-NA  
        class(glmm_fin)<-"try-error"
      }  
      Standard_errors<-matrix(NA,length(Delta_neu),length(Delta_neu))

      if(class(glmm_fin)!="try-error")
      {
        Delta_neu2<-Delta_neu
        Delta_neu2[aaa]<-glmm_fin$Delta
        Standard_errors[aaa,aaa]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
        Delta_neu<-Delta_neu2
        Eta_opt<-Z_alles%*%Delta_neu
        if(is.null(family$multivariate)){
          Mu_opt<-family$linkinv(Eta_opt)
        }else{
          Eta_cat <- matrix(Eta_opt, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          Mu_opt <- c(t(Mu_cat))
        }
      }else{
        glmm_fin<-list()
        complexity<-df
      }
      
      colnames(Standard_errors) <- rownames(Standard_errors) <- paste0("help.",1:nrow(Standard_errors))
      
      if(dim(X)[2]>0)
      {  
        names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
        colnames(Standard_errors)[1:dim(X)[2]]<-rownames(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      }
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        colnames(Standard_errors)[(dim(X)[2]+1):lin]<-rownames(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      colnames(Delta)<-final.names
      Delta_neu[1:lin] <- Delta_neu[1:lin][transf.names]
      Standard_errors[1:lin,1:lin] <-Standard_errors[1:lin,1:lin][transf.names,transf.names]
      names(Delta_neu)[1:lin] <- transf.names
      colnames(Standard_errors)[1:lin] <- rownames(Standard_errors)[1:lin] <- transf.names
      ## Transform the coefficients back to the original scale if the design
      ## matrix was standardized
      if(standardize){
        if(any.notpen)
        {
          Delta_neu[inotpen.which] <- (1 / scale.notpen) * Delta_neu[inotpen.which]
          if(class(glmm_fin)!="try-error")
            Standard_errors[inotpen.which] <- (1 / scale.notpen) * Standard_errors[inotpen.which]
        }
        ## For df > 1 we have to use a matrix inversion to go back to the
        ## original scale
        for(j in 1:length(ipen.which)){
          ind <- ipen.which[[j]]
          Delta_neu[ind] <- solve(scale.pen[[j]], Delta_neu[ind,drop = FALSE])
          if(class(glmm_fin)!="try-error")
          {
            Sc.help <- solve(scale.pen[[j]])
            Standard_errors[ind,ind] <- t(Sc.help)%*%(Standard_errors[ind,ind]%*%Sc.help)
          }
        }
      }
      
      Standard_errors <- sqrt(diag(Standard_errors))
      ## Need to adjust intercept if we have performed centering
      if(center){
        Delta_neu[intercept.which] <- Delta_neu[intercept.which] -
          sum(Delta_neu[1:lin][-intercept.which,drop = FALSE] * mu.x)   
      }
      
      
      aic<-NaN
      bic<-NaN
      
      if(is.element(family$family,c("gaussian", "binomial", "poisson","acat","cumulative"))) 
      {
        
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F, K=K)
        
        
        if(control$complexity!="hat.matrix")  
           complexity<-sum(Delta_neu[1:(lin)]!=0)
 
        aic<--2*loglik+2*complexity
        bic<--2*loglik+log(N)*complexity
      }else{
        warning("For the specified family (so far) no AIC and BIC are available!")  
      }

      ret.obj=list()
      ret.obj$aic<-aic
      ret.obj$bic<-bic
      ret.obj$Deltamatrix<-Delta
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$conv.step<-conv.step
      ret.obj$data<-data
      ret.obj$phi.med<-phi.med
      ret.obj$y <- y
      ret.obj$X <- cbind(X,U)
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
      return(ret.obj)
      ##############################################################  
      ######################## 2. Smooth ###########################  
      ##############################################################  
    }else{
      smooth<-control$smooth
      
      if(attr(terms(smooth$formula), "intercept")==0)
      {
        variables<-attr(terms(smooth$formula),"term.labels")
        smooth$formula<- "~ +1"
        for (ir in 1:length(variables))
          smooth$formula <- paste(smooth$formula, variables[ir], sep="+")
        smooth$formula <-formula(smooth$formula)
      }
      
      B <- model.matrix(smooth$formula, data)
      B.names<-attr(B,"dimnames")[[2]]
      
      B<-as.matrix(B[,-1])
      attr(B,"dimnames")[[2]]<-B.names[-1]
      
      nbasis<-smooth$nbasis
      diff.ord<-smooth$diff.ord
      spline.degree<-smooth$spline.degree
      knots.no<-nbasis-1
      if(spline.degree<3 && (spline.degree-diff.ord)<2)
        knots.no<-knots.no+1  
      penal<-smooth$penal
      
      if(!(diff.ord<spline.degree))
        stop("Order of differences must be lower than degree of B-spline polynomials!")
      
      m<-dim(B)[2]
      
      Phi<-numeric()
      
      for (r in 1:m)
      {
        Basis<-bs.design(B[,r],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
        Phi_temp<-cbind(Basis$X[,-1],Basis$Z)
        colnames(Phi_temp)<-paste(colnames(B)[r],rep(1:dim(Phi_temp)[2],each=1), sep=".")
        Phi<-cbind(Phi,Phi_temp)
      }
      
      dim.smooth<-dim(Phi)[2]
      
      if(is.null(smooth$start))
        smooth$start<-rep(0,dim.smooth)  
      
      smooth_null<-smooth$start
      
      if(lin>1)
      {
        Eta_start<-X%*%beta_null+Phi%*%smooth_null
      }else{
        Eta_start<-rep(beta_null,N)+Phi%*%smooth_null
      }
      
      if(is.null(family$multivariate)){
        D<-family$mu.eta(Eta_start)
        Mu<-family$linkinv(Eta_start)
        SigmaInv <- 1/family$variance(Mu)
      }else{
        Eta_cat <- matrix(Eta_start, byrow = TRUE, ncol = K)
        Mu_cat <- family$linkinv(Eta_cat)
        D <- family$deriv.mat(Mu_cat)
        SigmaInv <- family$SigmaInv(Mu_cat)
        Mu <- c(t(Mu_cat))
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,Phi)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,lin+dim.smooth)
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null
      active_old<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth))
      
      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      logLik.vec<-c()

      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      l=1
      
      if(diff.ord>1)
      {
        k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
      }else{
        k2<-rep(1,nbasis)
      }
      k22<-rep(k2,m)
      penal.vec<-penal*k22
      
          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)

          if(is.null(family$multivariate)){
            score_vec<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
          }else{
            score_vec<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
          }
          
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+dim.smooth)])
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:(lin+dim.smooth)],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf, K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
          Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*score_vec
          if(t_opt>t_edge)
            Delta[1,crit.obj$whichmin+q]<-0  
          Eta<-Z_alles%*%Delta[1,]

          if(is.null(family$multivariate)){
            D<-family$mu.eta(Eta)
            Mu<-family$linkinv(Eta)
            SigmaInv <- 1/family$variance(Mu)
          }else{
            Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
            Mu_cat <- family$linkinv(Eta_cat)
            D <- family$deriv.mat(Mu_cat)
            SigmaInv <- family$SigmaInv(Mu_cat)
            Mu <- c(t(Mu_cat))
          }
          
          logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F, K=K)
          
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
        betaact<-Delta[1,active]
        
        NRstep<-F; vorz <- T
        
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {

          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          

          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)
          
          if(is.null(family$multivariate)){
            score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l-1,]
          }else{
            score_old2<-score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l-1,]
          }
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda)
          }
          score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth)])
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin])
          t_edge<-crit.obj$min.rate
          
          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l-1,1:(lin+dim.smooth)],
                                             Grad=score_vec2,family=family,lower = 0, upper = Inf, K=K))
          
          t_opt<-optim.obj$par
          
          if(!NRstep && !(all(active_old==active)))
            NRstep <- T
          
          tryNR <- (t_opt<t_edge)  && NRstep  
          
          if(tryNR) 
          {
            lin_akt<-q+sum(!is.element(Delta[l-1,(q+1):lin],0))
            
            P_akt<-c(rep(0,lin_akt),penal.vec)
            if(is.null(family$multivariate)){
              D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+diag(P_akt)
            }else{
              F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_aktuell)))+diag(P_akt)
            }
            
            InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
            if(class(InvFisher)=="try-error")
              InvFisher<-try(solve(F_gross),silent=T)
            if(class(InvFisher)=="try-error")
            {
              tryNR <- FALSE
            }else{
              Delta.test<-Delta[l-1,active]+nue*InvFisher%*%score_old2[active]
              if(lin_akt>q)
              {
                vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
              }else{
                vorz<-T  
              }
            }
        }  

        score_old<-score_old2
        score_vec<-score_vec2
        
       if(!vorz)
        tryNR<-F
          
              if(tryNR)
              {
                Delta[l,active]<- Delta[l-1,active]+nue*InvFisher%*%score_old[active]
                NRstep<-T
#                print("NR-step!!!")
              }else{
                Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*score_vec
                if(t_opt>t_edge)
                  Delta[l,crit.obj$whichmin+q]<-0  
                NRstep<-F
              }
              Eta<-Z_alles%*%Delta[l,]

              if(is.null(family$multivariate)){
                D<-family$mu.eta(Eta)
                Mu<-family$linkinv(Eta)
                SigmaInv <- 1/family$variance(Mu)
              }else{
                Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
                Mu_cat <- family$linkinv(Eta_cat)
                D <- family$deriv.mat(Mu_cat)
                SigmaInv <- family$SigmaInv(Mu_cat)
                Mu <- c(t(Mu_cat))
              }
              
              logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F,K=K)
              
              active_old<-active
              active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth))
              Z_aktuell<-Z_alles[,active]
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              
            betaact<-Delta[l,active]

            Eta.ma[l+1,]<-Eta
            
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
      
      ######## Final calculation
      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }
      
      Delta_neu<-Delta[l,]
      Mu_opt<-Mu

      if(!final.re)
      {
      P_akt<-c(rep(0,lin_akt),penal.vec)
      if(is.null(family$multivariate)){
        D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
        W_opt <- D*SigmaInv*D
        F_gross <- t(Z_aktuell)%*%(Z_aktuell*W_opt)+diag(P_akt)
      }else{
        W_opt <- D%*%(SigmaInv%*%t(D))
        F_gross <- t(Z_aktuell)%*%(W_opt%*%Z_aktuell)+diag(P_akt)
        W_inv_t <- chol(W_opt)
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
      if(class(InvFisher2)!="try-error")
      {  
        if(is.null(family$multivariate)){
          FinalHat.df<-(Z_aktuell*sqrt(W_opt))%*%InvFisher2%*%t(Z_aktuell*sqrt(W_opt))
        }else{
          FinalHat.df<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
        }
        df<-sum(diag(FinalHat.df))
        if(control$overdispersion)
          phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-sum(diag(FinalHat.df)))
      }else{
        df <- NA
      }}
      
      aaa<-!is.element(Delta_neu[1:(lin)],0)
      if(final.re)
      {    
        ############ final re-estimation
          glmm_fin<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,K=K,
                                          Delta_start=Delta_neu[aaa],steps=control$maxIter,
                                          family=family,overdispersion=control$overdispersion,
                                          phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,K=K,
                                             Delta_start=Delta_start[aaa],steps=control$maxIter,
                                             family=family,overdispersion=control$overdispersion,
                                             phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }

        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {
          cat("Warning:\n")
          cat("Final Fisher scoring reestimation did not converge!\n")
        }
        
        #######
      }else{
        glmm_fin<-NA  
        class(glmm_fin)<-"try-error"
      }  
      Standard_errors<-matrix(NA,length(Delta_neu),length(Delta_neu))
      
      if(class(glmm_fin)!="try-error")
      {
        Delta_neu2<-Delta_neu
        EDF.matrix<-glmm_fin$EDF.matrix
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])  
        Delta_neu2[c(aaa,rep(T,dim.smooth))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,dim.smooth)),c(aaa,rep(T,dim.smooth))]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
        Delta_neu<-Delta_neu2
        Eta_opt<-Z_alles%*%Delta_neu
        if(is.null(family$multivariate)){
          Mu_opt<-family$linkinv(Eta_opt)
        }else{
          Eta_cat <- matrix(Eta_opt, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          Mu_opt <- c(t(Mu_cat))
        }
      }else{
        glmm_fin<-list()
        complexity<-df
        
        lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))

            P1<-c(rep(0,lin_akt),penal.vec)
            if(is.null(family$multivariate)){
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+diag(P1)
            }else{
              Fpart <- t(Z_aktuell)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_aktuell)))
              F_gross<-Fpart+diag(P1)
            }
            
        InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher)=="try-error")
          InvFisher<-try(solve(F_gross),silent=T)
        if(class(InvFisher)=="try-error")
        {
          warning("No EDF's for smooth functions available, as Fisher matrix not invertible!")
          complexity.smooth<-dim.smooth
        }else{  
          ###### EDF of spline; compare Wood's Book on page 167
          if(is.null(family$multivariate)){
            EDF.matrix<-InvFisher%*%(t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D))
          }else{
            EDF.matrix<-InvFisher%*%Fpart
          }
          complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])
        }
      }  
      
      if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
        complexity.smooth<-dim.smooth

      colnames(Standard_errors) <- rownames(Standard_errors) <- paste0("help.",1:nrow(Standard_errors))
      
      if(dim(X)[2]>0)
      {  
        names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
        colnames(Standard_errors)[1:dim(X)[2]]<-rownames(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      }
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        colnames(Standard_errors)[(dim(X)[2]+1):lin]<-rownames(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-rownames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Delta)<-c(old.names,colnames(Phi))
      
      Delta_neu[1:lin] <- Delta_neu[1:lin][transf.names]
      Standard_errors[1:lin,1:lin] <-Standard_errors[1:lin,1:lin][transf.names,transf.names]
      names(Delta_neu)[1:lin] <- transf.names
      colnames(Standard_errors)[1:lin] <- rownames(Standard_errors)[1:lin] <- transf.names
      ## Transform the coefficients back to the original scale if the design
      ## matrix was standardized
      if(standardize){
        if(any.notpen)
        {
          Delta_neu[inotpen.which] <- (1 / scale.notpen) * Delta_neu[inotpen.which]
          if(class(glmm_fin)!="try-error")
            Standard_errors[inotpen.which] <- (1 / scale.notpen) * Standard_errors[inotpen.which]
        }
        ## For df > 1 we have to use a matrix inversion to go back to the
        ## original scale
        for(j in 1:length(ipen.which)){
          ind <- ipen.which[[j]]
          Delta_neu[ind] <- solve(scale.pen[[j]], Delta_neu[ind,drop = FALSE])
          if(class(glmm_fin)!="try-error")
          {
            Sc.help <- solve(scale.pen[[j]])
            Standard_errors[ind,ind] <- t(Sc.help)%*%(Standard_errors[ind,ind]%*%Sc.help)
          }
        }
      }
      
      Standard_errors <- sqrt(diag(Standard_errors))
      ## Need to adjust intercept if we have performed centering
      if(center){
        Delta_neu[intercept.which] <- Delta_neu[intercept.which] -
          sum(Delta_neu[1:lin][-intercept.which,drop = FALSE] * mu.x)   
      }
      
      aic<-NaN
      bic<-NaN
      
      if(is.element(family$family,c("gaussian", "binomial", "poisson","acat","cumulative"))) 
      {
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F,K=K)
        
        if(control$complexity!="hat.matrix")  
        {  
          complexity<-sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
        }
        
        aic<--2*loglik+2*complexity
        bic<--2*loglik+log(N)*complexity
      }else{
        warning("For the specified family (so far) no AIC and BIC are available!")  
      }
      
      ret.obj=list()
      ret.obj$aic<-aic
      ret.obj$bic<-bic
      ret.obj$Deltamatrix<-Delta
      ret.obj$smooth<-Delta_neu[(lin+1):(lin+dim.smooth)]
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$data<-data
      ret.obj$B<-B
      ret.obj$nbasis<-nbasis
      ret.obj$spline.degree<-spline.degree
      ret.obj$diff.ord<-diff.ord
      ret.obj$knots.no<-knots.no
      ret.obj$conv.step<-conv.step
      ret.obj$phi.med<-phi.med
      ret.obj$complexity.smooth<-complexity.smooth
      ret.obj$y <- y
      ret.obj$X <- cbind(X,U)
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
      return(ret.obj)
    }  
    
    #######################################################################  
    ######################## no switch to Newton Raphson ###############
    #######################################################################  
  }else{
    #######################################################################  
    ###########################  1. No Smooth #############################  
    #######################################################################  
    if(is.null(control$smooth))
    {  
       if(lin>1)
      {
        Eta_start<-X%*%beta_null
      }else{
        Eta_start<-rep(beta_null,N)
      }
      
      if(is.null(family$multivariate)){
        D<-family$mu.eta(Eta_start)
        Mu<-family$linkinv(Eta_start)
        SigmaInv <- 1/family$variance(Mu)
      }else{
        Eta_cat <- matrix(Eta_start, byrow = TRUE, ncol = K)
        Mu_cat <- family$linkinv(Eta_cat)
        D <- family$deriv.mat(Mu_cat)
        SigmaInv <- family$SigmaInv(Mu_cat)
        Mu <- c(t(Mu_cat))
      }
      
      lin0<-sum(beta_null!=0)

      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U)
      
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,lin)
      Delta[1,1:lin]<-beta_null[final.names]

      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      logLik.vec<-c()

      control$epsilon<-control$epsilon*sqrt(ncol(Delta))
      
      Delta_start<-Delta[1,]
      
      l=1
      
      if(is.null(family$multivariate)){
        score_vec<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)
      }else{
        score_vec<-t(Z_alles)%*%D%*%SigmaInv%*%(y-Mu)
      }
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1)
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate

      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:lin],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf, K = K))
      
      t_opt<-optim.obj$par
      
       nue<-control$nue
      
        Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*score_vec
        if(t_opt>t_edge)# & half.index==0)
          Delta[1,crit.obj$whichmin+q]<-0  
        Eta<-Z_alles%*%Delta[1,]

        Eta.ma[2,]<-Eta
       
        if(is.null(family$multivariate)){
          D<-family$mu.eta(Eta)
          Mu<-family$linkinv(Eta)
          SigmaInv <- 1/family$variance(Mu)
        }else{
          Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          D <- family$deriv.mat(Mu_cat)
          SigmaInv <- family$SigmaInv(Mu_cat)
          Mu <- c(t(Mu_cat))
        }

        logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F, K=K)
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          
          if(control$print.iter)
            message("Iteration ",l)

          if(is.null(family$multivariate)){
            score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)
          }else{
            score_old2<-score_vec2<-t(Z_alles)%*%D%*%SigmaInv%*%(y-Mu)
          }
          
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda)
          }
          score_vec2<-c(score_vec2[1:q],grad.1)
          
          #F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin])
          t_edge<-crit.obj$min.rate
          #grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
          
          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l-1,1:lin],
                                             Grad=score_vec2,family=family,lower = 0, upper = Inf, K=K))
          
          t_opt<-optim.obj$par
          

          Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*score_vec2
            if(t_opt>t_edge)# & half.index==0)
              Delta[l,crit.obj$whichmin+q]<-0  
            
            Eta<-Z_alles%*%Delta[l,]


            if(is.null(family$multivariate)){
              D<-family$mu.eta(Eta)
              Mu<-family$linkinv(Eta)
              SigmaInv <- 1/family$variance(Mu)
            }else{
              Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
              Mu_cat <- family$linkinv(Eta_cat)
              D <- family$deriv.mat(Mu_cat)
              SigmaInv <- family$SigmaInv(Mu_cat)
              Mu <- c(t(Mu_cat))
            }

            logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F,K=K)
            
#            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
         Eta.ma[l+1,]<-Eta
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
        }}
      
      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }
      
      
      Delta_neu<-Delta[l,]
      Eta_opt<-Eta
      Mu_opt<-Mu
      
      SigmaInv_opt<-SigmaInv
      D_opt<-D

      if(!final.re)
      {
        if(is.null(family$multivariate)){
        D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
        W_opt <- D*SigmaInv*D
        F_gross <- t(Z_aktuell)%*%(Z_aktuell*W_opt)
        }else{
          W_opt <- D%*%(SigmaInv%*%t(D))
          F_gross <- t(Z_aktuell)%*%(W_opt%*%Z_aktuell)
          W_inv_t <- chol(W_opt)
        }
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)!="try-error")
        {  
          if(is.null(family$multivariate)){
            FinalHat.df<-(Z_aktuell*sqrt(W_opt))%*%InvFisher2%*%t(Z_aktuell*sqrt(W_opt))
          }else{
            FinalHat.df<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
          }
          
          df<-sum(diag(FinalHat.df))
          if(control$overdispersion)
            phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-sum(diag(FinalHat.df)))
          
        }else{
          df <- NA
        }}
      
      if(final.re)
      {    
        ############ final re-estimation
        aaa<-!is.element(Delta_neu[1:(lin)],0)
        

        glmm_fin<-try(glmm_final_noRE(y,Z_fastalles[,aaa],K=K,
                                      Delta_start=Delta_neu[aaa],steps=control$maxIter,
                                      family=family,overdispersion=control$overdispersion,
                                      phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = F)
 

        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_noRE(y,Z_fastalles[,aaa],K=K,
                                         Delta_start=Delta_start[aaa],steps=control$maxIter,
                                         family=family,overdispersion=control$overdispersion,
                                         phi=control$phi,print.iter.final=control$print.iter.final,
                                         eps.final=control$eps.final),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }
 
                  
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {
          cat("Warning:\n")
          cat("Final Fisher scoring reestimation did not converge!\n")
        }
        #######
      }else{
        glmm_fin<-NA  
        class(glmm_fin)<-"try-error"
      }  
      Standard_errors<-matrix(NA,length(Delta_neu),length(Delta_neu))
      
      if(class(glmm_fin)!="try-error")
      {
        Delta_neu2<-Delta_neu
        Delta_neu2[aaa]<-glmm_fin$Delta
        Standard_errors[aaa,aaa]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
        Delta_neu<-Delta_neu2
        Eta_opt<-Z_alles%*%Delta_neu
          if(is.null(family$multivariate)){
            Mu_opt<-family$linkinv(Eta_opt)
          }else{
            Eta_cat <- matrix(Eta_opt, byrow = TRUE, ncol = K)
            Mu_cat <- family$linkinv(Eta_cat)
            Mu_opt <- c(t(Mu_cat))
          }
        }else{
        glmm_fin<-list()
        complexity<-df
      }

      colnames(Standard_errors) <- rownames(Standard_errors) <- paste0("help.",1:nrow(Standard_errors))
      
      if(dim(X)[2]>0)
      {  
        names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
        colnames(Standard_errors)[1:dim(X)[2]]<-rownames(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      }
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        colnames(Standard_errors)[(dim(X)[2]+1):lin]<-rownames(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      colnames(Delta)<-final.names
      Delta_neu[1:lin] <- Delta_neu[1:lin][transf.names]
      Standard_errors[1:lin,1:lin] <-Standard_errors[1:lin,1:lin][transf.names,transf.names]
      names(Delta_neu)[1:lin] <- transf.names
      colnames(Standard_errors)[1:lin] <- rownames(Standard_errors)[1:lin] <- transf.names
      ## Transform the coefficients back to the original scale if the design
      ## matrix was standardized
      if(standardize){
        if(any.notpen)
        {
          Delta_neu[inotpen.which] <- (1 / scale.notpen) * Delta_neu[inotpen.which]
          if(class(glmm_fin)!="try-error")
            Standard_errors[inotpen.which] <- (1 / scale.notpen) * Standard_errors[inotpen.which]
        }
        ## For df > 1 we have to use a matrix inversion to go back to the
        ## original scale
        for(j in 1:length(ipen.which)){
          ind <- ipen.which[[j]]
          Delta_neu[ind] <- solve(scale.pen[[j]], Delta_neu[ind,drop = FALSE])
          if(class(glmm_fin)!="try-error")
          {
            Sc.help <- solve(scale.pen[[j]])
            Standard_errors[ind,ind] <- t(Sc.help)%*%(Standard_errors[ind,ind]%*%Sc.help)
          }
        }
      }
      
      Standard_errors <- sqrt(diag(Standard_errors))
      ## Need to adjust intercept if we have performed centering
      if(center){
        Delta_neu[intercept.which] <- Delta_neu[intercept.which] -
          sum(Delta_neu[1:lin][-intercept.which,drop = FALSE] * mu.x)   
      }
      
      aic<-NaN
      bic<-NaN
      
      if(is.element(family$family,c("gaussian", "binomial", "poisson", "acat" ,"cumulative"))) 
      {
        
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F, K=K)
        

        if(control$complexity!="hat.matrix")  
          complexity<-sum(Delta_neu[1:(lin)]!=0)

        aic<--2*loglik+2*complexity
        bic<--2*loglik+log(N)*complexity
      }else{
        warning("For the specified family (so far) no AIC and BIC are available!")  
      }
      
      ret.obj=list()
      ret.obj$aic<-aic
      ret.obj$bic<-bic
      ret.obj$Deltamatrix<-Delta
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$conv.step<-conv.step
      ret.obj$data<-data
      ret.obj$phi.med<-phi.med
      ret.obj$y <- y
      ret.obj$X <- cbind(X,U)
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
#      ret.obj$score_vec.unpen<- score_vec.unpen
      return(ret.obj)
      ##############################################################  
      ######################## 2. Smooth ###########################  
      ##############################################################  
    }else{
      smooth<-control$smooth
      
      if(attr(terms(smooth$formula), "intercept")==0)
      {
        variables<-attr(terms(smooth$formula),"term.labels")
        smooth$formula<- "~ +1"
        for (ir in 1:length(variables))
          smooth$formula <- paste(smooth$formula, variables[ir], sep="+")
        smooth$formula <-formula(smooth$formula)
      }
      
      B <- model.matrix(smooth$formula, data)
      B.names<-attr(B,"dimnames")[[2]]
      
      B<-as.matrix(B[,-1])
      attr(B,"dimnames")[[2]]<-B.names[-1]
      
      nbasis<-smooth$nbasis
      diff.ord<-smooth$diff.ord
      spline.degree<-smooth$spline.degree
      knots.no<-nbasis-1
      if(spline.degree<3 && (spline.degree-diff.ord)<2)
        knots.no<-knots.no+1  
      penal<-smooth$penal
      
      if(!(diff.ord<spline.degree))
        stop("Order of differences must be lower than degree of B-spline polynomials!")
      
      m<-dim(B)[2]
      
      Phi<-numeric()
      
      for (r in 1:m)
      {
        Basis<-bs.design(B[,r],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
        Phi_temp<-cbind(Basis$X[,-1],Basis$Z)
        colnames(Phi_temp)<-paste(colnames(B)[r],rep(1:dim(Phi_temp)[2],each=1), sep=".")
        Phi<-cbind(Phi,Phi_temp)
      }
      
      dim.smooth<-dim(Phi)[2]
      
      if(is.null(smooth$start))
        smooth$start<-rep(0,dim.smooth)  
      
      smooth_null<-smooth$start
      
      if(lin>1)
      {
        Eta_start<-X%*%beta_null+Phi%*%smooth_null
      }else{
        Eta_start<-rep(beta_null,N)+Phi%*%smooth_null
      }
      
      if(is.null(family$multivariate)){
        D<-family$mu.eta(Eta_start)
        Mu<-family$linkinv(Eta_start)
        SigmaInv <- 1/family$variance(Mu)
      }else{
        Eta_cat <- matrix(Eta_start, byrow = TRUE, ncol = K)
        Mu_cat <- family$linkinv(Eta_cat)
        D <- family$deriv.mat(Mu_cat)
        SigmaInv <- family$SigmaInv(Mu_cat)
        Mu <- c(t(Mu_cat))
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,Phi)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,(lin+dim.smooth))
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null

      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      logLik.vec<-c()
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      l=1
      
      if(diff.ord>1)
      {
        k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
      }else{
        k2<-rep(1,nbasis)
      }
      k22<-rep(k2,m)
      penal.vec<-penal*k22
      
          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)

          if(is.null(family$multivariate)){
            score_vec<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
          }else{
            score_vec<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
          }
          
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+dim.smooth)])
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:(lin+dim.smooth)],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf, K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
        Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*score_vec
        if(t_opt>t_edge)
          Delta[1,crit.obj$whichmin+q]<-0  
        Eta<-Z_alles%*%Delta[1,]

        if(is.null(family$multivariate)){
          D<-family$mu.eta(Eta)
          Mu<-family$linkinv(Eta)
          SigmaInv <- 1/family$variance(Mu)
        }else{
          Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          D <- family$deriv.mat(Mu_cat)
          SigmaInv <- family$SigmaInv(Mu_cat)
          Mu <- c(t(Mu_cat))
        }
        
        logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F,K=K)
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
      betaact<-Delta[1,active]
      
      vorz<-F
      
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)
          
          if(is.null(family$multivariate)){
            score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l-1,]
          }else{
            score_old2<-score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l-1,]
          }
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda)
          }
          score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth)])
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l-1,(q+1):lin])
          t_edge<-crit.obj$min.rate

          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l-1,1:(lin+dim.smooth)],
                                             Grad=score_vec2,family=family,lower = 0, upper = Inf,K=K))
          
          t_opt<-optim.obj$par
          
          score_vec<-score_vec2
          
            Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*score_vec
            if(t_opt>t_edge)
              Delta[l,crit.obj$whichmin+q]<-0
            Eta<-Z_alles%*%Delta[l,]
 
            if(is.null(family$multivariate)){
              D<-family$mu.eta(Eta)
              Mu<-family$linkinv(Eta)
              SigmaInv <- 1/family$variance(Mu)
            }else{
              Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
              Mu_cat <- family$linkinv(Eta_cat)
              D <- family$deriv.mat(Mu_cat)
              SigmaInv <- family$SigmaInv(Mu_cat)
              Mu <- c(t(Mu_cat))
            }
            
            logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F,K=K)
            
#            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
          betaact<-Delta[l,active]
                  
          Eta.ma[l+1,]<-Eta
          
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }
      
      Delta_neu<-Delta[l,]
      Mu_opt<-Mu

      if(!final.re)
      {
      P_akt<-c(rep(0,lin_akt),penal.vec)
      if(is.null(family$multivariate)){
        D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
        W_opt <- D*SigmaInv*D
        F_gross <- t(Z_aktuell)%*%(Z_aktuell*W_opt)+diag(P_akt)
      }else{
        W_opt <- D%*%(SigmaInv%*%t(D))
        F_gross <- t(Z_aktuell)%*%(W_opt%*%Z_aktuell)+diag(P_akt)
        W_inv_t <- chol(W_opt)
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
      if(class(InvFisher2)!="try-error")
      {  
        if(is.null(family$multivariate)){
          FinalHat.df<-(Z_aktuell*sqrt(W_opt))%*%InvFisher2%*%t(Z_aktuell*sqrt(W_opt))
        }else{
          FinalHat.df<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
        }
        df<-sum(diag(FinalHat.df))
        if(control$overdispersion)
          phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-sum(diag(FinalHat.df)))
      }else{
        df <- NA
      }}
      
      aaa<-!is.element(Delta_neu[1:(lin)],0)

      if(final.re)
      {    
        ############ final re-estimation
        
        glmm_fin<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,K=K,
                                             Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth))],steps=control$maxIter,
                                             family=family,overdispersion=control$overdispersion,
                                             phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,K=K,
                                                Delta_start=Delta_start[c(aaa,rep(T,dim.smooth))],steps=control$maxIter,
                                                family=family,overdispersion=control$overdispersion,
                                                phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
          if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }

        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {
          cat("Warning:\n")
          cat("Final Fisher scoring reestimation did not converge!\n")
        }
        
        #######
      }else{
        glmm_fin<-NA  
        class(glmm_fin)<-"try-error"
      }  
      
      Standard_errors<-matrix(NA,length(Delta_neu),length(Delta_neu))
      
      
      if(class(glmm_fin)!="try-error")
      {
        Delta_neu2<-Delta_neu
        EDF.matrix<-glmm_fin$EDF.matrix
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])  
        Delta_neu2[c(aaa,rep(T,dim.smooth))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,dim.smooth)),c(aaa,rep(T,dim.smooth))]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
        Delta_neu<-Delta_neu2
        Eta_opt<-Z_alles%*%Delta_neu
        if(is.null(family$multivariate)){
          Mu_opt<-family$linkinv(Eta_opt)
        }else{
          Eta_cat <- matrix(Eta_opt, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          Mu_opt <- c(t(Mu_cat))
        }
      }else{
        glmm_fin<-list()
        complexity<-df
        
        lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))

            P1<-c(rep(0,lin_akt),penal.vec)
            if(is.null(family$multivariate)){
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+diag(P1)
            }else{
              Fpart <- t(Z_aktuell)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_aktuell)))
              F_gross<-Fpart+diag(P1)
            }
            
        InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher)=="try-error")
          InvFisher<-try(solve(F_gross),silent=T)
        if(class(InvFisher)=="try-error")
        {
          warning("No EDF's for smooth functions available, as Fisher matrix not invertible!")
          complexity.smooth<-dim.smooth
        }else{  
          ###### EDF of spline; compare Wood's Book on page 167
          if(is.null(family$multivariate)){
            EDF.matrix<-InvFisher%*%(t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D))
          }else{
            EDF.matrix<-InvFisher%*%Fpart
          }
          complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])
        }
      }  
      
      if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
        complexity.smooth<-dim.smooth
      
      
      colnames(Standard_errors) <- rownames(Standard_errors) <- paste0("help.",1:nrow(Standard_errors))
      
      if(dim(X)[2]>0)
      {  
        names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
        colnames(Standard_errors)[1:dim(X)[2]]<-rownames(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      }
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        colnames(Standard_errors)[(dim(X)[2]+1):lin]<-rownames(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-rownames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Delta)<-c(old.names,colnames(Phi))
      
      Delta_neu[1:lin] <- Delta_neu[1:lin][transf.names]
      Standard_errors[1:lin,1:lin] <-Standard_errors[1:lin,1:lin][transf.names,transf.names]
      names(Delta_neu)[1:lin] <- transf.names
      colnames(Standard_errors)[1:lin] <- rownames(Standard_errors)[1:lin] <- transf.names
      ## Transform the coefficients back to the original scale if the design
      ## matrix was standardized
      if(standardize){
        if(any.notpen)
        {
          Delta_neu[inotpen.which] <- (1 / scale.notpen) * Delta_neu[inotpen.which]
          if(class(glmm_fin)!="try-error")
            Standard_errors[inotpen.which] <- (1 / scale.notpen) * Standard_errors[inotpen.which]
        }
        ## For df > 1 we have to use a matrix inversion to go back to the
        ## original scale
        for(j in 1:length(ipen.which)){
          ind <- ipen.which[[j]]
          Delta_neu[ind] <- solve(scale.pen[[j]], Delta_neu[ind,drop = FALSE])
          if(class(glmm_fin)!="try-error")
          {
            Sc.help <- solve(scale.pen[[j]])
            Standard_errors[ind,ind] <- t(Sc.help)%*%(Standard_errors[ind,ind]%*%Sc.help)
          }
        }
      }
      
      Standard_errors <- sqrt(diag(Standard_errors))
      ## Need to adjust intercept if we have performed centering
      if(center){
        Delta_neu[intercept.which] <- Delta_neu[intercept.which] -
          sum(Delta_neu[1:lin][-intercept.which,drop = FALSE] * mu.x)   
      }
      
      aic<-NaN
      bic<-NaN
      
      
      if(is.element(family$family,c("gaussian", "binomial", "poisson","acat","cumulative"))) 
      {
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F,K=K)
        
        if(control$complexity!="hat.matrix")  
        {  
           complexity<-sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
        }      
        aic<--2*loglik+2*complexity
        bic<--2*loglik+log(N)*complexity
      }else{
        warning("For the specified family (so far) no AIC and BIC are available!")  
      }
      
      ret.obj=list()
      ret.obj$aic<-aic
      ret.obj$bic<-bic
      ret.obj$Deltamatrix<-Delta
      ret.obj$smooth<-Delta_neu[(lin+1):(lin+dim.smooth)]
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$data<-data
      ret.obj$B<-B
      ret.obj$nbasis<-nbasis
      ret.obj$spline.degree<-spline.degree
      ret.obj$diff.ord<-diff.ord
      ret.obj$knots.no<-knots.no
      ret.obj$conv.step<-conv.step
      ret.obj$phi.med<-phi.med
      ret.obj$complexity.smooth<-complexity.smooth
      ret.obj$y <- y
      ret.obj$X <- cbind(X,U)
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
      return(ret.obj)
    }  
  }
  ##################################################################
  ##################################################################
}