est.glmmLasso.RE<-function(fix,rnd,data,lambda,family,final.re,switch.NR,control)
{  
  
  control<-do.call(glmmLassoControl, control)
  fix.old<-fix
  
  if(!is.null(fix))
  {  
  if(grepl("\\*", fix[3]))
    stop("Usage of '*' not allowed in formula! Please specify the corresponding variables separately.")  
  
  ic.dummy<-attr(terms(fix),"intercept")  
  y <- model.response(model.frame(fix, data))
  very.old.names<-attr(terms(fix),"term.labels")

  if(ic.dummy!=1 && sum(substr(very.old.names,1,9)=="as.factor")>0){
    fix.help <- update(fix,~ .+1)
    orig.names <- colnames(model.matrix(fix.help, data))[-1]
  }else{
    orig.names <- colnames(model.matrix(fix, data))
  }

  if(!is.null(control$index))
  {
    order.vec<-order(control$index)
    very.old.names<-very.old.names[order.vec]
    control$index<-control$index[order.vec]
  }else{
    control$index<-1:length(very.old.names)
  }
  
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
  
  ##############
  
  if(ncol(X)==1)
  {
    if(colnames(X)=="(Intercept)")
      stop("No terms to select! Use glmer, glmmPQL or glmmML!")
  }

  #old.names<-attr(X,"dimnames")[[2]]

  very.old.names<-very.old.names[!is.na(control$index)]
  ############
  ############
  ############
  }else{
    y <- control$y
    X <- control$X
    K <- control$K
    index.new <- control$index.new
    W.index <- control$W.index
    orig.names <- transf.names <- colnames(X)
    standardize <- center <- FALSE
  }
  
  
  if(control$print.iter)
    message("Iteration 1")

  if(is.list(rnd))
  {
    rnd.len<-length(rnd)
    
    if(rnd.len==1)
    {
      rndformula <- as.character(rnd)
      
      trmsrnd <- terms(rnd[[1]])
      
      if(!is.factor(data[,names(rnd)[1]]))
      {
        data[,names(rnd)[1]] <- as.factor(data[,names(rnd)[1]])
        warning("Cluster variable should be specified as a factor variable!")  
      }
      
      newrndfrml <- "~ -1"
      newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")
      
      if(length(attr(trmsrnd, "variables"))>1)
      {
        newrndfrml <- paste(newrndfrml,  
                            paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                              paste(lbl, names(rnd)[1], sep=":")
                            }), collapse=" + "), sep="+") 
      }
      
      
      W_start <- model.matrix(formula(newrndfrml), data)
      
      rnlabels<-terms(formula(newrndfrml))
      random.labels<-attr(rnlabels,"term.labels")
      k<-table(data[,colnames(data)==(names(rnd)[1])])   
      n<-length(k)
      s<-dim(W_start)[2]/n
      
      if(s>1)
      {
        W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
        for (i in 2:n)
          W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
      }else{
        W<-W_start
      }
    }else{     
      rndformula <- list()
      newrndfrml <- list()
      n <- numeric()
      s <- numeric()
      k<-NULL
      random.labels<-list()
      W<-NULL
      
      for (zu in 1:rnd.len)
      {
        rndformula[[zu]] <- as.character(rnd[[zu]])
        
        trmsrnd <- terms(rnd[[zu]])
        
        if(!is.factor(data[,names(rnd)[zu]]))
        {
          data[,names(rnd)[zu]] <- as.factor(data[,names(rnd)[zu]])
          warning("Cluster variable should be specified as a factor variable!")  
        }
        
        newrndfrml[[zu]] <- "~ -1"
        newrndfrml[[zu]] <- paste(newrndfrml[[zu]],  if(attr(trmsrnd, "intercept")) names(rnd)[zu] else "", sep=" + ")
        
        if(length(attr(trmsrnd, "variables"))>1)
        {
          newrndfrml[[zu]] <- paste(newrndfrml[[zu]],  
                                    paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                                      paste(lbl, names(rnd)[zu], sep=":")
                                    }), collapse=" + "), sep="+") }
        W_start <- model.matrix(formula(newrndfrml[[zu]]), data)
        
        
        rnlabels<-terms(formula(newrndfrml[[zu]]))
        random.labels[[zu]]<-attr(rnlabels,"term.labels")
        k1<-table(data[,colnames(data)==(names(rnd)[zu])])   
        n[zu]<-length(k1)
        s[zu]<-dim(W_start)[2]/n[zu]
        
        if(s[zu]>1)
        {
          W2<-W_start[,seq(from=1,to=1+(s[zu]-1)*n[zu],by=n[zu])]
          for (i in 2:n[zu])
            W2<-cbind(W2,W_start[,seq(from=i,to=i+(s[zu]-1)*n[zu],by=n[zu])])
        }else{
          W2<-W_start
        }
        W<-cbind(W,W2)
        k<-c(k,k1)
      }
    }
    subject.names<-names(rnd)

    if(!is.null(family$multivariate))
    {
      names.of.W <- colnames(W)
      W <- matrix(rep(W,each=K),ncol=ncol(W))
      colnames(W) <- names.of.W
    }
    
  }else{
    W<-rnd
    if(!is.null(control$W.index))
    {
      s<-length(unique(control$W.index))
      n<-ncol(W)/s
    }else{
      n<-ncol(W)
      s<-1
    }
    
    newrndfrml<-NULL  
    random.labels<-attr(W,"W.name")   
    k<-NULL
    rnd.len <- 1
    attr(rnd,"names")<-colnames(rnd)
    subject.names<-colnames(rnd)  
  }

  
  block<-as.numeric(table(index.new[!is.na(index.new)]))
  
  BLOCK<-FALSE
  if(!all(block==1))
    BLOCK<-TRUE
  
  lin<-ncol(X)
  
  if(is.null(control$start))
    control$start<-c(rep(0,(lin+n%*%s)))
  
  if(family$family=="cumulative" & all(control$start==0))
    control$start[1:K] <- ((1:K)-mean(1:K))/K
  
  if(is.null(control$q_start))
  {
    control$q_start<-rep(0.1,sum(s))
    if(sum(s)>1)
      control$q_start<-diag(control$q_start)
  }
  
  q_start<-control$q_start
  
  N<-length(y)
  
  beta_null<-control$start[1:lin]
  if(is.null(attr(beta_null,"names")))
    attr(beta_null,"names")<-orig.names
  beta_null<-beta_null[colnames(X)]
  
  ranef_null<-control$start[(lin+1):(lin+n%*%s)]
  
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
        Eta_start<-X%*%beta_null+W%*%ranef_null
      }else{
        Eta_start<-rep(beta_null,N)+W%*%ranef_null
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
      
      if(rnd.len==1)
      {
        lin0<-sum(beta_null!=0)
          Q_start<-q_start
      }else{
        lin0<-sum(beta_null!=0)
        if(all(s==1))
        {
          Q_start<-diag(diag(q_start),sum(s))
        }else{
          Q_start<-matrix(0,sum(s),sum(s))
          Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-ncol(X)

      Z_alles<-cbind(X,U,W)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,(lin+n%*%s))
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+n%*%s)]<-t(ranef_null)
      active_old<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      logLik.vec<-c()
      
      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      Q_inv<-NULL
      Q_inv.old.temp<-NULL
      Q_inv.start<-NULL

      Q<-list()
      Q[[1]]<-Q_start
      
      l=1
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin),rep(1/Q_start,n*s))
          P1<-diag(P1)
        }else{
          Q_inv.start<-chol2inv(chol(Q_start))
          P1<-matrix(0,lin+n%*%s,lin+n%*%s)
          for(j in 1:n)
            P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv.start
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin),rep(diag(Q_start)^(-1),n))
          P1<-diag(P1)
        }else{
          Q_inv.start<-list()
          Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
          P1<-matrix(0,lin+n%*%s,lin+n%*%s)
          for(jf in 1:n[1])
            P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv.start[[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            for(jf in 1:n[zu])
              P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
          }
        }
      }
      
      
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
      
      score_vec[(q+1):lin] <- grad.1
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
#      ranef.logLik<- -0.5*t(Delta_start[(lin+1):(lin+n%*%s)])%*%(P1[(lin+1):(lin+n%*%s),(lin+1):(lin+n%*%s)]%*%Delta_start[(lin+1):(lin+n%*%s)])

      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta_start[1:lin],ranef=Delta_start[(lin+1):(lin+n%*%s)],
                                         Grad=score_vec,family=family,P=diag(P1[(lin+1):(lin+n%*%s)]), 
                                         lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
      half.index<-0
      solve.test<-FALSE
      ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
      while(!solve.test)
      {  
        
        solve.test2<-FALSE  
        while(!solve.test2)
        {  
          
          if(half.index>50)
          {
            stop("Fisher matrix not invertible")
          }
          
          Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
          if(t_opt>t_edge & half.index==0)
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
          
          ranef.logLik<- -0.5*t(Delta[1,(lin+1):(lin+n%*%s)])%*%(P1[(lin+1):(lin+n%*%s),(lin+1):(lin+n%*%s)]%*%Delta[1,(lin+1):(lin+n%*%s)])
          
          logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
          
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
          if (control$method=="EM")
          {  
            if(rnd.len==1)
            {
              if(s==1)
              {
                P_akt<-diag(c(rep(0,lin_akt),rep(1/Q_start,n*s)))
                }else{
                P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                for(jf in 1:n)
                  P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.start
              }
            }else{
              if(all(s==1))
              {
                P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q_start)^(-1),n)))
              }else{
                P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                for(jf in 1:n[1])
                  P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.start[[1]]
                
                for (zu in 2:rnd.len)
                {
                  for(jf in 1:n[zu])
                    P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                          (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
                }
              }
            }
            if(is.null(family$multivariate)){
              D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
            }else{
              F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
            }
            
            ## include here              
            
            InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
            if(class(InvFisher2)=="try-error")
              InvFisher2<-try(solve(F_gross),silent=T)
            if(class(InvFisher2)=="try-error")
            {
              #stop("Fisher matrix not invertible")  
              half.index<-half.index+1  
            }else{
              solve.test2<-TRUE 
            }}else{
              solve.test2<-TRUE
            }
        }
        
        betaact<-Delta[1,active]
        
        if (control$method=="EM")
        {   
          ############################# Q update ################
          if(rnd.len==1)
          {
            Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
            for (i in 2:n)
              Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
            Q1<-1/n*Q1
          }else{
            Q1<-matrix(0,sum(s),sum(s))
            Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
            for (i in 2:n[1])
              Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
            Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
            
            for (zu in 2:rnd.len)
            {
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
              for (i in 2:n[zu])
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
            }
          }
        }else{
          if(is.null(family$multivariate)){
          Eta_tilde<-Eta+(y-Mu)*1/D
          }else{
            Eta_tilde<-Eta+solve(D)%*%(y-Mu)
          }
          Betadach<-Delta[1,1:(lin)]     
          aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
          X_aktuell<-Z_fastalles[,aktuell_vec]
          
          if(rnd.len==1)
          {
            
            if(s==1)
            {
              upp<-min(20,50*Q_start)
              low<-1e-14
              optim.obj<-nlminb(sqrt(Q_start),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
              Q1<-as.matrix(optim.obj$par)^2
            }else{
              q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
              up1<-min(20,50*max(q_start_vec))
              upp<-rep(up1,length(q_start_vec))
              low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
              optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
              Q1<-matrix(0,s,s)
              Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
              Q1<-Q1+t(Q1)
              diag(Q1)<-(optim.obj$par[1:s])
              
              #### Check for positive definitness ########
              for (ttt in 0:100)
              {
                Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                Q_solvetest<-try(solve(Q1))
                if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                  break
              }
            }   
          }else{
            if(all(s==1))
            {
              q_start_vec<-diag(q_start)
              upp<-rep(min(20,50*diag(q_start)),sum(s))
              low<-rep(1e-14,sum(s))
              optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
              Q1<-diag(optim.obj$par)^2
            }else{
              q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
              up1<-min(20,50*max(q_start_vec))
              low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))
              
              for (zu in 2:rnd.len)
              {
                q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                up1<-min(20,50*max(q_start_vec))
                low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
              }
              upp<-rep(up1,length(q_start_vec))
              optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
              optim.vec<-optim.obj$par
              
              
              Q1<-matrix(0,sum(s),sum(s))
              diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
              if(s[1]>1)
                Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
              optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
              
              for (zu in 2:rnd.len)
              {
                diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
                if(s[zu]>1)
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
                optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
              }
              
              #### Check for positive definitness ########
              for (ttt in 0:100)
              {
                Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                Q_solvetest<-try(solve(Q1))
                if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                  break
              }
            }
          }}
        
        Q[[2]]<-Q1
        
        NRstep<-F
        
        if(rnd.len==1)
        {
          if(s==1)
          {
            P1<-c(rep(0,lin),rep(1/Q1,n*s))
            P1<-diag(P1)
          }else{
            Q_inv<-solve(Q1)
            P1<-matrix(0,lin+n%*%s,lin+n%*%s)
            for(j in 1:n)
              P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
          }
        }else{
          if(all(s==1))
          {
            P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
            P1<-diag(P1)
          }else{
            Q_inv<-list()
            Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
            P1<-matrix(0,lin+n%*%s,lin+n%*%s)
            for(jf in 1:n[1])
              P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv[[1]]
            
            for (zu in 2:rnd.len)
            {
              Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
              for(jf in 1:n[zu])
                P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                   (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
            }
          }
        }
        
        if(is.null(family$multivariate)){
          score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
        }else{
          score_old2<- score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
        }  
        
        lambda.max<-max(abs(score_vec2[(q+1):lin]))
        
        if(BLOCK)
        {
          grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
        }else{
          grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
        }
        
        score_vec2[(q+1):lin] <- grad.1
        
        
        crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
        t_edge<-crit.obj$min.rate
        
        optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[1,1:lin],ranef=Delta[1,(lin+1):(lin+n%*%s)],
                                           Grad=score_vec2,family=family,P=diag(P1[(lin+1):(lin+n%*%s)]), 
                                           lower = 0, upper = Inf,K=K))
        
        t_opt<-optim.obj$par
        
        if(!NRstep && !(all(active_old==active)))
           NRstep <- T
           
           tryNR <- (t_opt<t_edge)  && NRstep  
        
        vorz <- T
        
        if(tryNR) 
        {
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          if(rnd.len==1)
          {
            if(s==1)
            {
              P_akt<-diag(c(rep(0,lin_akt),rep(1/Q1,n*s)))
            }else{
              P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
              for(jf in 1:n)
                P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv
            }
          }else{
            if(all(s==1))
            {
              P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q1)^(-1),n)))
            }else{
              P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
              for(jf in 1:n[1])
                P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv[[1]]
              
              for (zu in 2:rnd.len)
              {
                for(jf in 1:n[zu])
                  P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                        (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              }
            }
          }
          if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
          }else{
            F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
          }

          InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher)=="try-error")
            InvFisher<-try(solve(F_gross),silent=T)
          if(class(InvFisher)=="try-error")
          {
            half.index<-half.index+1  
          }else{
            solve.test<-TRUE 
            Delta.test<-Delta[1,active]+nue*InvFisher%*%score_old2[active]
            if(lin_akt>q)
            {
              vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
            }else{
              vorz<-T  
            }
          }
        }else{
          solve.test<-TRUE  
        }
      }   
      
      Eta.ma[2,]<-Eta
      
      score_old<-score_old2
      score_vec<-score_vec2
      Q1.old<-Q1
      Q_inv.old<-Q_inv    
      
      Q1.very.old<-Q_start
      Q_inv.very.old<-Q_inv.start
      
      active_old<-active
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)

          
          if(!vorz)
            tryNR<-F

          half.index<-0
          solve.test<-FALSE
          ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
          while(!solve.test)
          {  
            
            solve.test2<-FALSE  
            while(!solve.test2)
            {  
              if(half.index>50)
              {
                half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old
              }
              
              if(tryNR)
              {
                Delta[l,active_old]<-Delta[l-1,active_old]+nue*(0.5^half.index)*InvFisher%*%score_old[active_old]
                NRstep<-T
#                print("NR")
              }else{
                Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*nue*score_vec
                if(t_opt>t_edge & half.index==0)
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
              
              ranef.logLik<- -0.5*t(Delta[l,(lin+1):(lin+n%*%s)])%*%(P1[(lin+1):(lin+n%*%s),(lin+1):(lin+n%*%s)]%*%Delta[l,(lin+1):(lin+n%*%s)])
              
              logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
              
              active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,n%*%s))
              Z_aktuell<-Z_alles[,active]
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              
              if(control$method=="EM")
              {  
                if(rnd.len==1)
                {
                  if(s==1)
                  {
                    P_akt<-diag(c(rep(0,lin_akt),rep(1/Q1.old,n*s)))
                  }else{
                    P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                    for(jf in 1:n)
                      P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
                  }
                }else{
                  if(all(s==1))
                  {
                    P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n)))
                  }else{
                    P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                    for(jf in 1:n[1])
                      P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
                    
                    for (zu in 2:rnd.len)
                    {
                      for(jf in 1:n[zu])
                        P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                              (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
                    }
                  }
                }
                if(is.null(family$multivariate)){
                  D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
                }else{
                  F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
                }
                
                InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
                if(class(InvFisher2)=="try-error")
                  InvFisher2<-try(solve(F_gross),silent=T)
                if(class(InvFisher2)=="try-error")
                {
                  half.index<-half.index+1  
                }else{
                  solve.test2<-TRUE 
                }}else{
                  solve.test2<-TRUE 
                }
            }
            
            betaact<-Delta[l,active]
            
            if(control$method=="EM")
            {        
              ############################# Q update ################
              if(rnd.len==1)
              {
                Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[l,(lin+1):(lin+s)]%*%t(Delta[l,(lin+1):(lin+s)])
                for (i in 2:n)
                  Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])
                Q1<-1/n*Q1
              }else{
                Q1<-matrix(0,sum(s),sum(s))
                Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[l,(lin+1):(lin+s[1])]%*%t(Delta[l,(lin+1):(lin+s[1])])
                for (i in 2:n[1])
                  Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
                Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
                
                for (zu in 2:rnd.len)
                {
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
                  for (i in 2:n[zu])
                    Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
                }
              }  
            }else{
              if(is.null(family$multivariate)){
                Eta_tilde<-Eta+(y-Mu)*1/D
              }else{
                Eta_tilde<-Eta+solve(D)%*%(y-Mu)
              }
              Betadach<-Delta[l,1:(lin)]
              
              aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
              X_aktuell<-Z_fastalles[,aktuell_vec]
              
              
              if(rnd.len==1)
              {
                
                if(s==1)
                {
                  if(Q1<1e-14)
                    low<-0
                  
                  optim.obj<-nlminb(sqrt(Q1),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
                  Q1<-as.matrix(optim.obj$par)^2
                }else{
                  Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
                  optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
                  
                  Q1<-matrix(0,s,s)
                  Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
                  Q1<-Q1+t(Q1)
                  diag(Q1)<-(optim.obj$par[1:s])
                  
                  #### Check for positiv definitness ########
                  for (ttt in 0:100)
                  {
                    Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                    Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                    Q_solvetest<-try(solve(Q1))
                    if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                      break
                  }
                }   
              }else{
                if(all(s==1))
                {
                  Q1_vec<-diag(Q1)
                  optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                  Q1<-diag(optim.obj$par)^2
                }else{
                  Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
                  
                  for (zu in 2:rnd.len)
                    Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                  
                  optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                  optim.vec<-optim.obj$par
                  
                  Q1<-matrix(0,sum(s),sum(s))
                  diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
                  if(s[1]>1)
                    Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
                  optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
                  
                  for (zu in 2:rnd.len)
                  {
                    diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
                    if(s[zu]>1)
                      Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
                    optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
                  }
                  
                  #### Check for positive definitness ########
                  for (ttt in 0:100)
                  {
                    Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                    Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                    Q_solvetest<-try(solve(Q1))
                    if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                      break
                  }
                }
              }}
            
            Q[[l+1]]<-Q1
            
            if(rnd.len==1)
            {
              if(s==1)
              {
                P1<-c(rep(0,lin),rep(1/Q1,n*s))
                P1<-diag(P1)
              }else{
                Q_inv<-solve(Q1)
                P1<-matrix(0,lin+n%*%s,lin+n%*%s)
                for(j in 1:n)
                  P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
              }
            }else{
              if(all(s==1))
              {
                P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
                P1<-diag(P1)
              }else{
                Q_inv<-list()
                Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
                P1<-matrix(0,lin+n%*%s,lin+n%*%s)
                for(jf in 1:n[1])
                  P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv[[1]]
                
                for (zu in 2:rnd.len)
                {
                  Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
                  for(jf in 1:n[zu])
                    P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                       (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
                }
              }
            }
            
            if(is.null(family$multivariate)){
              score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l,]
            }else{
              score_old2<-score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l,]
            }  
            lambda.max<-max(abs(score_vec2[(q+1):lin]))
            
            if (BLOCK)
            {
              grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
            }else{
              grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
            }
            score_vec2[(q+1):lin] <- grad.1
            
            crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
            t_edge<-crit.obj$min.rate
            
            optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[l,1:lin],ranef=Delta[l,(lin+1):(lin+n%*%s)],
                                               Grad=score_vec2,family=family,P=diag(P1[(lin+1):(lin+n%*%s)]), 
                                               lower = 0, upper = Inf,K=K))
            
            t_opt<-optim.obj$par
            
            if(!NRstep && !(all(active_old==active)))
              NRstep <- T
            
            tryNR <- (t_opt<t_edge)  && NRstep  
            
            if(tryNR) 
            {
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              if(rnd.len==1)
              {
                if(s==1)
                {
                  P_akt<-diag(c(rep(0,lin_akt),rep((Q1^(-1)),n*s)))
                }else{
                  P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                  for(jf in 1:n)
                    P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv
                }
              }else{
                if(all(s==1))
                {
                  P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q1)^(-1),n)))
                }else{
                  P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                  for(jf in 1:n[1])
                    P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv[[1]]
                  
                  for (zu in 2:rnd.len)
                  {
                    for(jf in 1:n[zu])
                      P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                            (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
                  }
                }
              }
              if(is.null(family$multivariate)){
                D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
                F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
              }else{
                F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
              }

              InvFisher3<-try(chol2inv(chol(F_gross)),silent=T)
              if(class(InvFisher3)=="try-error")
                InvFisher3<-try(solve(F_gross),silent=T)
              if(class(InvFisher3)=="try-error")
              {
                half.index<-half.index+1  
              }else{
                
                solve.test<-TRUE 
                Delta.test<-Delta[l,active]+nue*InvFisher3%*%score_old2[active]
                if(lin_akt>q)
                {
                  vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
                }else{
                  vorz<-T  
                }
              }
            }else{
              solve.test<-TRUE  
            }
          }  
          
          if(tryNR) 
            InvFisher<-InvFisher3
          
          score_old<-score_old2
          score_vec<-score_vec2
          Q1.very.old<-Q1.old
          Q_inv.very.old<-Q_inv.old
          Q1.old<-Q1
          Q_inv.old<-Q_inv    
          
          Eta.ma[l+1,]<-Eta
          active_old<-active

          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
      if(control$method!="EM")
      {  
        if(rnd.len==1)
        {
          if(s==1)
          {
            P_akt<-diag(c(rep(0,lin_akt),rep(1/Q1.old,n*s)))
          }else{
            P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
            for(jf in 1:n)
              P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
          }
        }else{
          if(all(s==1))
          {
            P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n)))
          }else{
            P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
            for(jf in 1:n[1])
              P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
            
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
                P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                      (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
            }
          }
        }
        if(is.null(family$multivariate)){
          D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
        }else{
          F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
        }

        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
      }

        if(is.null(family$multivariate)){
          FinalHat<-(Z_aktuell*sqrt(D*SigmaInv*D))%*%(InvFisher2%*%t(Z_aktuell*sqrt(D*SigmaInv*D)))
        }else{
          W_inv_t <- chol(D%*%(SigmaInv%*%t(D)))
          FinalHat<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
        }
        df<-sum(diag(FinalHat))
      
      if(control$overdispersion)
        phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-df)

      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }
      
      Delta_neu<-Delta[l,]
      Mu_opt<-Mu
      Qfinal<-Q[[l+1]]
      
      aaa<-!is.element(Delta_neu[1:(lin)],0)

      if(final.re)
      {    
        ############ final re-estimation
        
        if(rnd.len==1 && s==1)
        {  
          Q.max<-max(sqrt(unlist(Q)))
          Q.min<-min(sqrt(unlist(Q)))
        }else{
          Q.max<-max(Qfinal)+1
          Q.min<-min(Qfinal)-1e-10
        }
        if(rnd.len==1)
        {
          glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=Qfinal,K=K,
                                   Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                                   family=family,method=control$method.final,overdispersion=control$overdispersion,
                                   phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                   Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=q_start,K=K,
                                      Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                                      family=family,method=control$method.final,overdispersion=control$overdispersion,
                                      phi=control$phi,print.iter.final=control$print.iter.final,
                                      eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }
        }else{
          glmm_fin<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,K=K,
                                                Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                                                family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          if(class(glmm_fin)=="try-error"|| glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=q_start,K=K,
                                                   Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                                                   family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                   phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                   eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)    
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
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
        Delta_neu2[c(aaa,rep(T,n%*%s))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,n%*%s)),c(aaa,rep(T,n%*%s))]<-glmm_fin$Standard_errors
        Qfinal<-glmm_fin$Q
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
        glmm_fin$ranef.logLik<-ranef.logLik
        complexity<-df
      }
      
      if(rnd.len==1)
      {
        if(s==1)
          Qfinal<-sqrt(Qfinal)
        
        if(!is.matrix(Qfinal))
          Qfinal<-as.matrix(Qfinal)
        colnames(Qfinal)<-random.labels
        rownames(Qfinal)<-random.labels
      }else{
        Qfinal_old<-Qfinal
        Qfinal<-list()
        Qfinal[[1]]<-as.matrix(Qfinal_old[1:s[1],1:s[1]])
        colnames(Qfinal[[1]])<-random.labels[[1]]
        rownames(Qfinal[[1]])<-random.labels[[1]]
        
        if(s[1]==1)
          Qfinal[[1]]<-sqrt(Qfinal[[1]])
        
        for (zu in 2:rnd.len)
        {
          Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
          colnames(Qfinal[[zu]])<-random.labels[[zu]]
          rownames(Qfinal[[zu]])<-random.labels[[zu]]
          
          if(s[zu]==1)
            Qfinal[[zu]]<-sqrt(Qfinal[[zu]])
          
        }
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
      names(Delta_neu)[(lin+1):(lin+n%*%s)]<-colnames(W)
      colnames(Standard_errors)[(lin+1):(lin+n%*%s)]<-rownames(Standard_errors)[(lin+1):(lin+n%*%s)]<-colnames(W)
      colnames(Delta)<-c(final.names,colnames(W))

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
        
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=glmm_fin$ranef.logLik,family=family,penal=T,K=K)
        
        
        if(control$complexity!="hat.matrix")  
        {  
          if(rnd.len==1)
          {
            complexity<-0.5*(s*(s+1))
          }else{
            complexity<-0.5*(s[1]*(s[1]+1))
            for(zu in 2:rnd.len)
              complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
          }
          complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)
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
      ret.obj$ranef<-Delta_neu[(lin+1):(lin+n%*%s)]
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$ranerror<-Standard_errors[(lin+1):(lin+n%*%s)]
      ret.obj$Q_long<-Q
      ret.obj$Q<-Qfinal
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$conv.step<-conv.step
      ret.obj$newrndfrml<-newrndfrml
      ret.obj$subject<-subject.names
      ret.obj$data<-data
      ret.obj$rnd.len<-rnd.len
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
        Eta_start<-X%*%beta_null+Phi%*%smooth_null+W%*%ranef_null
      }else{
        Eta_start<-rep(beta_null,N)+Phi%*%smooth_null+W%*%ranef_null
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
      
      if(rnd.len==1)
      {
        if(s==1)
          Q_start<-q_start
      }else{
        if(all(s==1))
        {
          Q_start<-diag(diag(q_start),sum(s))
        }else{
          Q_start<-matrix(0,sum(s),sum(s))
          Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,Phi,W)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,(lin+dim.smooth+n%*%s))
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null
      Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-t(ranef_null)
      active_old<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
      
      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      logLik.vec<-c()
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      Q_inv<-NULL
      Q_inv.old.temp<-NULL
      Q_inv.start<-NULL
      
      Q<-list()
      Q[[1]]<-Q_start
      
      l=1
      
      if(diff.ord>1)
      {
        k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
      }else{
        k2<-rep(1,nbasis)
      }
      k22<-rep(k2,m)
      penal.vec<-penal*k22
      
      Q_inv<-NULL
      Q_inv.old.temp<-NULL
      Q_inv.start<-NULL
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin),penal.vec,rep(1/Q_start,n*s))
          P1<-diag(P1)
        }else{
          Q_inv.start<-chol2inv(chol(Q_start))
          P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
          diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
          for(j in 1:n)
            P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin),penal.vec,rep(diag(Q_start)^(-1),n))
          P1<-diag(P1)
        }else{
          Q_inv.start<-list()
          Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
          P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
          diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
          for(jf in 1:n[1])
            P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            for(jf in 1:n[zu])
              P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
          }
        }
      }
      
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
      
      score_vec[(q+1):lin] <- grad.1
      
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      #ranef.logLik<- -0.5*t(Delta_start[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s),(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]%*%Delta_start[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta_start[1:(lin+dim.smooth)],ranef=Delta_start[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)],
                                         Grad=score_vec,family=family,P=diag(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]), 
                                         lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
      half.index<-0
      solve.test<-FALSE
      ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
      while(!solve.test)
      {  
        
        solve.test2<-FALSE  
        while(!solve.test2)
        {  
          
          if(half.index>50)
          {
            stop("Fisher matrix not invertible")
          }
          
          Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
          if(t_opt>t_edge & half.index==0)
            Delta[1,crit.obj$whichmin+q]<-0  
          Eta<-Z_alles%*%Delta[1,]

          if(is.null(family$multivariate)){
            D<-family$mu.eta(Eta)
            Mu<-family$linkinv(Eta)
            SigmaInv <- 1/family$variance(Mu)
          }else{
            Eta_cat <- matrix(Eta_start, byrow = TRUE, ncol = K)
            Mu_cat <- family$linkinv(Eta_cat)
            D <- family$deriv.mat(Mu_cat)
            SigmaInv <- family$SigmaInv(Mu_cat)
            Mu <- c(t(Mu_cat))
          }          
          ranef.logLik<- -0.5*t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s),(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]%*%Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])
          
          logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
          
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
          if (control$method=="EM")
          {  
            if(rnd.len==1)
            {
              if(s==1)
              {
                P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q_start,n*s)))
              }else{
                P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                for(jf in 1:n)
                  P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.start
              }
            }else{
              if(all(s==1))
              {
                P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q_start)^(-1),n)))
              }else{
                P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                for(jf in 1:n[1])
                  P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.start[[1]]
                
                for (zu in 2:rnd.len)
                {
                  for(jf in 1:n[zu])
                    P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                          (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
                }
              }
            }
            if(is.null(family$multivariate)){
              D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
            }else{
              F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
            }
            InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
            if(class(InvFisher2)=="try-error")
              InvFisher2<-try(solve(F_gross),silent=T)
            if(class(InvFisher2)=="try-error")
            {
              #stop("Fisher matrix not invertible")  
              half.index<-half.index+1  
            }else{
              solve.test2<-TRUE 
            }}else{
              solve.test2<-TRUE
            }
        }
        
        betaact<-Delta[1,active]
        
        if (control$method=="EM")
        {   
          ############################# Q update ################
          if(rnd.len==1)
          {
            Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)])
            for (i in 2:n)
              Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
            Q1<-1/n*Q1
          }else{
            Q1<-matrix(0,sum(s),sum(s))
            Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
            for (i in 2:n[1])
              Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
            Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
            
            for (zu in 2:rnd.len)
            {
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
              for (i in 2:n[zu])
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
            }
          }
        }else{
          if(is.null(family$multivariate)){
            Eta_tilde<-Eta+(y-Mu)*1/D
          }else{
            Eta_tilde<-Eta+solve(D)%*%(y-Mu)
          }
          Betadach<-Delta[1,1:(lin+dim.smooth)]     
          aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
          X_aktuell<-Z_fastalles[,aktuell_vec]
          
          if(rnd.len==1)
          {
            
            if(s==1)
            {
              upp<-min(20,50*Q_start)
              low<-1e-14
              optim.obj<-try(nlminb(sqrt(Q_start),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                    Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
              if(class(optim.obj)=="try-error")
                optim.obj<-try(bobyqa(sqrt(Q_start),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                      Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
              Q1<-as.matrix(optim.obj$par)^2
            }else{
              q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
              up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
              upp<-rep(up1,length(q_start_vec))
              low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
              #   kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
              optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                    Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
              Q1<-matrix(0,s,s)
              Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
              Q1<-Q1+t(Q1)
              diag(Q1)<-(optim.obj$par[1:s])
              
              #### Check for positive definitness ########
              for (ttt in 0:100)
              {
                Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                Q_solvetest<-try(solve(Q1))
                if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                  break
              }
            }   
          }else{
            if(all(s==1))
            {
              q_start_vec<-diag(q_start)
              upp<-rep(min(20,50*diag(q_start)),sum(s))
              low<-rep(1e-14,sum(s))
              optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
              Q1<-diag(optim.obj$par)^2
            }else{
              q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
              up1<-min(20,50*max(q_start_vec))
              low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))
              
              for (zu in 2:rnd.len)
              {
                q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                up1<-min(20,50*max(q_start_vec))
                low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
              }
              upp<-rep(up1,length(q_start_vec))
              optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
              optim.vec<-optim.obj$par
              
              Q1<-matrix(0,sum(s),sum(s))
              diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
              if(s[1]>1)
                Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
              optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
              
              for (zu in 2:rnd.len)
              {
                diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
                if(s[zu]>1)
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
                optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
              }
              
              #### Check for positive definitness ########
              for (ttt in 0:100)
              {
                Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                Q_solvetest<-try(solve(Q1))
                if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                  break
              }
            }
          }}
        
        Q[[2]]<-Q1

        NRstep<-F
        
        if(rnd.len==1)
        {
          if(s==1)
          {
            P1<-c(rep(0,lin),penal.vec,rep(1/Q1,n*s))
            P1<-diag(P1)
          }else{
            Q_inv<-chol2inv(chol(Q1))
            Q_inv.old.temp<-Q_inv
            P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
            diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
            for(j in 1:n)
              P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
          }
        }else{
          if(all(s==1))
          {
            P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
            P1<-diag(P1)
          }else{
            Q_inv<-list()
            Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
            P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
            diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
            for(jf in 1:n[1])
              P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
            
            for (zu in 2:rnd.len)
            {
              Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
              for(jf in 1:n[zu])
                P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                   (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
            }
          }
        }
        
        if(is.null(family$multivariate)){
          score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
        }else{
          score_old2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
        }  
        lambda.max<-max(abs(score_vec2[(q+1):lin]))
        
        if (BLOCK)
        {
          grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
        }else{
          grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
        }
        score_vec2[(q+1):lin] <- grad.1
        
        crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
        t_edge<-crit.obj$min.rate
        
        optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[1,1:(lin+dim.smooth)],ranef=Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)],
                                           Grad=score_vec2,family=family,P=diag(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]), 
                                           lower = 0, upper = Inf,K=K))
        
        t_opt<-optim.obj$par
        
        if(!NRstep && !(all(active_old==active)))
          NRstep <- T
        
        tryNR <- (t_opt<t_edge)  && NRstep  
        
        vorz <- T
        
        if(tryNR) 
        {
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          if(rnd.len==1)
          {
            if(s==1)
            {
              P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q1,n*s)))
            }else{
              P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
              diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
              for(jf in 1:n)
                P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                      (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
            }
          }else{
            if(all(s==1))
            {
              P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q1)^(-1),n)))
            }else{
              P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
              diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
              for(jf in 1:n[1])
                P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                      (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
              
              for (zu in 2:rnd.len)
              {
                for(jf in 1:n[zu])
                  P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                        (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              }
            }
          }
          if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
          }else{
            F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
          }
          InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher)=="try-error")
            InvFisher<-try(solve(F_gross),silent=T)
          if(class(InvFisher)=="try-error")
          {
            half.index<-half.index+1  
          }else{
            solve.test<-TRUE 
            Delta.test<-Delta[1,active]+nue*InvFisher%*%score_vec2[active]
            if(lin_akt>q)
            {
              vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
            }else{
              vorz<-T  
            }
          }
        }else{
          solve.test<-TRUE  
        }
      }   
      
      Eta.ma[2,]<-Eta
      
      score_old<-score_old2
      score_vec<-score_vec2
      Q1.old<-Q1
      Q_inv.old<-Q_inv    
      
      Q1.very.old<-Q_start
      Q_inv.very.old<-Q_inv.start

      active_old<-active
      
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {

          if(control$print.iter)
            message("Iteration ",l)

          if(!vorz)
            tryNR<-F
          
          half.index<-0
          solve.test<-FALSE
          ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
          while(!solve.test)
          {  
            
            solve.test2<-FALSE  
            while(!solve.test2)
            {  
              
              if(half.index>50)
              {
                half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old;
              }
              
              if(tryNR)
              {
                Delta[l,active_old]<- Delta[l-1,active_old]+nue*(0.5^half.index)*InvFisher%*%score_old[active_old]
                NRstep<-T
#                print("NR")
              }else{
                Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
                if(t_opt>t_edge & half.index==0)
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
              
              ranef.logLik<- -0.5*t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s),(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]%*%Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])
              
              logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
              
              active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
              Z_aktuell<-Z_alles[,active]
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              
              if (control$method=="EM")
              {  
                if(rnd.len==1)
                {
                  if(s==1)
                  {
                    P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q1.old,n*s)))
                  }else{
                    P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                    diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                    for(jf in 1:n)
                      P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                            (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.old
                  }
                }else{
                  if(all(s==1))
                  {
                    P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q1.old)^(-1),n)))
                  }else{
                    P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                    diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                    for(jf in 1:n[1])
                      P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                            (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.old[[1]]
                    
                    for (zu in 2:rnd.len)
                    {
                      for(jf in 1:n[zu])
                        P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                              (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
                    }
                  }
                 }
                if(is.null(family$multivariate)){
                  D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
                }else{
                  F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
                }
                InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
                if(class(InvFisher2)=="try-error")
                  InvFisher2<-try(solve(F_gross),silent=T)
                if(class(InvFisher2)=="try-error")
                {
                  half.index<-half.index+1  
                }else{
                  solve.test2<-TRUE 
                }}else{
                  solve.test2<-TRUE 
                }
            }
            
            betaact<-Delta[l,active]
            
            if (control$method=="EM")
            {        
              ############################# Q update ################
              if(rnd.len==1)
              {
                Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)])
                for (i in 2:n)
                  Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
                Q1<-1/n*Q1
              }else{
                Q1<-matrix(0,sum(s),sum(s))
                Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
                for (i in 2:n[1])
                  Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
                Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
                
                for (zu in 2:rnd.len)
                {
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
                  for (i in 2:n[zu])
                    Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
                }
              }  
            }else{
              if(is.null(family$multivariate)){
                Eta_tilde<-Eta+(y-Mu)*1/D
              }else{
                Eta_tilde<-Eta+solve(D)%*%(y-Mu)
              }
              
              Betadach<-Delta[l,1:(lin+dim.smooth)]
              
              aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
              X_aktuell<-Z_fastalles[,aktuell_vec]
              
              if(rnd.len==1)
              {
                
                if(s==1)
                {
                  if(Q1<1e-14)
                    low<-0
                  
                  optim.obj<-try(nlminb(sqrt(Q1),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
                  if(class(optim.obj)=="try-error")
                    optim.obj<-bobyqa(sqrt(Q1),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
                  
                  Q1<-as.matrix(optim.obj$par)^2
                }else{
                  Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
                  optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
                  
                  Q1<-matrix(0,s,s)
                  Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
                  Q1<-Q1+t(Q1)
                  diag(Q1)<-(optim.obj$par[1:s])
                  
                  #### Check for positiv definitness ########
                  for (ttt in 0:100)
                  {
                    Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                    Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                    Q_solvetest<-try(solve(Q1))
                    if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                      break
                  }
                }   
              }else{
                if(all(s==1))
                {
                  Q1_vec<-diag(Q1)
                  optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                  Q1<-diag(optim.obj$par)^2
                }else{
                  Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
                  
                  for (zu in 2:rnd.len)
                    Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                  
                  optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                  optim.vec<-optim.obj$par
                  
                  
                  Q1<-matrix(0,sum(s),sum(s))
                  diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
                  if(s[1]>1)
                    Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
                  optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
                  
                  for (zu in 2:rnd.len)
                  {
                    diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
                    if(s[zu]>1)
                      Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
                    optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
                  }
                  
                  #### Check for positive definitness ########
                  for (ttt in 0:100)
                  {
                    Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                    Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                    Q_solvetest<-try(solve(Q1))
                    if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                      break
                  }
                }
              }}

            Q[[l+1]]<-Q1
            
            if(rnd.len==1)
            {
              if(s==1)
              {
                P1<-c(rep(0,lin),penal.vec,rep(1/Q1,n*s))
                P1<-diag(P1)
              }else{
                Q_inv<-chol2inv(chol(Q1))
                P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
                diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
                for(j in 1:n)
                  P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
              }
            }else{
              if(all(s==1))
              {
                P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
                P1<-diag(P1)
              }else{
                Q_inv<-list()
                Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
                P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
                diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
                for(jf in 1:n[1])
                  P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
                
                for (zu in 2:rnd.len)
                {
                  Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
                  for(jf in 1:n[zu])
                    P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                       (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
                }
              }
            }
            
            if(is.null(family$multivariate)){
              score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l,]
            }else{
              score_old2<-score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l,]
            }  
            lambda.max<-max(abs(score_vec2[(q+1):lin]))
            
            if (BLOCK)
            {
              grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
            }else{
              grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
            }
            score_vec2[(q+1):lin] <- grad.1
            
            crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
            t_edge<-crit.obj$min.rate
            
            optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[l,1:(lin+dim.smooth)],ranef=Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)],
                                               Grad=score_vec2,family=family,P=diag(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]), 
                                               lower = 0, upper = Inf,K=K))
            
            t_opt<-optim.obj$par
            
            if(!NRstep && !(all(active_old==active)))
              NRstep <- T
            
            tryNR <- (t_opt<t_edge)  && NRstep  
            
            if(tryNR) 
            {
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              if(rnd.len==1)
              {
                if(s==1)
                {
                  P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q1,n*s)))
                }else{
                  P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                  diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                  for(jf in 1:n)
                    P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                          (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
                }
              }else{
                if(all(s==1))
                {
                  P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q1)^(-1),n)))
                }else{
                  P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                  diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                  for(jf in 1:n[1])
                    P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                          (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
                  
                  for (zu in 2:rnd.len)
                  {
                    for(jf in 1:n[zu])
                      P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                            (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
                  }
                }
              }
              if(is.null(family$multivariate)){
                D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
                F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
              }else{
                F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
              }
              InvFisher3<-try(chol2inv(chol(F_gross)),silent=T)
              if(class(InvFisher3)=="try-error")
                InvFisher3<-try(solve(F_gross),silent=T)
              if(class(InvFisher3)=="try-error")
              {
                half.index<-half.index+1  
              }else{
                
                
                solve.test<-TRUE 
                Delta.test<-Delta[l,active]+nue*InvFisher3%*%score_old2[active]
                if(lin_akt>q)
                {
                  vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
                }else{
                  vorz<-T  
                }
              }
            }else{
              solve.test<-TRUE  
            }
            
          }  
          
          if(tryNR) 
            InvFisher<-InvFisher3
          
          score_old<-score_old2
          score_vec<-score_vec2
          Q1.very.old<-Q1.old
          Q_inv.very.old<-Q_inv.old
          Q1.old<-Q1
          Q_inv.old<-Q_inv    
          
          Eta.ma[l+1,]<-Eta
          active_old<-active
          
          
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
      if(control$method!="EM")
      {  
        if(rnd.len==1)
        {
          if(s==1)
          {
            P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q1.old,n*s)))
          }else{
            P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
            for(jf in 1:n)
              P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                    (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.old
          }
        }else{
          if(all(s==1))
          {
            P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q1.old)^(-1),n)))
          }else{
            P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
            for(jf in 1:n[1])
              P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                    (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.old[[1]]
            
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
                P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                      (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
            }
          }
        }
        if(is.null(family$multivariate)){
          D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
        }else{
          F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
        }
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
      }          
      
      if(is.null(family$multivariate)){
        FinalHat<-(Z_aktuell*sqrt(D*SigmaInv*D))%*%(InvFisher2%*%t(Z_aktuell*sqrt(D*SigmaInv*D)))
      }else{
        W_inv_t <- chol(D%*%(SigmaInv%*%t(D)))
        FinalHat<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
      }
      df<-sum(diag(FinalHat))
      
      if(control$overdispersion)
        phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-df)
      
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
      Qfinal<-Q[[l+1]]
      
      aaa<-!is.element(Delta_neu[1:(lin)],0)
      
      if(final.re)
      {    
        ############ final re-estimation
        
        if(s==1)
        {  
          Q.max<-max(sqrt(unlist(Q)))
          Q.min<-min(sqrt(unlist(Q)))
        }else{
          Q.max<-max(Qfinal)+1
          Q.min<-min(Qfinal)-1e-10
        }
        
        
        if(rnd.len==1)
        {
          glmm_fin<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,n,penal.vec,q_start=Qfinal,K=K,
                                          Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                                          s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                          phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                          Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,n,penal.vec,q_start=q_start,K=K,
                                             Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                                             s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                             phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                             Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }
        }else{
          glmm_fin<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=Qfinal,K=K,
                                                       Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                                                       s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                       phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                                       Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=q_start,K=K,
                                                          Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                                                          s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                          phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                          eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
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
        EDF.matrix<-glmm_fin$EDF.matrix
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])  
        Delta_neu2<-Delta_neu
        Delta_neu2[c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,dim.smooth+n%*%s)),c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Standard_errors
        Qfinal<-glmm_fin$Q
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
        glmm_fin$ranef.logLik<-ranef.logLik
        complexity<-df
        
        lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))
        if(rnd.len==1)
        {
          if(s==1)
          {
            P1a<-diag(c(rep(0,lin_akt+dim.smooth),rep(1/Q1,n*s)))
          }else{
            P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            for(jf in 1:n)
            {
              P1a[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                  (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
            }  
          }
        }else{
          if(all(s==1))
          {
            P1a<-diag(c(rep(0,lin_akt+dim.smooth),rep(diag(Q1)^(-1),n)))
          }else{
            P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            for(jf in 1:n[1])
            {
              P1a[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                  (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
            }
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
              {
                P1a[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              }}
          }
        }
        if(class(InvFisher2)=="try-error")
        {
          warning("No EDF's for smooth functions available, as Fisher matrix not invertible!")
          complexity.smooth<-dim.smooth
        }else{  
          ###### EDF of spline; compare Wood's Book on page 167
          if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            EDF.matrix<-InvFisher2%*%(t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P1a)
          }else{
            EDF.matrix<-InvFisher2%*%(t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D)))%*%Z_aktuell+P1a)
          }
          complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])
        }
      }  
      
      if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
        complexity.smooth<-dim.smooth
      
      if(rnd.len==1)
      {
        if(s==1)
          Qfinal<-sqrt(Qfinal)
        
        if(!is.matrix(Qfinal))
          Qfinal<-as.matrix(Qfinal)
        colnames(Qfinal)<-random.labels
        rownames(Qfinal)<-random.labels
      }else{
        Qfinal_old<-Qfinal
        Qfinal<-list()
        Qfinal[[1]]<-as.matrix(Qfinal_old[1:s[1],1:s[1]])
        colnames(Qfinal[[1]])<-random.labels[[1]]
        rownames(Qfinal[[1]])<-random.labels[[1]]
        
        for (zu in 2:rnd.len)
        {
          Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
          colnames(Qfinal[[zu]])<-random.labels[[zu]]
          rownames(Qfinal[[zu]])<-random.labels[[zu]]
        }
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
      
      names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-rownames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      names(Delta_neu)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
      colnames(Standard_errors)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-rownames(Standard_errors)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
      colnames(Delta)<-c(final.names,colnames(Phi),colnames(W))

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
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=glmm_fin$ranef.logLik,family=family,penal=T,K=K)
        
        if(control$complexity!="hat.matrix")  
        {  
          if(rnd.len==1)
          {
            complexity<-0.5*(s*(s+1))
          }else{
            complexity<-0.5*(s[1]*(s[1]+1))
            for(zu in 2:rnd.len)
              complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
          }
          complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
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
      ret.obj$ranef<-Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$ranerror<-Standard_errors[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
      ret.obj$Q_long<-Q
      ret.obj$Q<-Qfinal
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$newrndfrml<-newrndfrml
      ret.obj$subject<-names(rnd)
      ret.obj$data<-data
      ret.obj$rnd.len<-rnd.len
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
        Eta_start<-X%*%beta_null+W%*%ranef_null
      }else{
        Eta_start<-rep(beta_null,N)+W%*%ranef_null
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
      
      
      if(rnd.len==1)
      {
        lin0<-sum(beta_null!=0)
        if(s==1)
        {
          Q_start<-diag(q_start,s)
        }else{
          Q_start<-q_start
        }
      }else{
        lin0<-sum(beta_null!=0)
        if(all(s==1))
        {
          Q_start<-diag(diag(q_start),sum(s))
        }else{
          Q_start<-matrix(0,sum(s),sum(s))
          Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,W)
      
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,(lin+n%*%s))
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+n%*%s)]<-t(ranef_null)

      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      logLik.vec<-c()

      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      Q_inv<-NULL
      Q_inv.old.temp<-NULL
      Q_inv.start<-NULL
      
      Q<-list()
      Q[[1]]<-Q_start
      
      l=1
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin),rep((Q_start^(-1)),n*s))
          P1<-diag(P1)
        }else{
          Q_inv.start<-chol2inv(chol(Q_start))
          P1<-matrix(0,lin+n%*%s,lin+n%*%s)
          for(j in 1:n)
            P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv.start
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin),rep(diag(Q_start)^(-1),n))
          P1<-diag(P1)
        }else{
          Q_inv.start<-list()
          Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
          P1<-matrix(0,lin+n%*%s,lin+n%*%s)
          for(jf in 1:n[1])
            P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv.start[[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            for(jf in 1:n[zu])
              P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
          }
        }
      }
      
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
      
      
      score_vec[(q+1):lin] <- grad.1
      
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
#      ranef.logLik<- -0.5*t(Delta_start[(lin+1):(lin+n%*%s)])%*%(P1[(lin+1):(lin+n%*%s),(lin+1):(lin+n%*%s)]%*%Delta_start[(lin+1):(lin+n%*%s)])
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta_start[1:lin],ranef=Delta_start[(lin+1):(lin+n%*%s)],
                                         Grad=score_vec,family=family,P=diag(P1[(lin+1):(lin+n%*%s)]), 
                                         lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
      half.index<-0
      
      solve.test2<-FALSE  
      while(!solve.test2)
      {  
        
        if(half.index>50)
        {
          stop("Fisher matrix not invertible")
        }
        
        Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
        if(t_opt>t_edge & half.index==0)
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
        
        ranef.logLik<- -0.5*t(Delta[1,(lin+1):(lin+n%*%s)])%*%(P1[(lin+1):(lin+n%*%s),(lin+1):(lin+n%*%s)]%*%Delta[1,(lin+1):(lin+n%*%s)])
        
        logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if (control$method=="EM")
        {  
          if(rnd.len==1)
          {
            if(s==1)
            {
              P_akt<-diag(c(rep(0,lin_akt),rep(1/Q_start,n*s)))
            }else{
              P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
              for(jf in 1:n)
                P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.start
            }
          }else{
            if(all(s==1))
            {
              P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q_start)^(-1),n)))
            }else{
              P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
              for(jf in 1:n[1])
                P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.start[[1]]
              
              for (zu in 2:rnd.len)
              {
                for(jf in 1:n[zu])
                  P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                        (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
              }
            }
          }
          if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
          }else{
            F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
          }
          InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher2)=="try-error")
            InvFisher2<-try(solve(F_gross),silent=T)
          if(class(InvFisher2)=="try-error")
          {
            #stop("Fisher matrix not invertible")  
            half.index<-half.index+1  
          }else{
            solve.test2<-TRUE 
          }}else{
            solve.test2<-TRUE
          }
      }
      
      betaact<-Delta[1,active]
      
      if (control$method=="EM")
      {   
        ############################# Q update ################
        if(rnd.len==1)
        {
          Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
          for (i in 2:n)
            Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
          Q1<-1/n*Q1
        }else{
          Q1<-matrix(0,sum(s),sum(s))
          Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
          for (i in 2:n[1])
            Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
          Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
            for (i in 2:n[zu])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }else{
        if(is.null(family$multivariate)){
          Eta_tilde<-Eta+(y-Mu)*1/D
        }else{
          Eta_tilde<-Eta+solve(D)%*%(y-Mu)
        }
        Betadach<-Delta[1,1:(lin)]     
        aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
        X_aktuell<-Z_fastalles[,aktuell_vec]
        
        if(rnd.len==1)
        {
          
          if(s==1)
          {
            upp<-min(20,50*Q_start)
            low<-1e-14
            optim.obj<-nlminb(sqrt(Q_start),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
            Q1<-as.matrix(optim.obj$par)^2
          }else{
            q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
            up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
            upp<-rep(up1,length(q_start_vec))
            low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
            optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
            Q1<-matrix(0,s,s)
            Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
            Q1<-Q1+t(Q1)
            diag(Q1)<-(optim.obj$par[1:s])
            
            #### Check for positive definitness ########
            for (ttt in 0:100)
            {
              Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
              Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
              Q_solvetest<-try(solve(Q1))
              if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                break
            }
          }   
        }else{
          if(all(s==1))
          {
            q_start_vec<-diag(q_start)
            upp<-rep(min(20,50*diag(q_start)),sum(s))
            low<-rep(1e-14,sum(s))
            optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
            Q1<-diag(optim.obj$par)^2
          }else{
            q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
            up1<-min(20,50*max(q_start_vec))
            low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))
            
            for (zu in 2:rnd.len)
            {
              q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
              up1<-min(20,50*max(q_start_vec))
              low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
            }
            upp<-rep(up1,length(q_start_vec))
            optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
            optim.vec<-optim.obj$par
            
            
            Q1<-matrix(0,sum(s),sum(s))
            diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
            if(s[1]>1)
              Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
            optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
            
            for (zu in 2:rnd.len)
            {
              diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
              if(s[zu]>1)
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
              optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
            }
            
            #### Check for positive definitness ########
            for (ttt in 0:100)
            {
              Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
              Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
              Q_solvetest<-try(solve(Q1))
              if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                break
            }
          }
        }}
      
      Q[[2]]<-Q1
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
          P1<-diag(P1)
        }else{
          Q_inv<-solve(Q1)
          P1<-matrix(0,lin+n%*%s,lin+n%*%s)
          for(j in 1:n)
            P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
          P1<-diag(P1)
        }else{
          Q_inv<-list()
          Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
          P1<-matrix(0,lin+n%*%s,lin+n%*%s)
          for(jf in 1:n[1])
            P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv[[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            for(jf in 1:n[zu])
              P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
          }
        }
      }
      
      if(is.null(family$multivariate)){
        score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
      }else{
        score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
      }  
      lambda.max<-max(abs(score_vec2[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      score_vec2[(q+1):lin] <- grad.1

      crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[1,1:lin],ranef=Delta[1,(lin+1):(lin+n%*%s)],
                                         Grad=score_vec2,family=family,P=diag(P1[(lin+1):(lin+n%*%s)]), 
                                         lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      Eta.ma[2,]<-Eta
      
      score_vec<-score_vec2
      Q1.old<-Q1
      Q_inv.old<-Q_inv    
      
      Q1.very.old<-Q_start
      Q_inv.very.old<-Q_inv.start
      
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
          half.index<-0
          
          solve.test2<-FALSE  
          while(!solve.test2)
          {  
            
            if(half.index>50)
            {
              half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old
            }
            
            Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
            if(t_opt>t_edge & half.index==0)
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

            ranef.logLik<- -0.5*t(Delta[l,(lin+1):(lin+n%*%s)])%*%(P1[(lin+1):(lin+n%*%s),(lin+1):(lin+n%*%s)]%*%Delta[l,(lin+1):(lin+n%*%s)])
            
            logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
            
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,n%*%s))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
            if(control$method=="EM")
            {  
              if(rnd.len==1)
              {
                if(s==1)
                {
                  P_akt<-diag(c(rep(0,lin_akt),rep(1/Q1.old,n*s)))
                }else{
                  P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                  for(jf in 1:n)
                    P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
                }
              }else{
                if(all(s==1))
                {
                  P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n)))
                }else{
                  P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                  for(jf in 1:n[1])
                    P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
                  
                  for (zu in 2:rnd.len)
                  {
                    for(jf in 1:n[zu])
                      P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                            (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
                  }
                }
              }
              if(is.null(family$multivariate)){
                D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
                F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
              }else{
                F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
              }
              InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
              if(class(InvFisher2)=="try-error")
                InvFisher2<-try(solve(F_gross),silent=T)
              if(class(InvFisher2)=="try-error")
              {
                half.index<-half.index+1  
              }else{
                solve.test2<-TRUE 
              }}else{
                solve.test2<-TRUE 
              }
          }
          
          betaact<-Delta[l,active]
          
          if (control$method=="EM")
          {        
            ############################# Q update ################
            if(rnd.len==1)
            {
              Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[l,(lin+1):(lin+s)]%*%t(Delta[l,(lin+1):(lin+s)])
              for (i in 2:n)
                Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])
              Q1<-1/n*Q1
            }else{
              Q1<-matrix(0,sum(s),sum(s))
              Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[l,(lin+1):(lin+s[1])]%*%t(Delta[l,(lin+1):(lin+s[1])])
              for (i in 2:n[1])
                Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
              Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
              
              for (zu in 2:rnd.len)
              {
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
                for (i in 2:n[zu])
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
              }
            }  
          }else{
            if(is.null(family$multivariate)){
              Eta_tilde<-Eta+(y-Mu)*1/D
            }else{
              Eta_tilde<-Eta+solve(D)%*%(y-Mu)
            }
            
            Betadach<-Delta[l,1:(lin)]
            
            aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
            X_aktuell<-Z_fastalles[,aktuell_vec]
            
            
            if(rnd.len==1)
            {
              
              if(s==1)
              {
                if(Q1<1e-14)
                  low<-0
                
                optim.obj<-nlminb(sqrt(Q1),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
                Q1<-as.matrix(optim.obj$par)^2
              }else{
                Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
                optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
                
                Q1<-matrix(0,s,s)
                Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
                Q1<-Q1+t(Q1)
                diag(Q1)<-(optim.obj$par[1:s])
                
                #### Check for positiv definitness ########
                for (ttt in 0:100)
                {
                  Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                  Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                  Q_solvetest<-try(solve(Q1))
                  if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                    break
                }
              }   
            }else{
              if(all(s==1))
              {
                Q1_vec<-diag(Q1)
                optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                Q1<-diag(optim.obj$par)^2
              }else{
                Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
                
                for (zu in 2:rnd.len)
                  Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                
                optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                optim.vec<-optim.obj$par
                
                Q1<-matrix(0,sum(s),sum(s))
                diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
                if(s[1]>1)
                  Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
                optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
                
                for (zu in 2:rnd.len)
                {
                  diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
                  if(s[zu]>1)
                    Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
                  optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
                }
                
                #### Check for positive definitness ########
                for (ttt in 0:100)
                {
                  Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                  Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                  Q_solvetest<-try(solve(Q1))
                  if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                    break
                }
              }
            }}
          
          Q[[l+1]]<-Q1
          
          if(rnd.len==1)
          {
            if(s==1)
            {
              P1<-c(rep(0,lin),rep(1/Q1,n*s))
              P1<-diag(P1)
            }else{
              Q_inv<-solve(Q1)
              P1<-matrix(0,lin+n%*%s,lin+n%*%s)
              for(j in 1:n)
                P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
            }
          }else{
            if(all(s==1))
            {
              P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
              P1<-diag(P1)
            }else{
              Q_inv<-list()
              Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
              P1<-matrix(0,lin+n%*%s,lin+n%*%s)
              for(jf in 1:n[1])
                P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv[[1]]
              
              for (zu in 2:rnd.len)
              {
                Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
                for(jf in 1:n[zu])
                  P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                     (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              }
            }
          }
          
          if(is.null(family$multivariate)){
            score_old2<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l,]
          }else{
            score_old2<-score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l,]
          }  
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
          }
          
          score_vec2[(q+1):lin] <- grad.1
          
          t_edge<-crit.obj$min.rate

          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[l,1:lin],ranef=Delta[l,(lin+1):(lin+n%*%s)],
                                             Grad=score_vec2,family=family,P=diag(P1[(lin+1):(lin+n%*%s)]), 
                                             lower = 0, upper = Inf,K=K))
          
          t_opt<-optim.obj$par
          
          score_vec<-score_vec2
          Q1.very.old<-Q1.old
          Q_inv.very.old<-Q_inv.old
          Q1.old<-Q1
          Q_inv.old<-Q_inv    
          
          Eta.ma[l+1,]<-Eta
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
      if(control$method!="EM")
      {  
        if(rnd.len==1)
        {
          if(s==1)
          {
            P_akt<-diag(c(rep(0,lin_akt),rep(1/Q1.old,n*s)))
          }else{
            P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
            for(jf in 1:n)
              P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
          }
        }else{
          if(all(s==1))
          {
            P_akt<-diag(c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n)))
          }else{
            P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
            for(jf in 1:n[1])
              P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
            
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
                P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                      (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
            }
          }
        }
        if(is.null(family$multivariate)){
          D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
        }else{
          F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
        }
        
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
      }
      
      if(is.null(family$multivariate)){
        FinalHat<-(Z_aktuell*sqrt(D*SigmaInv*D))%*%(InvFisher2%*%t(Z_aktuell*sqrt(D*SigmaInv*D)))
      }else{
        W_inv_t <- chol(D%*%(SigmaInv%*%t(D)))
        FinalHat<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
      }
      df<-sum(diag(FinalHat))
      
      if(control$overdispersion)
        phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-df)
      
      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }

      Delta_neu<-Delta[l,]
      Mu_opt<-Mu
      Qfinal<-Q[[l+1]]
      
      aaa<-!is.element(Delta_neu[1:(lin)],0)
      
      if(final.re)
      {    
        ############ final re-estimation
 
        if(rnd.len==1 && s==1)
        {  
          Q.max<-max(sqrt(unlist(Q)))
          Q.min<-min(sqrt(unlist(Q)))
        }else{
          Q.max<-max(Qfinal)+1
          Q.min<-min(Qfinal)-1e-10
        }
        
        if(rnd.len==1)
        {
          glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=Qfinal,K=K,
                                   Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                                   family=family,method=control$method.final,overdispersion=control$overdispersion,
                                   phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                   Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=q_start,K=K,
                                      Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                                      family=family,method=control$method.final,overdispersion=control$overdispersion,
                                      phi=control$phi,print.iter.final=control$print.iter.final,
                                      eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }
        }else{
          glmm_fin<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,K=K,
                                                Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                                                family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          if(class(glmm_fin)=="try-error"|| glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=q_start,K=K,
                                                   Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                                                   family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                   phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                   eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)    
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }
        }
        
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {
          cat("Warning:\n")
          cat("Final Fisher scoring reestimation did not converge!\n")
        }
        #######
        ranef.logLik<-glmm_fin$ranef.logLik
      }else{
        glmm_fin<-NA  
        class(glmm_fin)<-"try-error"
      }  
      Standard_errors<-matrix(NA,length(Delta_neu),length(Delta_neu))
      
      if(class(glmm_fin)!="try-error")
      {
        Delta_neu2<-Delta_neu
        Delta_neu2[c(aaa,rep(T,n%*%s))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,n%*%s)),c(aaa,rep(T,n%*%s))]<-glmm_fin$Standard_errors
        Qfinal<-glmm_fin$Q
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
        glmm_fin$ranef.logLik<-ranef.logLik
        complexity<-df
      }

      if(rnd.len==1)
      {
        if(s==1)
          Qfinal<-sqrt(Qfinal)
        
        if(!is.matrix(Qfinal))
          Qfinal<-as.matrix(Qfinal)
        colnames(Qfinal)<-random.labels
        rownames(Qfinal)<-random.labels
      }else{
        Qfinal_old<-Qfinal
        Qfinal<-list()
        Qfinal[[1]]<-as.matrix(Qfinal_old[1:s[1],1:s[1]])
        colnames(Qfinal[[1]])<-random.labels[[1]]
        rownames(Qfinal[[1]])<-random.labels[[1]]
        
        if(s[1]==1)
          Qfinal[[1]]<-sqrt(Qfinal[[1]])
        
        for (zu in 2:rnd.len)
        {
          Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
          colnames(Qfinal[[zu]])<-random.labels[[zu]]
          rownames(Qfinal[[zu]])<-random.labels[[zu]]
          
          if(s[zu]==1)
            Qfinal[[zu]]<-sqrt(Qfinal[[zu]])
          
        }
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
      
      names(Delta_neu)[(lin+1):(lin+n%*%s)]<-colnames(W)
      colnames(Standard_errors)[(lin+1):(lin+n%*%s)]<-rownames(Standard_errors)[(lin+1):(lin+n%*%s)]<-colnames(W)
      colnames(Delta)<-c(final.names,colnames(W))
      
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
        
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
        
        if(control$complexity!="hat.matrix")  
        {  
          if(rnd.len==1)
          {
            complexity<-0.5*(s*(s+1))
          }else{
            complexity<-0.5*(s[1]*(s[1]+1))
            for(zu in 2:rnd.len)
              complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
          }
          complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)
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
      ret.obj$ranef<-Delta_neu[(lin+1):(lin+n%*%s)]
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$ranerror<-Standard_errors[(lin+1):(lin+n%*%s)]
      ret.obj$Q_long<-Q
      ret.obj$Q<-Qfinal
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$conv.step<-conv.step
      ret.obj$newrndfrml<-newrndfrml
      ret.obj$subject<-subject.names
      ret.obj$data<-data
      ret.obj$rnd.len<-rnd.len
      ret.obj$phi.med<-phi.med
      ret.obj$y <- y
      ret.obj$X <- cbind(X,U)
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
      #    ret.obj$logLik.test<-logLik.test
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
        Eta_start<-X%*%beta_null+Phi%*%smooth_null+W%*%ranef_null
      }else{
        Eta_start<-rep(beta_null,N)+Phi%*%smooth_null+W%*%ranef_null
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
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          Q_start<-diag(q_start,s)
        }else{
          Q_start<-q_start
        }
      }else{
        if(all(s==1))
        {
          Q_start<-diag(diag(q_start),sum(s))
        }else{
          Q_start<-matrix(0,sum(s),sum(s))
          Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,Phi,W)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,(lin+dim.smooth+n%*%s))
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null
      Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-t(ranef_null)

      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      active_old<-!is.element(Delta[1,],0)
      logLik.vec<-c()
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      Q_inv<-NULL
      Q_inv.old.temp<-NULL
      Q_inv.start<-NULL
      
      Q<-list()
      Q[[1]]<-Q_start
      
      l=1
      
      if(diff.ord>1)
      {
        k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
      }else{
        k2<-rep(1,nbasis)
      }
      k22<-rep(k2,m)
      penal.vec<-penal*k22
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin),penal.vec,rep(1/Q_start,n*s))
          P1<-diag(P1)
        }else{
          Q_inv.start<-chol2inv(chol(Q_start))
          P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
          diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
          for(j in 1:n)
            P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin),penal.vec,rep(diag(Q_start)^(-1),n))
          P1<-diag(P1)
        }else{
          Q_inv.start<-list()
          Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
          P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
          diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
          for(jf in 1:n[1])
            P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            for(jf in 1:n[zu])
              P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
          }
        }
      }
      
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
      
      score_vec[(q+1):lin] <- grad.1
      
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
#      ranef.logLik<- -0.5*t(Delta_start[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s),(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]%*%Delta_start[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta_start[1:(lin+dim.smooth)],ranef=Delta_start[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)],
                                         Grad=score_vec,family=family,P=diag(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]), 
                                         lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      nue<-control$nue
      
      half.index<-0
      
      solve.test2<-FALSE  
      while(!solve.test2)
      {  
        
        if(half.index>50)
        {
          stop("Fisher matrix not invertible")
        }
        
        Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
        if(t_opt>t_edge & half.index==0)
          Delta[1,crit.obj$whichmin+q]<-0  
        Eta<-Z_alles%*%Delta[1,]
        
        if(is.null(family$multivariate)){
          D<-family$mu.eta(Eta)
          Mu<-family$linkinv(Eta)
          SigmaInv <- 1/family$variance(Mu)
        }else{
          Eta_cat <- matrix(Eta_start, byrow = TRUE, ncol = K)
          Mu_cat <- family$linkinv(Eta_cat)
          D <- family$deriv.mat(Mu_cat)
          SigmaInv <- family$SigmaInv(Mu_cat)
          Mu <- c(t(Mu_cat))
        }          
        ranef.logLik<- -0.5*t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s),(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]%*%Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])
        
        logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if (control$method=="EM")
        {  
          if(rnd.len==1)
          {
            if(s==1)
            {
              P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q_start,n*s)))
            }else{
              P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
              diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
              for(jf in 1:n)
                P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.start
            }
          }else{
            if(all(s==1))
            {
              P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q_start)^(-1),n)))
            }else{
              P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
              diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
              for(jf in 1:n[1])
                P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.start[[1]]
              
              for (zu in 2:rnd.len)
              {
                for(jf in 1:n[zu])
                  P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                        (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
              }
            }
          }
          if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
          }else{
            F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
          }
          InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher2)=="try-error")
            InvFisher2<-try(solve(F_gross),silent=T)
          if(class(InvFisher2)=="try-error")
          {
            #stop("Fisher matrix not invertible")  
            half.index<-half.index+1  
          }else{
            solve.test2<-TRUE 
          }}else{
            solve.test2<-TRUE
          }
      }
      
      betaact<-Delta[1,active]
      
      if (control$method=="EM")
      {   
        ############################# Q update ################
        if(rnd.len==1)
        {
          Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)])
          for (i in 2:n)
            Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
          Q1<-1/n*Q1
        }else{
          Q1<-matrix(0,sum(s),sum(s))
          Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
          for (i in 2:n[1])
            Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
          Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
            for (i in 2:n[zu])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }else{
        if(is.null(family$multivariate)){
          Eta_tilde<-Eta+(y-Mu)*1/D
        }else{
          Eta_tilde<-Eta+solve(D)%*%(y-Mu)
        }
        Betadach<-Delta[1,1:(lin+dim.smooth)]     
        aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
        X_aktuell<-Z_fastalles[,aktuell_vec]
        
        if(rnd.len==1)
        {
          
          if(s==1)
          {
            upp<-min(20,50*Q_start)
            low<-1e-14
            optim.obj<-try(nlminb(sqrt(Q_start),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                  Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
            if(class(optim.obj)=="try-error")
              optim.obj<-try(bobyqa(sqrt(Q_start),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                    Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
            Q1<-as.matrix(optim.obj$par)^2
          }else{
            q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
            up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
            upp<-rep(up1,length(q_start_vec))
            low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
            #   kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
            optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                  Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
            Q1<-matrix(0,s,s)
            Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
            Q1<-Q1+t(Q1)
            diag(Q1)<-(optim.obj$par[1:s])
            
            #### Check for positive definitness ########
            for (ttt in 0:100)
            {
              Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
              Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
              Q_solvetest<-try(solve(Q1))
              if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                break
            }
          }   
        }else{
          if(all(s==1))
          {
            q_start_vec<-diag(q_start)
            upp<-rep(min(20,50*diag(q_start)),sum(s))
            low<-rep(1e-14,sum(s))
            optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
            Q1<-diag(optim.obj$par)^2
          }else{
            q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
            up1<-min(20,50*max(q_start_vec))
            low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))
            
            for (zu in 2:rnd.len)
            {
              q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
              up1<-min(20,50*max(q_start_vec))
              low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
            }
            upp<-rep(up1,length(q_start_vec))
            optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
            optim.vec<-optim.obj$par
            
            Q1<-matrix(0,sum(s),sum(s))
            diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
            if(s[1]>1)
              Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
            optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
            
            for (zu in 2:rnd.len)
            {
              diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
              if(s[zu]>1)
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
              optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
            }
            
            #### Check for positive definitness ########
            for (ttt in 0:100)
            {
              Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
              Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
              Q_solvetest<-try(solve(Q1))
              if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                break
            }
          }
        }}
      
      Q[[2]]<-Q1
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin),penal.vec,rep(1/Q1,n*s))
          P1<-diag(P1)
        }else{
          Q_inv<-chol2inv(chol(Q1))
          Q_inv.old.temp<-Q_inv
          P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
          diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
          for(j in 1:n)
            P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
          P1<-diag(P1)
        }else{
          Q_inv<-list()
          Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
          P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
          diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
          for(jf in 1:n[1])
            P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
          
          for (zu in 2:rnd.len)
          {
            Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            for(jf in 1:n[zu])
              P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
          }
        }
      }
      
      if(is.null(family$multivariate)){
        score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
      }else{
        score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
      }  
      lambda.max<-max(abs(score_vec2[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      score_vec2[(q+1):lin] <- grad.1

      crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate

      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[1,1:(lin+dim.smooth)],ranef=Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)],
                                         Grad=score_vec2,family=family,P=diag(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]), 
                                         lower = 0, upper = Inf,K=K))
      
      t_opt<-optim.obj$par
      
      Eta.ma[2,]<-Eta
      
      score_vec<-score_vec2
      Q1.old<-Q1
      Q_inv.old<-Q_inv    
      
      Q1.very.old<-Q_start
      Q_inv.very.old<-Q_inv.start
      
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
          half.index<-0
          
          solve.test2<-FALSE  
          while(!solve.test2)
          {  
            
            if(half.index>50)
            {
              half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old;
            }
            
            Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
            if(t_opt>t_edge & half.index==0)
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
            
            ranef.logLik<- -0.5*t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s),(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]%*%Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])
            
            logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=ranef.logLik,family=family,penal=T,K=K)
            
            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
            if(control$method=="EM")
            {  
              if(rnd.len==1)
              {
                if(s==1)
                {
                  P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q1.old,n*s)))
                }else{
                  P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                  diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                  for(jf in 1:n)
                    P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                          (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.old
                }
              }else{
                if(all(s==1))
                {
                  P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q1.old)^(-1),n)))
                }else{
                  P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                  diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                  for(jf in 1:n[1])
                    P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                          (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.old[[1]]
                  
                  for (zu in 2:rnd.len)
                  {
                    for(jf in 1:n[zu])
                      P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                            (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
                  }
                }
              }
              if(is.null(family$multivariate)){
                D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
                F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
              }else{
                F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
              }
              InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
              if(class(InvFisher2)=="try-error")
                InvFisher2<-try(solve(F_gross),silent=T)
              if(class(InvFisher2)=="try-error")
              {
                half.index<-half.index+1  
              }else{
                solve.test2<-TRUE 
              }}else{
                solve.test2<-TRUE 
              }
          }
          
          betaact<-Delta[l,active]
          
          if (control$method=="EM")
          {        
            ############################# Q update ################
            if(rnd.len==1)
            {
              Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)])
              for (i in 2:n)
                Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
              Q1<-1/n*Q1
            }else{
              Q1<-matrix(0,sum(s),sum(s))
              Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
              for (i in 2:n[1])
                Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
              Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
              
              for (zu in 2:rnd.len)
              {
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
                for (i in 2:n[zu])
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
              }
            }  
          }else{
            if(is.null(family$multivariate)){
              Eta_tilde<-Eta+(y-Mu)*1/D
            }else{
              Eta_tilde<-Eta+solve(D)%*%(y-Mu)
            }
            
            Betadach<-Delta[l,1:(lin+dim.smooth)]
            
            aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
            X_aktuell<-Z_fastalles[,aktuell_vec]
            
            if(rnd.len==1)
            {
              
              if(s==1)
              {
                if(Q1<1e-14)
                  low<-0
                
                optim.obj<-try(nlminb(sqrt(Q1),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
                if(class(optim.obj)=="try-error")
                  optim.obj<-bobyqa(sqrt(Q1),likelihood_nlminb,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
                
                Q1<-as.matrix(optim.obj$par)^2
              }else{
                Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
                optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
                
                Q1<-matrix(0,s,s)
                Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
                Q1<-Q1+t(Q1)
                diag(Q1)<-(optim.obj$par[1:s])
                
                #### Check for positiv definitness ########
                for (ttt in 0:100)
                {
                  Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                  Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                  Q_solvetest<-try(solve(Q1))
                  if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                    break
                }
              }   
            }else{
              if(all(s==1))
              {
                Q1_vec<-diag(Q1)
                optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                Q1<-diag(optim.obj$par)^2
              }else{
                Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
                
                for (zu in 2:rnd.len)
                  Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                
                optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,SigmaInv=SigmaInv,family=family,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                optim.vec<-optim.obj$par
                
                Q1<-matrix(0,sum(s),sum(s))
                diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
                if(s[1]>1)
                  Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
                optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
                
                for (zu in 2:rnd.len)
                {
                  diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
                  if(s[zu]>1)
                    Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
                  optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
                }
                
                #### Check for positive definitness ########
                for (ttt in 0:100)
                {
                  Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
                  Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
                  Q_solvetest<-try(solve(Q1))
                  if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
                    break
                }
              }
            }}
          
          Q[[l+1]]<-Q1
          
          if(rnd.len==1)
          {
            if(s==1)
            {
              P1<-c(rep(0,lin),penal.vec,rep(1/Q1,n*s))
              P1<-diag(P1)
            }else{
              Q_inv<-chol2inv(chol(Q1))
              P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
              diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
              for(j in 1:n)
                P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
            }
          }else{
            if(all(s==1))
            {
              P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
              P1<-diag(P1)
            }else{
              Q_inv<-list()
              Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
              P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
              diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
              for(jf in 1:n[1])
                P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
              
              for (zu in 2:rnd.len)
              {
                Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
                for(jf in 1:n[zu])
                  P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                     (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              }
            }
          }
          
          if(is.null(family$multivariate)){
            score_vec2<-t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l,]
          }else{
            score_vec2<-t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l,]
          }  
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
          }
          score_vec2[(q+1):lin] <- grad.1
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
          t_edge<-crit.obj$min.rate

          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt,y=y,X=Z_alles,fixef=Delta[l,1:(lin+dim.smooth)],ranef=Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)],
                                             Grad=score_vec2,family=family,P=diag(P1[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]), 
                                             lower = 0, upper = Inf,K=K))
          
          t_opt<-optim.obj$par
          
          score_vec<-score_vec2
          Q1.very.old<-Q1.old
          Q_inv.very.old<-Q_inv.old
          Q1.old<-Q1
          Q_inv.old<-Q_inv    
          
          Eta.ma[l+1,]<-Eta
          
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
        }}
      
      ######## Final calculation
      if(control$method!="EM")
      {  
        if(rnd.len==1)
        {
          if(s==1)
          {
            P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(1/Q1.old,n*s)))
          }else{
            P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
            for(jf in 1:n)
              P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                    (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.old
          }
        }else{
          if(all(s==1))
          {
            P_akt<-diag(c(rep(0,lin_akt),penal.vec,rep(diag(Q1.old)^(-1),n)))
          }else{
            P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
            for(jf in 1:n[1])
              P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                    (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.old[[1]]
            
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
                P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                      (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
            }
          }
        }
        if(is.null(family$multivariate)){
          D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P_akt
        }else{
          F_gross<-t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D))%*%Z_aktuell)+P_akt
        }
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
      }          
      
      if(is.null(family$multivariate)){
        FinalHat<-(Z_aktuell*sqrt(D*SigmaInv*D))%*%(InvFisher2%*%t(Z_aktuell*sqrt(D*SigmaInv*D)))
      }else{
        W_inv_t <- chol(D%*%(SigmaInv%*%t(D)))
        FinalHat<-W_inv_t%*%(Z_aktuell%*%(InvFisher2%*%(t(Z_aktuell)%*%t(W_inv_t))))
      }
      df<-sum(diag(FinalHat))
      
      if(control$overdispersion)
        phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-df)
      
      conv.step<-l
      phi.med<-phi
      
      if(conv.step==control$steps)
      {
        cat("Warning:\n")
        cat("Algorithm did not converge!\n")
      }
      
      Delta_neu<-Delta[l,]
      Eta_opt<-Z_alles%*%Delta_neu
      Mu_opt<-Mu
      Qfinal<-Q[[l+1]]
      
      aaa<-!is.element(Delta_neu[1:(lin)],0)
      
      if(final.re)
      {    
        ############ final re-estimation
        
        if(s==1)
        {  
          Q.max<-max(sqrt(unlist(Q)))
          Q.min<-min(sqrt(unlist(Q)))
        }else{
          Q.max<-max(Qfinal)+1
          Q.min<-min(Qfinal)-1e-10
        }
        
        if(rnd.len==1)
        {
          glmm_fin<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,n,penal.vec,q_start=Qfinal,K=K,
                                          Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                                          s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                          phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                          Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,n,penal.vec,q_start=q_start,K=K,
                                             Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                                             s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                             phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                             Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
          }
        }else{
          glmm_fin<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=Qfinal,K=K,
                                                       Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                                                       s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                       phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                                       Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=q_start,K=K,
                                                          Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                                                          s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                          phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                          eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
            }
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
        EDF.matrix<-glmm_fin$EDF.matrix
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])  
        Delta_neu2<-Delta_neu
        Delta_neu2[c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,dim.smooth+n%*%s)),c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Standard_errors
        Qfinal<-glmm_fin$Q
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
        glmm_fin$ranef.logLik<-ranef.logLik
        complexity<-df
        
        lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))
        if(rnd.len==1)
        {
          if(s==1)
          {
            P1a<-diag(c(rep(0,lin_akt+dim.smooth),rep(1/Q1,n*s)))
          }else{
            P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            for(jf in 1:n)
            {
              P1a[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                  (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
            }  
          }
        }else{
          if(all(s==1))
          {
            P1a<-diag(c(rep(0,lin_akt+dim.smooth),rep(diag(Q1)^(-1),n)))
          }else{
            P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            for(jf in 1:n[1])
            {
              P1a[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                  (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
            }
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
              {
                P1a[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              }}
          }
        }
        if(class(InvFisher2)=="try-error")
        {
          warning("No EDF's for smooth functions available, as Fisher matrix not invertible!")
          complexity.smooth<-dim.smooth
        }else{  
          ###### EDF of spline; compare Wood's Book on page 167
          if(is.null(family$multivariate)){
            D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
            EDF.matrix<-InvFisher2%*%(t(Z_aktuell)%*%(Z_aktuell*D*SigmaInv*D)+P1a)
          }else{
            EDF.matrix<-InvFisher2%*%(t(Z_aktuell)%*%(D%*%(SigmaInv%*%t(D)))%*%Z_aktuell+P1a)
          }
          complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])
        }
      }  
      
      if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
        complexity.smooth<-dim.smooth
 
      if(rnd.len==1)
      {
        if(s==1)
          Qfinal<-sqrt(Qfinal)
        
        if(!is.matrix(Qfinal))
          Qfinal<-as.matrix(Qfinal)
        colnames(Qfinal)<-random.labels
        rownames(Qfinal)<-random.labels
      }else{
        Qfinal_old<-Qfinal
        Qfinal<-list()
        Qfinal[[1]]<-as.matrix(Qfinal_old[1:s[1],1:s[1]])
        colnames(Qfinal[[1]])<-random.labels[[1]]
        rownames(Qfinal[[1]])<-random.labels[[1]]
        
        for (zu in 2:rnd.len)
        {
          Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
          colnames(Qfinal[[zu]])<-random.labels[[zu]]
          rownames(Qfinal[[zu]])<-random.labels[[zu]]
        }
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
      
      names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-rownames(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      names(Delta_neu)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
      colnames(Standard_errors)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(Standard_errors)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
      colnames(Delta)<-c(final.names,colnames(Phi),colnames(W))
      
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
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=glmm_fin$ranef.logLik,family=family,penal=T,K=K)
        
        if(control$complexity!="hat.matrix")  
        {  
          if(rnd.len==1)
          {
            complexity<-0.5*(s*(s+1))
          }else{
            complexity<-0.5*(s[1]*(s[1]+1))
            for(zu in 2:rnd.len)
              complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
          }
          complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
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
      ret.obj$ranef<-Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$ranerror<-Standard_errors[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
      ret.obj$Q_long<-Q
      ret.obj$Q<-Qfinal
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$newrndfrml<-newrndfrml
      ret.obj$subject<-names(rnd)
      ret.obj$data<-data
      ret.obj$rnd.len<-rnd.len
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
      ret.obj$rnd <- rnd
      return(ret.obj)
    }  
  }
  ##################################################################
  ##################################################################
}
