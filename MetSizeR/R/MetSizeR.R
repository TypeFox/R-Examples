
MetSizeR <-
function()
{
  options(guiToolkit="RGtk2") 
  #################################################################################### 
  # function to scale covariates to lie between 0 and 1 for stability                                                               #
  ####################################################################################
  standardize<-function(C)
  {
    for(i in 1:ncol(C))
    {
      rg<-range(C[,i])
      C[,i]<-(C[,i]-min(C[,i]))/(rg[2] - rg[1])
    } # end for
    C
  }
  
  #################################################################################### 
  # function for simulating data using pilot data
  ####################################################################################
  sim.pilot.data<-function(n1,n2,p,Zerop,Ip,Zeroq,Iq,ZeroL,IL,ZeroL1,IL1,eta.sd,eta_sc.sd,sig,W,Alpha,mu,mod,model)
  {
    n<-n1+n2
    if((mod[1]==model)|(mod[2]==model))
    {
      if(mod[1]==model)
      { 
        u<-rmvnorm(n,Zeroq,Iq)
      }else{
        C<-rmvnorm(n,ZeroL,IL)
        C<-standardize(C)  		        ## Standardize covariates for stability
        C<-rbind(rep(1,n), t(C))
        u<-rmvnorm(n,Zeroq,Iq)+t(Alpha%*%C)
      }#ifppca
      x<-rmvnorm(n,Zerop,sig*Ip)+tcrossprod(u,W)+ matrix(mu, n, p, byrow=TRUE)
    }else{
      ## SV model on the errors     
      eta.true<-rnorm(1,0,eta.sd)
      
      ## SV model on the scores
      eta_sc.true<-c(rmvnorm(1, Zeroq, eta_sc.sd))
      
      ## DPPCA model
      u<-rmvnorm(n,Zeroq,exp(eta_sc.true)*Iq)
      x<-rmvnorm(n,Zerop,exp(eta.true)*Ip)+tcrossprod(u,W)
    }#long
    return(x)
  }
  
  #################################################################################### 
  # function for simulating data without pilot data
  ####################################################################################
  sim.pilot<-function(n1,n2,p,Zerop,Ip,q,Zeroq,Iq,ZeroL,IL,ZeroL1,IL1,eta.sd,eta_sc.sd,alpha.sigma,beta.sigma,mod,model,ao1,bo)
  {
    n<-n1+n2
    if((mod[1]==model)|(mod[2]==model))
    {
      sig<-1/rgamma(1,alpha.sigma,beta.sigma)
      if(mod[1]==model)
      {
        u<-rmvnorm(n,Zeroq,Iq)
      }else{
        Alpha<-rmvnorm(q,ZeroL1,3*IL1)
        C<-rmvnorm(n,ZeroL,IL)
        C<-standardize(C)			        ## Standardize covariates for stability
        C<-rbind(rep(1, n), t(C))
        u<-rmvnorm(n,Zeroq,Iq)+t(Alpha%*%C)
      }#ifppca
      v<-1/rgamma(q,ao1,bo)
      W<-rmvnorm(p,Zeroq,v*Iq)
      x<-rmvnorm(n,Zerop,sig*Ip)+tcrossprod(u,W)
    }else{
      
      ## SV model on the errors
      eta.true<-rnorm(1,0,eta.sd)
      
      ## SV model on the scores
      eta_sc.true<-c(rmvnorm(1,Zeroq,eta_sc.sd))         
      
      ## DPPCA model
      v<-1/rgamma(q,ao1,bo)
      u<-rmvnorm(n,Zeroq,exp(eta_sc.true)*Iq)
      W<-rmvnorm(p,Zeroq,v*Iq)
      x<-rmvnorm(n,Zerop,exp(eta.true)*Ip)+tcrossprod(u,W)
    }#long
    return(x)
  }
  
  #################################################################################### 
  # function for sampling from the null distribution 
  ####################################################################################
  samp.dist<-function(T,S,TS,x,y,n11,n22,in1n2,nn2,cpcf)
  {
    for(t in 1:T)
    {
      yperm<-sample(y, replace=FALSE)
      x1<-x[yperm==1,]
      x2<-x[yperm==2,]
      Sj<-sqrt(in1n2*(n11*diag(var(x1))+n22*diag(var(x2)))/nn2)  # FASTER!!!
      S[,t]<-Sj+sort(Sj)[cpcf]                     ## corrected standard deviation
      TS[,t]<-(colMeans(x1)-colMeans(x2))/S[,t]
    }#t
    return(list(S=S, TS=TS))
  }
  
  #################################################################################### 
  # function to see if a number is an integer
  #################################################################################### 
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  #################################################################################### 
  # function for estimating sample size 
  ####################################################################################
  metsize<-function(pilot=pilot, n1=4, n2=4, p=200, prop=0.25, covars=covars, ncovar=0, model="PPCA", plot.prop=FALSE, target.fdr= 0.05, Targeted=FALSE)
  { 
    mod<-c("PPCA", "PPCCA", "DPPCA")
    if(prop>0.5) 
    {
      cat("Warning! The proportion of significant metabolites is usually expected to be less than 0.5.")
    }
    if(length(pilot)!=0)
    {
      if(mod[2]==model)
      {
        ncovar<-ncol(covars)
        resppcca<-ppcca.metabol(pilot, covars, 2, 2, scale="unit")     
        W<-resppcca$loadings
        sig<-resppcca$sig
        Alpha<-resppcca$coefficients 
        mu<-colMeans(pilot)
        p<-ncol(pilot)
      }else{
        pilot<-sweep(pilot, 2, colMeans(pilot), "-")
        resppca<-ppca.metabol(pilot, 2, 2, scale="unit") 
        W<-resppca$loadings
        sig<-resppca$sig
        mu<-colMeans(pilot)
        p<-ncol(pilot)
      }#ifcovars 
    }#ifpilot
    if(plot.prop==FALSE)
    {
      mprop<-c(prop)                       ## proportion of significant metabolites
      mrange<-ceiling(mprop*p)                                 ## different number significant metabolites
      mtry <- length(mrange)                         ## number of different m values
      nprop.increment <- 1		           ## Sample size increments
      sfactors<-(n1+n1)*nprop.increment  
    }else{
      mprop<-c(prop,seq(0.1,0.5,0.1))                       ## proportion of significant metabolites
      mrange<-ceiling(mprop*p)                                 ## different number significant metabo
      mtry <- length(mrange)                         ## number of different m values  
      nprop.increment <- c(1,2,4,7)		           ## Sample size increments
      sfactors<-(n1+n2)*nprop.increment                       ## four different sample sizes to be considered                        
    }
    
    ## Setting up the initial values and the prior parameters
    T<-20                                         ## number of permutations to estimate the sampling distribution.
    Sim<-20  				              ## Number of pilot data sets simulated.
    L<-ncovar                                     ## number of covariates
    n<-n1+n2                                      ## sample size for pilot data
    q<-2                                          ## dimension of the latent space
    alpha.sigma<-ao<-5                                ## the scale parameter of the prior distribution of the variance.
    beta.sigma<-2*(alpha.sigma-1)                 ## the shape parameter of the prior distribution of the variance. 
    bo<-c(0.5*(ao-1),0.25*(ao-1))
    ao1<-rep(ao,q)
    n.increment = c(1:8)		              ## Sample size increments considered.
    ntry = length(n.increment)                    ## number of different samples sizes to be considered.
    
    ## Storing statistics for each sample size
    fdr50_sim = matrix(NA, ntry, 1)                    ## FDR values for each sample simulated from the PPCCA model
    fdr90_sim = matrix(NA, ntry, 1)
    fdr10_sim = matrix(NA, ntry, 1)
    fdr_sim = matrix(NA, ntry, Sim)  
    
    ## Storing statistics for each m values
    fdr50_prop = matrix(NA, mtry-1, length(sfactors))
    fdr90_prop = matrix(NA, mtry-1, length(sfactors))
    fdr10_prop = matrix(NA, mtry-1, length(sfactors))
    Add<-array(NA, c(p,T,ntry))
    TSstore<-array(NA, c(p,T,ntry))
    
    delta<-qnorm(0.99)   # Based on variance of underlying simulation model
    if(length(pilot)!=0){delta<-qnorm(0.89)}
    cf<-0.05
    Zeroq<-rep(0,q) 
    Iq<-diag(q)
    Zerop<-rep(0,p)
    Ip<-diag(p)
    if(mod[2]==model){ZeroL<-rep(0,L); IL<-diag(L); ZeroL1<-rep(0,(L+1)); IL1<-diag(L+1)}
    if(mod[3]==model){ phi.true<-0.8; v2.true<-0.1; eta.sd<-sqrt(v2.true/(1-phi.true^2))
                       phi_sc.true<-rep(0.8,q); v2_sc.true<-rep(0.1,q)
                       eta_sc.sd<-diag(sqrt(v2_sc.true/(1-phi_sc.true^2)))}
    TS<-S<-matrix(NA,p,T) 
    cpcf<-ceiling(p*cf)
    for(m in 1:mtry)
    {
      pstat<-1-(mrange[m]/p)
      ind <- matrix(FALSE, p, T)
      pos <- sample(1:p, size=mrange[m])         ## sampling mrange[m] metabolites
      ind[pos,] <- TRUE               ## matrix indicating truly significant and non significant metabolites
      i<-0
      for(k in 1:ntry)
      {
        n1star<-n1*n.increment[k]    ## sample size for treatment group 1
        n2star<-n2*n.increment[k]    ## sample size for treatment group 2 
        nstar<-n1star+n2star
        in1n2<-1/n1star+1/n2star; n11<-n1star-1; n22<-n2star-1; nn2<-n1star+n2star-2
        Add.sd<-sqrt(in1n2)
        
        y<-c(rep(1,n1star),rep(2,n2star))
        if(m==1)               ## If influence of different proportion of significant metabolites is not of interest...  
        {
          for(s in 1:Sim)            ## assessing the effect of repeated simulations from the underlying model.
          {
            ## Simulating the pilot data
            if(length(pilot)!=0)
            {
              x<-sim.pilot.data(n1star,n2star,p,Zerop,Ip,Zeroq,Iq,ZeroL,IL,ZeroL1,IL1,eta.sd,eta_sc.sd,sig,W,Alpha,mu,mod,model)
            }else{
              x<-sim.pilot(n1star,n2star,p,Zerop,Ip,q,Zeroq,Iq,ZeroL,IL,ZeroL1,IL1,eta.sd,eta_sc.sd,alpha.sigma,beta.sigma,mod,model,ao1,bo)
            }
            
            ## Estimating the sampling distribution of the test statistic using permutations.
            res.sampdist<-samp.dist(T,S,TS,x,y,n11,n22,in1n2,nn2,cpcf)
            TS<-res.sampdist$TS
            S<-res.sampdist$S
            
            ## Store the test statistics for the current sample size.
            TSstore[,,k]<-TS
            
            ## calculating the shift in metabolites in grp 2
            vars<-S/Add.sd
            Add[,,k]<-delta/(vars*Add.sd)
            
            ## estimating the FDR
            tsB<-TS
            tsB[pos,] <- tsB[pos,] + Add[pos,,k]            ## Add an increment factor to the t.statistic values of the truly significant metabolites 
            atsB <- abs(tsB) 
            crit <- quantile(atsB, pstat)         ## identifying a cut-off point (m-th largest absolute value of the p TsB values)     
            errors <- (colSums(atsB > crit & !ind))/(colSums(atsB > crit))        ## number of false positives 
            fdr_sim[k,s] <- quantile(errors[!is.na(errors)], 0.5)      ## median FDR of the Tstar permutations
          }#s
          
          ## assessing the effect of repeated simulations from the underlying model.
          emp<-quantile(fdr_sim[k,], c(0.1,0.5,0.9))
          fdr10_sim[k] <- emp[1]
          fdr50_sim[k] <- emp[2]
          fdr90_sim[k] <- emp[3]
        }else{
          ## assessing the effect of varying the proportion of truly significant metabolites (increasing m values) on four different sample sizes.
          if(any(nprop.increment==n.increment[k]))
          {
            i<-i+1
            ## estimating the FDR
            tsB<-TSstore[,,k]
            tsB[pos,] <- tsB[pos,] + Add[pos,,k]            ## adding a shift component to truly significant metabolites
            atsB <- abs(tsB)
            crit <- quantile(atsB, 1 - (mrange[m]/p))
            errors <- (colSums(atsB > crit & !ind ))/(colSums(atsB > crit)) ## (number of false positives for Tstar permutations)/(number of metabolites declared significant for Tstar permutations)
            emp<-quantile(errors[!is.na(errors)], c(0.1,0.5,0.9))
            fdr10_prop[m-1,i] <- emp[1]
            fdr50_prop[m-1,i] <- emp[2]          ## median FDR for Tstar permutations permutations
            fdr90_prop[m-1,i] = emp[3]
          }#if(any)     
        }#if
      }#k
    }#m
    
    results_sim <- cbind(n.increment*n, fdr50_sim, fdr90_sim, fdr10_sim)
    results_prop <- cbind(mprop[-1], fdr50_prop, fdr10_prop, fdr90_prop)
    colnames(results_sim) <- c("sample size","fdr50_sim", "fdr90_sim", "fdr10_sim")
    colnames(results_prop) <- c("prop", paste("fdr50_prop", sfactors, sep = "_"), paste("fdr10_prop", sfactors, sep = "_"), paste("fdr90_prop", sfactors, sep = "_"))
    results_prop <- list(results_prop, sfactors)
    names(results_prop) <- c("results_prop", "sample_sizes")
    
    ######################################### Determine the sample size at which the FDR line is equal to 0.05
    opty <- rep(0, 2)
    ind1 <- min(c(1:nrow(results_sim))[results_sim[,2]<0.05])
    ind2 <- max(c(1:ind1)[results_sim[1:ind1,2]>0.05])
    opty <- c(results_sim[ind1,2], results_sim[ind2,2])
    optx <- c(results_sim[ind1,1], results_sim[ind2,1])
    optres<-lm(opty~optx)
    nhat <- round((target.fdr - optres$coef[1])/optres$coef[2])
    if(is.na(nhat)){print("Sorry: an error has occurred. Please rerun the function."); stop()}
    if(n1==n2)
    {
      if(is.wholenumber(nhat/2)){n1<-n2<-nhat/2}else{nhat<- nhat +1; n1<-n2<-nhat/2}
    }else{
      if(is.wholenumber(n1*nhat/(n1+n2))){n1 <- n1*nhat/(n1+n2); n2 <- nhat-n1
      }else{
        n1.user <- n1
        n2.user <- n2
        n1 <- ceiling(n1.user*nhat/(n1.user+n2.user)) 
        n2 <- ceiling(n2.user*nhat/(n1.user+n2.user))
        nhat <- n1 + n2
      }
    }
    est<-c(nhat, n1, n2)
    names(est)<-c("n","n1","n2")
    
    ############################################## Plotting the sample size estimation results
    if(plot.prop == FALSE)
    {
      par(mfrow=c(1,1))
      ## Plot to illustrate variability due to simulation of the pilot data
      plot(results_sim[, "sample size"], results_sim[, "fdr50_sim"], xlab="Sample size", ylab="FDR", pch=16, col=2, ylim = c(0, 1))
      lines(results_sim[, "sample size"], results_sim[, "fdr50_sim"], col = 2, lwd=2)
      lines(results_sim[, "sample size"], results_sim[, "fdr90_sim"], col = 2, lty=2, lwd=2)
      lines(results_sim[, "sample size"], results_sim[, "fdr10_sim"], col = 2, lty=2, lwd=2)
      abline(h = target.fdr, lty = 3, col=1)
      legend("topleft", bty="n", paste("FDR = ", target.fdr), col=1, lty=3)
      arrows(nhat, 0.8, nhat,-0.03, length = 0.1,lwd=2, col=3)
      title(paste("Sample size estimation"), font=1)
      text(nhat+8, 0.89, labels=substitute(hat(n)==nhat, list(nhat=nhat)), col=3, cex=1.2)
      text(nhat+8, 0.85, labels=substitute((list(n[1]==n1,n[2]==n2)), list(n1=n1,n2=n2)), col=3, cex=1.2)
    }else{
      par(mfrow=c(2,2), mar=c(4,3,1.5,0.5), oma=c(0,0,2,0),mgp=c(2,1,0)) 
      for(k in 1:length(sfactors)) 
      {
        ## variability due to different number of k values
        if(Targeted){plot(results_prop$results_prop[, "prop"], results_prop$results_prop[, paste("fdr50_prop",sfactors[k],sep="_")], xlab="Proportion of significant metabolites", ylab="FDR", pch=16, col=2, ylim = c(0, 1))}
        else{plot(results_prop$results_prop[, "prop"], results_prop$results_prop[, paste("fdr50_prop",sfactors[k],sep="_")], xlab="Proportion of significant bins", ylab="FDR", pch=16, col=2, ylim = c(0, 1))}
        lines(results_prop$results_prop[, "prop"], results_prop$results_prop[, paste("fdr50_prop",sfactors[k],sep="_")], col = 2, lwd=2)
        lines(results_prop$results_prop[, "prop"], results_prop$results_prop[, paste("fdr90_prop",sfactors[k],sep="_")], col = 2, lty=2, lwd=2)
        lines(results_prop$results_prop[, "prop"], results_prop$results_prop[, paste("fdr10_prop",sfactors[k],sep="_")], col = 2, lty=2, lwd=2)
        abline(h = target.fdr, lty = 3, col=1)
        legend("topleft", bty="n", paste("FDR = ", target.fdr), col=1, lty=3)
        title(paste("Sample size=", sfactors[k]), cex = 0.7)
      }
      if(Targeted){title("Varying the proportion of significant metabolites", outer = T)}
      else{title("Varying the proportion of significant bins", outer = T)}
    }
    return(list(nhat=est, results_sim = results_sim, results_prop = results_prop, p=p, prop=prop, ncovar=ncovar, model=model, n1=n1, n2=n2, target.fdr = target.fdr))
  }#End of metsize function
  
  
  #################################################################################### 
  # function to find data and covariates files, read in and quit                                    #
  ####################################################################################
  pilot<- NULL                  # Pilot dataset
  covars<- NULL                 # covariates
  pilot.open<- function(h,...)
  {
    datapath<- gfile("Select the data file", filter = list("txt" = list(patterns = c("*.txt","*.TXT"))), type="open")
    pilot<<- read.table(datapath,sep="\t",header=TRUE)
  }
  cov.open<- function(h,...)
  {
    covpath<- gfile("Select the covariates file",filter = list("txt" = list(patterns = c("*.txt","*.TXT"))),type="open")
    covars<<- as.matrix(read.table(covpath,sep="\t",header=TRUE))
  }
  demo.pilot.open<- function(h,...)
  {
    pilot<<- read.table(paste(file.path(path.package(package="MetSizeR")[1]),"/extdata/nmr_spectra.txt",sep=""),sep="\t",header=TRUE)[,-1]
    covars<<- as.matrix(read.table(paste(file.path(path.package(package="MetSizeR")[1]),"/extdata/nmr_spectra.txt",sep=""),sep="\t",header=TRUE)[,1])
  }

  ####################################################################################
  # GUI for sample size calculation with/without pilot dataset                       #
  ####################################################################################
  ssize.est<- function(inputData, inputCov, Targeted=FALSE)
  {
    pilot<-inputData
    covars<-inputCov
    
    # GUI handler
    up.size<- function(h,...)
    {
      mod=svalue(models)
      if(((length(pilot)!=0)&(length(covars)==0)&(mod=="PPCCA")))
      {
        gmessage("Please upload the covariates of the pilot data first",icon="info")
      }else{
        
        n1=as.numeric(svalue(ng1)); n2=as.numeric(svalue(ng2))
        if(is.na(n1)|is.na(n2)){n1=n2=4}
        if((n1>15)|(n2>15)|(n1<3)|(n2<3))
        {  
          if((n1<3)|(n2<3))
          {
            gmessage("n1 or n2 values can not be less than three",icon="info")
          }else{
            gmessage("n1 or n2 values can not be greater than 15",icon="info")
          }
        }else{
          enabled(propsig)<-enabled(Mod)<-enabled(tfdr)<-enabled(Sz)<-enabled(gch)<-enabled(cal)<-FALSE
          if(length(pilot)!=0){ncov=0}else{ncov=svalue(ncov); enabled(spec)<-FALSE}
          enabled(fin)<-TRUE
          svalue(sampleSize)<- "is running please wait..."
          p=svalue(nspectralBins)
          prop=svalue(sigMet)
          t.fdr=svalue(targetfdr)
          if(svalue(dis)){ ggraphics(dpi = 100, ps = 9)}
          if(Targeted){Targ=TRUE}else{Targ=FALSE}
          Sav<-FALSE 
          if(Sav<-svalue(dis))
          {
            pdf(paste(getwd(),"/MetSizeRplot.pdf",sep="")) 
            metsize(pilot=pilot, n1=n1, n2=n2, p=p, prop=prop, covars=covars, ncovar=ncov, model=mod, plot.prop=svalue(gch), target.fdr= t.fdr, Targeted=Targ)
            dev.off() 
          }else{metsize(pilot=pilot, n1=n1, n2=n2, p=p, prop=prop, covars=covars, ncovar=ncov, model=mod, plot.prop=svalue(gch), target.fdr= t.fdr, Targeted=Targ)}
          enabled(fin)<-FALSE
          if(length(pilot)==0){enabled(spec)<-TRUE}
          enabled(propsig)<-enabled(Mod)<-enabled(tfdr)<-enabled(Sz)<-enabled(gch)<-enabled(cal)<-TRUE
          svalue(sampleSize)<- "has finished!!!"
        }#if n
      }#if mod
    }# end up.size function
    
    # GUI  components
    availModels <- c("PPCA", "PPCCA", "DPPCA")
    nspectralBins <- gradio(c(50,100,200,300,500,1000), selected=3)
    nMet <- gradio(c(20,30,50,70,100,200), selected=3)
    sigMet <- gspinbutton(from=0.1, to = 1, by = 0.1, value=0.2)
    models <- gcombobox(availModels,selected=1)
    targetfdr <- gslider(from=0.01,to=0.2,by=.01, value=0.05)
    
    # GUI layout MetSizeR without pilot data
    if(length(pilot)==0)
    {
      if(Targeted){ss.est<<- gwindow("MetSizeR without pilot targeted data")}
      else{ss.est<<- gwindow("MetSizeR without pilot NMR data")}
      pg <- gpanedgroup(container = ss.est, horizontal=TRUE)
      nb <- ggroup(horizontal = FALSE, container = pg)
      if(Targeted)
      {spec<-gframe("Metabolites", container=nb )
       add(spec, nMet )
      }else{spec<-gframe("Spectral bins", container=nb )
            add(spec, nspectralBins)}
      if(Targeted){propsig <- gframe("Proportion of significant metabolites", container=nb )}
      else{ propsig <- gframe("Proportion of significant bins", container=nb )}
      add(propsig, sigMet, expand=TRUE)
      Mod <- gframe("Model", container=nb)
      Mod.pos<- glayout(container=Mod)
      Mod.pos[1,1]<-models<-gcombobox(availModels,handler=function(h,...)
      {
        if (svalue(models)=="PPCA"){enabled(ncov)<- FALSE}
        if (svalue(models)=="DPPCA"){enabled(ncov)<- FALSE}
        if (svalue(models)=="PPCCA"){enabled(ncov)<- TRUE}
      },selected=1) 
      Mod.pos[1,2]<- "ncovars="
      Mod.pos[1,3]<-ncov<-gspinbutton(from=1, to = 10, by = 1, value=1)
      enabled(ncov)<- FALSE
      tfdr <- gframe("Target FDR", container=nb )
      add(tfdr,targetfdr , expand=TRUE)
      Sz<- gframe("Sample size per group", container=nb )
      Sz.pos<- glayout(container=Sz)
      Sz.pos[1,1]<- glabel("n1=")
      Sz.pos[1,2]<-ng1<-gedit("4")
      Sz.pos[2,1]<- glabel("n2=")
      Sz.pos[2,2]<-ng2<-gedit("4")
      cal<-gbutton("calculate", container=nb, handler = up.size)
      if(Targeted){gch<- gcheckbox("Plot proportion of significant metabolites", container= nb, handler = up.size)}
      else{gch<- gcheckbox("Plot proportion of significant bins", container= nb, handler = up.size)}
      dis<- gcheckbox("Save results in R directory", container= nb, handler=up.size)
      gseparator()
      fin<- gframe("", container=nb )
      fin.stat<- glayout(container=fin)
      fin.stat[10,1]<-glabel("MetSizeR status:", handler=function(h,...){enabled(sampleSize)<-FALSE})
      fin.stat[11,1,anchor=c(-1,-1)]<-sampleSize<- gedit("idle...")  
      enabled(fin)<-FALSE
      add(pg, ggraphics(dpi = 100, ps = 9))
    }else{
      
      # GUI layout MetSizeR with pilot data
      if(Targeted){ss.est<<- gwindow("MetSizeR with pilot targetd data")}
      else{ss.est<<- gwindow("MetSizeR with pilot NMR data")}
      pg <- gpanedgroup(container = ss.est, horizontal=TRUE)
      nb <- ggroup(horizontal = FALSE, container = pg)
      if(Targeted)
      {propsig <- gframe("Proportion of significant metabolites", container=nb )}
      else{propsig <- gframe("Proportion of significant bins", container=nb )}
      add(propsig,sigMet, expand=TRUE)
      Mod <- gframe("Model", container=nb)
      Mod.pos<- glayout(container=Mod)
      Mod.pos[1,1]<-models<-gcombobox(availModels,selected=1) 
      tfdr <- gframe("Target FDR", container=nb )
      add(tfdr,targetfdr, expand=TRUE)
      Sz<- gframe("Sample size per group", container=nb )
      Sz.pos<- glayout(container=Sz)
      Sz.pos[1,1]<- glabel("n1=")
      Sz.pos[1,2]<-ng1<-gedit("4")
      Sz.pos[2,1]<- glabel("n2=")
      Sz.pos[2,2]<-ng2<-gedit("4")
      cal<-gbutton("Calculate", container=nb, handler = up.size)
      if(Targeted){gch<- gcheckbox("Plot proportion of significant metabolites", container= nb, handler = up.size)}
      else{gch<- gcheckbox("Plot proportion of significant bins", container= nb, handler = up.size)}
      dis<- gcheckbox("Save results in R directory", container= nb, handler=up.size)
      gseparator()
      fin<- gframe("", container=nb )
      fin.stat<- glayout(container=fin)
      fin.stat[24,1]<-glabel("MetSizeR status:", handler=function(h,...){enabled(sampleSize)<-FALSE})
      fin.stat[25,1,anchor=c(-1,-1)]<-sampleSize<- gedit("is idle...")  
      enabled(fin)<-FALSE
      add(pg, ggraphics(dpi = 100, ps = 9))
    }
  }## end ssize.est function
  
  ####################################################################################
  # MetSizeR GUI component                                                              # 
  ####################################################################################
  cat("Welcome to MetSizeR!! \n")
  
  MetSizeR.menu<- gwindow("MetSizeR",container=TRUE)
  gimage(paste(file.path(path.package(package="MetSizeR")[1]),"/extdata/cover.jpg",sep=""),container= MetSizeR.menu)
  #addHandlerDestroy(MetSizeR.menu, handler=function(h,...){dispose(ssize.est)})
  ss.est<- gwindow(visible=FALSE)
  list_all = list(
    File=list(
      open = list(handler=pilot.open,icon="open"), 
      covariates = list(handler=cov.open,icon="plot"),
      demo_nmr_pilot_data = list(handler=demo.pilot.open, icon="newplot"),
      quit = list(handler=function(h,..){
        cat("you are quitting MetSizeR... \n")
        cat("byebye \n")
        dispose(MetSizeR.menu)
      }, icon="close")), 
    Sample_size = list(
      pilot_data = list(
        NMR_data = list(handler=function(h,..)
        {
          if(length(pilot)==0)
          {
            gmessage("Please upload the data first",icon="info")
          }else{
            ssize.est(inputData=pilot, inputCov=covars)
          }
        }),         
        Targeted_data = list(handler=function(h,..)
        {
          if(length(pilot)==0)
          {
            gmessage("Please upload the data first",icon="info")
          }else{
            #pilot<-log(pilot)
            ssize.est(inputData=pilot, inputCov=covars, Targeted=TRUE)
          }
        })),
      no_pilot_data  = list(
        NMR_analysis = list(handler=function(h,..)
        {
          pilot=NULL; covars=NULL
          ssize.est(inputData=pilot, inputCov=covars)
        }),
        Targeted_analysis = list(handler=function(h,..)
        {
          pilot=NULL; covars=NULL
          ssize.est(inputData=pilot, inputCov=covars, Targeted=TRUE)
        }))),
    Help = list(
      manual = list(handler=function(h,..)
      { print(help(MetSizeR))},icon="info")
      ))
  menu<- gmenu(list_all, container = MetSizeR.menu)
}
