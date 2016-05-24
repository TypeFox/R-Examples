  .packageName<-'JOP'                                                                                   
    ########################################################
    ####### Function that counts zeros after comma   #######
    ########################################################
    countdig<-function(x)
    {
      count<-0
      if(abs(x)<1e-8){
      return(0)}
      while(abs(x)<1)
      {
        x<-x*10
        count<-count+1
      }
      return(count)
    }
    ########################################################
    ####### Transformation to polar coordinates ############
    ########################################################
    trafopar<-function(x,optreg)
    {
        n   <- length(x)
        if (n==1)
        {
        return(cat("Error: Function needs more dimensions!\n"))
        }
        if(optreg=="box")
        {
          return(x)
        }
        y   <- numeric(n)
        r   <- x[1]
        phi <- numeric(length(n-1))
        
        for (i in 1:(n-1))
        {
            phi[i] <- x[i+1]
        }
        
        y[n] <- r * sin(phi[n-1])
        y[1] <- r * cos(phi[1])
        y[2] <- r * sin(phi[1])
        
          if (n==2)
        {
        return(y)
        }
        
        for (i in 3:n)
        {
            for (j in 1:(i-1))
            {
            y[j] <- y[j]*cos(phi[i-1])      
            }
        
        y[i] <- r * sin(phi[i-1])  
        }
        
        if(optreg=="sphere")
        {
          return(y)
        }
    }
 
      is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
     
         ###############################################
         ###############################################
         #########  Main Function 'JOP'       ##########
         ###############################################
         ###############################################
      
    JOP<-function(datax,datay,tau="min",Wstart=-5,Wend=5,numbW=10,d=NULL,optreg="sphere",Domain=NULL,form.mean=NULL,form.disp=NULL,family.mean=gaussian(),dlink="log",mean.model=NULL,var.model=NULL,joplot=FALSE,solver="solnp")
    {
    ######### Checks if the input is correctly passed
    ######### START CHECK
      if(is.null(datax) && is.null(datay))
      {
        return(cat("\n Error: Experimental plan is needed!\n"))
      }
      
      if(dlink=="log")
      {
        dlinkvec<-vector("list",dim(datay)[2])
        for(i in 1:dim(datay)[2])
        {
          dlinkvec[[i]]<-dlink
        }
      }
      if(length(family.mean)!=dim(datay)[2])
      {
         family.mean<-vector("list",dim(datay)[2])
         for(i in 1:dim(datay)[2])
         {
            family.mean[[i]]<-gaussian()
         }
      }
      
      if(is.data.frame(datax)==FALSE||is.data.frame(datay)==FALSE)
      {
        return(cat("\n Error: datax and datay have to be data frames\n"))
      }
      
      nx<-dim(datax)[2]
      ny<-dim(datay)[2]
      

      if(is.wholenumber(numbW)==FALSE)
      {
        return(cat("\n Error: numbW has to be a whole number!\n"))
      }
      
      if(!is.null(d))
      {
        if(!is.double(d))
        {
          return(cat("\n Error: d has to be numeric!\n"))
        }
      }
      if(is.null(d))
      {
        d<-c(1,rep(0,ny-1))
      }    
  
      if(numbW<2)
      {
        return(cat("\n Error: numbW has to be a positive whole number greater than 2!\n"))
      }
      
      if(Wstart>=Wend)
      {
        return(cat("\n Error: Wend have to be greater than Wstart\n"))
      }
      
      if(nx<2)
      {
        return(cat("\n Error: More then 1 variable required!!\n"))
      }
      
      if(ny<2)
      {
        return(cat("\n Error: More then 1 response required!!\n"))
      }
      if(length(d)!=ny)
      {
        return(cat("\n Error: Dimension of d do not match to number of responses!\n"))
      }
      if(!is.list(tau))
      {
        taucor<-0
        if(tau=="min" || tau=="max")
        {
          tau1<-vector("list",ny)
          for(i in 1:ny)
          {
            tau1[[i]]<-tau
          }
          tau<-tau1  
          taucor<-1
        }
        if(taucor==0)
        {
          return(cat("\n Error: tau is not correctly specified!\n")) 
        }
      }
      if(length(tau)!=ny)
      {
        return(cat("\n Error: Dimension of tau do not match to number of responses\n"))
      }
      for(i in 1:length(tau))
      {
        if(!is.double(tau[[i]])&&tau[[i]]!="max"&&tau[[i]]!="min")
        {
          return(cat("\n Error: tau has to be either numeric or max or min!\n"))
        }
      }
      tau1<-tau
      if(!is.null(mean.model) && !is.null(mean.model))
      { 
        if(is.list(mean.model)==FALSE || is.list(var.model)==FALSE)
        {
          return(cat("\n Error: mean.model and var.model have to be of type 'list'\n"))
        }
        for(o in 1:length(mean.model))
        {
          if(is.function(mean.model[[o]])==FALSE || is.function(var.model[[o]])==FALSE)
          {
            return(cat("\n Error: mean.model and var.model have to contain functions!\n"))  
          }
          else
          {
            testmean<-try(mean.model[[o]](as.numeric(datax[1,])),silent=TRUE)
            testvar<-try(var.model[[o]](as.numeric(datax[1,])),silent=TRUE)
            if(is.character(testmean) || is.character(testvar))
            {
              return(cat(paste("\n Error: Check your functions for mean and dispersion for response",names(datay)[o],"!\n",sep="")))
            }
          }
        }
      }
      if(!is.null(Domain))
      {
        if(dim(Domain)[1]!=nx && dim(Domain)[2]!=2)
        {
          return(cat("Check Dimensions of Domain!\n"))
        }
        if(!is.numeric(Domain))
        {
          return(cat("Domain entries have to be numeric!\n"))
        }
        optreg<-"box"
      }
    ######### END CHECK
    
  ##################################################
  ##################################################
  ##                   
  ## Automatic model building
  ##
  ##### START MODEL BUILDING
    shifter<-0
    if(is.null(mean.model) || is.null(var.model))
    {
      shifter<-1
      cat("Automatic Modeling starts...\n")
  
      xnames<-names(datax)
      resp<-names(datay)
      quadeff<-NULL
      for(i in 1:length(xnames))
      {
        quadeff[i]<-paste("I(",xnames[i],"^2)",sep="")
      }
      rights<-paste(paste("(",paste(xnames,collapse="+"),")^2+",sep=""),paste(quadeff,collapse="+"),sep="")
          
      #################
      ## Mean-Modell ##
      #################
      flist<-vector("list",ny)
      names(flist)<-names(datay)  
      for(i in 1:ny)
      {
        flist[[i]]<-as.formula(paste(paste(resp[i],"~",sep=""),rights,sep=""))
      }
  
      #######################
      ## Dispersion-Modell ##
      #######################
      dispf<-vector("list",ny)
      names(dispf)<-names(datay) 
       
      for(i in 1:ny)
      {
        dispf[[i]]<-as.formula(paste("~",paste(xnames,collapse="+"),sep=""))#as.formula(paste("~",rights,sep=""))#
      }
      
      outmod<-vector("list",ny)
      names(outmod)<-resp
      qval<-0.1
      xnames<-names(datax)
      
      for(i in 1:ny)
      {
          k<-1
          if(!is.null(form.mean[[i]]) && !is.null(form.disp[[i]]))
          {
            dataset<-data.frame(datay[i],datax)
            dlink<-dlinkvec[[i]]
            invisible(capture.output(outmod[[i]]<-try(dglm(form.mean[[i]],form.disp[[i]],data=dataset,family=family.mean[[i]],dlink=dlink,method="reml"),silent=TRUE)))
            if(is.character(outmod[[i]]))
            {
              return(cat(paste("\n Error: Model building failed for ",names(datay)[i],"\n\n","Check the distribution assumption or the link function\n",sep="")))
            }  
            k<-0
          }
          formm<-0
          if(!is.null(form.mean[[i]]) && is.null(form.disp[[i]]))
          {
            flist[[i]]<-form.mean[[i]]
            k<-1
            formm<-1
          }
          formd<-0
          if(is.null(form.mean[[i]]) && !is.null(form.disp[[i]]))
          {
            dispf[[i]]<-form.disp[[i]]
            k<-1
            formd<-1
          }
          if((is.null(form.mean[[i]])||formm==1) && (is.null(form.disp[[i]])||formd==1))
          {
            while(k==1)
            {
              dataset<-data.frame(datay[i],datax)
              dlink<-dlinkvec[[i]]
              invisible(capture.output(outmod[[i]]<-try(dglm(flist[[i]],dispf[[i]],data=dataset,family=family.mean[[i]],dlink=dlink,method="reml"),silent=TRUE)))
              if(is.character(outmod[[i]]))
              {
                return(cat(paste("\n Error: Model building failed for ",names(datay)[i],"\n\n","Check the distribution assumption or the link function\n",sep="")))
              }
              outm<-summary(outmod[[i]])$coefficients
              outd<-summary(outmod[[i]])$dispersion.summary$coefficients
              helpa<-outm[,4]
              helpb<-outd[,4]
              mina<-ifelse(length(helpa)>1,max(helpa[-1]),0.5*qval)
              mina<-ifelse(formm==1,0.5*qval,mina)
              minb<-ifelse(length(helpb)>1,max(helpb[-1]),0.5*qval)
              minb<-ifelse(formd==1,0.5*qval,minb)
              if(minb<qval && mina<qval)
              {
                k<-0
                coeff<-outmod[[i]]$dispersion.fit$coefficients
                ### The following routine checks if all main effects of higher order terms are included
                ## START
                addcoeff<-NULL
                iter<-1
                if(length(coeff)>1)
                {
                  for(ii in 2:length(coeff))
                  {
                     for(jj in 1:length(xnames))
                     {
                       if(names(coeff)[ii]!=xnames[jj] && names(coeff)[ii]==paste("I(",xnames[jj],"^2)",sep=""))
                       {
                          addcoeff[iter]<-xnames[jj]
                          iter<-iter+1   
                       }
                       for(kk in jj:length(xnames))
                       {
                          if(names(coeff)[ii]!=xnames[kk] && kk!=jj && names(coeff)[ii]==paste(xnames[jj],":",xnames[kk],sep="") || coeff[ii]==paste(xnames[kk],":",xnames[jj],sep=""))
                          {
                            addcoeff[iter]<-xnames[kk]
                            iter<-iter+1
                          }
                       }
                     }
                  }
                }
                ## end
                updisp<-0
                if(iter==1)
                {           
                  dispf[[i]]<-formula(dispf[[i]])
                }
                if(iter>1)
                {
                  updisp<-1
                  dispf[[i]]<-update.formula(formula(dispf[[i]]),paste("~.",paste(addcoeff,collapse="+"),sep="+"))
                }
                
                coeff<-outmod[[i]]$coefficients
                ### The following routine checks if all main effects of higher order terms are included
                ## START
                addcoeff<-NULL
                iter<-1
                if(length(coeff)>1)
                {
                  for(ii in 2:length(coeff))
                  {
                     for(jj in 1:length(xnames))
                     {
                       if(names(coeff)[ii]!=xnames[jj] && names(coeff)[ii]==paste("I(",xnames[jj],"^2)",sep=""))
                       {
                          addcoeff[iter]<-xnames[jj]
                          iter<-iter+1   
                       }
                       for(kk in jj:length(xnames))
                       {
                          if(names(coeff)[ii]!=xnames[kk] && kk!=jj && names(coeff)[ii]==paste(xnames[jj],":",xnames[kk],sep="") || coeff[ii]==paste(xnames[kk],":",xnames[jj],sep=""))
                          {
                            addcoeff[iter]<-xnames[kk]
                            iter<-iter+1
                          }
                       }
                     }
                  }
                }
                ## end
                if(iter==1)
                {           
                  flist[[i]]<-formula(flist[[i]])
                }
                if(iter>1 || updisp==1)
                {
                  flist[[i]]<-update.formula(formula(flist[[i]]),paste(resp[i],"~.+",paste(addcoeff,collapse="+"),sep=""))
                  dlink<-dlinkvec[[i]]
                  invisible(capture.output(outmod[[i]]<-try(dglm(flist[[i]],dispf[[i]],data=dataset,family=family.mean[[i]],dlink=dlink,method="reml"),silent=TRUE)))
                  if(is.character(outmod[[i]]))
                  {
                    return(cat(paste("\n Error: Model building failed for ",names(datay)[i],"\n\n","Check the distribution assumption or the link function\n",sep="")))
                  }
                }

              }
              if(minb>qval && minb>=mina)
              { 
                  dispf[[i]]<-as.formula(paste("~",ifelse(length(helpb)==2,"1",paste(names(helpb[-1])[-which(helpb[-1]==minb)],collapse="+")),sep=""))
              }
              if(mina>qval && mina>minb)
              { 
                  flist[[i]]<-as.formula(paste(resp[i],"~",ifelse(length(helpa)==2,"1",paste(names(helpa[-1])[-which(helpa[-1]==mina)],collapse="+")),sep=""))
              }
            }
          }
      }
      mean.model<-vector("list",ny)
      var.model<-vector("list",ny)
      meanmd<-function(x,p)
      {
        allx<-NULL
        xnames<-names(datax)
        coeff<-outmod[[p]]$coefficients
        value<-coeff[1]
        if(length(coeff)>1)
        {
          for(i in 2:length(coeff))
          {
             for(j in 1:length(xnames))
             {
               if(names(coeff)[i]==xnames[j])
               {
                  value<-value+coeff[i]*x[j]
               }
               if(names(coeff)[i]==paste("I(",xnames[j],"^2)",sep=""))
               {
                  value<-value+coeff[i]*x[j]^2
               }
               for(k in j:length(xnames))
               {
                  if(k!=j && names(coeff)[i]==paste(xnames[j],":",xnames[k],sep="") || names(coeff)[i]==paste(xnames[k],":",xnames[j],sep=""))
                  {
                    value<-value+coeff[i]*x[j]*x[k]
                  }
               }
             }
          }
        }
        names(value)<-NULL
        return(outmod[[p]]$family$linkinv(value))
      }
      varmd<-function(x,p)
      {
        allx<-NULL
        xnames<-names(datax)
        coeff<-outmod[[p]]$dispersion.fit$coefficients
        value<-coeff[1]
        if(length(coeff)>1)
        {
          for(i in 2:length(coeff))
          {
           for(j in 1:length(xnames))
           {
             if(names(coeff)[i]==xnames[j])
             {
                value<-value+coeff[i]*x[j]
             }
             if(names(coeff)[i]==paste("I(",xnames[j],"^2)",sep=""))
             {
                value<-value+coeff[i]*x[j]^2
             }
             for(k in j:length(xnames))
             {
                if(k!=j && names(coeff)[i]==paste(xnames[j],":",xnames[k],sep="") || names(coeff)[i]==paste(xnames[k],":",xnames[j],sep=""))
                {
                  value<-value+coeff[i]*x[j]*x[k]
                }
             }
           }
          }
        }
        names(value)<-NULL
        return(outmod[[p]]$family$variance(meanmd(x,p))*outmod[[p]]$dispersion.fit$family$linkinv(value))
      }

#      meanmd<-function(x,p)
#      {
#        x<-as.data.frame(t(x))
#        names(x)<-names(datax)
#        value<-predict(outmod[[p]],x,type="response")
#        return(value)
#      }
#      varmd<-function(x,p)
#      {
#        x<-as.data.frame(t(x))
#        names(x)<-names(datax)
#        value<-predict(outmod[[p]]$dispersion.fit,x,type="response")
#        return(value)
#      }

      for(ii in 1:ny)
      { 
        ## Final models:
        mean.model[[ii]]<-meanmd
        var.model[[ii]]<-varmd
      }
      cat("\n")
      cat("Model building finished ....\n ")
    }
    ##### END MODEL BUILDING
    
    ####
    #  Model building finished
    #
    
    ## if some of the target values are min or max, then 
    ## JOP calculates the maximum or minimum and then takes
    ## the calculated values as targets.
    
    ### BEGIN TARGET VALUES
      if(optreg=="sphere")
      {
        lower<-c(0,rep(-pi,nx-1))
        normvec<-NULL
        for(i in 1:dim(datax)[1])
        {
          normvec[i]<-sqrt(t(as.numeric(datax[i,]))%*%as.numeric(datax[i,]))
        }
        upper<-c(max(normvec),rep(pi,nx-1))
      }
      if(optreg=="box")
      {                
        lower<-NULL
        upper<-NULL
        for(i in 1:nx)
        {
          lower[i]<-min(datax[,i])
          upper[i]<-max(datax[,i])
        }
      }
      if(!is.null(Domain))
      {
        lower<-Domain[,1]
        upper<-Domain[,2]
      }
      for(i in 1:ny)
      {  
        if(tau[[i]]=="max")
        {
          helpfun<-function(x)
          {
            return(-ifelse(shifter==1,mean.model[[i]](trafopar(x,optreg),i),mean.model[[i]](trafopar(x,optreg))))
          }
          tauval<-gosolnp(fun=helpfun,LB=lower,UB=upper,control=list(trace=0))
          tau[[i]]<--tauval$values[length(tauval$values)]
        }
        if(tau[[i]]=="min")
        {
          helpfun<-function(x)
          {
            return(ifelse(shifter==1,mean.model[[i]](trafopar(x,optreg),i),mean.model[[i]](trafopar(x,optreg))))
          }
          tauval<-gosolnp(fun=helpfun,LB=lower,UB=upper,control=list(trace=0))
          tau[[i]]<-tauval$values[length(tauval$values)]      
        }
      }
      tau<-as.numeric(as.vector(tau))
      for(i in 1:ny)
      {
        tau[i]<-round(tau[i],digits=2+countdig(tau[i]))
      }
    ###### END TARGET VALUES 
    
    
    ################################################## 
    ####       
    #  Weight matrices
    #
    W<-vector("list",numbW)
    logat<-seq(Wstart,Wend,(Wend-Wstart)/(numbW-1))
    for(i in 1:(numbW))
    {
      W[[i]]<-diag(exp(d*logat[i])) 
    }
    ####
    #  Standardization matrix
    #
    zw<-NULL
    diagA<-NULL
    for(i in 1:ny)
    {
      for(j in 1:dim(datax)[1])
      {
        if(shifter==1)
        {
          zw[j]<-var.model[[i]](as.numeric(datax[j,]),i) 
        }
        else
        {
          zw[j]<-var.model[[i]](as.numeric(datax[j,])) 
        }
      }
      diagA[i]<-1/sqrt(mean(zw))
      zw<-NULL 
    }  
    A<-diag(diagA)
    ####
    #   Cost Matrix C
    for(i in 1:length(W))
    {
      W[[i]]<-t(A)%*%W[[i]]%*%A
    }       
    cat("\n")
    cat("Cost matrices calculated ....\n ")       
           #################################
           #################################
           ######### Optimization ##########
           #################################
           #################################
    ####
    #  Matrix with optimal Parameter Settings
    #
    optmatrix<-matrix(NaN,ncol=nx,nrow=numbW)
      
    ####    
    #  Matrix with optimal Responses
    #
    reoptmatrix<-matrix(NaN,ncol=ny,nrow=numbW)
  
  
    ####
    #  specifying constraints constraints
    #
    if(optreg=="sphere")
    {
      start1<-c(max(datax)/2,rep(-pi/2,nx-1))
      start2<-c(max(datax)/2,rep(0,nx-1))
      start3<-c(max(datax)/2,rep(pi/2,nx-1))          
    }
    if(optreg=="box")
    {
      start1<-rep(max(datax)/2,nx)
      start2<-rep(-max(datax)/2,nx)
      start3<-rep(0,nx)                 
    }
    deviation<-matrix(NaN,ncol=ny,nrow=numbW)
    optval<-NULL
    cat("\n")
    cat("Optimization starts ....\n ")
    pb <- txtProgressBar(min = 0, max = numbW, style = 3)
    for(i in 1:numbW)
    {
    ####
    #  Building risc function for every cost matrix
    #
      riscfun<-function(x)
      {
        varval<-NULL
        meanval<-NULL
        for(j in 1:ny)
        {
          if(shifter==1)
          {
            varval[j]<-as.numeric(var.model[[j]](trafopar(x,optreg),j))
            meanval[j]<-as.numeric(mean.model[[j]](trafopar(x,optreg),j))    
          }
          else
          {
            varval[j]<-as.numeric(var.model[[j]](trafopar(x,optreg)))
            meanval[j]<-as.numeric(mean.model[[j]](trafopar(x,optreg)))              
          }
        }
        return(sum(diag(W[[i]]%*%diag(varval)))+sum(t(meanval-tau)%*%W[[i]]%*%(meanval-tau)))  
      }
      ####
      #  Solver starts with three initial vectors 
      #  Then take best solution
      #
    if(solver=="gosolnp")
    {
      opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,control=list(trace=0),rseed=100+i)
    }
    if(solver=="solnp")
    {
      opt1<-solnp(start1,fun=riscfun,LB=lower,UB=upper,control=list(trace=0))
      opt<-opt1
      opt2<-solnp(start2,fun=riscfun,LB=lower,UB=upper,control=list(trace=0))
      if(opt$values[length(opt$values)]>opt1$values[length(opt1$values)])
      {
        opt<-opt2
      }
      opt3<-solnp(start3,fun=riscfun,LB=lower,UB=upper,control=list(trace=0))    
      if(opt$values[length(opt$values)]>opt3$values[length(opt3$values)])
      {
        opt<-opt3
      }
    }
      optval[i]<-opt$values[length(opt$values)]
      optmatrix[i,]<-trafopar(opt$pars,optreg)
      for(k in 1:ny)
      {
        if(shifter==1)
        {
          reoptmatrix[i,k]<-as.numeric(mean.model[[k]](trafopar(opt$pars,optreg),k))
          deviation[i,k]<-as.numeric(var.model[[k]](trafopar(opt$pars,optreg),k))
        }
        else
        {
          reoptmatrix[i,k]<-as.numeric(mean.model[[k]](trafopar(opt$pars,optreg)))
          deviation[i,k]<-as.numeric(var.model[[k]](trafopar(opt$pars,optreg)))        
        }
      }
    
      ####
      #   Solver finished
      #
      setTxtProgressBar(pb, i)                                             
      }
    close(pb)
    cat("Optimization finished ....\n ")
    cat("\n")
    ####
    #  Optimization finished
    #  
    
    ####
    #  Create Output-List
    #
     
    variance<-deviation
    deviation<-sqrt(deviation)
    sirisk<-matrix(NA,nrow=numbW,ncol=ny)
    sirisknorm<-NULL
    for(i in 1:numbW)
    {
      for(j in 1:ny)
      {
        sirisk[i,j]<-(variance[i,j]+(reoptmatrix[i,j]-tau[j])^2)/diagA[j]^(-2)
      }
      sirisknorm[i]<-sqrt(sum(sirisk[i,]^2))
    }
    xlu<-which(sirisknorm==min(sirisknorm))[1]
    xaxis1<-1:numbW
    xaxis1names<-paste("W",xaxis1,sep="")
    
    dimnames(optmatrix)<-list(xaxis1names,names(datax))
    dimnames(reoptmatrix)<-list(xaxis1names,names(datay))
    dimnames(deviation)<-list(xaxis1names,names(datay))
    optval<-as.matrix(optval)
    dimnames(optval)<-list(xaxis1names,c(""))
    
    if(shifter==1)
    {
    ####
    #   Output DGLM if no models specified by user
    #
      optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau,tau1,outmod,optmatrix[xlu,],reoptmatrix[xlu,],c(Wstart,Wend),d,numbW)
      names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValueJOP","TargetValueUSER","DGLM","RiskminimalParameters", "RiskminimalResponses","ValW","d","numbW") 
    }
    else
    {
      optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau,tau1,optmatrix[xlu,],reoptmatrix[xlu,],c(Wstart,Wend),d,numbW)
      names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValueJOP","TargetValueUSER","RiskminimalParameters", "RiskminimalResponses","ValW","d","numbW") 
    }
    ### Output is of Classe "JOP"
    class(optimres)<-"JOP"  
    ####
    #  Joint optimization Plot 
    #  if joplot=T
    if(joplot==TRUE)plot(optimres)
    
    return(optimres)
  }     

