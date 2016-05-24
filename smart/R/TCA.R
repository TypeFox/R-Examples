TCA <-
function(x,K,para,method="kendall", algorithm="tp", max.iter=200,verbose=TRUE,eps.conv=1e-3)
{
    if(length(para)!=K){ 
        cat("length of para is not equal to K \n")
        return(NULL)
    }
    if(!method%in%c("pearson","kendall","spearman","npn","ns")){ 
        cat("method not provided correctly, should be one of (pearson,kendall,spearman,npn,ns) \n")
        return(NULL)
    }

  	gcinfo(FALSE)
	  n = nrow(x);
	  d = ncol(x);
	  namecol=colnames(x)
	  fit = list()
	  fit$cov.input = isSymmetric(x);
	  if(fit$cov.input){
		   if(verbose) cat("The input is identified as the covriance matrix.\n")
		   S = cov2cor(x);
	  }
	  if(!fit$cov.input)
	  {
	     tmp = smart.npn(x,npn.func=method,verbose=verbose)
		   S = tmp$cov
		   y = tmp$scaled
		   rm(tmp)
	  }
	
	  xdat=x
	  rm(x)
	  gc()
     
     
    ### SPCA Algorithm 
    if(algorithm == "spca"){
   		 if(verbose){
   		   	cat("Conducting SPCA....")
        	flush.console()
   		 }        
       output = spca(S,K,para,type="Gram",sparse="varnum",max.iter=max.iter,trace=verbose,eps.conv=eps.conv)
       fit$loadings = M = output$loadings; pev = NULL 
       for(j in 1:K){
           Vk=M[,c(1:j)]; tmp=Vk%*%solve(t(Vk)%*%Vk)%*%t(Vk); 
           tmp=tmp%*%S%*%tmp;tmp=sum(diag(tmp))
           pev=c(pev,tmp/d)           
       }
       fit$pev=pev
     	 if(verbose){
   		    cat("done\n")
          flush.console()
       }
       rm(output,M,Vk,tmp,pev)
    }
    
    ### PMD Algorithm
    if(algorithm == "pmd"){
   		 if(verbose){
   		   	cat("Conducting PMD....")
        	flush.console()
   		 }
       output = SPC(S,sumabsv=para,K=K,niter=max.iter,trace=verbose)
       M = output$v; rownames(M)=namecol; colnames(M)=paste("PC",1:K,sep="")
       fit$loadings = M
       fit$pev = output$prop.var.explained
     	 if(verbose){
   		    cat("done\n")
          flush.console()
       }
       rm(output, M)
    }
    
    ### Truncated Power Algorithm
    if(algorithm == "tp"){
   		 if(verbose){
   		   	cat("Conducting TP....")
        	flush.console()
   		 }        
       
       M=NULL;sim.num=1;pev=NULL;A=S 
       while(sim.num <= K){ 
           x0 = SPC(A,sumabsv=sqrt(d)/2,K=1,trace=F)$v
           x = x0 
           tmp = A%*%x0; tmp = tmp/sqrt(sum(tmp^2)); x = x0
           trh = sort(abs(tmp),decreasing=T)[para[sim.num]]
           xt = tmp; xt[abs(xt)<trh] = 0
           xt = xt/sqrt(sum(xt^2))
           sim = 0
           while(sqrt(sum((xt-x)^2))>eps.conv & sim<max.iter){
               tmp = A%*%xt; tmp = tmp/sqrt(sum(tmp^2)); x = xt
               trh = sort(abs(tmp),decreasing=T)[para[sim.num]]
               xt = tmp; xt[abs(xt)<trh] = 0
               xt = xt/sqrt(sum(xt^2))
               sim = sim+1
               if(verbose) cat(sim)
           }
           if(verbose) cat("\n")
           M = cbind(M,xt) 
           A=(diag(dim(A)[2])-xt%*%t(xt))%*%A%*%(diag(dim(A)[2])-xt%*%t(xt))     
           sim.num = sim.num+1              
       }         
       for(j in 1:K){
           Vk=M[,c(1:j)]; tmp=Vk%*%solve(t(Vk)%*%Vk)%*%t(Vk); 
           tmp=tmp%*%S%*%tmp;tmp=sum(diag(tmp))
           pev=c(pev,tmp/d)           
       }  
       rownames(M)=namecol;colnames(M)=paste("PC",1:K,sep="")
       fit$loadings=M
       fit$pev=pev      
    	 if(verbose){
   		    cat("done\n")
          flush.console()
       }  
       rm(M,Vk,tmp,pev,trh,xt,x,sim,x0) 
    }    
    
    if(!fit$cov.input){
       if(method!= "pearson"){
           fit$PC = y%*%fit$loadings
           rm(y)
       }else{
           fit$PC = xdat%*%fit$loadings
           rm(xdat)
       }
    }   
    
    fit$method = method
    fit$algorithm = algorithm
    fit$K = K
    
    class(fit) = "TCA"
    return(fit)
}
