
pred.y  <- function (y, x=NULL, B, sigma, lambda, w, model = "LN", t.outl=0.5) 
{
  # model - valori ammissibili: "N", "LN" 

   vars <- check.vars(y,x,model,parent="pred.y")

  if (vars$ret == -9) {
     stop(vars$msg.err) 
  }
  if (vars$ret != 0) {
     warning(vars$msg.err) 
  } 
  y <- as.matrix(vars$y)
  x <- as.matrix(vars$x) 
  nomicol.y <- colnames(y)
  id=1:nrow(y)
   ncoly <- ncol(y)
   n <- nrow(y)

   ncolx <- ncol(x)
   nomiy <- paste("y",1:ncoly,sep="")
   nomix <- paste("x",1:ncolx,sep="")
   mat.pat <- (is.na (as.matrix(y))) | ( is.infinite(y) )   
   mat.pat <- 1-mat.pat
   
   colnames(x) <- nomix
   colnames(y) <- nomiy
   colnames(mat.pat) <- paste(nomiy, ".pat", sep="")

   pattern <- apply(mat.pat,1,function(k) { appo<- ""; for (i in 1:length(k) ) { appo<-paste(appo,k[i],sep="") };appo  })
#   ll <- levels(as.factor(pattern))
   appo.df <- data.frame(x, y, 1:n, mat.pat, w,  pattern)

   colnames(appo.df) <- c(nomix, nomiy, "id", colnames(mat.pat), "w","pattern" )
   appo <- split (appo.df,pattern)
   stime<- NULL
   tau2 <- NULL
   for ( i in 1:length(appo))  {  
      pattern.corr <- appo[[i]]    
      y.corr <- as.matrix(pattern.corr[,nomiy])
      x.corr <- as.matrix(pattern.corr[,nomix])
      ind.obs <- as.logical(pattern.corr[1,colnames(mat.pat)])
      ind.mis <- !ind.obs


      mu <- as.matrix(x.corr %*% B )
      mu.o <- mu[,ind.obs,drop=FALSE]
      mu.m <- mu[,ind.mis,drop=FALSE]
      sigma.oo <- sigma [ind.obs, ind.obs, drop=FALSE]
      sigma.mm <- sigma [ind.mis, ind.mis, drop=FALSE]
      sigma.om <- sigma [ind.obs, ind.mis, drop=FALSE]    
        
      sigma.bar.oo <- (lambda * sigma.oo) / (lambda + 1)    
      sigma.bar.om <- (lambda * sigma.om) / (lambda + 1)      # Non si usa
      
          
    ############   attenzione: modifiche fatte da ugo: !!!!!  
      if (sum(ind.obs) == 0)  {
         
         mu.bar1 <-  mu.bar2 <- mu.m
         diag.sigma.bar1 <- diag.sigma.bar2 <- diag(sigma.mm)
       
         tau1 <- 1-pattern.corr[,"w"]   ### tau1 può essere qualsiasi numero perchè non ha più senso !!!
            
         
      } else {
  
        sigma.res.mo <-  sigma.mm - (t(sigma.om) %*% solve(sigma.oo) %*% sigma.om)
  
        sigma.bar.mm <- sigma.mm - ((t(sigma.om) %*% solve(sigma.oo) %*% sigma.om ) / (lambda + 1))   
        
        tau1 <-  post.prob(y.corr[,ind.obs,drop=FALSE], x.corr, B[,ind.obs,drop=FALSE],
        sigma.oo, 1-w, lambda)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        
        mu.tilde.obs <- (y.corr[,ind.obs,drop=FALSE] + (lambda * mu.o)) / (lambda + 1) 
        
        mu.bar1 <- cbind(y.corr[,ind.obs,drop=FALSE], mu.m + (y.corr[,ind.obs,drop=FALSE] - mu.o) %*% solve(sigma.oo) %*% sigma.om)
        mu.bar2 <- cbind(mu.tilde.obs, mu.m + (mu.tilde.obs - mu.o) %*% solve(sigma.oo) %*% sigma.om)   

        diag.sigma.bar1 <- c(rep(0,ncol(mu.o)), as.numeric(diag(sigma.res.mo)))
        diag.sigma.bar2 <- c(diag(sigma.bar.oo), as.numeric(diag(sigma.bar.mm)))  

           }
           
############   attenzione: modifiche fatte da ugo: !!!!!        
        if (model == "LN") {
        
        # diag.sigma.bar1 <- c(rep(0,ncol(mu.o)),  diag(sigma.res.mo))
        # diag.sigma.bar2 <- c(diag(sigma.bar.oo), diag(sigma.bar.mm))     
                                                                                                   
         mu.bar <- as.vector(tau1) * (t(exp (t(mu.bar1) + (0.5 * diag.sigma.bar1)))) + 
          as.vector(1 - tau1) * (t(exp (t(mu.bar2) + (0.5 * diag.sigma.bar2))))
        }
        else
         mu.bar <- as.vector(tau1) * mu.bar1 +  as.vector(1 - tau1) * mu.bar2    # previsione modello Normale
        
   
      colnames(mu.bar) <- c(nomiy[ind.obs],nomiy[ind.mis])     
      stima.corr <- matrix (NA, nrow(pattern.corr), ncoly+1)

    
      stima.corr <- cbind(mu.bar [,nomiy,drop=FALSE], 1-tau1)

      stime <- rbind(stime,cbind(stima.corr,pattern.corr[,c( "pattern", "id")]))

   }
   stime <- stime[order(stime[,"id"]),]
   ypred <- as.matrix(stime[,1:ncoly, drop=FALSE])
   if (is.null(nomicol.y))
       colnames(ypred)<-paste("ypred",1:ncoly,sep="") 
   else       
       colnames(ypred) <- paste(nomicol.y,"p",sep=".")
   tau<-stime[,ncoly+1]
   outlier<-(stime[,ncoly+1] > t.outl)+0
   pattern<- stime[,"pattern"]
   
#   list(ypred = as.matrix(stime[,1:ncoly, drop=FALSE]), 
#        tau=stime[,ncoly+1], 
#        outlier=(stime[,ncoly+1] > t.outl)+0, 
#        pattern = stime[,"pattern"]
#        )    

    data.frame(cbind(ypred,tau,outlier),pattern)             
} 

  
