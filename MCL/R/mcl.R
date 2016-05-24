mcl <-
function(x, addLoops = NULL, expansion = 2, inflation = 2, allow1 = FALSE, max.iter = 100, ESM = FALSE ){
    
    
    if(is.null(addLoops)){
      stop("addLoops has to be TRUE or FALSE")
    }
    
    
    if (addLoops) diag(x) <- 1    
    
    
    adj.norm <- apply(x[,], MARGIN=2, FUN=function(Var){
      Var/sum(Var)
    })
    
    
    a <- 1
    
    
    repeat{
      
      expans <- adj.norm %^% expansion
      
      
      infl <- expans ^ inflation
            
      
      infl.norm <- apply(infl[,], MARGIN=2, FUN=function(Var){
        Var/sum(Var)
        
      })
      
      
      if(identical(infl.norm,adj.norm)) {
        ident <- TRUE
        break
      }
      

      if(a==max.iter) {
        ident <- FALSE
        a <- a+1
        break
      }
      
      
      adj.norm <- infl.norm
      a<- a+1
    }
    
    
    
    if(!is.na(infl.norm[1,1]) & ident){
      
      count <- 0 
      for(i in 1:ncol(infl.norm)){
        if(sum(abs(infl.norm[i,])) != 0) {
          count <- count+1
        }
      }
      
      neu <- matrix(nrow=count, ncol=ncol(infl.norm)) 
      
      zeile <- 1
      for(i in 1:nrow(infl.norm)){
        if(sum(infl.norm[i,]) != 0) {
          for(j in 1:ncol(infl.norm)) {
            
            neu[zeile,j]<-infl.norm[i,j]
          }
          zeile <- zeile+1
        }
      }
      
         
      for(i in 1:nrow(neu)){
        for(j in 1:ncol(neu)) {
          if((neu[i,j] < 1) & (neu[i,j] > 0)){
            neu[,j] <- 0
            neu[i,j] <- 1
          }
        }
      }
      
      
      for(i in 1:nrow(neu)){
        for (j in 1:ncol(neu)){
          if(neu[i,j] != 0){
            neu[i,j] <- i
          }
        }
      }
      
      ClusterNummern <- sum(neu[,1])
      for(j in 2:ncol(neu)){
        ClusterNummern <- c(ClusterNummern,sum(neu[,j]))
      }
      
      
    } 
    
    ifelse(!(!is.na(infl.norm[1,1]) & ident), output <- paste("An Error occurred at iteration", a-1),
      {
      if(!allow1){
        dub <- duplicated(ClusterNummern) + duplicated(ClusterNummern,fromLast = T)
        for(i in 1:length(dub)){
          if(dub[[i]]==0) ClusterNummern[[i]]<-0
        }
      }
      
      #### dimnames for infl.norm
      dimnames(infl.norm) <- list(1:nrow(infl.norm), 1:ncol(infl.norm))
      
      output <- list()
      output[[1]] <- length(table(ClusterNummern))
      output[[2]] <- a-1 
      output[[3]] <- ClusterNummern
      output[[4]] <- infl.norm
      
      names(output) <-c("K", "n.iterations","Cluster",
                        "Equilibrium.state.matrix")
    }
    )
  ifelse(ESM==TRUE,return(output),return(output[-4]))
}
