# grDM alternative

grDMb <-
  function(aDD,gr,design,TYPE=TYPE)
{
    #ctrl
    
    if(any(is.na(gr))){stop("Missing values inside the grouping vector are not allowed")}
    
    if(TYPE=="NRM")
    {
      
      if(nlevels(gr) == 1) # if no groups
      {
        cats <- sapply(aDD,function(x)x$anz_cat)
        nro  <- sum(cats*2)
        cateach <- rep(cats,each=2)
        partlist <- lapply(cateach, function(cei)
            {
              Qelement <- matrix(1/cei, cei-1, cei)
              diag(Qelement[,-1]) <- 1/cei - 1
              t(Qelement)
            })
        
       Q  <- diagblock(partlist)
        
        #naming
        prae <- rep(paste("I",1:length(cats),sep=""),cats*2)
        app1 <- unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(1:AA,2),sep="")))
        rwn  <- paste(prae,app1,sep="")
        rownames(Q) <- rwn
        colnames(Q) <- paste("eta",1:ncol(Q),sep="")
        # DONE
        
      } else if(nlevels(gr) > 1) # if there are any groups (> than 1 group)!
      {
        numbcat <- length(levels(gr))

        cats <- sapply(aDD,function(x)x$anz_cat)
        nro  <- sum(cats*2)
        cateach <- rep(cats,each=2)
        partlist <- lapply(cateach, function(cei)
        {
          Qelement <- matrix(1/cei, cei-1, cei)
          diag(Qelement[,-1]) <- 1/cei - 1
          t(Qelement)
        })
        
        Qbb     <- diagblock(partlist)
        mult1   <- matrix(0,numbcat,numbcat)
        
        # names with multiple groups
        prae1 <- rep(rep(paste("I",1:length(cats),sep=""),cats*2),numbcat)
        prae <- paste(rep(paste("G",1:numbcat,sep=""),each=nro),prae1,sep="")
        app1 <- rep(unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(1:AA,2),sep=""))),numbcat)
        rwn  <- paste(prae,app1,sep="")
        
        ##################
        #  DIF ? #########
        ##################
        
        if(all(design == "nodif"))
        {
          mult1[,1] <- 1
          Q <- mult1 %x% Qbb
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
          
        } else if(all(design == "dif1"))
        {
          # design with big identity matrix
          diag(mult1) <- 1
          Q <- mult1 %x% Qbb
          
          
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
          
        } else if(all(design == "dif2"))
        {
          # design where zetas are estimated within each group - lambdas remaining the same for each group
          mult1[lower.tri(mult1,diag=TRUE)] <- 1
          Q <- mult1 %x% Qbb
          
          whzeta <- grep("zeta",rwn)
          whzeta1 <- whzeta[whzeta > nrow(Qbb)]
          
          Q[whzeta1,1:ncol(Qbb)] <- 0 
          
          whlam <- grep("lam",rwn)
          whlam1 <- whlam[whlam > nro]
          Q[whlam1,-(1:ncol(Qbb))] <- 0
          
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
          
        } else if(all(design == "dif3"))
        {
          # design where lambda is estimated for each group while the zetas remain the same for all different groups
          mult1[lower.tri(mult1,diag=TRUE)] <- 1
          Q <- mult1 %x% Qbb
          
          whzeta <- grep("lam",rwn)
          whzeta1 <- whzeta[whzeta > nrow(Qbb)]
          
          Q[whzeta1,1:ncol(Qbb)] <- 0 
          
          whlam <- grep("zeta",rwn)
          whlam1 <- whlam[whlam > nro]
          Q[whlam1,-(1:ncol(Qbb))] <- 0
          
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
        } else if(is.list(design))
        {

          
          # hier jetzt falls ein eigenes design angegeben wird.
          #mult1[,1] <- 1
          mult1[lower.tri(mult1,diag=TRUE)] <- 1
          
          Q <- mult1 %x% Qbb
          rownames(Q) <- rwn
          
          
          spaltBEG <- c(0,rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1 
          spaltEND <- c(rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))
          #aDD
          zeilBEG <- c(0,rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1
          zeilEND <- c(rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))
          
          prmA <- c("zeta","lam")
          

          zerg <- mapply(function(begz,endz,galaZ) # speichert alle zeilen
          {

            prot <- mapply(function(begsp,endsp,grunr) # speichert je eine vollständige zeile
            {
              
              Qtemp <- Q[begz:endz,begsp:endsp] #auswahl des quadranten
              
              
              for(EACH in 1:length(design))
              {
                gala  <- design[[EACH]]
                loe   <- which(gala[galaZ,] != grunr)
                prm   <- prmA[EACH]
                
                if(length(loe) == 0)
                {
                  #Qtemp
                  next
                  
                } else {
                  
                  iaus <- paste("^.+(",paste("I",loe,collapse="|",sep=""),")",prm,sep="")
                  wo <- grep(iaus,rownames(Qtemp),value=F,perl=TRUE)
                  Qtemp[wo,] <- 0
                  
                }
              }
              Qtemp
            },begsp=spaltBEG, endsp=spaltEND, grunr=1:nlevels(gr),SIMPLIFY=FALSE)
            
            
            
          },begz=zeilBEG, endz=zeilEND,galaZ=1:nlevels(gr),SIMPLIFY=FALSE)
          
          
          zwQ <- lapply(zerg,function(XX)do.call(cbind,XX))
          Q  <- do.call(rbind,zwQ)
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
        } else {stop("Check your input for argument: 'design' ")}
        
      } 
      
    } 
    
    # -----------------------------------------------#
    ############## ---->  NLM <---- ##################
    ##################################################
    
    if(TYPE=="NLM")
    {
      
      if(nlevels(gr) == 1) # if no groups
      {
        
        cats <- sapply(aDD,function(x)x$anz_cat)
        cateach <- rep(cats,each=2)
        partlist <- lapply(cateach, function(cei)
        {
          Qelement <- matrix(0,cei, cei-1)
          Qelement[1,1] <- 1
          ceiNLM <- cei -1 
          Qelement1 <- matrix(1/ceiNLM, ceiNLM-1, ceiNLM)
          diag(Qelement1[,-1]) <- 1/ceiNLM - 1
          Qelement[-1,-1] <- t(Qelement1)
          Qelement
        })
        
        Q  <- diagblock(partlist)
        
        #naming
        prae <- rep(paste("I",1:length(cats),sep=""),cats*2)
        app1 <- unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(0:(AA-1),2),sep=""),simplify=FALSE ))
        app1[grep("zeta0",app1)] <- "beta"
        app1[grep("lam0",app1)]  <- "alpha"
        
        rwn  <- paste(prae,app1,sep="")
        rownames(Q) <- rwn
        colnames(Q) <- paste("eta",1:ncol(Q),sep="")
        # DONE
        
      } else if(nlevels(gr) > 1)
      {
        
        numbcat <- nlevels(gr)

        cats <- sapply(aDD,function(x)x$anz_cat)
        nro  <- sum(cats*2)
        cateach <- rep(cats,each=2)
        partlist <- lapply(cateach, function(cei)
        {
          Qelement <- matrix(0,cei, cei-1)
          Qelement[1,1] <- 1
          ceiNLM <- cei -1 
          Qelement1 <- matrix(1/ceiNLM, ceiNLM-1, ceiNLM)
          diag(Qelement1[,-1]) <- 1/ceiNLM - 1
          Qelement[-1,-1] <- t(Qelement1)
          Qelement
        })
        
        Qbb  <- diagblock(partlist)
        
        
        mult1   <- matrix(0,numbcat,numbcat)
        
        # names with multiple groups
        prae1 <- rep(rep(paste("I",1:length(cats),sep=""),cats*2),numbcat)
        prae <- paste(rep(paste("G",1:numbcat,sep=""),each=nro),prae1,sep="")
        
        app1 <- rep(unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(0:(AA-1),2),sep=""))),numbcat)
        app1[grep("zeta0",app1)] <- "beta"
        app1[grep("lam0",app1)]  <- "alpha"
        
        rwn  <- paste(prae,app1,sep="")
        
        
        ##################
        #  DIF ? #########
        ##################   
        
        if(all(design == "nodif"))
        {
          mult1[,1] <- 1
          Q <- mult1 %x% Qbb
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
          
        } else if(is.list(design))
        {
          
          mult1[lower.tri(mult1,diag=TRUE)] <- 1
          
          Q <- mult1 %x% Qbb
          rownames(Q) <- rwn
          
          
          spaltBEG <- c(0,rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1 
          spaltEND <- c(rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))
          #aDD
          zeilBEG <- c(0,rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1
          zeilEND <- c(rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))
          
          
          prmA <- c("alpha","beta","zeta","lam")

          zerg <- mapply(function(begz,endz,galaZ) # speichert alle zeilen
          {

            prot <- mapply(function(begsp,endsp,grunr) # speichert je eine vollständige zeile
            {
              
              Qtemp <- Q[begz:endz,begsp:endsp] #auswahl des quadranten
              
              
              for(EACH in 1:length(design))
              {
                gala  <- design[[EACH]]
                loe   <- which(gala[galaZ,] != grunr)
                prm   <- prmA[EACH]
                
                if(length(loe) == 0)
                {
                  #Qtemp
                  next
                  
                } else {
                  
                  iaus <- paste("^.+(",paste("I",loe,collapse="|",sep=""),")",prm,sep="")
                  wo <- grep(iaus,rownames(Qtemp),value=F,perl=TRUE)
                  Qtemp[wo,] <- 0
                  
                }
              }
              Qtemp
            },begsp=spaltBEG, endsp=spaltEND, grunr=1:nlevels(gr),SIMPLIFY=FALSE)
            
            
            
          },begz=zeilBEG, endz=zeilEND,galaZ=1:nlevels(gr),SIMPLIFY=FALSE)
          
          
          zwQ <- lapply(zerg,function(XX)do.call(cbind,XX))
          Q  <- do.call(rbind,zwQ)
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")

        } else {stop("Check your input for argument: 'design' ")}
        
      }
      
      
      
    }
    
    return(Q)
    
  }

















