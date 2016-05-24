#rm(list=ls()); i <- 1
#library(lavaan)
#load("~/lavori/proveR/influence.SEM/data/jeff.rda")
#model <- model00
#data <- X
#fit0 <- sem(model, data=data)
#inspect(fit0,"rsquare")
#parm <- PAR

parinfluence <-
function(parm,model,data,cook=FALSE,...) {
  fit0 <- sem(model, data, ...)
    
  (E <- parameterEstimates(fit0))
  (E$parm <- paste(E$lhs,E$op,E$rhs,sep=""))
  (th0 <- E$est[E$parm%in%parm])
  
  ## controllo dei parametri che non ci sono nella VCOV
  if (sum(colnames(vcov(fit0))%in%parm)!=length(parm)) {
    warning(paste("Dparm estimates not available for the following parameters:",
                  paste(parm[!(parm%in%(colnames(vcov(fit0))))],collapse=", ")))
    parmS <- parm[(colnames(vcov(fit0))%in%parm)]
    ### ristrutturo S se parmS!=parm
    S2 <- diag(NA,length(parm))
    colnames(S2) <- rownames(S2) <- parm  
  } else {
    parmS <- parm
  }
  
  Dparm <- NULL; THi <- NULL
  ## tolgo le soglie e le intercette 
  ## creano problemi con i dati ordinali nel 
  ## test successivo any(fit@Fit@est[var.idx]<0))
  #(LPT <- LPT[which((LPT$op!="|")&(LPT$op!="~1")),])
  (var.idx <- which(E$op=="~~" & E$lhs==E$rhs))
  
  has.tcltk <- requireNamespace("tcltk")
  if (has.tcltk) 
    pb <- tkProgressBar("parinfluence", "Inspecting case ", 0, nrow(data))
    
  for (i in 1:nrow(data)) {
    
    if (has.tcltk) 
      setTkProgressBar(pb, i, label = sprintf(paste("Inspecting case", i,"of",nrow(data))))
    
    fit <- try(sem(model,data[-i,],...),TRUE)
    #fit <- try(sem(model,data=data[-i,]),TRUE)
    
    if (class(fit)=="try-error") {
      Dparm <- rbind(Dparm,rep(NA,length(parm)))
      THi <- rbind(THi,rep(NA,length(parm)))
    } else {
      
      (LPT <- parameterEstimates(fit))
      LPT$parm <- paste(LPT$lhs,LPT$op,LPT$rhs,sep="")
      
      if (length(which(is.na(LPT$est[var.idx])))>0) {
        LPT$est[var.idx] <- ifelse(is.na(LPT$est[var.idx]),-1,LPT$est[var.idx])
      }
      
      if ((!fit@Fit@converged)|(length(var.idx)>0L && any(LPT$est[var.idx]<0))) {
        Dparm <- rbind(Dparm,rep(NA,length(parm)))
        THi <- rbind(THi,LPT$est[LPT$parm%in%parm])
      } else {
        thi <- LPT$est[LPT$parm%in%parm]; THi <- rbind(THi,thi)
        S <- try(vcov(fit)[parmS,parmS],TRUE)
        
        if (class(S)=="try-error") {
          Dparm <- rbind(Dparm,rep(NA,length(parm)))
        } else {
          ## gestisce il caso parmS!=parm
          if (exists("S2")) {
            S2[parmS,parmS] <- S
            S <- S2            
          }
          if (length(parm)>1) {
            (S <- sqrt(diag(S)))  
          } else {
            (S <- sqrt(S))  
          }
          Di <- (th0-thi)/S
          Dparm <- rbind(Dparm,Di)                
        }
      }
    }
  }
  
  if (has.tcltk) close(pb)
  
  colnames(THi) <- colnames(Dparm) <- parm
  
  if (cook) {
    gCD <- Dparm^2
    return(list(gCD=gCD,Dparm=Dparm,THi=THi))
  } else {
    return(list(Dparm=Dparm,THi=THi))
  }
  
}

#data <- data.frame(lapply(data[,c("m1","m2","m3")],as.numeric))
#(TH <- parinfluence(parm,model,data))
#Mpar <- parinfluence(PAR,model,X)
