additional_compPV <-
function(by, svydat, pvcat)  
{
  
  #########################################################
  ###########        N cases       ######################
  #########################################################   
  tabnamsplit <- all.vars(by)
  
  NCply <- lapply(pvcat, function(p)
            {
              #FORM <- paste("~",p,collapse=" ")
              tabnamsplitA <- c(tabnamsplit,p)
              aggcom       <- paste("list(",paste0("svydat$variables$",tabnamsplitA,collapse=","),")",sep="")
              Ncases       <- aggregate(svydat$variables[,1], eval(parse(text=aggcom)), FUN=length)
              colnames(Ncases)[length(tabnamsplitA)+1] <- p
              return(Ncases)
            })
            
  

  anzvar <- length(tabnamsplit)+1
  
  
  gemerg   <- mergeALL(NCply)
  gemerg1n <- sapply(gemerg[,-c(1:anzvar)], function(Y) ifelse(is.na(Y),0,Y)) # nimmt nur die anzahlen raus
  roM      <- rowMeans(gemerg1n)
  NCASALL  <- data.frame(gemerg[,c(1:anzvar)],"Ncases" = round(roM))
  
  
  #########################################################
  ###########     Sum of weights    ####################
  #########################################################   
  
  
  SWply <- lapply(pvcat, function(p)
            {
              #FORM <- paste("~",p,collapse=" ")
              tabnamsplitA <- c(tabnamsplit,p)
              aggcom       <- paste("list(",paste0("svydat$variables$",tabnamsplitA,collapse=","),")",sep="")
              Sumweights      <- aggregate(svydat$pweights, eval(parse(text=aggcom)), FUN=sum)
              colnames(Sumweights)[length(tabnamsplitA)+1] <- p
              return(Sumweights)
            })
            
  gemSW     <- mergeALL(SWply)
  gemSW1n   <- sapply(gemSW[,-c(1:anzvar)], function(Y) ifelse(is.na(Y),0,Y)) # nimmt nur die anzahlen raus
  roMW      <- rowMeans(gemSW1n)
  SWALL     <- data.frame(gemSW[,c(1:anzvar)],"Ncases" = round(roMW,1))

  
  return(list(Ncases=NCASALL, Sumweights=SWALL))
  
}
