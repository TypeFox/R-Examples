pv_by <-
function(by, svydat, pvs, FUNC)
{
  
  avar <- c(pvs, all.vars(by))
  
  dum1 <- apply(svydat$variables[,avar], 1, function(x) any(is.na(x)))
  svydat.sub <- subset(svydat, dum1 == FALSE)
  #svydat.sub <- svydat
  
  # um nachher die richtige anzahl an spalten rausnehmen zu können (ergmat)
  howmanyby <- 1:length(all.vars(by))
  
  ## Hier werden die PV durchgerattert
  rawerg <- lapply(pvs, function(p)
      {
      FORM <- paste("~",p,collapse=" ")
      
      erg1 <- svyby(formula = eval(parse(text=FORM)), 
                    by= by, 
                    svydat, FUN=eval(parse(text=FUNC)),
                    na.rm=TRUE, na.rm.by=TRUE)
      return(erg1)
      })

  # nimm nur die numeric variables raus
  ergmat <- lapply(rawerg, function(x)
              {
                u1 <- as.matrix(x[,-howmanyby]) 
                zw1 <- u1[,"se"]^2
                cbind(u1, varvse = zw1)
              })
              

  # mach ein 3d array mit den numerischen variablen
  dreid <- simplify2array(ergmat)
  
  
  # summiere ueber die 3.dimension also ?ber die replizierten 10 PVs
  meanit <- apply(dreid,1:2,function(x)
              {
               mean(x)  
              })
              
  # bilde die var() ?ber die PVs replizierten statistiken
  sdit <- apply(dreid,1:2,var)[,1]
  
  # dranhaengen um nachher geschickt transform machen zu k?nnen
  meanit2 <- cbind(meanit,varrep=sdit)
  #cat("hier bin ich")
  # addiere jetzt auf geeignete weise die imputationsvarianz
  se_end <- sqrt(meanit2[,"varvse"] + ((1 + 1/length(pvs)) *sdit))

  meanit3 <- cbind(meanit2,se_end)
  endergv  <- meanit3[,c(1,ncol(meanit3))]
  
  if(dim(dreid)[1] == 1)
        {
        names(endergv) <- c("mean", "SE")  
        } else {
               colnames(endergv) <- c("mean", "SE")
               }
  
 
  #########################################################
  ########### additional information ######################
  ######################################################### 
  # country ID # kommt erst in der user-main-function dazu!
  # N of cases
  # sum of SPFWT0
  # nimm: svydat.sub als grundlage, sonst gibts probleme
  # wenn fehlende werte in den NA's sind
  
  
  ADC <- additional_comp(by=by,svydat=svydat.sub)
  
  
  if(dim(dreid)[1] == 1)
  {
  
    
  endergvM <- matrix(endergv,ncol=2,dimnames=list(1,names(endergv)))
  enderg   <- data.frame("Group1"=unique(svydat$variables[,all.vars(by)]),endergvM)  
    
    
    
    
  } else {
    

    # necessary because the order of by-variable combination is different between the functions
    # ------------------------------------------------------------------
    splitnames <- strsplit(rownames(endergv),"\\.")
    splitndf   <- data.frame(matrix(unlist(splitnames),ncol=length(splitnames[[1]]), byrow=TRUE))
    colnames(splitndf) <- paste0("Group",1:length(all.vars(by)))
    enderg             <- data.frame(splitndf,endergv)
    rownames(enderg)   <- NULL    
 
  }
  
  colnames(ADC$Ncases)     <- c(paste0("Group",1:length(all.vars(by))),"Number.of.cases")
  colnames(ADC$Sumweights) <- c(paste0("Group",1:length(all.vars(by))),"Sum.of.weights")
  
  
  
  return(list(enderg=enderg, Ncases=ADC$Ncases, Sumweights=ADC$Sumweights))
  
}
