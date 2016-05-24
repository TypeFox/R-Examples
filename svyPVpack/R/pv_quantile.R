pv_quantile <-
function(by, svydat, pvs, quantile, ...)
{
  
  
  if(length(quantile) < 2) stop("Please submit a vector of quantiles with >=2 elements!")
  
  
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
                  svydat, FUN=svyquantile,
                  na.rm=TRUE, na.rm.by=TRUE,quantile=quantile,...)
    return(erg1)
  })
  
  # nimm nur die numeric variables raus
  ergmat <- lapply(rawerg, function(x)
  {
    u1 <- as.matrix(x[,-howmanyby])
    
    wose <- grep("^se\\d*$",colnames(u1),perl=TRUE)
    # cbind garantiert mir, dass auch bei nur einer kategorie das ganze ein zeilenvector ist!
    zw1  <- rbind(u1[,wose]^2) # SE wird quadriert  ##!
    colnames(zw1) <- paste0("varvse",1:ncol(zw1))
    cbind(u1, zw1)
    
  })
  
  
  # mach ein 3d array mit den numerischen variablen
  dreid <- simplify2array(ergmat)
  
  
  # summiere ueber die 3.dimension also ?ber die replizierten 10 PVs
  meanit <- apply(dreid,1:2,function(x)
  {
    mean(x)  
  })
  
  
  
  # bilde die var() ?ber die PVs replizierten statistiken
  sdita <- apply(dreid,1:2,var)
  wose  <- grep("^(se|varvse)\\d*$",colnames(sdita),perl=TRUE)
  sdit  <- rbind(sdita[,-wose]) # nimmt nur die means raus (weil mich nur die VAR der SE interessiert) ##!
  colnames(sdit) <- paste0("varrep",1:ncol(sdit))
  
  
  # nimmt hier jetzt nur die mittlere varianz raus
  VARV <- rbind(meanit[,grep("^varvse\\d*$",colnames(meanit))]) ##!
  # endberechnung
  SEend <- sqrt(VARV + ((1 + 1/length(pvs)) * sdit))
  # beschriftung
  endergt1  <- rbind(meanit[,-wose])
  
  
  
  quN <- paste0("q",quantile)
  
  colnames(endergt1) <- quN
  
  colnames(SEend) <- paste0("se:",quN)
  
  # fertig
  endergv  <- cbind(endergt1,SEend)
  
  
  
  
  if(dim(dreid)[1] == 1)
  {
    grnam    <- unique(svydat$variables[,all.vars(by)])
    splitndf <- data.frame("Group1"=grnam)
    enderg   <- data.frame(splitndf,endergv)
    colnames(enderg)   <- c("Group1",colnames(endergv))
    rownames(enderg)   <- NULL
    
  } else {
    # necessary because the order of by-variable combination is different between the functions
    # ------------------------------------------------------------------
    splitnames <- strsplit(rownames(endergv),"\\.")
    splitndf   <- data.frame(matrix(unlist(splitnames),ncol=length(splitnames[[1]]), byrow=TRUE))
    # colnames(splitndf) <- paste0("Group",1:length(all.vars(by)))
    enderg             <- data.frame(splitndf,endergv)
    colnames(enderg)   <- c(paste0("Group",1:length(all.vars(by))),colnames(endergv))
    rownames(enderg)   <- NULL    
  }
  
  
  #########################################################
  ########### additional information ######################
  ######################################################### 
  # country ID # kommt erst in der user-main-function dazu!
  # N of cases
  # sum of SPFWT0
  # nimm: svydat.sub als grundlage, sonst gibts probleme
  # wenn fehlende werte in den NA's sind  
  
  ADC <- additional_comp(by=by, svydat=svydat.sub)
  
  colnames(ADC$Ncases)     <- c(paste0("Group",1:(length(howmanyby))),"Number.of.cases")
  colnames(ADC$Sumweights) <- c(paste0("Group",1:(length(howmanyby))),"Sum.of.weights")
  
  
  return(list(Ncases=ADC$Ncases, Sumweights=ADC$Sumweights,enderg=enderg))
  
 
  
}
