pv_perc <-
function(by, svydat, pvcat)
{
# pvcat  =  die namen der variablen der kategorisierten PVs als vector ?bergeben
# by     =  wonach kategorisiert werden soll
# svydat =  svyrepdesign

  
  #################################
  #### check ######################

  # hier wird ?berpr?ft ob eh alle PV-variablen 'factors' sind.
  ISFACTOR<- sapply(svydat$variables[,names(svydat$variables) %in% pvcat],is.factor)
  if(!all(ISFACTOR)){stop("Check whether all the variables provided with the 'by' statement are of class factor.")}
  
  
  ################################
  
  howmanyby <- 1:length(all.vars(by))
  
  # hier werden die zellen% berechnet inkl. SE.
  ## das funktioniert indem man svymean als funktion nimmt, und als formula einen factor!
  
  rawerg <- lapply(pvcat, function(p)
                {
                  FORM <- paste("~",p,collapse=" ")
                  erg1 <- svyby(formula = eval(parse(text=FORM)), 
                                by= by, 
                                design=svydat, FUN=svymean, na.rm=TRUE)  
                  return(erg1)
                })
  
  
  # nimm nur die numeric variables raus
  ergmat <- lapply(rawerg, function(x)
              {
                #u1   <- as.matrix(x[,sapply(x,is.numeric)]) 
                u1   <- as.matrix(x[,-howmanyby])
                wose <- grep("^se\\d*$",colnames(u1),perl=TRUE)
                # rbind garantiert mir, dass auch bei nur einer kategorie das ganze ein zeilenvector ist!
                zw1  <- rbind(u1[,wose]^2) # SE wird quadriert
                colnames(zw1) <- paste0("varvse",1:ncol(zw1))
                cbind(u1, zw1)
              })
  
  
  # mach ein 3d array mit den numerischen variablen
  dreid <- simplify2array(ergmat)  
  
  
  # bilde den Mittelwert ?ber die 3.dimension also ?ber die replizierten 10 PVs
  meanit <- apply(dreid,1:2,function(x)
              {
                mean(x)  
              })
  
  
  # bilde die var() ?ber die PVs replizierten statistiken
  sdita <- apply(dreid,1:2,var)
  wose  <- grep("^(se|varvse)\\d*$",colnames(sdita),perl=TRUE)
  sdit  <- rbind(sdita[,-wose]) # nimmt nur die means raus (weil mich nur die VAR der SE interessiert)
  colnames(sdit) <- paste0("varrep",1:ncol(sdit))
              
  
  # nimmt hier jetzt nur die mittlere varianz raus
  VARV <- rbind(meanit[,grep("^varvse\\d*$",colnames(meanit))])
  # endberechnung
  SEend <- sqrt(VARV + ((1 + 1/length(pvcat)) * sdit))
  # beschriftung
  endergt1  <- rbind(meanit[,-wose])
  #colnames(endergt1) <- gsub(pvcat[1],"Percent ",colnames(endergt1))
  colnames(endergt1) <- gsub(pvcat[1],"",colnames(endergt1))
  
  colnames(SEend) <- gsub("Percent ","SE ", colnames(endergt1))
  
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
            #splitnames <- as.list((rownames(endergv))
            splitndf   <- data.frame(matrix(unlist(splitnames),ncol=length(splitnames[[1]]), byrow=TRUE))
            # colnames(splitndf) <- paste0("Group",1:length(all.vars(by)))
            enderg             <- data.frame(splitndf,endergv)
            colnames(enderg)   <- c(paste0("Group",1:length(all.vars(by))),colnames(endergv))
            rownames(enderg)   <- NULL    
            }
  

  
  endnum <- enderg[ ,-(1:length(all.vars(by)))]
  numbCAT <- ncol(endnum)/2 # number of categories
  endnam <- enderg[rep(1:nrow(enderg),each=numbCAT) , (1:length(all.vars(by)))]
  
  eins <- apply(endnum,1,function(x) list(mymat=matrix(x,nrow=numbCAT)))
  zwei <- lapply(eins,function(A)A[[1]])
  drei <- do.call("rbind",zwei)
  colnames(drei) <- c("Proportion","Proportion:SE")
  
  CN <- colnames(enderg)[-howmanyby][1:numbCAT]
  
  endergF <- data.frame(endnam, CN , drei)
  colnames(endergF)[1:(length(howmanyby)+1)] <- paste0("Group",1:(length(howmanyby)+1))
  rownames(endergF) <- NULL
  
  
  #########################################################
  ########### additional information ######################
  ######################################################### 
  # country ID # kommt erst in der user-main-function dazu!
  # N of cases
  # sum of SPFWT0
  # nimm: svydat.sub als grundlage, sonst gibts probleme
  # wenn fehlende werte in den NA's sind
  
  
  # hier sollte man sich überlegen ob man die zusatzinfos wirklich will, und wenn ja wie.
  # denn die ANZAHL DER PERSONEN JE GRUPPE ist nichtmehr so einfach zu bestimmen, weil die Gruppen ja 
  # AUCH plausible values sind!!!
  
  
  ADC <- additional_compPV(by=by, svydat=svydat, pvcat=pvcat)
  
  colnames(ADC$Ncases)     <- c(paste0("Group",1:(length(howmanyby)+1)),"Number.of.cases")
  colnames(ADC$Sumweights) <- c(paste0("Group",1:(length(howmanyby)+1)),"Sum.of.weights")
  

  return(list(endergF=endergF, Ncases=ADC$Ncases, Sumweights=ADC$Sumweights))

}
