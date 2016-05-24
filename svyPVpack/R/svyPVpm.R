svyPVpm <-
function(by, svydat, pvs, colN=FALSE)  
#function(by, svydat, pvs, addcountry=TRUE)
{

# check input

checkds <- data.frame(svydat$variables[,all.vars(by)])
only1gr <- all(sapply(checkds,function(A)length(unique(A))) == 1)




if(only1gr)
  { # fuer den Fall dass es nur eine Gruppe gibt
    
  # do what you have to do
  meanerg <- pv_by(by=by, svydat=svydat, pvs=pvs,FUNC="svymean")
  sderg   <- pv_by(by=by, svydat=svydat, pvs=pvs,FUNC="svyby_SD")  
  
  no.vars <- length(all.vars(by))
  
  nams       <- paste(unique(gsub("\\d*","",pvs,perl=TRUE)),collapse="|")
  meannam    <- colnames(meanerg$enderg)
  meannam[no.vars + 2] <- paste0("mean:",meannam[no.vars + 2])
  colnames(meanerg$enderg)[(no.vars+1):(no.vars+2)] <- paste0(nams,"_",meannam[(no.vars+1):(no.vars+2)])
  
  
  sdnam <- c("stddev","stddev:SE")
  colnames(sderg$enderg)[(no.vars+1):(no.vars+2)] <- paste0(nams,"_",sdnam)
  
  cellp <- data.frame("Group1"=meanerg$enderg[1,1],"Percent"=100, "SE.Percent"=0)
  
  
  ### M E R G E
  
  mergelist <- list(meanerg$Ncases,meanerg$Sumweights,cellp,meanerg$enderg,sderg$enderg)
  pm        <- mergeALL(mergelist)


    
  } else { # fuer den Fall von mehreren Gruppen
    
        # do what you have to do
        meanerg <- pv_by(by=by, svydat=svydat, pvs=pvs,FUNC="svymean")
        sderg   <- pv_by(by=by, svydat=svydat, pvs=pvs,FUNC="svyby_SD") ##! problem
        cellp   <- opv_perc(by=by, svydat=svydat, pvs=pvs)  
        
        no.vars <- length(all.vars(by))
        
        nams       <- paste(unique(gsub("\\d*","",pvs,perl=TRUE)),collapse="|")
        meannam    <- colnames(meanerg$enderg)
        meannam[no.vars + 2] <- paste0("mean:",meannam[no.vars + 2])
        colnames(meanerg$enderg)[(no.vars+1):(no.vars+2)] <- paste0(nams,"_",meannam[(no.vars+1):(no.vars+2)])
        
        sdnam <- c("stddev","stddev:SE")
        colnames(sderg$enderg)[(no.vars+1):(no.vars+2)] <- paste0(nams,"_",sdnam)
        
        
        ### M E R G E
        
        mergelist <- list(meanerg$Ncases,meanerg$Sumweights,cellp,meanerg$enderg,sderg$enderg)
        pm        <- mergeALL(mergelist)        
        

        }
 

### um die ordnung der factors gleich zu lassen (vor allem wichtig bezogen auf grafiken) wird hier nochmal umgeordnet so wie es im datensatz ?blich ist
mybys <- all.vars(by)

#### NEU um keine NAs zu erzeugen ########
# facordall <- mapply(function(x,number){
#   
#   if(is.factor(pm[,number]))
#     {
#     factor(pm[,number], levels=levels(svydat$variables[[x]]))
#     } else 
#       {
#       factor(pm[,number], levels=levels(as.factor(svydat$variables[[x]])))
#       }
# }, x=mybys, number=1:length(mybys),SIMPLIFY=FALSE)


#facordallDF <- data.frame(facordall)
facordallDF <- fALL(mybys,pm, svydat)

pm[,1:length(mybys)] <- facordallDF


if(colN)
{
  colnames(pm)[1:length(mybys)] <- c(mybys) 
  
}

# if(addcountry)
# {  
# pm  <- data.frame("Country"=unique(svydat$variables$CNTRYID), pm) 
# }


return(pm)
}
