getXList <- function(result, group=vector(), file="") {
  modtype <- result$currModel@modellist[[1]]@mod_type
  resultlist <- result$currModel@fit@resultlist
  m <- result$currModel@modellist
  t <- result$currTheta
  XList <- vector("list", length=length(m))
  tauList <- muList <- list()

  if(modtype == "kin") {
    groups <- result$currModel@groups
    multimodel <- result$currModel
    f1<-function(x){x[[1]]}   
    f2<-function(x){x[[2]]}
    ## will plot the first concentration from each dataset
    grtoplot <- vector("list", length(m))
    
    for(i in 1:length(m)) {
      if(length(group)==0) {
        cnt <- 1
        notfound <- TRUE 
        while(notfound) {
          for(j in 1:length(groups[[cnt]])) {
            if(groups[[cnt]][[j]][2] == i) {
              grtoplot[[i]]<-list(groups[[cnt]],j)
              notfound<-FALSE
            }
          }
          cnt<-cnt+1
        }	
      }
      else{
        grtoplot[[i]] <- list(groups[[group]], 1)
      }
    }
    for(i in 1:length(m)) {
      group <- grtoplot[[i]][[1]]
      place <-  grtoplot[[i]][[2]]
      dset <- group[[place]][2]
      irfmu <- unlist(lapply(resultlist[[i]]@irfvec, f1))
      irftau <- unlist(lapply(resultlist[[i]]@irfvec, f2))
      muList[[i]] <- irfmu 
      tauList[[i]] <- irftau
      XList[[i]] <- getConToPlot(getKinConcen(
      group, multimodel, t, oneDS = place), m[[i]]@cohspec, m[[i]]@cohcol)
    }
  }
  if(modtype=="spec") {
    for(i in 1:length(m)) {
      if(m[[i]]@timedep)
        specpar <- specparF(t[[i]]@specpar, m[[i]]@x[group], 
                            group, m[[i]]@specref, m[[i]]@specdispindex, 
                            t[[i]]@specdisppar, parmufunc = m[[i]]@parmufunc)
              else 
                specpar <- t[[i]]@specpar 
      
      XList[[i]] <- doClpConstr(specModel(specpar, m[[i]]),
                                clp_ind = group, 
                                clpCon = m[[i]]@clpCon, clpequ = t[[i]]@clpequ, 
                                num_clpequ = length(m[[i]]@clpequspec), 
                                usecompnames0 = m[[i]]@usecompnames0, 
                                usecompnamesequ = m[[i]]@usecompnamesequ)
      
    }
  }
  if(modtype=="mass") {
    for (i in 1:length(m)) {
     XList[[i]] <- compModelMass(theta = t[[i]], model = m[[i]])
   }
  }
  for(i in 1:length(XList)) {
    xdim <- dim(XList[[i]])
    attributes(XList[[i]])<-NULL
    dim(XList[[i]])<- xdim
    if(file!="") 
      for(i in 1:length(XList))
        write.table(XList[[i]], file=paste(file,
                                  "_concen_dataset_", i, ".txt",
                                  sep=""), quote = FALSE,
                    row.names = m[[i]]@x)
  }
  
  XList
}

getX <- function(result, group=vector(), dataset=1, file="",lreturnA=FALSE,lreturnC=FALSE) {
  # minimal model of getXList for one dataset case 
  modtype <- result$currModel@modellist[[1]]@mod_type
  resultlist <- result$currModel@fit@resultlist
  m <- result$currModel@modellist
  t <- result$currTheta
  tauList <- muList <- list()
  if(modtype == "kin") {
    groups <- result$currModel@groups
    multimodel <- result$currModel
    f1<-function(x){x[[1]]}   
    f2<-function(x){x[[2]]}
    ## will plot the first concentration from each dataset
    grtoplot <- list()
    if(length(group)==0) {
      cnt <- 1
      notfound <- TRUE 
      while(notfound) {
        for(j in 1:length(groups[[cnt]])) {
          if(groups[[cnt]][[j]][2] == dataset) {
            grtoplot[[dataset]]<-list(groups[[cnt]],j)
            notfound<-FALSE
          }
        }
        cnt<-cnt+1
      }
    }
    else{
      grtoplot[[dataset]] <- list(groups[[group]], 1)
    }
    group <- grtoplot[[dataset]][[1]]
    place <-  grtoplot[[dataset]][[2]]
    dset <- group[[place]][2]
    irfmu <- unlist(lapply(resultlist[[dataset]]@irfvec, f1))
    irftau <- unlist(lapply(resultlist[[dataset]]@irfvec, f2))
      muList[[dataset]] <- irfmu 
    tauList[[dataset]] <- irftau
    if (lreturnA||lreturnC)
       { group <- list(c(1,dataset))
        XList <- getKinConcen(group, multimodel, t, oneDS = 1,lreturnA=lreturnA,lreturnC=lreturnC)
        }
#        XList <- getKinConcen(group, multimodel, t, oneDS = place,lreturnA=lreturnA,lreturnC=lreturnC)
   else
    XList <- getConToPlot(getKinConcen(group, multimodel, t, oneDS = place),
                          m[[dataset]]@cohspec, m[[dataset]]@cohcol)
  }
  if(modtype=="spec") {
    if(m[[dataset]]@timedep)
      specpar <- specparF(t[[dataset]]@specpar, m[[dataset]]@x[group], 
                    group, m[[dataset]]@specref, m[[dataset]]@specdispindex, 
                    t[[dataset]]@specdisppar,
                          parmufunc = m[[dataset]]@parmufunc)
    else 
      specpar <- t[[dataset]]@specpar 
    
    XList <- doClpConstr(specModel(specpar, m[[dataset]]),
                         clp_ind = group, 
                         clpCon = m[[dataset]]@clpCon,
                         clpequ = t[[dataset]]@clpequ, 
                         num_clpequ = length(m[[dataset]]@clpequspec), 
                         usecompnames0 = m[[dataset]]@usecompnames0, 
                         usecompnamesequ = m[[dataset]]@usecompnamesequ)
    
  }
  if(modtype=="mass") 
    XList <- compModelMass(theta = t[[dataset]], model = m[[dataset]])
  xdim <- dim(XList)
  attributes(XList)<-NULL
  dim(XList)<- xdim
  if(file!="") 
    write.table(XList, file=paste(file,
                         "_concen_dataset_", dataset, ".txt",
                         sep=""), quote = FALSE,
                row.names = m[[dataset]]@x)
  XList
}
