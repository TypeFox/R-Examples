add.center <- function(design, ncenter, distribute=NULL, ...){
    if (!"design" %in% class(design)) stop("design must be of class design")
    di <- design.info(design)
    if (missing(ncenter)) stop("ncenter must be specified")
    if (!(substr(di$type,1,4)=="FrF2" | substr(di$type,1,2)=="pb" | (substr(di$type,1,14)=="full factorial" & all(di$nlevels==2))))
       stop("center points only available for FrF2 and pb type designs")
    if (di$type=="FrF2.splitplot") 
       stop("currently, center points for splitplot designs are not supported")
    if (di$type=="FrF2.param") 
       stop("currently, center points for long version parameter designs are not supported")
    if (length(grep("center",di$type))>0)
       stop("design has center points already")
    if (!is.numeric(ncenter)) stop("ncenter must be a number")
    if (!length(ncenter)==1) stop("ncenter must be a number")
    if (!ncenter==floor(ncenter)) stop("ncenter must be an integer number")
    if (is.null(distribute)){
      if (ncenter==0 | !di$randomize) distribute <- 1
         else distribute <- min(ncenter, 3)}
    if (!is.numeric(distribute)) stop("distribute must be a number")
    if (!distribute==floor(distribute)) stop("distribute must be an integer number")
    if (distribute < 1 | distribute > max(1,min(di$nruns+1, ncenter)))
       stop("distribute must be at least 1 and at most min(ncenter, nruns+1 of the design)")
    if (di$randomize & distribute==1 & ncenter > 1) warning("running all center point runs together is usually not a good idea.")
    
    if (any(is.na(sapply(factor.names(design),"is.numeric"))))
       stop("Center points are implemented for experiments with all factors quantitative only.")
       

    design <- qua.design(design, quantitative="all")
    di <- design.info(design)
    clevels <- sapply(factor.names(design), "mean")
    di$ncube <- di$nruns
    di$ncenter <- ncenter
    if (!is.null(di$blocksize)) di$nruns <- di$nruns + di$nblocks*ncenter else
    di$nruns <- di$nruns+ncenter
    
    desnum <- desnum(design)
    ro <- run.order(design)
    if (!all(ro$run.no == sort(ro$run.no)) & (di$replications>0 | !is.null(di$blocksize))) 
        stop("center points cannot be added to designs that are not in original order")
    hilf <- undesign(design)

    ### 15 Feb 2013: removed restriction to bbreps==1 (as now always factor)
    if (!is.null(di$blocksize)) block.contrasts <- contrasts(hilf[,di$block.name])
    
    ## positions for adding center points (added after these positions)
    ## !!! blocks, replications and repeated measurements to be treated correctly
    if (!is.null(di$blocksize)){
         ## distribute within blocks
        if (distribute==1)
            addpos <- di$blocksize 
        else
            addpos <- c(0,round((1:(distribute-1))*di$blocksize/(distribute-1)))
    }
    else{
        if (distribute==1) 
             addpos <- di$ncube 
        else
          addpos <- c(0,round((1:(distribute-1))*di$ncube/(distribute-1)))
    }
    ## determine numbers of center points to be added at each position
        nrest <- ncenter%%distribute
        if (nrest==0) nadd <- rep(round(ncenter/distribute),distribute)
        if (nrest==1) nadd <- c(1,rep(0,distribute-1))+rep(floor(ncenter/distribute),distribute)
        if (nrest==2) nadd <- c(1,rep(0,distribute-2),1)+rep(floor(ncenter/distribute),distribute)
        if (nrest>2) nadd <- rep(floor(ncenter/distribute),distribute) + sample(c(rep(1,nrest),rep(0,distribute-nrest)))
        ## repeated measurements also for center points (assuming they are done for measurement accuracy)
        if (di$repeat.only) if (!is.null(di$blocksize)) nadd <- nadd*di$wbreps else nadd <- nadd*di$replications

    ## split vector for design and run order data frame
        if (!is.null(di$blocksize)){
          if (di$wbreps==1 | di$repeat.only)
             together <- paste(rep(1:(di$nblocks*di$bbreps),each=di$blocksize*di$wbreps),
                          rep(cut(1:di$blocksize,unique(c(0,addpos))),each=di$wbreps,di$nblocks*di$bbreps),sep=".")
          else
             together <- paste(rep(1:(di$nblocks*di$bbreps*di$wbreps),each=di$blocksize),
                          rep(cut(1:di$blocksize,unique(c(0,addpos))),di$nblocks*di$bbreps*di$wbreps),sep=".")
             }
        else {if (di$replications==1 | (di$replications > 1 & di$repeat.only)) 
                together <- cut(rep(1:di$ncube,each=di$replications), unique(c(0,addpos)))
              else together <- paste(rep(1:di$replications,each=di$ncube),rep(cut(1:di$ncube,unique(c(0,addpos))),di$replications),sep=".")
                }
        ## the following command prevents split from reordering the data
        together <- factor(together, levels=unique(together))
   ## split the data frame
        getrennt <- split(hilf,together)
          ## added 15 Feb 2013
          ## obtaine levels for the new factors including the center points
          rostdlev <- levels(ro$run.no.in.std.order)
          rostdrplev <- levels(ro$run.no.std.rp)
          levgleich <- identical(ro$run.no.in.std.order,ro$run.no.std.rp)
          
          ## extend levels to include all center point related ones
          ## for run.no.in.std.order
          if (all(sapply(rostdlev, function(obj) length(grep(".", obj, fixed=TRUE)) > 0 ) )){
              punkt1 <- sapply(rostdlev, function(obj) regexpr(".",obj,fixed=TRUE))
              if (all(sapply(substr(rostdlev,punkt1+1,999), function(obj) length(grep(".", obj, fixed=TRUE))>0)))
              punkt2 <- sapply(substr(rostdlev, punkt1+1, 999), function(obj) regexpr(".", obj, fixed=TRUE))+punkt1
              else punkt2 <- 1000  ## repeated measurements or replications only
              if (punkt2[1] < 1000) rostdlev <- c(unique(paste("0", substr(rostdlev, punkt1+1, punkt2-1), "0", sep=".")),
                  rostdlev) else rostdlev <- c(unique(paste("0", substr(rostdlev, punkt1+1, punkt2-1), sep=".")),
                  rostdlev)         
          } 
          else rostdlev <- c("0", rostdlev)
          ## for run.no.std.rp
          if (levgleich) rostdrplev <- rostdlev else{
          if (all(sapply(rostdrplev, function(obj) length(grep(".", obj, fixed=TRUE))>0))){
              punkt1 <- sapply(rostdrplev, function(obj) regexpr(".",obj,fixed=TRUE))
              if (all(sapply(substr(rostdrplev,punkt1+1,999), function(obj) length(grep(".", obj, fixed=TRUE))>0))){
              punkt2 <- sapply(substr(rostdrplev, punkt1+1, 999), function(obj) regexpr(".", obj, fixed=TRUE))+punkt1
              punkt3 <- sapply(substr(rostdrplev, punkt2+1, 999), function(obj) regexpr(".", obj, fixed=TRUE))+punkt2
              if (all(sapply(substr(rostdrplev,punkt3+1,999), function(obj) length(grep(".", obj, fixed=TRUE))>0)))
                 punkt4 <- sapply(substr(rostdrplev, punkt3+1, 999), function(obj) regexpr(".", obj, fixed=TRUE))+punkt3 else
                 punkt4 <- 1000
              if (punkt4[1] == 1000) rostdrplev <- c(unique(paste("0", 
                   substr(rostdrplev, punkt1+1, punkt2-1), "0", 
                   substr(rostdrplev, punkt3+1, punkt4-1), sep=".")), rostdrplev)
                 else {if (di$repeat.only) rostdrplev <- c(unique(paste("0", 
                   substr(rostdrplev, punkt1+1, punkt2-1), "0", 
                   substr(rostdrplev, punkt3+1, punkt4-1), 1:di$wbreps, sep=".")), rostdrplev)
                   else rostdrplev <- c(unique(paste("0", 
                   substr(rostdrplev, punkt1+1, punkt2-1), "0", 
                   substr(rostdrplev, punkt3+1, punkt4-1), substr(rostdrplev, punkt4+1, 999), sep=".")), rostdrplev)
              }
              }
              else {punkt2 <- 1000 ## repeated measurements or replications only
              rostdrplev <- c(unique(paste("0", substr(rostdrplev, punkt1+1, punkt2-1), 
                     sep=".")), rostdrplev)
              }            
          } 
          else rostdrplev <- c("0", rostdrplev)
          }
          ## end of added 15 Feb 2013
        ro.getrennt <- split(ro,together)
        ## wird das überhaupt gebraucht???
        ro$run.no.in.std.order <- as.character(ro$run.no.in.std.order)
        ro$run.no.std.rp <- as.character(ro$run.no.std.rp)
        
        blockid <- rep("1", length(getrennt))

        von <- cumsum(c(1,lapply(getrennt, "nrow")))[-(length(getrennt)+1)]
        bis <- cumsum(lapply(getrennt, "nrow"))
        
        ## added 15 Feb 2013
        std <- lapply(ro.getrennt, function(obj) {
           hilf <- t(as.matrix(sapply(obj$run.no.in.std.order, function(obj2) unlist(strsplit(as.character(obj2),".",fixed=TRUE)))))
           if (nrow(hilf)==1) hilf <- t(hilf)
           hilf
           })
        stdrp <- lapply(ro.getrennt, function(obj) {
            hilf <- t(as.matrix(sapply(obj$run.no.std.rp, function(obj2) unlist(strsplit(as.character(obj2),".",fixed=TRUE)))))
            if (nrow(hilf)==1) hilf <- t(hilf)
            hilf
            })
        ## get rid of rp-entry
        if (is.null(di$blocksize) & di$replications > 1 & di$repeat.only) 
            stdrp <- lapply(stdrp, function(obj) obj[,1:(ncol(obj)-1)])   
        if (!is.null(di$blocksize)){ 
            blockid <- sapply(getrennt, function(obj) as.character(obj[1,di$block.name]))
            ## get rid of rp-entry in stdrp for blocked designs
            if (di$wbreps > 1 & di$repeat.only) 
                stdrp <- lapply(stdrp, function(obj) obj[,1:(ncol(obj)-1)])
            }
        ## end of added 15 Feb 2013
        if (!is.null(di$blocksize)) 
            blockid <- sapply(getrennt, function(obj) as.character(obj[1, 
               di$block.name]))
        if (di$replications > 1 & is.null(di$blocksize) & !di$repeat.only){ 
            blockid <- rep(1:di$replications, each=max(1,distribute-1))
            }
        if (!is.null(di$blocksize)) if (di$wbreps > 1 & !di$repeat.only){
            blockid <- paste(blockid, rep(1:di$wbreps, 
                    each=max(1,distribute-1), 
                    times=di$nblocks*di$bbreps),sep=".")
            }
    ## create center point data frames to be interleaved
    #    if (!is.null(di$nblocks)) no.center.groups <- distribute*di$nblocks*di$bbreps else 
    #        no.center.groups <- distribute*di$replications
        centers <- vector("list")
        ros <- lapply(1:length(std), function(obj) vector("list"))
    
        ## centers and ros for WITHIN each block
        more <- setdiff(colnames(hilf), c(di$block.name, names(di$factor.names)))
            ## columns other than design factors, e.g. responses
      if (ncenter>0){
        for (i in 1:distribute){
              cnext <- data.frame(matrix(clevels, ncol=di$nfactors, nrow=nadd[i], 
                                                 byrow=TRUE, dimnames=list(NULL, names(di$factor.names))))
              if (length(more)>0) cnext <- cbind(cnext, matrix(NA,nrow=nadd[i],ncol=length(more),dimnames=list(NULL,more)))
              centers <- c(centers, list(cnext))
              
              ## simultaneous creation of inserts for all blocks via lapply
              if (!is.null(di$blocksize)){
                 stdnext <- lapply(1:length(std), function(obj) rep(paste(c(0,std[[obj]][1,2],0),collapse="."), nadd[i]))
                 stdrpnext <- stdnext
                 if (di$bbreps>1 | (di$wbreps>1 & !di$repeat.only)) stdrpnext <- 
                       mapply(paste, stdrpnext, lapply(stdrp, function(obj) paste(obj[1,4], collapse=".")), 
                       sep=".", SIMPLIFY=FALSE)
                 if (di$bbreps>1 & di$wbreps>1 & !di$repeat.only) stdrpnext <- 
                       mapply(paste, stdrpnext, lapply(stdrp, function(obj) paste(obj[1,5], collapse=".")), 
                       sep=".", SIMPLIFY=FALSE)
                 if (di$wbreps>1 & di$repeat.only) 
                       stdrpnext <- lapply(stdrpnext, function(obj) paste(obj, rep(1:di$wbreps, round(nadd[i]/di$wbreps)),sep="."))
                 }
              else {
                stdnext <- stdrpnext <- lapply(1:length(std), function(obj) rep(0, nadd[i]) )
                if (di$replications > 1){ 
                if (!di$repeat.only) 
                  stdrpnext <- mapply(paste, stdrpnext, lapply(stdrp, function(obj) obj[1,2]),sep=".", SIMPLIFY=FALSE)
                  else stdrpnext <- lapply(stdrpnext, function(obj) paste(obj, rep(1:di$replications,round(nadd[i]/di$replications)),sep="."))
                }
                }
                stdnext <- lapply(stdnext, function(obj) factor(obj, levels=rostdlev))
                stdrpnext <- lapply(stdrpnext, function(obj) factor(obj, levels=rostdrplev))
                bothnext <- mapply(function(...) list(data.frame(...)), stdnext, lapply(1:length(std), function(obj) rep(0,nadd[i])),
                    stdrpnext, SIMPLIFY=FALSE)
                bothnext <- lapply(bothnext, function(obj){colnames(obj[[1]]) <- c("run.no.in.std.order","run.no","run.no.std.rp"); obj})
                ros <- mapply(c,ros, bothnext, SIMPLIFY=FALSE)
              }
    ## interleave
         ## first center point entry of first (replication) block
          if (!is.null(di$blocksize)) 
              new <- cbind(matrix(hilf[1,di$block.name],nrow=nadd[1],ncol=1, dimnames=list(NULL, di$block.name)),centers[[1]])
          else new <- centers[[1]]
          ronew <- ros[[1]][[1]] 
          blockids <- unique(blockid)
          if (distribute==1){
            ## center points appended at end of each (replication) block
            if (length(blockids)==1){
               ## one (replication) block only
               new <- rbind(hilf, new)
               ronew <- rbind(ro, ronew)
               }
            else{
               ## more than one (replication) block
               ## as distribute equals 1, blockids == blockid
               new <- rbind(getrennt[[1]], new)
               ronew <- rbind(ro.getrennt[[1]], ronew)
               for (i in 2:length(blockids)){
                  if (!is.null(di$blocksize)) 
                  new <- rbind(new, getrennt[[i]], 
                      cbind(matrix(getrennt[[i]][1,di$block.name],nrow=nadd[1],ncol=1, dimnames=list(NULL, di$block.name)),centers[[1]]))
                  else new <- rbind(new, getrennt[[i]], centers[[1]])
                  ronew <- rbind(ronew, ro.getrennt[[i]], ros[[i]][[1]])
               }
              }
            }
            else{
              ### distribute > 1
              for (i in 1:length(blockids)){
              iblock <- which(blockid==blockids[i])
              if (i>1){ 
                  ## append first center points series of new block
                  ## for first block done already
                  if (!is.null(di$blocksize)) 
                  new <- rbind(new, 
                      cbind(matrix(getrennt[[iblock[1]]][1,di$block.name],nrow=nadd[1],ncol=1, dimnames=list(NULL, di$block.name)),centers[[1]]))
                  else new <- rbind(new, centers[[1]])
                  ronew <- rbind(ronew, ros[[iblock[1]]][[1]])
                  }
              for (j in 1:length(iblock)){
                  ## append next data with subsequent center points for the current (replication) block 
                  if (!is.null(di$blocksize)) 
                  new <- rbind(new, getrennt[[iblock[j]]],
                      cbind(matrix(getrennt[[iblock[j]]][1,di$block.name],nrow=nadd[j+1],ncol=1, dimnames=list(NULL, di$block.name)),centers[[j+1]]))
                  else new <- rbind(new, getrennt[[iblock[j]]], centers[[j+1]])
                  ronew <- rbind(ronew, ro.getrennt[[iblock[j]]], ros[[iblock[1]]][[j+1]])
                  }
              }
              }
        if (!is.null(di$blocksize)){
          new[[di$block.name]] <- factor(new[[di$block.name]], levels(design[[di$block.name]]))
          if (exists("block.contrasts")) contrasts(new[,di$block.name]) <- block.contrasts
          desnum <- model.matrix(formula(paste("~",paste(c(di$block.name,names(di$factor.names)),collapse="+"))),data=new)[,-1]
        }
        else
        desnum <- model.matrix(formula(paste("~",paste(names(di$factor.names),collapse="+"))),data=new)[,-1]
        if (length(setdiff(colnames(new),c(di$block.name, names(di$factor.names))))>0){ 
               anhaeng <- as.matrix(new[,setdiff(colnames(new),c(di$block.name, names(di$factor.names))),drop=FALSE])
               storage.mode(anhaeng) <- "numeric"
               desnum <- cbind(desnum, anhaeng)
               }
     
     rownames(new) <- rownames(desnum) <- ronew$run.no <- 1:nrow(new)
     class(new) <- c("design","data.frame")
     desnum(new) <- desnum
     
     run.order(new) <- ronew
     di$type <- paste(di$type,"center",sep=".")
     ## added in order to support usage of function code.design and steepest ascent analysis
     di$coding <- make.formulas(paste("x", 1:length(di$factor.names), sep=""), di$factor.names)
     design.info(new) <- di
     }
     else {
     ### ??? was wird noch benötigt für ccd.augment???
        new <- design
        ronew <- run.order(design)
        if (!any(is.na(as.numeric(as.character(ronew$run.no.in.std.order)))))
            ronew$run.no.in.std.order <- as.numeric(as.character(ronew$run.no.in.std.order))
        design.info(new) <- c(design.info(new), ncenter=0, ncube=nrow(new))
        run.order(new) <- ronew
        }
     new
}