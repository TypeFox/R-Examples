"[.design" <- function(x,i,j,drop.attr=TRUE,drop=FALSE){
    #   if (missing(i)) i <- 1:nrow(x)
    #   if (missing(j)) j <- 1:ncol(x)
    creator <- sys.call()
    
    if (missing(j)){
       if (identical(sort(i),1:nrow(x)) | identical(i,rep(TRUE, nrow(x)))) {
           ## only rearrange rows
           class(x) <- "data.frame"
           x<-x[i,]
           attr(x,"run.order") <- attr(x,"run.order")[i,]
           attr(x,"desnum") <- attr(x,"desnum")[i,]
           di <- attr(x,"design.info")
           di$creator <- list(original=di$creator, modify=creator)
           attr(x,"design.info") <- di
           class(x) <- c("design","data.frame")
           return(x)}
       else{
           if (is.logical(i)){
              if (!length(i) == nrow(x)) stop("i has wrong length")
              i <- which(i)
              } 
           nnew <- length(i)
           repl <- nnew/nrow(x)
           di <- design.info(x)
           ro <- run.order(x)
           dn <- desnum(x)
           tab <- table(i)
           ## replicate with reshuffle
           if (all(1:nrow(x) %in% as.numeric(names(tab))) & 
                 all(as.numeric(names(tab)) %in% 1:nrow(x)) & min(tab)==max(tab)){
               ## full replication of original experiment
               ## allow proper repeated measurements and proper blocked replication only
               ## otherwise loose attributes
               replicated <- FALSE
               if (di$replications>1) replicated <- TRUE
               if (!is.null(di$wbreps)) if (di$wbreps>1 | di$bbreps>1) replicated <- TRUE
               if (!replicated){
                   ## repeat.only
                   if (all(rep(i[nnew%%repl==1],each=repl)==i)){
                   class(x) <- "data.frame"
                   ro <- ro[rep(1:nrow(x),each=repl),]
                   ro$run.no <- 1:nnew
                   ro$run.no.std.rp <- paste(sapply(strsplit(as.character(ro$run.no.std.rp),".",fixed=TRUE),
                             function(obj) {if (length(obj)==1) obj else paste(obj[-length(obj)],collapse=".")}), 
                                   rep(1:(repl*di$replications),di$nruns),sep=".")
                   x <- x[i,]
                   rownames(x) <- 1:nnew
                   rownames(ro) <- 1:nnew
                   if (!di$type=="FrF2.blocked") di$replications <- repl
                       else di$wbreps <- repl
                       di$repeat.only <- TRUE
                   attr(x,"design.info") <- di
                   attr(x,"run.order") <- ro
                   dn <- dn[i,]
                   rownames(dn) <- 1:nnew
                   attr(x,"desnum") <- dn
                   class(x) <- c("design","data.frame")
                   return(x)
                   }
                   else{
                       ## not repeat.only
                       ##cannot handle blocked designs
                       if (di$type == "FrF2.blocked") warning("design was reduced to data.frame without any attributes")
                       else{
                           proper <- TRUE
                           for (a in 1:repl)
                              if (!all(1:nrow(x) %in% i[(a-1)*repl+(1:nrow(x))])) proper <- FALSE
                           if (proper){
                             ## repl proper replications
                             class(x) <- "data.frame"
                             ro <- ro[i,]
                             ro$run.no <- 1:nnew
                             ro$run.no.std.rp <- paste(sapply(strsplit(as.character(ro$run.no.std.rp),".",fixed=TRUE),
                             function(obj) {if (length(obj)==1) obj else paste(obj[-length(obj)],collapse=".")}), 
                                   rep(1:(repl*di$replications),di$nruns),sep=".")
                             x <- x[i,]
                             rownames(x) <- 1:nnew
                             rownames(ro) <- 1:nnew
                             di$replications <- repl
                             di$repeat.only <- FALSE
                             di$creator <- list(original=di$creator, modify=creator)
                             attr(x,"design.info") <- di
                             attr(x,"run.order") <- ro
                             dn <- dn[i,]
                             rownames(dn) <- 1:nnew
                             attr(x,"desnum") <- dn
                             class(x) <- c("design","data.frame")
                             return(x)
                           }
                       }
                   }}
               ## only rearrange rows with replication
               ## both repeat.only
               if (replicated & di$repeat.only & all(rep(i[nnew%%repl==1],each=repl)==i)){ 
                   ## can handle blocked designs
                   class(x) <- "data.frame"
                   ro <- ro[rep(1:nrow(x),each=repl),]
                   ro$run.no <- 1:nnew
                   ro$run.no.std.rp <- paste(sapply(strsplit(as.character(ro$run.no.std.rp),".",fixed=TRUE),
                             function(obj) {if (length(obj)==1) obj else paste(obj[-length(obj)],collapse=".")}), 
                                   rep(1:(repl*di$replications),di$nruns),sep=".")
                   x <- x[i,]
                   rownames(x) <- 1:nnew
                   rownames(ro) <- 1:nnew
                   di$replications <- di$replications*repl
                   di$creator <- list(original=di$creator, modify=creator)
                   attr(x,"design.info") <- di
                   attr(x,"run.order") <- ro
                   dn <- dn[i,]
                   rownames(dn) <- 1:nnew
                   attr(x,"desnum") <- dn
                   class(x) <- c("design","data.frame")
                   return(x)
                   }
               ## both not repeat.only
               if (replicated & !di$repeat.only){
                   ##cannot handle blocked designs
                      proper <- TRUE
                      for (a in 1:repl)
                         if (!all(1:nrow(x) %in% i[(a-1)*repl+(1:nrow(x))])) proper <- FALSE
                      if (proper){
                        ## repl proper replications
                        class(x) <- "data.frame"
                        ro <- ro[i,]
                        ro$run.no <- 1:nnew
                        ro$run.no.std.rp <- paste(sapply(strsplit(as.character(ro$run.no.std.rp),".",fixed=TRUE),
                             function(obj) {if (length(obj)==1) obj else paste(obj[-length(obj)],collapse=".")}), 
                                   rep(1:(repl*di$replications),di$nruns),sep=".")
                        x <- x[i,]
                        rownames(x) <- 1:nnew
                        rownames(ro) <- 1:nnew
                        di$replications <- repl
                        di$repeat.only <- FALSE
                        di$creator <- list(original=di$creator, modify=creator)
                        attr(x,"design.info") <- di
                        attr(x,"run.order") <- ro
                        dn <- dn[i,]
                        rownames(dn) <- 1:nnew
                        attr(x,"desnum") <- dn
                        class(x) <- c("design","data.frame")
                        return(x)
                     }
                  }
                  }
       ## subset rows
       if (!drop.attr){
         class(x) <- "data.frame"
           aus <- x[i, ,drop=drop]
           class(aus) <- c("design","data.frame")
           attr(aus, "desnum") <- dn[i, ,drop=drop]
           attr(aus, "run.order") <- ro[i,,drop=FALSE]
           attr(aus, "design.info") <- list(type="subset of design", 
               subset.rows=i, nruns=nnew, nfactors=di$nfactors, factor.names=di$factor.names,
               replications=1,repeat.only=di$repeat.only, seed=di$seed, 
               randomize=di$randomize, creator=list(original=di$creator, modify=creator), 
               orig.design.info = di)
       }
       else{
           attr(x, "desnum") <- NULL
           attr(x, "run.order") <- NULL
           attr(x, "design.info") <- NULL
           class(x) <- "data.frame"
           aus <- x[i,,drop=drop]
         }
         }
         ## next brace is end of missing j
       }
    else { class(x) <- "data.frame"
           aus <- x[i,j,drop=drop]}
          
       aus
       }

desnum <- function(design){
     if (!"design" %in% class(design)) stop("desnum is applicable for class design only.")
     else attr(design,"desnum")
 }
`desnum<-` <- function(design, value){
     if (!"design" %in% class(design)) stop("desnum<- is applicable for class design only.")
     if (!is.matrix(value)) stop("value for desnum must be a matrix")
     if (!nrow(value)==nrow(design)) 
         stop("mismatch between numbers of rows for value and design")
     if (!ncol(value)>=ncol(design)) 
         stop("value does not contain enough columns")
     attr(design,"desnum") <- value
     design
 }
 
run.order <- function(design){
     if (!"design" %in% class(design)) stop("run.order is applicable for class design only.")
     else attr(design,"run.order")
 }
`run.order<-` <- function(design, value){
     if (!"design" %in% class(design)) stop("run.order<- is applicable for class design only.")
     if (!is.data.frame(value)) stop("value for run.order must be a data frame")
     if (!nrow(value)==nrow(design)) stop("value and design must have the same number of rows")
     if (!(all(c("run.no.in.std.order","run.no","run.no.std.rp") %in% colnames(value))
          | all(c("run.no.in.std.order","run.no.1","run.no.std.rp.1") %in% colnames(value)))) 
         stop("value does not contain all necessary columns")
         ## covers the long version and the wide reshape with standard settings
     attr(design,"run.order") <- value
     design
 }
 
design.info <- function(design){
     if (!"design" %in% class(design)) stop("design.info is applicable for class design only.")
     else attr(design,"design.info")
 }

`design.info<-` <- function(design, value){
     if (!"design" %in% class(design)) stop("design.info<- is applicable for class design only.")
     if (!is.list(value)) stop("value for design.info must be a list")
     if (!value$nruns*value$replications==nrow(design)){
         if (is.null(value$wbreps)) stop("mismatch between content of value and number of rows in design")
         else if(!value$nruns*value$bbreps*value$wbreps==nrow(design)) 
               stop("mismatch between content of value and number of rows in design")}
     if (!all(c("type","nruns","nfactors","factor.names","replications", "randomize","seed", "repeat.only", "creator") %in% names(value))) 
         stop("value does not contain all necessary elements")
     attr(design,"design.info") <- value
     design
 }

factor.names <- function(design){
     if (!"design" %in% class(design)) stop("design.info is applicable for class design only.")
     else attr(design,"design.info")$factor.names
 }

fnmap <- function(design){
  ## auxiliary function for function factor.names<-
  ## function to map each factor level of the R factor 
  ## to the respective position of the factor level 
  ## in the factor.names element of the 
  ## design.info attribute
  fn <- factor.names(design)
  nlevels <- sapply(fn, length)
  nf <- length(fn)
  facs <- sapply(design, is.factor)
  maps <- mapply(function(obj) 1:obj, nlevels, SIMPLIFY=FALSE)
  for (i in 1:nf) 
      if (facs[i]) maps[[i]] <- sapply(levels(design[[i]]), function(obj) which(fn[[i]]==obj))
  maps
}

`factor.names<-` <- function(design, contr.modify=TRUE, levordold=FALSE, value){
   di <- design.info(design)
   if (!(is.list(value) | is.character(value))) stop("value must be a list or a character vector")
   fnold <- factor.names(design)
   fnnold <- names(fnold)
   fnmap <- fnmap(design)

   if (is.character(value)){ 
      if (!length(value)==length(fnold)) stop("value has the wrong length")
      names(fnold) <- value
      value <- fnold
      ## now value is a list
      }

   if (!length(unique(names(value)))==length(value)) 
      stop("factor names are not unique")
   
   for (i in 1:length(value)) if (identical(value[[i]],"")) value[[i]] <- fnold[[i]][fnmap[[i]]]
   
   if (!length(value)==length(fnold)) stop("value has wrong length")
   ## fac <- sapply(design, "is.factor")
   if (any(sapply(value,function(obj) !length(obj)==length(unique(obj)))))
      stop("duplicate factor levels")
   nlevelsold <- sapply(fnold,"length")
   nlevelsnew <- sapply(value, "length")
   if (!all(nlevelsold==nlevelsnew)) 
        stop("some elements of factor.names do not have the right length")
   
   for (i in 1:length(value)){
        ersetze <- FALSE
        if (nlevelsnew[i]>2 | is.null(di$quantitative[i])) ersetze <- TRUE
        if (nlevelsnew[i]==2 & !is.null(di$quantitative)){
            if (!di$quantitative[i]) ersetze <- TRUE
        }
        if (ersetze){
            if (!is.factor(design[[fnnold[i]]]))
               design[[fnnold[i]]] <- as.factor(design[[fnnold[i]]])
            if (is.factor(design[[fnnold[i]]])) {
                #lev <- as.list(fnold[[i]])
                #names(lev) <- value[[i]] 
                if (levordold)
                  levels(design[[fnnold[i]]]) <- value[[i]]
                else
                  levels(design[[fnnold[i]]]) <- value[[i]][fnmap[[i]]]
                  }
 #           else design[[fnnold[i]]] <- 
 #                factor(design[[fnnold[i]]],levels=fnold[[i]],labels=value[[i]])
           if (contr.modify){
            if (nlevelsnew[i]==2) contrasts(design[[i]]) <- contr.FrF2(2)
            else {
              if (is.numeric(value[[i]])) contrasts(design[[i]]) <- contr.poly(nlevelsnew[i],scores=value[[i]])
                else contrasts(design[[i]]) <- contr.treatment(nlevelsnew[i])}
          }
        }
        else{
        design[[fnnold[i]]] <- (design[[fnnold[i]]] - mean(fnold[[i]]))/(max(fnold[[i]])-min(fnold[[i]])) * 
               (max(value[[i]])-min(value[[i]]))+mean(value[[i]])}
   }
   
   ## must be outside loop, because otherwise problems with names occurring in both versions
   colnames(design)[sapply(fnnold, function(obj) which(colnames(design)==obj))] <- names(value)
   if (!is.null(di$aliased)){
        dial <- strsplit(di$aliased$legend,"=")
        newnames <- names(value)
        dial <- lapply(dial, function(obj){ obj[2] <- newnames[which(fnnold==obj[2])]
                                           obj})
        di$aliased$legend <- sapply(dial, function(obj) paste(obj,collapse="="))
   }
   ## how about blocks (block variable is currently a character variable; why ?
   di$factor.names <- value
   attr(design, "design.info") <- di
   if (!di$type %in% c("lhs","ccd","bbd","bbd.blocked")) 
        desnum(design) <- model.matrix(~.,design)[,-1,drop=FALSE] 
   else {hilf <- desnum(design)
          colnames(hilf)[sapply(fnnold, function(obj) which(colnames(hilf)==obj))] <- names(value)
       desnum(design) <- hilf
       }
   design
}

response.names <- function(design){
     if (!"design" %in% class(design)) stop("design.info is applicable for class design only.")
     else attr(design,"design.info")$response.names
 }

`response.names<-` <- function(design, remove=FALSE, value){
   di <- design.info(design)
   if (!(is.character(value) | is.null(value))) 
       stop("value must be a character vector or NULL")
   if (!length(unique(value))==length(value)) 
      stop("response.names elements are not unique")
   ## dont know why changing response.names was suppressed earlier for wide designs
   ## deactivated that 29 Jan 2011; changed it to impossibility to remove response
   ##      in order to be on the safe side
#   if (!is.null(design.info(design)$responselist))
#      stop("this is a design in wide format, for which function response.names currently does not work")
   if (!is.null(design.info(design)$responselist) && remove){
      warning("this is a design in wide format, for which function response.names does not remove responses")
      remove <- FALSE
      }
   rnold <- di$response.names
   newresp <- setdiff(value, rnold)
   dropresp <- setdiff(rnold, value)

   newrespdrop <- character(0)
   
   if (length(newresp)>0){ 
      for (i in 1:length(newresp)) 
         if (!newresp[i] %in% colnames(design)){
           design[,newresp[i]] <- rep(NA,nrow(design))
           hilf <- desnum(design)
           hilf <- cbind(hilf, as.matrix(design[,newresp[i]]))
           attr(design, "desnum") <- hilf
           }
         else if (!is.numeric(design[,newresp[i],drop=TRUE])){
            newrespdrop <- c(newrespdrop, newresp[i])
         } 
         }
   if (length(newrespdrop>0)) 
       warning("non-numeric response not permitted, responses ", paste(newrespdrop, collapse=","), " not valid")

   di$response.names <- setdiff(value, newrespdrop)
   attr(design,"design.info") <- di
   if (length(dropresp)>0){
      if (remove) {
          hilf <- desnum(design)
          for (i in 1:length(dropresp)){ 
              design[dropresp[i]] <- NULL
              hilf <- hilf[,setdiff(1:ncol(hilf),which(colnames(hilf)==dropresp[i]))]
              }
          attr(design, "desnum") <- hilf
          message("previous responses ", paste(dropresp, collapse=","), " have been removed from the design")
      }
      else
      message("previous responses ", paste(dropresp, collapse=","), " are not considered responses any longer")
   }
   design
}

undesign <- function(design){
   if (!"design" %in% class(design)) stop("design must be of class design")
   ## make design loose its class design and all related attributes
   attr(design,"desnum") <- NULL
   attr(design,"run.order") <- NULL
   attr(design,"design.info") <- NULL
   class(design) <- setdiff(class(design),"design")
   design
}

redesign <- function(design, undesigned){
   if (!"design" %in% class(design)) stop("design must be of class design")
   if (!is.data.frame(undesigned)) stop("undesigned must be a data frame")
   if (!nrow(undesigned) == nrow(design)) stop("design and undesigned must have the same number of rows")
   if (!all(undesigned[,colnames(design)]==design))
       stop("undesign must contain all data columns from design with identical content" )
   class(undesigned) <- c("design",class(undesigned))
   desnum(undesigned) <- desnum(design)
       newcols <- setdiff(colnames(undesigned),colnames(design))
       if (length(newcols)>0){ 
              newdat <- undesigned[,newcols,drop=FALSE]
              for (i in 1:length(newcols))
                  newdat[,newcols[i]] <- as.numeric(newdat[,newcols[i]])
              newdat <- as.matrix(newdat)
              desnum(undesigned) <- cbind(desnum(undesigned), newdat)
              }
   run.order(undesigned) <- run.order(design)
   design.info(undesigned) <- design.info(design)
   undesigned
}

col.remove <- function(design, colnames){
    if (!"design" %in% class(design)) stop("design must be of class design")
    di <- design.info(design)
    if (!is.character(colnames)) stop("colnames must be character")
    if (any(colnames %in% names(di$factor.names))) 
       stop("design factors cannot be removed")
    if (!is.null(di$block.name)) 
      if (di$block.name %in% colnames) 
         stop("the block factor cannot be removed")
    if (length(loeschresp <- intersect(colnames, di$response.names)) > 0){
         loeschrest <- setdiff(colnames, loeschresp)
         if (length(loeschrest)>0){
          hilf <- desnum(design)
          for (i in 1:length(loeschrest)){ 
              design[loeschrest[i]] <- NULL
              hilf <- hilf[,setdiff(1:ncol(hilf),which(colnames(hilf)==loeschrest[i]))]
              }
              attr(design, "desnum") <- hilf
          }
          response.names(design, remove=TRUE) <- setdiff(di$response.names, loeschresp)
    } 
    else {
          hilf <- desnum(design)
          for (i in 1:length(colnames)){ 
              design[colnames[i]] <- NULL
              hilf <- hilf[,setdiff(1:ncol(hilf),which(colnames(hilf)==colnames[i]))]
              }
              attr(design, "desnum") <- hilf
    }
    design
}