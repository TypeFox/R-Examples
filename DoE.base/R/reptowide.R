 reptowide <- function(design, constant=NULL, ...){
    if (!"design" %in% class(design)) stop("design must be of class design")
    di <- design.info(design)
    if (!di$repeat.only){
         warning("no reshaping has been done; reshaping works for designs with repeated measurements only")
         return(invisible(design))
         } 
    if (di$replications==1){
      if (is.null(di$wbreps)){
         warning("no reshaping has been done; reshaping works for designs with repeated measurements only")
         return(invisible(design))
         }
      else if (di$wbreps==1){
         warning("no reshaping has been done; for blocked designs, reshaping only works for repeated measurements within blocks")
         return(invisible(design))
         }
         }
    ## now all cases without repeated measurements have been stopped with a warning
    
    
    hilf <- design
    ro <- run.order(hilf)
    ro$replct <- sapply(lapply(strsplit(as.character(ro$run.no.std.rp),".",fixed=TRUE),"rev"),function(obj) obj[1])
    nrepeat.only <- max(as.numeric(names(table(ro$replct))))
    hilf$run.no.in.std.order <- ro$run.no.in.std.order
    after <- rbind(ro$run.no, rownames(hilf))
    hilf$replct <- ro$replct
    if (is.null(response.names(hilf)))
       response.names(hilf) <- "y"  ## new variable with missings
    fn <- names(factor.names(hilf))
    if (di$type == "FrF2.blocked") fn <- c(di$block.name,fn)
    desnum <- desnum(hilf)
    di <- design.info(hilf)  ## redo in case response names have been updated
    ## extracted all design info now
    hilf <- undesign(hilf) 
    aus <- reshape(hilf, v.names=c(setdiff(colnames(hilf), c("replct","run.no.in.std.order",fn,constant))), 
           idvar=c("run.no.in.std.order",fn),timevar="replct",direction="wide")
    ro <- reshape(ro, v.names=c("run.no","run.no.std.rp"), ids = "run.no", 
             idvar="run.no.in.std.order",timevar="replct",direction="wide")
    ## bring back into class design 
    ## adjust replication info and response names
    ## provide additional design.info for reshaping back to long
        if (!is.null(di$wbreps)) di$wpreps <- round(di$wbreps/nrepeat.only)
           else di$replications <- round(di$replications/nrepeat.only)
        rnlong <- di$response.names
        restlong <- setdiff(colnames(hilf),c(fn, rnlong, constant, "run.no.in.std.order", "replct"))
        di$response.names <- c(t(outer(rnlong, attr(aus,"reshapeWide")$times, paste, sep=".")))
        di$format <- "repeatedMeasuresWide"
        di$responselist <- as.data.frame(t(outer(rnlong, attr(aus,"reshapeWide")$times, paste, sep=".")),
               stringsAsFactors = FALSE)
        colnames(di$responselist) <- rnlong
        if (length(restlong)>0){
           di$restlist <- as.data.frame(t(outer(restlong, attr(aus,"reshapeWide")$times, paste, sep=".")),
               stringsAsFactors = FALSE)
           colnames(di$restlist) <- restlong
        }
        if (!is.null(di$restlist)) restlnames <- c(as.matrix(di$restlist)) else restlnames <- NULL
        cn <- colnames(aus)
        vor <- cn[1:(min(which(cn %in% c(di$response.names,restlnames)))-1)]
        rest <- setdiff(cn, c(di$response.names, restlnames, vor))
        aus <- aus[,c(vor,di$response.names,restlnames, rest)]
    ### ??? switch off warning here for NA from transformation ???
        desnum <- cbind(model.matrix(as.formula(paste("~",paste(fn,collapse="+"))),data=aus)[,-1], 
            as.matrix(as.data.frame(lapply(aus[,setdiff(colnames(aus),c(fn, "run.no.in.std.order"))], "as.numeric"))))
        aus <- aus[,setdiff(colnames(aus),"run.no.in.std.order")]
        ## bring column names of desnum into line with those of design for FrF2 and pb designs
        ##    model.matrix has messed them up for the factors
        if (substr(di$type,1,4)%in%c("FrF2","pb") & !di$type=="FrF2.blocked") colnames(desnum) <- colnames(aus)
        if (di$type=="FrF2.blocked") {colnames(desnum)[-(1:(di$nblocks-1))] <- colnames(aus)[-1]}
        ## remove reshape wide information
        ## leave run.order attribute in the changed form, if not too annoying to distinguish these cases
        attr(aus, "reshapeWide") <- NULL
        attr(ro, "reshapeWide") <- NULL
    ## bring ro back to normal column content (additional columns do not add info)
    ro$run.no <- rank(ro$run.no.1)
    ro$run.no.std.rp <- sapply(lapply(strsplit(as.character(ro$run.no.std.rp.1),".",fixed=TRUE), 
                      function(obj) obj[-length(obj)]),
                     "paste",collapse=".")
    ro <- ro[,c("run.no.in.std.order","run.no","run.no.std.rp")]
    rownames(aus) <- rownames(desnum) <- rownames(ro) <- ro$run.no
    ## attach attributes to design again
    class(aus) <- c("design", class(aus))
    desnum(aus) <- desnum
    run.order(aus) <- ro
    di$repeat.only <- FALSE
    design.info(aus) <- di
    aus
 }
