paramtowide <- function(design, constant=NULL, ...){
    if (!"design" %in% class(design)) stop("design must be of class design")
    di <- design.info(design)
    if (!length(grep("param",di$type))>0)
         stop("this function works for parameter designs only (inner / outer array)")
    if (length(grep("wide",di$type))>0)
         stop("this design is already in wide format")
    
    hilf <- design
    ro <- run.order(hilf)
    ro$run.no.outer <- sapply(lapply(strsplit(as.character(ro$run.no),"_",fixed=TRUE),"rev"),function(obj) obj[1])
    ro$run.no <- sapply(strsplit(as.character(ro$run.no),"_",fixed=TRUE),function(obj) obj[1])
    ro$run.no.std.rp <- sapply(strsplit(as.character(ro$run.no.std.rp),"_",fixed=TRUE),function(obj) obj[1])
    ro$run.no.in.std.order <- sapply(strsplit(as.character(ro$run.no.in.std.order),"_",fixed=TRUE),function(obj) obj[1])
    nouter <- length(table(ro$run.no.outer))
    hilf$run.no <- ro$run.no
    hilf$run.no.outer <- ro$run.no.outer
    if (is.null(response.names(hilf)))
       response.names(hilf) <- "y"  ## new variable with missings
    fn <- names(factor.names(hilf))
    desnum <- desnum(hilf)
    di <- design.info(hilf)  ## redo in case response names have been updated
    fninner <- names(di$inner)
    fnouter <- names(di$outer)
    des.outer <- design[,fnouter][1:nouter,]
    ## extracted all design info now
    hilf <- undesign(hilf) 
    constant <- c(constant)
    hilf <- hilf[,setdiff(colnames(hilf),fnouter)]
    aus <- reshape(hilf, v.names=c(setdiff(colnames(hilf), c("run.no.outer","run.no",fn,constant))), 
           idvar=c("run.no",fninner),timevar="run.no.outer", direction="wide")
    ### March 7 2016: corrected wrong "run.no.in.standard.order" to 
    ###     "run.no.in.std.order" (thanks to Bill Dunlap for spotting this)
    ro <- reshape(ro, v.names=c("run.no.in.standard.order"), 
             idvar="run.no",timevar="run.no.outer",direction="wide")
    ## bring back into class design 
    ## adjust replication info and response names
    ## provide additional design.info for reshaping back to long
    ## for the moment decided to disable switching back to long
        rnlong <- di$response.names
        restlong <- setdiff(colnames(hilf),c(fn, rnlong, constant, "run.no", "run.no.outer"))
        di$response.names <- c(t(outer(rnlong, attr(aus,"reshapeWide")$times, paste, sep=".")))
        di$format <- "innerouterWide"
        di$type <- paste(di$cross.types[[1]],".paramwide",sep="")
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
        if (min(which(cn %in% c(di$response.names,restlnames))) < Inf)
             vor <- cn[1:(min(which(cn %in% c(di$response.names,restlnames)))-1)]
             else vor <- character(0)
        rest <- setdiff(cn, c(di$response.names, restlnames, vor))
        aus <- aus[,c(vor,di$response.names,restlnames, rest)]
    ### ??? switch off warning here for NA from transformation ???
        desnum <- cbind(model.matrix(as.formula(paste("~",paste(fninner,collapse="+"))),data=aus)[,-1], 
            as.matrix(as.data.frame(lapply(aus[,setdiff(colnames(aus),c(fninner, "run.no"))], "as.numeric"))))
        aus <- aus[,setdiff(colnames(aus),"run.no")]
        ## bring column names of desnum into line with those of design for FrF2 and pb designs
        ##    model.matrix has messed them up for the factors
        if (substr(di$type,1,4)%in%c("FrF2","pb")) colnames(desnum) <- colnames(aus)
        ## remove reshape wide information
        ## leave run.order attribute in the changed form, if not too annoying to distinguish these cases
        attr(aus, "reshapeWide") <- NULL
        attr(ro, "reshapeWide") <- NULL
    ## bring ro back to normal column content (additional columns do not add info)
    ro$run.no.std.rp <- ro$run.no.in.std.order
    ro <- ro[,c("run.no.in.std.order","run.no","run.no.std.rp")]
    rownames(aus) <- rownames(desnum) <- rownames(ro) <- ro$run.no
    ## attach attributes to design again
    class(aus) <- c("design", class(aus))
    desnum(aus) <- desnum
    run.order(aus) <- ro
    
    ## treat design info

    ## mandatory elements
    di$nruns <- di$cross.nruns[1]
    di$cross.nruns <- NULL
    di$nfactors <- di$cross.nfactors[1]
    di$factor.names <- di$factor.names[1:di$nfactors]
    di$replications <- di$cross.replications[1]
    di$repeat.only <- di$cross.repeat.only[1]
    di$randomize <- di$cross.randomize[1]
    di$seed <- di$cross.seed[1]
    di$outer <- des.outer
    di$inner <- NULL
    di$cross.randomize <- NULL
    di$cross.seed <- NULL
    di$cross.repeat.only <- NULL
    di$cross.replications <- NULL
    ## elements not always present
    di$aliased <- di$aliased[[1]]
    di$generators <- di$generators[[1]]
    di$catlg.entry <- di$catlg.entry[[1]]
    di$map <- di$cross.map[[1]]
    di$cross.map <- NULL
    di$clear <- di$clear[1]
      if (!is.null(di$clear)) if (is.na(di$clear)) di$clear <- NULL
    di$res3 <- di$res3[1]
      if (!is.null(di$res3)) if (is.na(di$res3)) di$res3 <- NULL
    di$nlevels <- di$cross.nlevels[[1]]
    di$cross.nlevels <- NULL
    di$selected.columns <- di$selected.columns[[1]]
    di$selected.columns <- NULL
    di$generating.oa <- di$generating.oa[1]
       if (identical(di$generating.oa,"")) di$generating.oa <- NULL
    di$origin <- di$origin[1]
       if (identical(di$origin,"")) di$origin <- NULL
    di$comment <- di$comment[1]
       if (identical(di$comment,"")) di$comment <- NULL
    di$residual.df <- di$residual.df[1]
        if (!is.null(di$residual.df)) if (is.na(di$residual.df)) di$residual.df <- NULL
    di$quantitative <- di$quantitative[1:di$nfactors]
    di$ncube <- di$ncube[1]
        if (!is.null(di$ncube)) if (is.na(di$ncube)) di$ncube <- NULL
    di$ncenter <- di$ncenter[1]
        if (!is.null(di$ncenter)) if (is.na(di$ncenter)) di$ncenter <- NULL
    design.info(aus) <- di
    aus
 }
