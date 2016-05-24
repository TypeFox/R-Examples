
direct.sum <- function (D1, D2, ..., tiebreak = letters)
{
    ## taken from package conf.design by Bill Venables
    ## commented out class "design"
    ## added code for keeping contrasts
    l <- list(...)
    if (length(l))
        return(Recall(D1, Recall(D2, ..., tiebreak = tiebreak[-1]),
            tiebreak = tiebreak))
    f1 <- lapply(D1, "is.factor")
    cc <- function(x) {if (is.factor(x)) contrasts(x) else NULL}
    c1 <- lapply(D1, "cc")
    E1 <- lapply(D1, function(x, n2) rep(x, rep(n2, length(x))),
        nrow(D2))
    f2 <- lapply(D2, "is.factor")
    c2 <- lapply(D2, "cc")
    E2 <- lapply(D2, function(x, n1) rep(x, n1), nrow(D1))
    D <- c(E1, E2)
    if (any(i <- duplicated(names(D))))
        names(D)[i] <- paste(names(D)[i], tiebreak[1], sep = "")
    ## reinstate contrasts for factors to previous state
    ff <- function(x, f, c){
         if (f) contrasts(x) <- c
         x}
    D <- mapply("ff", D, f=c(f1,f2), c=c(c1,c2),SIMPLIFY=FALSE)
    D <- as.data.frame(D)
 #   class(D) <- c("design", class(D))
    D
}



cross.design <- function (design1, design2, ..., randomize=TRUE, seed=NULL)
{
    creator <- sys.call()
    workhorse <- function(design1, design2, ...){
    ## taken from package conf.design by Bill Venables
    ## needed as inner function, since postprocessing is not to be racalled
    nam <- deparse(substitute(design2))  ## just in case, for one-column design2
    l <- list(...)
    if (length(l))
        return(Recall(design1, Recall(design2, ...)))
    if (!"design" %in% class(design1)){ 
            warning("design1 is not of class design, all design information is lost")
            cn <- NULL
            if (is.numeric(design2)) cn <- nam
            if (!is.data.frame(design2)) design2 <- data.frame(design2)
            if (is.null(cn)) colnames(design2) <- nam
            aus <- direct.sum(design1, design2, ..., tiebreak = letters)
            if (randomize){ 
               if (!is.null(seed)) set.seed(seed)
               aus <- aus[sample(nrow(aus)),]
               }
            return(aus)
        }
    else{
        di1 <- design.info(design1)
        if (di1$type %in% c("FrF2.blocked", "FrF2.splitplot", "crossed", "FrF2.param", "param") | length(grep("folded",di1$type))>0) 
              stop("crossing blocked designs, splitplot designs, folded designs, crossed designs or parameter designs is not supported")
        if (!is.null(di1$format)) stop("crossing wide format designs is not supported")
        if (!is.null(di1$response.names)) stop("crossing designs with responses is not supported")
        ro1 <- run.order(design1)
        if (di1$repeat.only) stop("only last design can have repeat.only replications")
        des1 <- undesign(design1)
        desn1 <- desnum(design1)
    if ("design" %in% class(design2)){
        di2 <- design.info(design2)
        if (any(di2$type %in% c("FrF2.blocked", "FrF2.splitplot", "crossed", "FrF2.param", "param")) | length(grep("folded",di2$type))>0) 
              stop("crossing blocked designs, splitplot designs, folded designs, crossed designs or parameter designs is not supported")
        if (!is.null(di2$format)) stop("crossing wide format designs is not supported")
        if (!is.null(di2$response.names)) stop("crossing designs with responses is not supported")
        ro2 <- run.order(design2)
        des2 <- undesign(design2)
        desn2 <- desnum(design2)
    }else
        if (is.data.frame(design2) | is.matrix(design2) | is.list(design2) | is.array(design2))
             stop("design2 must be a vector or a data frame of class design")
        else {  ##redo vector into design object
            tab <- table(design2)
            repl <- 1
                if (!min(tab)==max(tab)) type2 <- "vector.unbalanced"
                     else if (max(tab)>1) {type2 <- "vector.replicated"
                                           repl <- max(tab)
                                           }
                          else type2 <- "vector"
            nruns <- length(design2)
            if (type2=="vector.replicated") nruns <- round(nruns/repl)
            nfactors <- 1
            factor.names <- list(names(tab))
            names(factor.names) <- nam 
            ro2 <- data.frame(run.no.in.std.order = match(design2, names(tab)),
                         run.no = 1:length(design2))
            ro2$run.no.std.rp <- ro2$run.no.in.std.order
               if (!type2=="vector") for (i in 1:length(design2))
                    ro2$run.no.std.rp[i] <- paste(ro2$run.no.std.rp[i],
                         cumsum(ro2$run.no.in.std.order==ro2$run.no.std.rp[i])[i],".")
            desn2 <- NULL
            if (is.numeric(design2))
                         desn2 <- matrix(design2,ncol=1, dimnames=list(NULL, nam))
            if (is.character(design2)) design2 <- factor(design2, levels=unique(design2))
            des2 <- design2 <- as.data.frame(design2)
            if (!is.null(nam)) colnames(des2) <- colnames(design2) <- nam
            class(design2) <- c("design","data.frame")
            if (is.null(desn2))
               desnum(design2) <- desn2 <- model.matrix(as.formula(paste("~",nam)),design2)[,-1,drop=FALSE]
               else desnum(design2) <- desn2
               run.order(design2) <- ro2
               design.info(design2) <- di2 <- list(type=type2, nruns=nruns,
                      nfactors=nfactors, factor.names=factor.names,
                      replications=repl, repeat.only=NULL,
                      randomize=NULL, seed=NULL, creator=nam)
            }
    if (any(duplicated(c(colnames(design1),colnames(design2)))))
       stop ("duplicated factor names are not permitted when crossing designs")
    D <- as.data.frame(direct.sum(des1,des2))
    class (D) <- c("design", "data.frame")
    Dn <- as.matrix(direct.sum(as.data.frame(desn1), as.data.frame(desn2)))
    rownames(Dn) <- rownames(D)
    desnum(D) <- Dn
    touter <- function(obj1,obj2,FUN,...) t(outer(obj1,obj2,FUN=FUN,...))
    ## March 7 2016: removed superfluous "t" in the line below 
    ## (thanks to Bill Dunlap for spotting this bug)
    ro <- as.data.frame(mapply("touter",ro1, ro2, "paste", sep="_"))
    run.order(D) <- ro
    
    ## create reasonable content for design.info
    ## accomodate randomize and replications
    ## function for combining design info from two designs
    cc <- function(d1,d2){
         if (is.list(d1) & is.list(d2)) return(list(d1,d2))
         if (is.list(d1) & !is.list(d2)) return(c(d1,list(d2)))
         if (is.list(d2) & !is.list(d1)) return(c(list(d1),d2))
         if (!(is.list(d1) | is.list(d2))) return(c(list(d1),list(d2)))
         }
    cc.alias <- function(d1,d2){
         if (is.list(d1) & is.list(d2)) return(list(d1,d2))
         if (is.list(d1) & !is.list(d2)) return(list(d1,list(d2)))
         if (is.list(d2) & !is.list(d1)) return(list(list(d1),d2))
         if (!(is.list(d1) | is.list(d2))) return(list(list(d1),list(d2)))
         }
    cc.quan <- function(d1,d2){
         if (is.null(d1)) d1 <- rep(NA, di1$nfactors)
         if (is.null(d2)) d2 <- rep(NA, di2$nfactors)
         return(c(d1, d2))
         }
    di <- vector("list") ## empty list
    
    for (nn in union(names(di1), names(di2))){
         if (nn=="quantitative") di[[nn]] <- cc.quan(di1[[nn]],di2[[nn]])
         if (nn=="aliased") di[[nn]] <- cc.alias(di1[[nn]],di2[[nn]])
         if (!nn %in% c("quantitative","aliased"))
            di[[nn]] <- cc(di1[[nn]],di2[[nn]])
      }
    
    ## manually combine infos that is requested to be correct for interim 
    ## processing of recursive procedure (otherwise assignment of design.info throws error)
    di$factor.names <- c(factor.names(design1),factor.names(design2))
    if (!is.list(di$selected.columns)) 
        di$selected.columns <- list(di1$selected.columns ,di2$selected.columns )
    if (is.null(di$cross.nruns)) di$cross.nruns <- unlist(di$nruns)
       else di$cross.nruns <- c(di1$nruns, di$cross.nruns)
    di$nruns <- prod(unlist(di$nruns))
    if (is.null(di$cross.replications)) di$cross.replications <- unlist(di$replications)
       else di$cross.replications <- c(di1$replications, di$cross.replications)
    di$replications <- prod(unlist(di$replications))

    design.info(D) <- di
    if (any(unlist(di$repeat.only)) & di$replications>di2$replications) 
        warning("repeat.only replications and proper replications mixed in one crossed design, this does not work with any post-processing!")
       }
    D
    }
    D <- workhorse(design1, design2, ...)
    
    ## no postprocessing, if design1 was not a design
    if (!"design" %in% class(D)) return(D) 
    
    ## postprocessing for designs
    di <- design.info(D)
       ## has valid lists for many entries
       ## has list with always only 2 elements for nruns and replications
    
    ## modify design info
    
    ## mandatory entries
        ## nruns and replications treated inside the workhorse function already
    di$cross.nruns <- unlist(di$cross.nruns)
    di$cross.replications <- unlist(di$cross.replications)
    di$cross.nfactors <- unlist(di$nfactors)
    di$nfactors <- sum(di$cross.nfactors)
    di$cross.types <- unlist(di$type)
    di$type <- "crossed"
    di$cross.randomize <- unlist(di$randomize)
    di$cross.seed <- di$seed
    di$randomize <- randomize
    di$seed <- seed
    if (is.null(di$seed)) di["seed"] <- list(NULL)
    di$creator <- list(original=di$creator, modify=creator)
    di$cross.repeat.only <- unlist(di$repeat.only)
    di$repeat.only <- any(di$cross.repeat.only)

    ## character string entries that may or may not be filled
    unliststr <- function(obj) if (is.null(obj)) "" else obj
    di$origin <- sapply(di$origin, "unliststr")
    if (all(di$origin=="")) di$origin <- NULL
    di$comment <- sapply(di$comment, "unliststr")
    if (all(di$comment=="")) di$comment <- NULL
    di$generating.oa <- sapply(di$generating.oa, "unliststr")
    if (all(di$generating.oa=="")) di$generating.oa <- NULL
    di$format <- sapply(di$format, "unliststr")
    if (all(di$format=="")) di$format <- NULL
    
    ## integer or logical entries that may or may not be filled
    unlistnumlog <- function(obj, replace=NA) if (is.null(obj)) replace else obj
    if (is.null(unlist(di$clear))) di$clear <- NULL else 
          di$clear <- sapply(di$clear, "unlistnumlog")
    if (is.null(unlist(di$res3))) di$res3 <- NULL else 
          di$res3 <- sapply(di$res3, "unlistnumlog")
    if (is.null(unlist(di$ncube))) di$ncube <- NULL else 
          di$ncube <- sapply(di$ncube, "unlistnumlog")
    if (is.null(unlist(di$ncenter))) di$ncenter <- NULL else 
          di$ncenter <- sapply(di$ncenter, "unlistnumlog")
    if (is.null(unlist(di$residual.df))) di$residual.df <- NULL else 
          di$residual.df <- sapply(di$residual.df, "unlistnumlog")
       

    ## vector or list entries that may or may not be filled
        ## selected.columns treated inside function already    
    if (is.null(unlist(di$selected.columns))) di$selected.columns <- NULL
        di$cross.selected.columns <- di$selected.columns
        di$selected.columns <- NULL
    if (is.null(unlist(di$map))) di$map <- NULL
        di$cross.map <- di$map
        di$map <- NULL
    if (is.null(unlist(di$nlevels))) di$nlevels <- NULL 
       di$cross.nlevels <- di$nlevels
       if (!any(sapply(di$nlevels,"is.null"))) di$nlevels <- unlist(di$nlevels)
          else di$nlevels <- NULL
    
    if (is.null(unlist(di$aliased))) di$aliased <- NULL
    if (is.null(unlist(di$generators))) di$generators <- NULL
    if (is.null(unlist(di$catlg.entry))) di$catlg.entry <- NULL

    
    ## postprocess alias list, if exists
    if (!is.null(di$aliased)){
          dia <- di$aliased
          ## postprocess generators, if exist
              dig <- di$generators   ## may be NULL
          if (di$nfactors<=50){ 
                    nam <- Letters[1:di$nfactors]
                    sepchar <- ""
                    }
                    else{ 
                    nam <- paste("F",1:di$nfactors,sep="")
                    sepchar <- ":"
                    }
          for (i in 1:length(di$aliased)){
             if (!is.null(dia[[i]])){
             if (di$cross.nfactors[i] <=50){ 
                    namalt <- Letters[1:di$cross.nfactors[i]]
                    sepcharalt <- ""
                    }
                    else{ 
                    namalt <- paste("F",1:di$cross.nfactors[i],sep="")
                    sepcharalt <- ":"
                    }
             if (i>1)
             nami <- nam[(1:di$cross.nfactors[i])+sum(di$cross.nfactors[1:(i-1)])]
             else nami <- nam[1:di$cross.nfactors[1]]
             if (!all(namalt==nami)){
             ## replace backward so that overlap between namalt and nam 
             ## is not problematic (namalt>=nam)
             for (j in length(namalt):1){ 
                 ## identical sepchar
                 if (sepcharalt==":" | sepchar==""){
                    ## the following code relies on the fact tha both sepchars are equal 
                    ## (which they should be, since >50 factors cannot occur in individual design if not also in crossed)
                    
                    dia[[i]][["legend"]] <- sub(paste("^",namalt[j],"=",sep=""),paste(nami[j],"=",sep=""),dia[[i]][["legend"]])
                    ## main
                    ## factor before equal or next factor
                    dia[[i]][["main"]] <- gsub(paste(namalt[j],"([=",sepcharalt,"[:alpha:]]{1})",sep=""),
                            paste(nami[j],"\\1",sep=""),dia[[i]][["main"]])
                    ## last factor
                    dia[[i]][["main"]] <- sub(paste(namalt[j],"$",sep=""), nami[j], dia[[i]][["main"]])
                    
                    ## now fi2
                    dia[[i]][["fi2"]] <- gsub(paste(namalt[j],"([=",sepcharalt,"[:alpha:]]{1})",sep=""),
                            paste(nami[j],"\\1",sep=""),dia[[i]][["fi2"]])
                    dia[[i]][["fi2"]] <- gsub(paste(namalt[j],"$",sep=""), nami[j],dia[[i]][["fi2"]])
                    
                    ## finally fi3
                    if (!is.null(dia[[i]][["fi3"]])){
                        dia[[i]][["fi3"]] <- gsub(paste(namalt[j],"([=",sepcharalt,"[:alpha:]]{1})",sep=""),
                            paste(nami[j],"\\1",sep=""),dia[[i]][["fi3"]])
                        dia[[i]][["fi3"]] <- gsub(paste(namalt[j],"$",sep=""), nami[j],dia[[i]][["fi3"]])
                        }
                    ## dig, if not null
                    if (!is.null(dig[[i]])){
                      ## factor before equal or next factor
                      dig[[i]] <- gsub(paste(namalt[j],"([=",sepcharalt,"[:alpha:]]{1})",sep=""),
                              paste(nami[j],"\\1",sep=""),dig[[i]])
                      ## last factor
                      dig[[i]] <- sub(paste(namalt[j],"$",sep=""), nami[j], dig[[i]])
                    }
                    }
                 else{
                    ## sepcharalt "" and sepchar ":"
                    ## also means that alt is an individual letter, while new is not
                    ## only letter that could be problematic: F
                    dia[[i]][["legend"]] <- sub(paste("^",namalt[j],sep=""),paste(nami[j],sep=""),dia[[i]][["legend"]])
                    dia[[i]][["main"]] <- sub(paste("^",namalt[j],sep=""),paste(nami[j],sep=""),dia[[i]][["main"]])
                    dia[[i]][["main"]] <- gsub(paste("=",namalt[j],sep=""),paste("=",nami[j],sepchar,sep=""),dia[[i]][["main"]])
                    dia[[i]][["main"]] <- gsub(paste(namalt[j],"=",sep=""),paste(nami[j],"=",sep=""),dia[[i]][["main"]])
                    dia[[i]][["main"]] <- gsub(paste(namalt[j],"([[:alpha:]{1}])",sep=""),paste(nami[j],sepchar,"\\1",sep=""),dia[[i]][["main"]])
                    dia[[i]][["main"]] <- gsub(paste(namalt[j],"$",sep=""),nami[j],dia[[i]][["main"]])
                    dia[[i]][["fi2"]] <- sub(paste("^",namalt[j],sep=""),paste(nami[j],sepchar,sep=""),dia[[i]][["fi2"]])
                    dia[[i]][["fi2"]] <- sub(paste("=",namalt[j],sep=""),paste("=",nami[j],sepchar,sep=""),dia[[i]][["fi2"]])
                    dia[[i]][["fi2"]] <- sub(paste(namalt[j],"=",sep=""),paste(nami[j],"=",sep=""),dia[[i]][["fi2"]])
                    dia[[i]][["fi2"]] <- gsub(paste(namalt[j],"([[:alpha:]{1}])",sep=""),paste(nami[j],sepchar,"\\1",sep=""),dia[[i]][["fi2"]])
                    dia[[i]][["fi2"]] <- sub(paste(namalt[j],"$",sep=""),nami[j],dia[[i]][["fi2"]])
                    if (!is.null(dia[[i]][["fi3"]])){
                    dia[[i]][["fi3"]] <- sub(paste("^",namalt[j],sep=""),paste(nami[j],sepchar,sep=""),dia[[i]][["fi3"]])
                    dia[[i]][["fi3"]] <- sub(paste("=",namalt[j],sep=""),paste("=",nami[j],sepchar,sep=""),dia[[i]][["fi3"]])
                    dia[[i]][["fi3"]] <- sub(paste(namalt[j],"=",sep=""),paste(nami[j],"=",sep=""),dia[[i]][["fi3"]])
                    dia[[i]][["fi3"]] <- gsub(paste(namalt[j],"([[:alpha:]{1}])",sep=""),paste(nami[j],sepchar,"\\1",sep=""),dia[[i]][["fi3"]])
                    dia[[i]][["fi3"]] <- sub(paste(namalt[j],"$",sep=""),nami[j],dia[[i]][["fi3"]])
                    }
                    if (!is.null(dig[[i]])){
                      ## factor before equal or next factor
                      dig[[i]] <- sub(paste("^",namalt[j],sep=""),paste(nami[j],sep=""),dig[[i]])
                      ## last factor
                      dig[[i]] <- gsub(paste("=",namalt[j],sep=""),paste("=",nami[j],sepchar,sep=""),dig[[i]])
                      dig[[i]] <- gsub(paste(namalt[j],"=",sep=""),paste(nami[j],"=",sep=""),dig[[i]])
                      dig[[i]] <- gsub(paste(namalt[j],"([[:alpha:]{1}])",sep=""),paste(nami[j],sepchar,"\\1",sep=""),dig[[i]])
                      dig[[i]] <- gsub(paste(namalt[j],"$",sep=""),nami[j],dig[[i]])
                    }
                    }
                    
         }
        }}}
         di$aliased <- dia 
         di$generators <- dig
        }
    
    
    design.info(D) <- di
    ## now randomize if requested
    if (randomize){ 
        if (!is.null(seed)) set.seed(seed)
        if (!di$repeat.only) D <- D[sample(nrow(D)),]
        else {repl.repeatonly <- di$cross.replications[length(di$cross.replications)]
               D <- D[rep(repl.repeatonly*(sample(1:round(nrow(D)/repl.repeatonly))-1),each=repl.repeatonly) + 
                    rep(1:repl.repeatonly,round(nrow(D)/repl.repeatonly))]}
        }
    D
}

