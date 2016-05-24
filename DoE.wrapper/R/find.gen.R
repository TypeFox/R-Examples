find.gen <- function(design,...){
    if (!"design" %in% class(design))
        stop("find.gen works on class design objects only")
    di <- design.info(design)
    if ((length(grep("full factorial",di$type))>0 & all(di$nlevels==2)))
         return("full factorial")
    if (!(substr(di$type,1,4)=="FrF2"))
        stop("this is not a regular (fractional) factorial 2-level design")
  #  if (length(grep("blocked", di$type))>0 | length(grep("splitplot", di$type))>0)
  #      stop("find.gen does not work for blocked or splitplot designs")
    if (length(grep("splitplot", di$type))>0)
        stop("find.gen does not work for splitplot designs")
   gen <- NULL

   ## extract and remove old generator information
   ## also in form catlg.entry (or base.design)
   if (!is.null(di$catlg.entry)){
      gen <- di$catlg.entry[[1]]$gen
      if (is.null(di$base.design) & !is.null(di$catlg.entry)) di$base.design <- names(di$catlg.entry)
      di$catlg.entry <- NULL
   }
   else if (!is.null(di$generators))
        gen <- sapply(strsplit(di$generators,"="),
        function(obj) if (substr(obj[2],1,1)=="-") -1*which(names(Yates)==substring(obj[2],2))
                      else which(names(Yates) == obj[2]))
   if (is.null(gen) & !is.null(di$base.design)){ 
   if (substring(di$base.design,1,18)=="generator columns:")
      gen <- eval(parse(text=paste("c(",gsub("generator columns:", "",di$base.design),")")))
   else if (di$base.design %in% names(catlg)) gen <- catlg[di$base.design][[1]]$gen
   }
  ## now gen is a vector of signed columns, if it has been present before or NULL otherwise
  gen
}

generators.from.design <- function(design, ...){
## function to find generators from a design 
## that need not necessarily refer to the standard base factors
    if (!"design" %in% class(design))
        stop("generators.from.design works on class design objects only")
    di <- design.info(design)
    gen <- find.gen(design)
    
    ## relies on the fact that the first k factors are base factors
    ### must take care of map in estimable designs!!!
    ### try with FrF2(32,9, estimable=c("AC","BC","AB")) for reshuffled,
    ###          FrF2(32,9, estimable="CJ") for basic
    if (is.null(gen)){
           if (di$nfactors > 50) stop("generators.from.design does not work for more than 50 factors (and will presumably break down much earlier)")
           if (!is.null(di$ncube)) k <- round(log2(di$ncube)) else k <- round(log2(di$nruns))
           fn <- names(di$factor.names)
           g <- di$nfactors - k
           if (is.null(di$ncube))
              linmod <- lm(formula(paste("I(1:",(di$nruns*di$replications),")~(.)^",k+1,sep="")), design[,names(di$factor.names)])
           else
              linmod <- lm(formula(paste("I(1:",(di$ncube*di$replications),")~(.)^",k+1,sep="")), design[iscube(design),names(di$factor.names)])
           al <- aliases(linmod)
           sel <- al[[2]][1:di$nfactors]   ## aliases of base and generated factors

             sel <- sel[sapply(sel, function(obj) !obj[1] %in% fn[di$map[[1]][1:k]])]
                 ## generated factors
             hilf <- lapply(sel, function(obj){
                    for (i in 1:k){
                        ## complete
                        obj <- gsub(paste("^",fn[di$map[[1]]][i],"$",sep=""), "", obj)
                        ## first
                        obj <- gsub(paste("^",fn[di$map[[1]]][i],":",sep=""), ":", obj)
                        ## middle
                        obj <- gsub(paste(":",fn[di$map[[1]]][i],":",sep=""), ":", obj)
                        ## last
                        obj <- gsub(paste(":",fn[di$map[[1]]][i],"$",sep=""), ":", obj)
                      }
                      obj <- gsub(":","",obj)
                    })
           hilf <- sapply(hilf, function(obj) which(obj==""))  ## the one that contains base factors only
           gen <- sapply(sel, function(obj) paste(Letters[which(fn[di$map[[1]]] == obj[1])],"~",sep="")) ## beginning of generator equation
           if (!length(gen)==g) stop("a problem in determining the generator!")
           for (i in 1:g){ 
                gen[i] <- paste(gen[i],
                                    paste(Letters[which(fn[di$map[[1]][1:k]] %in% unlist(strsplit(sel[[i]][hilf[i]],":",fixed=TRUE)))],collapse="*"),
                                    sep="")
                                    }
    }
    else {
           if (!is.null(di$ncube)) k <- round(log2(di$ncube)) else k <- round(log2(di$nruns))
           gen <- paste(Letters[(k+1):di$nfactors],"~",sapply(Yates[gen], function(obj) paste(Letters[obj],collapse="*")),sep="")
    }
    gen <- lapply(gen, function(obj){ 
        for (i in 1:length(di$factor.names)) obj <- gsub(Letters[i], paste("x",i,sep=""), obj)
        as.formula(obj)
        })
    gen
}