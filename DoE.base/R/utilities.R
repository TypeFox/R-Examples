make.generators <- function(name,liste){
       g <- length(liste)
       k <- length(name)-g
       generators <- rep(list(""),g)
       for (i in 1:g) generators[[i]] <- formula(paste(name[k+i],"~",paste(name[liste[[i]]],collapse="*"),sep=""))
       generators
   }
   
make.formulas <- function(orignames, factor.names){
   ## function for creating coding information
   ## for response surface analysis
   aus <- eval(parse(text=paste("list(",
                paste(paste(orignames,"~(",names(factor.names),"-",sapply(factor.names,"mean"),")/",
                       sapply(factor.names,function(obj) (obj[2]-obj[1])/2), sep=""),collapse=","),
                ")",sep="")))
   names(aus) <- orignames
   aus
}

ord <- function(matrix, decreasing=FALSE){
    ## determines ordering vector that orders matrix 
    ##      w.r.t. first columns, second column etc.
    text <- "order(matrix[,1]"
    if (ncol(matrix)>1)
    for (i in 2:ncol(matrix)) text <- paste(text,",matrix[,",i,"]",sep="")
    text <- paste(text,", decreasing=",decreasing,")")
    eval(parse(text=text))
}


## new function des.recode
des.recode <- function (var, recodes, as.factor.result, char) 
{
## fix level order, changed with version 0.27 
## also simplified the way the function operates
## recodes is a long string with assignments separated by ;  
## intention: leave level orders unchanged for factors var
##            choose level orders according to recodes order for non-factor var
##            whenever all levels are given in recodes
##            otherwise have as.factor determine the level order (like so far)
## I got rid of the rev, but don't know what the rev was for so far!
##   
    is.fac <- is.factor(var)
    if (missing(as.factor.result)) 
        as.factor.result <- is.fac
    if (missing(char)) char <- FALSE

    recode.list <- strsplit(recodes, ";")[[1]]
    recode.list <- strsplit(recode.list, "=")
    oldcodes <- sapply(recode.list, function(obj) obj[1])
    if (char)
    newcodes <- sapply(recode.list, function(obj) obj[2])
    else newcodes <- sapply(recode.list, function(obj) eval(parse(text=obj[2]), 
                      envir = parent.frame(), enclos = sys.frame(0)))

    result <- var
    if (is.fac) 
        result <- as.character(result)
        
    for (i in 1:length(oldcodes)){
            if (is.na(oldcodes[i])) 
                result[is.na(var)] <- newcodes[i]
            else result[var == oldcodes[i]] <- newcodes[i]
        }
    
## fix level order, changed with version 0.27 
    if (as.factor.result) {
        if (is.fac){ 
           levnew <- levels(var)
           for (i in 1:length(oldcodes)){
               if (is.na(oldcodes[i])) 
                  levnew[is.na(levnew)] <- newcodes[i]
               else levnew[levnew == oldcodes[i]] <- newcodes[i]
               }
           result <- factor(result, levels=levnew)
        }
        else{ 
        if (length(newcodes)==length(unique(result)) & length(setdiff(newcodes, result))==0)    
           result <- factor(result, levels=newcodes)
        else result <- as.factor(result)
        }
        }   
    result
}
### needed for old level order in oa.design
des.recode.old <- function (var, recodes, as.factor.result, char) 
{
    recode.list <- rev(strsplit(recodes, ";")[[1]])
    is.fac <- is.factor(var)
    if (missing(as.factor.result)) 
        as.factor.result <- is.fac
    if (missing(char)) char <- FALSE
    result <- var
    if (is.fac) 
        result <- as.character(result)
        
    for (term in recode.list){
        set <- eval(parse(text = strsplit(term, "=")[[1]][1]))
        if (!char)
        target <- eval(parse(text = strsplit(term, "=")[[1]][2]), 
            envir = parent.frame(), enclos = sys.frame(0))
        else 
        target <- strsplit(term, "=")[[1]][2]
        for (val in set){
            if (is.na(val)) 
                result[is.na(var)] <- target
            else result[var == val] <- target
        }
    }
    if (as.factor.result) 
        result <- as.factor(result)
    result
}

## accessor functions for attributes of orthogonal arrays
origin<-function(ID) attr(ID,"origin")
## not needed, as it is available in package base
#comment<-function(ID) attr(ID,"comment")

generators <- function(design, ...){
   UseMethod("generators")
}
generators.design <- function(design, ...){
    ## extract generating contrasts for all FrF2 designs
    ## special care is needed for splitplot, hard and blocked designs
    ## and also estimable
    aus <- NULL
    di <- design.info(design)
    
    if (di$type=="planor") aus <- list(generators=di$generator)
    else{
    
    ### make sure that all functions use the correct catlg (from catlg.name entry)
    catlg.name <- di$catlg.name
    if (is.null(catlg.name)) catlg.name <- "catlg"
    if (!exists(catlg.name)) {
        catlg <- try(eval(parse(text=catlg.name)), silent=TRUE)
        if ("try-error" %in% class(catlg))
        stop("alias information can only be provided, if ", catlg.name, " is available")
        }
    else catlg <- get(catlg.name)   ## within this function, default catlg to the current catalogue
    if (!"catlg" %in% class(catlg)) stop("alias information can not be provided, \nbecause ", catlg.name, " is not a valid design catalogue")
    
    if (length(grep("FrF2", di$type)) == 0) 
          stop("generators are only determined for regular fractional factorial 2-level designs.")
    else{
        k <- round(log2(di$nruns))
        ## prevent execution of function for blocked or splitplot designs generated 
        ## with versions of FrF2 before 1.1
        neuver <- FALSE
        if (!is.null(di$FrF2.version))
            if (compareVersion(di$FrF2.version, "1.1") >= 0) neuver <- TRUE
        if ((length(grep("splitplot",di$type)) > 0 | length(grep("blocked",di$type))>0) & !neuver) 
              stop("generators cannot be extracted from blocked or splitplot designs created with FrF2 versions before 1.1.")
  #      if (length(grep("blocked",di$type))>0 & nfac.catlg(get(di$catlg)[di$base.design]) > di$nfactors) 
  #            stop("generators cannot be extracted from blocked designs created with blockpick.big or with user-specified individual blocking factors")
        if (length(grep("param",di$type)) > 0 | length(grep("folded",di$type)) > 0)
              stop("generators cannot be calculated for folded or parameter designs.")
        if (!is.null(di$catlg.entry)){
             ## catalogue entries for block or split-plot designs are in base.design
             gen <- di$catlg.entry[[1]]$gen
             if(di$nfactors <= 50)
             aus <- list("generators"=paste(Letters[(round(log2(di$catlg.entry[[1]]$nruns),0)+1) : 
                      di$catlg.entry[[1]]$nfac], 
                       unlist(names(Yates)[gen]),sep="="))
             else aus <- list("generators"=paste(paste("F",(round(log2(di$catlg.entry[[1]]$nruns),0)+1) : 
                      di$catlg.entry[[1]]$nfac,sep=""), unlist(names(Yates)[gen]),
                      sep="="))
             #else aus <- paste(paste("F",(round(log2(di$catlg.entry[[1]]$nruns),0)+1) : 
             #         di$catlg.entry[[1]]$nfac,sep=""), sapply(Yates[gen], function(obj2) paste(paste("F", obj2, sep=""),collapse=":")),
             #         sep="=")
          }
        else{ if (!is.null(di$generators))
               aus <- list("generators"=di$generators)
        else{## estimable has map only, named by name of the base design
             ## must be treated separately, because map can refer to all factors
             ##      not only base factors
            if (is.null(di$base.design) & !is.null(di$map)){ 
                 ## determine unmapped generators
                 hilf <- generators(names(di$map), select.catlg=catlg)[[1]]
                 aus <- list("generators"=sort(chartr(paste(Letters[di$map[[1]]],collapse=""),paste(Letters[1:di$nfactors],collapse=""),hilf)))
                 }
            }
        if (!is.null(di$base.design)){ 
               ### can happen for blocked (both versions) and splitplot
               ### can be character string starting with "generator columns:"
               ###   or a design name from catlg (the latter is resolved with catlg loaded only
               
               ### complications arise from blockpick.big (because of additional block generators)
               ###   and from block factors specified individually
               ### and from splitplot because of generated columns potentially moving up to the front
               ### (depending on resolution of whole plot portion of the design)
               hilf.name <- character(0)
               if (length(grep("generator columns:", di$base.design))>0){
                    di$base.design <- as.numeric(unlist(strsplit(unlist(strsplit(di$base.design,c(" "))),","))[-c(1,2)])
                    if (length(di$base.design)==0) return(list(generators="full factorial"))
                    }
               else {
                  if (exists("catlg")){
                    if (di$base.design %in% names(catlg)){
                       hilf.name <- di$base.design
                       di$base.design <- catlg[[di$base.design]]$gen
                    }
                    else stop("For generator information, you need to load the catalogue ", catlg.name, ".")
                  }
               }
               ## Now di$base.design is numeric vector of column numbers
               if (is.null(di$map)) di$map <- 1:k
               if (is.null(di$orig.fac.order)) di$orig.fac.order <- 1:di$nfactors
               mLetters <- Letters[invperm(di$orig.fac.order)]  
                  ## letters that only refer to experimental factors
                  ## appropriately re-arranged for split-plot designs
                   di$base.design <- Yates[di$base.design]
                   di$base.design <- lapply(di$base.design, function(obj) sort(invperm(di$map)[obj]))
               if (!length(di$base.design) + k > di$nfactors){
                   ## no extra block factors apart from experimental factors
                   di$base.design <- paste(mLetters[(k+1):(di$nfactors)],
                                              sapply(di$base.design, function(obj) paste(sort(mLetters[obj]),collapse="")),sep="=")
               }
               k.block.add <- 0
              if (length(di$base.design) + k > di$nfactors){
               #    di$base.design <- paste(Letters[(k+1):(di$nfactors)],
               #                               sapply(di$base.design, function(obj) paste(sort(Letters[obj]),collapse="")),sep="=")
               ### for blockpick.big - generated blocked designs and added factors for blocking:
               ### returned wrong results for blockpick.big and did not work for added factors for blocking 
               ###      before version 0.23-2
                         k.block <- round(log2(di$nblocks)) 
                         k.block.add <- length(di$base.design) + k - di$nfactors
                  ## k.block.add > 0, but potentially < k.block
                       ## base factors as generators or not ???
                       if (!identical(names(di$base.design),names(Yates[di$block.gen]))) 
                          di$base.design <- paste(Letters[(k-k.block.add+1):(di$nfactors)],
                                          sapply(di$base.design, function(obj) 
                                          paste(c(paste("b",1:k.block.add,sep=""),Letters)[obj],collapse="")),sep="=")
                          else di$base.design <- "full factorial"
              }
               if (!is.null(di$block.gen)){ 
                    hilf <- di$block.gen
                    if (k.block.add > 0) {
                        ## special treatment of blocked full factorials
                        if (k == di$nfactors) return(list("generators for design itself"=di$base.design, "block generators"=names(Yates[hilf])))
                        ## other blocked designs
                         hilf <- paste("block generators", paste(paste("b",1:k.block.add,sep=""), collapse=" ")) 
                         if (k.block > k.block.add) hilf <- paste(hilf, names(Yates)[di$block.gen[-(1:k.block.add)]])
                         hilf <- rbind(hilf, 
                                       paste("from Yates matrix columns", paste(di$block.gen, collapse=" ")), 
                                       paste("of base design", hilf.name, "in catalogue", di$catlg))
                         rownames(hilf) <- rep("",3)
                         colnames(hilf) <- ""
                         if (!is.null(di$map)){
                             if (!identical(di$map, 1:k)){
                             hilf <- rbind(hilf, paste("base factors remapped as", paste(di$map, collapse=" ")))
                             rownames(hilf) <- rep("", 4) }
                             }
                         }
                     else{
                      if (is.list(hilf)) 
                       hilf <- names(Yates)[sapply(hilf, function(obj) which(names(Yates)==paste(Letters[sort(obj)],collapse="")))]
                       else hilf <- names(Yates)[hilf]
                      }
                       aus <- c(list(di$base.design), list(hilf))
                      names(aus) <- c(paste("generators for design itself"),
                          "block generators")
                      }
               else {aus <- list(di$base.design)
                      names(aus) <- paste("generators")
               }
              }
        }  
    }
    } ## end FrF2
    aus
}

gen.fun <- function(obj,num=FALSE){
       ## obj must be a single catlg entry
             gen <- obj$gen
             if(obj$nfac <= 50)
             aus <- paste(Letters[(round(log2(obj$nruns),0)+1) : 
                      obj$nfac], 
                       unlist(names(Yates)[gen]),sep="=")
           #  else aus <- paste(paste("F",(round(log2(obj$nruns),0)+1) : 
           #           obj$nfac,sep=""), sapply(Yates[gen], function(obj2) paste(paste("F", obj2, sep=""),collapse=":")),
           #           sep="=")
             else aus <- paste(paste("F",(round(log2(obj$nruns),0)+1) : 
                      obj$nfac,sep=""), unlist(names(Yates)[gen]),
                      sep="=")
             if (num) gen else aus
        }


generators.catlg <- function(design, ...){
    ## design is a list of class catlg
    if (!"catlg" %in% class(design))
       stop("This function works on class catlg objects only.")
    lapply(design, gen.fun)
}

generators.character <- function (design, select.catlg=catlg, ...) 
{
     catlg.name <- deparse(substitute(select.catlg))
     ## select.catlg is used for looking up the design name
    if (!is.character(design)) 
        stop("This function works on character strings only.")
    if (!exists(catlg.name)) {
        catlg <- try(eval(parse(text=catlg.name)), silent=TRUE)
        if ("try-error" %in% class(catlg))
        stop("alias information can only be provided, if ", catlg.name, " is available")
        }
    else catlg <- get(catlg.name)
    if (!"catlg" %in% class(catlg)) stop(catlg.name, " is not a valid design catalogue")
    if (!all(design %in% names(catlg))) 
        stop("character string design contains invalid elements")
    generators(select.catlg[design])
}


invperm <- function (perm) 
{
    sort(perm, index.return = TRUE)$ix
}

PFTs.from.variants <- function(array, variants, R=3, rela=TRUE){
  ## function to calculate a list of (relative) projection frequency tables from 
  ##    an array and a matrix of column numbers 
  
  ## array is an orthogonal array
  ## variants is a matrix the rows of which contain distinct column numbers
  ##    pertaining to array
  ## R is the resolution of the design (3 or 4)
  if (max(variants)>ncol(array)) stop("invalid variants for array")
  if (R==3) PFTs <- lapply(1:nrow(variants), function(obj) P3.3(array[,variants[obj,]], rela=rela))
  else PFTs <- lapply(1:nrow(variants), function(obj) P4.4(array[,variants[obj,]], rela=rela))
  PFTs
}

matrix.fromPFTs <- function(PFTs){
   ## function to bring a list of (R)PFTs with possibly different entries into matrix form
   ## for easy comparison
    zeilen <- sort(unique(unlist(lapply(PFTs, function(obj) obj[,1]))))
    spalten <- 1:length(PFTs)
    pfts <- matrix(0, nrow=length(zeilen), ncol=length(PFTs), dimnames=list(zeilen, spalten))
    for (i in 1:length(PFTs))
       pfts[as.character(PFTs[[i]][,1]),i] <- PFTs[[i]][,2]
    pfts
   }

rankPFT <- function(pfts){
   ## input is an output object of matrix.from.PFTs
   hilf <- t(pfts[nrow(pfts):1,])
   invperm(ord(hilf))
}

bestPFT <- function(pfts){
   ## input is an output object of matrix.from.PFTs
   hilf <- t(pfts[nrow(pfts):1,])
   resort <- ord(hilf)
   hilf <- hilf[resort,,drop=FALSE]
   best <- TRUE
   zeile <- 1
   hb <- hilf[1,]
   nbest <- 1
   while (best & zeile < nrow(hilf)){
        zeile <- zeile + 1
            if (all(hilf[zeile,]==hb)) nbest <- zeile
            else best <- FALSE
                }
   invperm(resort) <= nbest
}