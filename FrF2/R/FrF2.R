FrF2 <- function(nruns=NULL, nfactors=NULL, 
                 factor.names = if(!is.null(nfactors)) {if(nfactors<=50) Letters[1:nfactors] 
                                else paste("F",1:nfactors,sep="")} else NULL, 
                 default.levels = c(-1,1), ncenter=0, center.distribute=NULL,
                 generators=NULL, design=NULL, resolution=NULL, select.catlg=catlg, 
                 estimable=NULL, clear=TRUE, method="VF2", sort="natural", 
                 res3=FALSE, max.time=60, 
                 perm.start=NULL, perms=NULL, MaxC2=FALSE, 
                 replications=1, repeat.only=FALSE, 
                 randomize=TRUE, seed=NULL, alias.info=2, 
                 blocks=1, block.name="Blocks", bbreps=replications, wbreps=1, alias.block.2fis = FALSE,
                 hard=NULL, check.hard=10, WPs=1, nfac.WP=0, WPfacs=NULL, check.WPs=10, ...){
## Version 1.7: added the LAD method and made it the default for clear designs
## method option allows to change it to VF2
creator <- sys.call()
catlg.name <- deparse(substitute(select.catlg))
    nichtda <- "try-error" %in% class(try(eval(parse(text=paste(catlg.name,"[1]",sep=""))),silent=TRUE))
    if (nichtda){
      ## provide valid catalogue names from FrF2.catlg128 (added for version 1.6-6
      
      catlgs128 <- c("catlg128.8to15","catlg128.26to33",paste("catlg128",16:25,sep="."))
      
      if (catlg.name %in% catlgs128){
             if (!require("FrF2.catlg128", quietly=TRUE, character.only=TRUE)) 
                  stop("Package FrF2.catlg128 is not available")
             if (packageVersion("FrF2.catlg128") < numeric_version(1.2)){ 
                  if (catlg.name %in% catlgs128[c(1,3:11)])
                  stop("For this version of package FrF2.catlg128,\n",
                       "load ", catlg.name, " with the command data(", catlg.name,")\n",
                       "and then rerun the FrF2 command.\n",
                       "Alternatively, install the latest version of package FrF2.catlg128.")
                  else stop("You need to get the latest version of package FrF2.catlg128 for using ", catlg.name)
                       }
       }
     else stop(catlg.name, " not available")
     }
    if (!"catlg" %in% class(select.catlg)) stop("invalid choice for select.catlg")
## check validity of center point options
    if (!is.numeric(ncenter)) stop("ncenter must be a number")
    if (!length(ncenter)==1) stop("ncenter must be a number")
    if (!ncenter==floor(ncenter)) stop("ncenter must be an integer number")
    if (is.null(center.distribute)){
      if (!randomize) center.distribute <- min(ncenter, 1)
         else center.distribute <- min(ncenter, 3)}
    if (!is.numeric(center.distribute)) stop("center.distribute must be a number")
    if (!center.distribute==floor(center.distribute)) stop("center.distribute must be an integer number")
    if (center.distribute > ncenter)
       stop("center.distribute can be at most ncenter")
    if (randomize & center.distribute==1) warning("running all center point runs together is usually not a good idea.")
    block.name <- make.names(block.name) ## make block.name a valid R name

## check that no incompatible options are used together
if (ncenter>0 & !identical(WPs,1)) stop("center points for split plot designs are not supported")
if (!(is.null(generators) | (identical(WPs,1) | !is.null(WPfacs)))) stop("generators can only be used with split-plot designs, if WPfacs are specified.")
if (!is.null(nruns)) if (ncenter>0) if (center.distribute > nruns + 1) 
    stop("center.distribute must not be larger than nruns+1")
if (!(is.null(design) | is.null(estimable))) stop("design and estimable must not be specified together.")
if (!(is.null(generators) | is.null(design))) stop("generators and design must not be specified together.")
if (is.null(nruns) & !(is.null(generators))) stop("If generators is specified, nruns must be given.")
if (!(is.null(generators) | is.null(estimable))) stop("generators and estimable must not be specified together.")
if (!(identical(blocks,1) | is.null(estimable))) stop("blocks and estimable must not be specified together.")
    ## 1.7 it might be desirable to alleviate this constraint, and it might also be possible!
if (!(identical(WPs,1) | is.null(estimable))) stop("WPs and estimable must not be specified together.")
if (!(is.null(hard) | is.null(estimable))) stop("hard and estimable must not be specified together.")
if (!(identical(blocks,1) | identical(WPs,1))) stop("blocks and WPs must not be specified together.")
if (!(identical(blocks,1) | is.null(hard))) stop("blocks and hard must not be specified together.")
if (!(identical(WPs,1) | is.null(hard))) stop("WPs and hard must not be specified together.")
if (identical(blocks,1) & !identical(wbreps,1)) stop("wbreps must not differ from 1, if blocks = 1.")
if (!(is.null(WPfacs) | identical(WPs,1)) & is.null(design) & is.null(generators))
          stop("WPfacs requires explicit definition of a design via design or generators.")
if (identical(nfac.WP,0) & is.null(WPfacs) & !identical(WPs,1)) 
          stop("WPs whole plots require specification of whole plot factors through nfac.WP or WPfacs!")
    ## thus, wbreps = 1 can always be assumed for non-split-plot situations,
    ## which reduces necessary case distinctions
#if ((!identical(WPs,1)) & is.null(nruns)) stop("WPs only works if nruns is specified.")
    ### ??? or is that needed for the automatic version only ???
    if (!(is.null(resolution) | is.null(estimable))) stop("You can only specify resolution OR estimable.")
    if (!(is.null(resolution) | is.null(nruns))) warning("resolution is ignored, if nruns is given.")

## simple checks for individual option entries
    if (default.levels[1]==default.levels[2]) stop("Both default levels are identical.")
    if (!(is.logical(clear) & is.logical(res3) & is.logical(MaxC2) & is.logical(repeat.only) 
         & is.logical(randomize) & is.logical(alias.block.2fis) ))
         stop("clear, res3, MaxC2, repeat.only, randomize, and alias.block.2fis must be logicals (TRUE or FALSE).")
    if (!is.numeric(max.time)) 
         stop("max.time must be a positive maximum run time for searching a design with estimable given and clear=FALSE.")
    if (!is.numeric(check.hard)) stop("check.hard must be an integer number.")
    if (!is.numeric(check.WPs)) stop("check.WPs must be an integer number.")
    check.hard <- floor(check.hard)  ## avoid stupid errors
    check.WPs <- floor(check.WPs)  ## avoid stupid errors
    if (!is.numeric(bbreps)) stop("bbreps must be an integer number.")
    if (!is.numeric(wbreps)) stop("wbreps must be an integer number.")
    if (!is.numeric(replications)) stop("replications must be an integer number.")
    if (bbreps > 1  & identical(blocks,1) & !replications > 1) 
        stop("Use replications, not bbreps, for specifying replications for unblocked designs.")

    if (!alias.info %in% c(2,3)) 
         stop("alias.info can be 2 or 3 only.")
    if (!(is.numeric(default.levels) | is.character(default.levels))) 
         stop("default.levels must be a numeric or character vector of length 2")
    if (!length(default.levels) ==2) 
         stop("default.levels must be a numeric or character vector of length 2")
    if (!(is.null(hard) | is.numeric(hard))) 
         stop("hard must be numeric.")
    if (!(is.null(resolution) | is.numeric(resolution))) 
         stop("resolution must be numeric.")
    if (is.numeric(resolution)) if(!(resolution == floor(resolution) & resolution>=3)) 
         stop("resolution must be an integer number (at least 3), if specified.")
    res.WP <- NULL   ## for simple way of checking whether split-plot design later on

## treatment of one specifically-chosen design entry
    if (!is.null(design)){ 
        ## design is looked for in the catalogue specified in select.catlg
        if (!is.character(design)) stop("design must be a character string.")
        if (!length(design)==1) stop("design must be one character string.")
           if (design %in% names(select.catlg)){ 
                cand <- select.catlg[design]
                if (!is.null(nruns)) {if (!nruns==cand[[1]]$nruns) 
                          stop("selected design does not have the desired number of runs.")}
                       else nruns <- cand[[1]]$nruns
                if (!is.null(factor.names)) {if (!length(factor.names)==cand[[1]]$nfac) 
                          stop("selected design does not have the number of factors specified in factor.names.")}
                if (!is.null(nfactors)) {if (!nfactors==cand[[1]]$nfac) 
                          stop("selected design does not have the number of factors specified in nfactors.")}
                       else nfactors <- cand[[1]]$nfac
                }
                else stop("invalid entry for design")
    }
    
    ## should occur after special treatment for specific design
    if (!is.null(nruns)){
       k <- round(log2(nruns))
       if (!2^k==nruns) stop("nruns must be a power of 2.")
       if (nruns < 4 | nruns > 4096) stop("less than 4 or more than 4096 runs are not covered by FrF2.")
       }

    ## check factor specifications
    if (is.null(factor.names) & is.null(nfactors) & (is.null(nruns) | is.null(generators)) & is.null(estimable)) 
         stop("The number of factors must be specified via nfactors, via factor.names, via estimable, through selecting 
         one specific catalogued design or via nruns together with generators.")
    if (!is.null(factor.names) & !(is.character(factor.names) | is.list(factor.names)) ) 
         stop("factor.names must be a character vector or a list.")
    if (is.null(nfactors)) {if (!is.null(factor.names)) nfactors <- length(factor.names)
                          else if (!is.null(generators)) nfactors <- length(generators)+k 
                          }
    if (!is.null(estimable)) {
              if (!is.character(sort)) stop("option sort must be a character string")
              if (!is.character(method)) stop("option method must be a character string")
              if (!sort %in% c("natural","high","low")) stop("invalid choice for option sort")
              if (clear && !method %in% c("LAD","VF2")) stop("invalid choice for option method")
              estimable <- estimable.check(estimable, nfactors, factor.names)
              if (is.null(nfactors)) nfactors <- estimable$nfac
              estimable <- estimable$estimable
              if (is.null(nruns)) {
                 ## make even
                 nruns <- nfactors+ncol(estimable)+1 + (nfactors+ncol(estimable)+1)%%2
                 if (!isTRUE(all.equal(log2(nruns) %% 1,0))) nruns <- 2^(floor(log2(nruns))+1)
                 k <- round(log2(nruns))
                 if (k<3) stop("Please specify nruns and/or nfactors. Calculated values are unreasonable.")
           }
                if (is.null(perm.start)) perm.start <- 1:nfactors
                    else if (!is.numeric(perm.start)) 
                         stop ("perm.start must be NULL or a numeric permutation vector of length nfactors.")
                if (!all(sort(perm.start)==1:nfactors)) 
                         stop ("perm.start must be NULL or a numeric permutation vector of length nfactors.")
                if (!is.null(perms)) {
                   if (!is.matrix(perms) | !is.numeric(perms)) stop("perms must be a numeric matrix.")
                   if (!ncol(perms)==nfactors) stop ("matrix perms must have nfactors columns.")
                   if (any(apply(perms,1,function(obj) any(!sort(obj)==1:nfactors)))) 
                        stop("Each row of perms must be a permutation of 1:nfactors.")
                   }
    }
          ## from here on, estimable is a numeric matrix with two rows
          ## and nfactors is known
    if (!nfactors==floor(nfactors)) 
         stop("nfactors must be an integer number.")
    if (!is.null(factor.names) & !length(factor.names)==nfactors) 
         stop("There must be nfactors factor names, if any.")
    if (is.null(factor.names)) 
         if(nfactors<=50) factor.names <- Letters[1:nfactors] else factor.names <- paste("F",1:nfactors,sep="")
    
    ## check factor level specifications
    if (!((is.character(default.levels) | is.numeric(default.levels)) & length(default.levels)==2) ) 
                 stop("default.levels must be a vector of 2 levels.")
    if (is.list(factor.names)){ 
        if (is.null(names(factor.names))){
            if (nfactors<=50) names(factor.names) <- Letters[1:nfactors]
            else names(factor.names) <- paste("F", 1:nfactors, sep="")
        } 
        if (any(factor.names=="")) factor.names[which(factor.names=="")] <- list(default.levels)}
        else {hilf <- vector("list",nfactors)
              names(hilf) <- factor.names
              hilf[1:nfactors]<-list(default.levels)
              factor.names <- hilf}
    ## from here on, factor.names is a named list
    ## make all names valid R names
      names(factor.names) <- make.names(names(factor.names), unique=TRUE)
    
    if (ncenter > 0) if(any(is.na(sapply(factor.names,"is.numeric"))))
       stop("Center points are implemented for experiments with all factors quantitative only.")

    ### prepare generators case
    if (!is.null(generators)){
               ### nruns and k are always given for this situation
              generators <- gen.check(k, generators)
                 g <- nfactors - k
                 if (!length(generators)== g) 
                 stop("This design in ", nruns, " runs with ", nfactors," factors requires ", g, " generators.")
              
               res <- NA; nclear.2fis<-NA; clear.2fis<-NULL;all.2fis.clear<-NA
               if (g<10) wl <- words.all(k, generators,max.length=6)
               else if (g<15) wl <- words.all(k, generators,max.length=5)
               else if (g<20) wl <- words.all(k, generators,max.length=4)
               else if (g>=20) wl <- alias3fi(k, generators, order=2)  ## not a word list
               WLP <- NULL
               if (g < 20){
               WLP <- wl$WLP
                 res <- min(as.numeric(names(WLP)[which(WLP>0)]))
                 if (res==Inf) {if (g<10) res="7+"
                   else if (g<15) res="6+"
                   else if (g<20) res="5+" }
                 }
               else{
                   if (!is.list(wl)) res="5+"
                   else{
                   if (length(wl$"main")>0) res="3"
                      else if (length(wl$"fi2")>0) res="4"
                           else res="5+"}
               }

               gen <- sapply(generators,function(obj) which(sapply(Yates[1:(nruns-1)],
                             function(obj2) isTRUE(all.equal(sort(abs(obj)),obj2)))))
               gen <- gen*sapply(generators, function(obj) sign(obj[1]))

               cand <- list(custom=list(res=res, nfac=nfactors, nruns=nruns, 
                    gen=gen, 
                    WLP=WLP, nclear.2fis=nclear.2fis, clear.2fis=clear.2fis, all.2fis.clear=all.2fis.clear))
                   ## needs to be list of list, because access later is always via cand[1]
               class(cand) <- c("catlg","list")
      }
    ## prepare blocks
    ## 1.8: if blocks and estimable go together, does this preparation need to accomodate estimable ?
    ## 1.8: so far, estimable is a two-row matrix, and nfactor and nruns are known (lines 135 to 170)
    ## 1.8: provisionally, it looks as though no adjustments are needed
    ## 1.8: only perhaps checks against using block generators with estimable (???)
    if (!identical(blocks,1)) {
             blocks <- block.check(k, blocks, nfactors, factor.names)  
             if (is.list(blocks)) k.block <- length(blocks)
             block.auto=FALSE   ## default for block.auto, set to TRUE below, if TRUE
             map <- NULL   ## initialize, will be filled for blockpick.big
             }
             ## blocks is now list of generators or number of blocks
    if (!is.list(blocks)){
    if (blocks>1){
       block.auto=TRUE
       if (is.null(nruns)) 
            stop("blocks>1 only works if nruns is specified.")
       k.block <- round(log2(blocks))
       if (blocks > nruns/2) 
            stop("There cannot be more blocks than half the run size.")
       if (nfactors+blocks-1>=nruns) 
            stop(paste(nfactors, "factors cannot be accomodated in", nruns, "runs with", blocks, "blocks."))
       ntreat <- nfactors
#       nfactors <- nfactors + k.block     #### omitted for new version, that alone is probably not yet the solution
    }
    }
    ### treat hard before WPs, so that it can be handled by split-plot approach
    if (!is.null(hard)){
      if (is.null(generators)){ 
           ## otherwise, cand has already been specified
           if (is.null(nruns)){
              ## resolution requested, k not yet defined
                 cand <- select.catlg[which(res.catlg(select.catlg)>=resolution & nfac.catlg(select.catlg)==nfactors)]
                 if (length(cand)==0) {
                    message("full factorial design needed for achieving requested resolution")
                    k <- nfactors
                    nruns <- 2^k
                    cand <- list(list(gen=numeric(0)))
                    }
                 else {
                     nruns <- min(nruns.catlg(cand))
                     k <- round(log2(nruns))
                     cand <- cand[which(nruns.catlg(cand)==nruns)]
                    }
           }
           else {
               ## nruns specified
               if (nfactors > k)
               cand <- select.catlg[which(nfac.catlg(select.catlg)==nfactors & nruns.catlg(select.catlg)==nruns)]
               else cand <- list(list(gen=numeric(0)))
           }
      }
      if (hard == nfactors) stop("It does not make sense to choose hard equal to nfactors.")
      if (hard >= nruns/2) 
          warning ("Do you really need to declare so many factors as hard-to-change ?")
      nfac.WP <- hard
         ## determine WPs (implying k.WP)
         if (hard < nruns/2){
         WPs <- NA
         ## for full factorials, WPs <- 2^hard, 
         ## otherwise achievable WPs is determined using leftadjust 
         if (length(cand[[1]]$gen) > 0)
         for (i in 1:min(length(cand),check.hard)){
             leftadjust.out <- leftadjust(k,cand[[i]]$gen, early=hard, show=1)
             if (is.na(WPs) | WPs > 2^leftadjust.out$k.early) 
                 WPs <- 2^leftadjust.out$k.early
             }
             else WPs <- 2^hard
         }
      if (hard>=nruns/2 | WPs==nruns) {
           warning("There are so many hard-to-change factors that no good special design could be found.
                Randomization has been switched off.")
           randomize <- FALSE
           WPs <- 1
           ## desmat will be the special slow-change matrix
           leftadjust.out <- leftadjust(k,cand[[1]]$gen,early=hard,show=1)
           generators <- leftadjust.out$gen
           ## with this generator and the hard-to-change factors in the earliest columns
      }
    }
    ## prepare WPs
    if (!identical(WPs,1)) {
             if (is.null(nruns)) stop("WPs>1 only works if nruns is specified.")
             if (WPs > nruns/2) stop("There cannot be more whole plots (WPs) than half the run size.")
             k.WP <- round(log2(WPs))
             if (!WPs == 2^k.WP) stop("WPs must be a power of 2.")
             if (!is.null(WPfacs) & nfac.WP==0){
                    nfac.WP <- length(WPfacs)
                    if (nfac.WP < k.WP) stop("WPfacs must specify at least log2(WPs) whole plot factors.")
                 }
             if (nfac.WP==0) stop("If WPs > 1, a positive nfac.WP or WPfacs must also be given.")
             if (nfac.WP < k.WP) {
                 add <- k.WP - nfac.WP
                 names.add <- rep(list(default.levels),add)
                 names(names.add) <- paste("WP",(nfac.WP+1):(nfac.WP+add),sep="")
                 nfactors <- nfactors + add
                 factor.names <- c(factor.names[1:nfac.WP],names.add,factor.names[-(1:nfac.WP)])
                 nfac.WP <- k.WP
                 warning("There are fewer factors than needed for a full factorial whole plot design. ", add, 
                   " dummy splitting factor(s) have been introduced.")
              }
             if (!is.null(WPfacs)) {
                WPfacs <- WP.check(k, WPfacs, nfac.WP, nfactors, factor.names)  
                ## WPfacs is now character vector of the form Fno or single number of whole plots
                WPsnum <- as.numeric(chartr("F", " ", WPfacs))
                ## the factor positions within the factor names list from which the design is generated
                WPsorig <- WPsnum 
             }
             else {WPsorig <- WPsnum <- 1:nfac.WP }  ## first nfac.WP factors are whole plot factors
             
             ## here, the original WPsorig has been documented
             ## this must not be overwritten later!
    }
    
    if (!is.null(nruns)){
       ## find FrF2 for given number of runs 
       ## or one specific design (because nruns has already been calculated if necessary)
       ## or estimable given (because nruns has already been calculated if necessary)
       if (nfactors<=k & identical(blocks,1) & identical(WPs,1)) {
                    ### full factorial without blocks or WPs finally treated here
                    ### result immediately returned, i.e. function ends here in this case
                    if (nfactors==k) aus <- fac.design(2, k, factor.names=factor.names, 
                             replications=replications, repeat.only=repeat.only, 
                             randomize=randomize, seed=seed)
                    else aus <- fac.design(2, nfactors, factor.names=factor.names, 
                             replications=replications*2^(k-nfactors), repeat.only=repeat.only, 
                             randomize=randomize, seed=seed)
                    aus <- qua.design(aus, quantitative="none", contrasts=rep("contr.FrF2",nfactors))
                    ## add center points, if requested
                    if (ncenter>0) aus <- add.center(aus, ncenter, distribute=center.distribute)
                    di <- design.info(aus)
                    di$creator <- creator
                    di <- c(di, list(FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version))
                    design.info(aus) <- di
                    return(aus)
              }
       else {
       if (nfactors < k) stop("A full factorial for nfactors factors requires fewer than nruns runs. 
           Please reduce the number of runs and introduce replications instead.")
       if (nfactors == k) {
           generators <- as.list(numeric(0))
               cand <- list(custom=list(res=Inf, nfac=nfactors, nruns=nruns, 
                    gen=numeric(0), 
                    WLP=c(0,0,0,0), nclear.2fis=choose(k,2), clear.2fis=combn(k,2), all.2fis.clear="all"))
                   ## needs to be list of list, because access later is always via cand[1]
               class(cand) <- c("catlg","list")
           }
             ### fractional factorial
         if (nfactors > nruns - 1) 
               stop("You can accomodate at most ",nruns-1," factors in a FrF2 design with ",nruns," runs." )
         g <- nfactors - k  ## number of generators needed
            if (!is.null(estimable)) {
                ## determine design that satisfies estimability requests
                      desmat <- estimable(estimable, nfactors, nruns, 
                           clear=clear, res3=res3, max.time=max.time, select.catlg=select.catlg, 
                           method=method,sort=sort,
                           perm.start=perm.start, perms=perms, order=alias.info)
                      design.info <- list(type="FrF2.estimable", 
                             nruns=nruns, nfactors=nfactors, factor.names=factor.names, catlg.name = catlg.name,
                             map=desmat$map, aliased=desmat$aliased, clear=clear, res3=res3, 
                             FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version)
                      desmat <- desmat$design
                      desmat <- as.matrix(sapply(desmat,function(obj) as.numeric(as.character(obj))))
                      rownames(desmat) <- 1:nrow(desmat)
                  }
            else if (is.null(generators) & is.null(design)) 
                 cand <- select.catlg[nruns.catlg(select.catlg)==nruns & nfac.catlg(select.catlg)==nfactors]
               ## candidate designs for further steps
               ## resolution not here, because ignored for given nruns

            ## treating blocked designs
            ## 1.7 here is where the main changes are needed for allowing to block 
            ## 1.7    designs with estimable 2fis;
            ## 1.7    probably restrict choose(nruns - 1 - nfactors - length(estimable), k.block)
            ## 1.7    also generally probably choose(nruns - 1 - nfactors - ncol2fis, k.block)
            ## 1.7     if alias.block.2fis = FALSE; but do I get ncol2fis?
            block.gen <- NULL    ## for later checks
            if (!is.list(blocks)){
                if (blocks > 1) {
                    ### small case
                    if (g==0 | choose(nruns - 1 - nfactors, k.block) < 100000){
                    for (i in 1:length(cand)){
                      ### loop through possible generator designs from best to worse overall
                      ### break stops the loop as soon as design has been found
                      if (g==0) {blockpick.out <- try(blockpick(k, gen=0, 
                            k.block=k.block, show=1, alias.block.2fis = alias.block.2fis),TRUE)
                            }
                      else {
                      if (is.null(generators))
                      blockpick.out <- try(blockpick(k, design=names(cand[i]), 
                            k.block=k.block, show=1, alias.block.2fis = alias.block.2fis),TRUE)
                      else  blockpick.out <- try(blockpick(k, gen=cand[[i]]$gen, 
                            k.block=k.block, show=1, alias.block.2fis = alias.block.2fis),TRUE)
                      }

                      if (!"try-error" %in% class(blockpick.out)) {
                         blocks <- blockpick.out$blockcols   ## column numbers in Yates matrix
                         block.gen <- blocks ## for design.info
                         cand <- cand[i]
                         cand[[1]]$gen <- c(cand[[1]]$gen,blocks)
                         blocks <- nfactors+(1:k.block) ## now in terms of factors
                      #### nfactors changed, will be reduced again later
                         nfactors <- nfactors+k.block
                         g <- g+k.block
                      ### adjust factor.names
                         hilf <- factor.names
                         factor.names <- vector("list", nfactors)
                         factor.names[-blocks] <- hilf
                         factor.names[blocks] <- list(default.levels)
                         names(factor.names) <- c(names(hilf),paste("b",1:k.block,sep=""))
                         blocks <- as.list(blocks)   ### rest is treated with the manual routine
                         break
                      }
                    }## else treats big case
                    }
                else{
                    nfactors <- nfactors+k.block
                    g <- g+k.block
                    hilf <- factor.names
                    factor.names <- vector("list", nfactors)
                    factor.names[(k.block+1):nfactors] <- hilf
                    factor.names[1:k.block] <- list(default.levels)
                    names(factor.names) <- c(paste("b",1:k.block,sep=""),paste(names(hilf)))
                    cand <- select.catlg[nfac.catlg(select.catlg)==nfactors & nruns.catlg(select.catlg)==nruns]
                    if (!(is.null(generators) & is.null(design))) 
                         warning("For this request, generator or design specifications have been ignored, 
                              because the block allocation procedure for big problems was used.")
                    for (i in 1:length(cand)){
                      blockpick.out <- try(blockpick.big(k, gen=cand[[i]]$gen, 
                            k.block=k.block, show=1, alias.block.2fis = alias.block.2fis),TRUE)
                      if (!"try-error" %in% class(blockpick.out)) {
                         cand <- cand[i]
                         cand[[1]]$gen <- blockpick.out$gen[1,]  ## includes block generators
                         map <- blockpick.out$perms[1,]
                         blocks <- as.list(1:k.block) ## first k.block generator columns
                                                      ### rest is treated with the manual routine
                         block.gen <- 2^(0:(k.block-1))
                         break
                      }
                    }
                    } ## end else (for big case)
                if (alias.block.2fis & !is.list(blocks)) 
                         stop("no adequate block design found")
                if ((!alias.block.2fis) & !is.list(blocks)) 
                         stop("no adequate block design found with 2fis unconfounded with blocks")
                }
                }
            if (is.list(blocks)) {
                # can be a pre-treated automatic situation or a manual specification
                hilf.gen <- c(2^(0:(k-1)), cand[[1]]$gen)    ### cand[[1]] exists ?
                hilf.block.gen <- sapply(blocks, function(obj) 
                       as.intBase(paste(rowSums(do.call(cbind,lapply(obj, function(obj2) digitsBase(hilf.gen[obj2],2,k))))%%2,collapse="")))
                k.block.add <- length(intersect(hilf.block.gen, hilf.gen))
                if (is.null(block.gen)) {
                   ## manual allocation 
                   ntreat <- nfactors - k.block.add
                   block.gen <- hilf.block.gen
                   }
                if (k.block > 1){
                ### check that manual choices for block entries contain independent entries only
                hilf <- hilf.block.gen
                for (i in 2:k.block){
                       sel <- combn(k.block,i)
                       for (j in 1:ncol(sel)){
                          neu <- as.intBase(paste(rowSums(do.call(cbind,lapply(sel[,j], function(obj) digitsBase(hilf.block.gen[obj],2,k))))%%2,collapse=""))
                          if (neu %in% hilf) stop("specified blocks is invalid (dependencies)")
                          else hilf <- c(hilf, neu)
                          }
                   }
                rm(hilf)
                }
            }
            ## automatic treatment of split-plot designs
            if (WPs > 1){ 
               WP.auto <- FALSE
               map <- 1:k
               orignew <- WPsorig  ## orignew will be changed according to map later
                                   ## for cases without WPfacs
                                   ## this is for orig.fac.order element of design.info, 
                                   ## not for re-arranging columns in design
               if (is.null(WPfacs)){
                    WP.auto <- TRUE
                    ## max.res.5 is the max number of factors for resolution 5, 
                    ##    if k (or k.WP is of size ...
                    max.res.5 <- c(1,2,3, 5, 6, 8, 11, 17)
                    for (i in 1:length(cand)){
                      ## prevent wasting search time on impossible setups
                      if (is.null(generators)){
                        if (cand[[i]]$res>=5 & nfac.WP > max.res.5[k.WP]) next
                        if (cand[[i]]$res>=4 & nfac.WP > WPs/2) next
                      }
                      if (nfac.WP > WPs/2 | nfac.WP <= k.WP)
                      splitpick.out <- try(splitpick(k, cand[[i]]$gen, k.WP=k.WP, nfac.WP=nfac.WP, show=1),TRUE)
                      else  splitpick.out <- try(splitpick(k, cand[[i]]$gen, k.WP=k.WP, nfac.WP=nfac.WP, 
                             show=check.WPs),TRUE)
                      if (!"try-error" %in% class(splitpick.out)) {
                         WPfacs <- 1:k.WP    
                         if (nfac.WP > k.WP) WPfacs <- c(WPfacs, (k + 1):(k+nfac.WP-k.WP))
                         cand <- cand[i]
                         cand[[1]]$gen <- splitpick.out$gen[1,]
                         res.WP <- splitpick.out$res.WP[1]
                         map <- splitpick.out$perms[1,]
                         break
                      }
                    }
                if (is.null(res.WP)){
 #                     if (nruns >= 2^nfactors) stop("currently, splitplot designs for full factorials are not covered.")
                      if (nruns >= 2^nfactors) {
                          res.WP <- Inf
                          WP.auto <- TRUE
                          WPfacs <- 1:k.WP
                          if (!k.WP == nfac.WP) stop(nfac.WP, " whole plot factors cannot be accomodated in ", (2^k.WP), 
                              " whole plots for a full factorial. Please request smaller design with replication instead.")
                          cand <- list(list(gen=numeric(0)))  ## in order to be treatable below
                      }
                      else stop("no adequate splitplot design found")
                }
  #             WPsorig <- WPsnum <- WPfacs    ## this was a mistake, leading to other than the first nfac.WP factors being WP factors
                orignew <- WPsnum <- WPfacs
                WPfacs <- paste("F",WPfacs,sep="")  ## treat with manual below
                }
              }
              ## next closing brace for !(nfactors<=k & identical(blocks,1) & identical(WPs,1))
              }
              ### next closing brace for !is.null(nruns)
              } 
else {
       ## nruns not given and not already assigned for estimable
       ## find smallest FrF2 that fulfills the requirements regarding resolution/estimability
       if (is.null(resolution) & is.null(estimable)) 
                 stop("At least one of nruns or resolution or estimable must be given.")
       if (!is.null(resolution)){ 
                 cand <- select.catlg[which(res.catlg(select.catlg)>=resolution & nfac.catlg(select.catlg)==nfactors)]
                 if (length(cand)==0) {
                    message("full factorial design needed")
                    aus <- fac.design(2, nfactors, factor.names=factor.names, 
                             replications=replications, repeat.only=repeat.only, 
                             randomize=randomize, seed=seed)
                    for (i in 1:nfactors) 
                           if (is.factor(aus[[i]])) contrasts(aus[[i]]) <- contr.FrF2(2)
                    ## add center points, if requested
                    if (ncenter>0) aus <- add.center(aus, ncenter, distribute=center.distribute)
                    return(aus)
                    }
            }
       else{ cand <- select.catlg[which(nfac.catlg(select.catlg)==nfactors)]   
           ## no prior restriction by resolution
         if (length(cand)==0) 
             stop("No design listed in catalogue ", deparse(substitute(select.catlg)), " fulfills all requirements.")
           }
    }
    ## standard FrF2 situations
    ## maximize clear 2fis among maximum resolution designs (res3=FALSE)
    ##     or among all designs (res3=TRUE)
    ##   changed 28 01 2011 from always maximum resolution
    if (MaxC2 & is.null(estimable) & is.null(generators) ) {
           if (!res3)
           cand <- cand[which.max(sapply(cand[which(sapply(cand, 
                   function(obj) obj$res)==max(sapply(cand, function(obj) obj$res)))], 
                   function(obj) obj$nclear.2fis))]
           else
           cand <- cand[which.max(sapply(cand, 
                   function(obj) obj$nclear.2fis))]
            }
    ### creation of actual design 
    if (is.null(nruns)) {nruns <- cand[[1]]$nruns 
                         k <- round(log2(nruns))
                         g <- nfactors - k}
    if (is.null(estimable)){
       destxt <- "expand.grid(c(-1,1)"
       for (i in 2:k) destxt <- paste(destxt,",c(-1,1)",sep="")
       destxt <- paste("as.matrix(",destxt,"))",sep="")
       desmat <- eval(parse(text=destxt))
       if (is.character(WPfacs) | is.list(blocks)) {
          desmat <- desmat[,k:1]  ## slow first rather than fast first
          ##rownames(desmat) <- slowfast(k)  ## 17 Feb 2013; this would make the 
          ##   run numbers refer to the Yates matrix row order with A changing fastest
          ##   but would also change row order of the resulting design
          ##   therefore not done; 
          }
       if (!is.null(hard)) {
             desmat <- rep(c(-1,1),each=nruns/2)
             for (i in 2:k) desmat <- cbind(desmat,rep(c(1,-1,-1,1),times=(2^i)/4,each=nruns/(2^i)))
             }
       if (g>0)
       for (i in 1:g) 
       desmat <- cbind(desmat, sign(cand[[1]]$gen[i][1])*apply(desmat[,unlist(Yates[abs(cand[[1]]$gen[i])])],1,prod))
       if (WPs > 1) {
          if (!WP.auto) {
            hilf <- apply(desmat[,WPsorig,drop=FALSE],1,paste,collapse="")
            if (!length(table(hilf))==WPs) 
            stop("The specified design creates ",
                 length(table(hilf)), " and not ", WPs, " whole plots.")
                 for (j in setdiff(1:nfactors,WPsorig))
                     if (!length(table(paste(hilf,desmat[,j],sep="")))>WPs) 
                         stop("Factor ", names(factor.names)[j], " is also a whole plot factor.")
            ## added 19 Feb 2013
            if (nfac.WP<3) res.WP <- Inf
            else res.WP <- GR((3-desmat[,WPsorig,drop=FALSE])%/%2)$GR%/%1           
          }
          }
       
       ## slow changing order, if required
      # if ((!is.character(WPfacs)) & !is.null(hard) ) 
      #     desmat <- desmat[,order(c(2^(0:(k-1)),abs(cand[[1]]$gen)))]

      ## make sure matrix desmat has rownames (14 Feb 2013)
      if (is.null(rownames(desmat))) rownames(desmat) <- 1:nruns

       if (is.list(blocks)) {
           ## manually blocked designs and continuation of automatically blocked designs
            if (is.null(block.gen)) block.gen <- blocks
            hilf <- blocks
            for (i in 1:k.block) hilf[[i]] <- apply(desmat[,hilf[[i]],drop=FALSE],1,prod)
            Blocks <- factor(as.numeric(factor(apply(matrix(as.character(unlist(hilf)),ncol=k.block),
                        1,paste,collapse=""))))
            hilf <- order(Blocks, as.numeric(rownames(desmat)))   ## rownames(desmat) added 17 Feb 2013
            desmat <- desmat[hilf,]
            Blocks <- Blocks[hilf]
            contrasts(Blocks) <- contr.FrF2(levels(Blocks))
            nblocks <- 2^length(blocks)
            blocksize <- nruns / nblocks
            block.no <- paste(Blocks,rep(1:blocksize,nblocks),sep=".")
            }
       if (is.character(WPfacs)) {
            ## make WP factor columns the first ones
            ## make factor.names match this order 
            #### for automatic selection the first WPs factors are WP factors !!!
            #if (nfac.WP > k.WP){
            
            ## are WPsnum is the current factor numbering 
            ## WPsorig is the original numbering before column permutations 
               desmat <- desmat[,c(WPsnum,setdiff(1:nfactors,WPsnum))]
               factor.names <- factor.names[c(WPsorig,setdiff(1:nfactors,WPsorig))]

               plotsize <- round(nruns/WPs)
               
               if (is.null(hard)){
                 hilf <- ord(cbind(desmat[,1:nfac.WP],as.numeric(rownames(desmat))))  ## as.numeric(rownames) added 17 Feb 2013
                 desmat <- desmat[hilf,]
               }
               wp.no <- paste(rep(1:WPs,each=plotsize),rep(1:plotsize,WPs),sep=".")
               
            }
    }

    ## until Feb 14 2013, rownames(desmat) were defined here, 
    ##    which did not make sense, especially since they should be available 
    ##    in already resorted form above for block and wp
    ## 14 Feb 2013: auskommentiert
    ## rownames(desmat) <- 1:nruns

    ## handle randomization and replication
    if (randomize & !is.null(seed)) set.seed(seed)
    if (!(is.list(blocks) | WPs > 1)){
      ## simple randomization situations
      rand.ord <- rep(1:nruns,replications)
      if (replications > 1 & repeat.only) rand.ord <- rep(1:nruns,each=replications)
      if (randomize & !repeat.only) for (i in 1:replications) 
                  rand.ord[((i-1)*nruns+1):(i*nruns)] <- sample(nruns)
      if (randomize & repeat.only) rand.ord <- rep(sample(1:nruns), each=replications)
      }
      else {   
          if (is.list(blocks)){
          #### implement randomization and replication for blocks
          ## bbreps=replications, wbreps=1
          ## if ((!repeat.only) & !randomize)
              rand.ord <- rep(1:nruns, bbreps * wbreps)
              if ((!repeat.only) & !randomize)
                 for (i in 0:(nblocks-1))
                    for (j in 1:wbreps)
                    rand.ord[(i*blocksize*wbreps+(j-1)*blocksize+1):((i+1)*blocksize*wbreps+j*blocksize)] <- 
                         (i*blocksize+1):((i+1)*blocksize)
                    rand.ord <- rep(rand.ord[1:(nruns*wbreps)],bbreps)
              if (repeat.only & !randomize)
                    for (j in 1:wbreps)
                    rand.ord[(i*blocksize*wbreps + (j-1)*blocksize + 1) : 
                          (i*blocksize*wbreps + j*blocksize)] <- 
                               sample((i%%nblocks*blocksize+1):(i%%nblocks*blocksize+blocksize))
              if (wbreps > 1 & repeat.only) rand.ord <- rep(1:nruns,bbreps, each=wbreps)
              if ((!repeat.only) & randomize)
                 for (i in 0:(nblocks*bbreps-1))
                    for (j in 1:wbreps)
                    rand.ord[(i*blocksize*wbreps + (j-1)*blocksize + 1) : 
                          (i*blocksize*wbreps + j*blocksize)] <- 
                               sample((i%%nblocks*blocksize+1):(i%%nblocks*blocksize+blocksize))
                               
              if (repeat.only & randomize)
                 for (i in 0:(nblocks*bbreps-1))
                    rand.ord[(i*blocksize*wbreps + 1) : 
                          ((i+1)*blocksize*wbreps)] <- rep(sample((((i%%nblocks)*blocksize)+1):
                                        ((i%%nblocks+1)*blocksize)),each=wbreps)
              }
          else {
          #### implement randomization and replication for split-plots
          ####   standard approach; replications replicate complete plots
          ####   repeat.only generate repeated measurements within plots
                rand.ord <- rep(1:nruns,replications)
                if (replications > 1 & repeat.only) 
                                  rand.ord <- rep(1:nruns,each=replications)
                if ((!repeat.only) & randomize){
                    ## whole plots are replicated
                    ## might be more wise to use larger design

                    ## randomization on two levels
                    ## first: randomization within each whole plot
                    for (i in 1:(WPs*replications)) 
                        rand.ord[((i-1)*plotsize+1):(i*plotsize)] <- 
                            sample(rand.ord[((i-1)*plotsize+1):(i*plotsize)])
                    ## second: randomization of whole plots (in blocks per replication)
                    for (i in 1:replications){
                    ## only if not hard to change
                    if (is.null(hard)){
                        WPsamp <- sample(WPs) 
                        WPsamp <- (rep(WPsamp,each=plotsize)-1)*plotsize + rep(1:plotsize,WPs)
                        rand.ord[((i-1)*plotsize*WPs+1):(i*plotsize*WPs)] <- 
                            rand.ord[(i-1)*plotsize*WPs + WPsamp]
                        }
                    }
                    }
                if (repeat.only & randomize){
                    ## repeated measurements within whole plots
                    ## repetition of complete whole plots is not covered, 
                    ##     as I don't think that this makes sense

                    ## first: randomization within each whole plot
                    for (i in 1:WPs) 
                        rand.ord[((i-1)*plotsize*replications+1):(i*plotsize*replications)] <- 
                            rand.ord[(i-1)*plotsize*replications + 
                               rep(replications*(sample(plotsize)-1),each=replications) + 
                               rep(1:2,each=plotsize)]
                    ## second: randomization of whole plots
                    ## only if not hard to change
                    if (is.null(hard)){
                        WPsamp <- sample(WPs) 
                        WPsamp <- (rep(WPsamp,each=plotsize*replications)-1)*plotsize*replications + 
                             rep(1:(plotsize*replications),WPs)
                        rand.ord <- rand.ord[WPsamp]
                    }
                    }
          }
      }          
    orig.no <- rownames(desmat)
    orig.no <- orig.no[rand.ord]
    ## row added 27 01 2011 (for proper ordering of factor levels)
    orig.no.levord <- sort(as.numeric(orig.no),index=TRUE)$ix

    rownames(desmat) <- NULL
    desmat <- desmat[rand.ord,]
        
    if (is.list(blocks)) {
              Blocks <- Blocks[rand.ord]
              block.no <- block.no[rand.ord]
              }
    if (WPs > 1) wp.no <- wp.no[rand.ord]
      colnames(desmat) <- names(factor.names)
      ## adapt original number to replications
      if (is.list(blocks)) orig.no <- paste(orig.no,block.no,sep=".")
      if (WPs > 1) orig.no <- paste(orig.no,wp.no,sep=".")
      orig.no.rp <- orig.no
      if (bbreps * wbreps > 1){
           if (bbreps > 1) {
                ## !repeat.only covers all blocked cases and the repeat.only standard cases
                ## since bbreps stands in for replications
                if (repeat.only & !is.list(blocks))
                orig.no.rp <- paste(orig.no.rp, rep(1:bbreps,nruns),sep=".")
                else
                orig.no.rp <- paste(orig.no.rp, rep(1:bbreps,each=nruns*wbreps),sep=".")
           }
           if (wbreps > 1){
                ## blocked with within block replications
                if (repeat.only) 
                     orig.no.rp <- paste(orig.no.rp, rep(1:wbreps,nruns*bbreps),sep=".")
                else orig.no.rp <- paste(orig.no.rp, rep(1:wbreps,each=blocksize,nblocks*bbreps),sep=".")
           }
    }

    desdf <- data.frame(desmat)
    quant <- rep(FALSE,nfactors)
    for (i in 1:nfactors) {
        ## make all variables into factors
        desdf[,i] <- des.recode(desdf[,i],"-1=factor.names[[i]][1];1=factor.names[[i]][2]") 
        quant[i] <- is.numeric(desdf[,i])
        desdf[,i] <- factor(desdf[,i],levels=factor.names[[i]]) 
        contrasts(desdf[,i]) <- contr.FrF2(2)
        }

    if (is.list(blocks)) {
        desdf <- cbind(Blocks, desdf)
        ## reassign levels for Blocks to concatenation of manual build factors
        ## if they are built from only single factors
        ## !!! it is important to maintain the original order of numbered blocks factor
        ##     for more than four blocks (because e.g. grouping 1,2,3,7 vs. 4,5,6,8 not
        ##     orthogonal to everything else !
        hilf <- blocks
        if (all(sapply(hilf,length)==1) & !block.auto) {
            hilf <- as.numeric(hilf) + 1
            levels(desdf$Blocks) <- unique(apply(desdf[,hilf,drop=FALSE],1,paste,collapse=""))
            }
        colnames(desdf)[1] <- block.name
        ## delete any single block generator columns
        hilf <- blocks
        hilf <- as.numeric(hilf[which(sapply(hilf,length)==1)])
        if (length(hilf)>0) {
            desdf <- desdf[,-(hilf+1)]
            desmat <- desmat[,-hilf]
            factor.names <- factor.names[-hilf]
        }
        ## prepend block model matrix to desmat (not only generators but all block main effects)  
        desmat <- cbind(model.matrix(~desdf[,1])[,-1],desmat)
        colnames(desmat)[1:(2^k.block-1)] <- paste(block.name,1:(2^k.block-1),sep="")
        if (alias.info==3)
        hilf <- aliases(lm((1:nrow(desmat))~(.)^3,data=data.frame(desmat)))
        else 
        hilf <- aliases(lm((1:nrow(desmat))~(.)^2,data=data.frame(desmat)))
        ## check against blockpick.out here, if possible?
        ## blockpick.out$alias.2fis.block
        
        ## delete all aliases of Block factors themselves
        ## with or without "-"
        ## as these are not interesting
        if (length(hilf$aliases) > 0)
        for (i in 1:length(hilf$aliases)) {
             txt <- hilf$aliases[[i]]
             if (length(grep(paste("^-?",block.name,sep=""),txt)) > 0)
                  txt <- txt[-grep(paste("^-?",block.name,sep=""),txt)] 
                  if (length(txt)==length(grep("^-",txt))) txt <- sub("-","",txt)
                  hilf$aliases[[i]] <- txt
             }
        ## determine treatment effects that are aliased with blocks
        aliased.with.blocks <- hilf$aliases[1:(2^k.block-1)]
           aliased.with.blocks <- unlist(aliased.with.blocks)
        ## interim result for updating alias information
        if (nfactors<=50) leg <- paste(Letters[1:ntreat],names(factor.names),sep="=")
        else leg <- paste(paste("F",1:ntreat,sep=""),names(factor.names),sep="=")
        if (length(aliased.with.blocks)==0) aliased.with.blocks <- "none"
           else {
             aliased.with.blocks <- recalc.alias.block(aliased.with.blocks, leg)
             aliased.with.blocks <- aliased.with.blocks[ord(data.frame(nchar(aliased.with.blocks),aliased.with.blocks))]
             }
        ## determine treatment effects that are aliased with each other
        aliased <- hilf$aliases[-(1:(2^k.block-1))]
        aliased <- aliased[which(sapply(aliased,length)>1)]
        ## update: same format like aliased element of unblocked designs
        if (length(aliased)>0) aliased <- struc.aliased(recalc.alias.block(aliased, leg), nfactors, alias.info)

           
        ## prepare design info for blocked designs
        ntreat <- ncol(desdf) - 1
        if (block.auto) factor.names <- factor.names[1:ntreat]
 #       if (block.auto) factor.names <- factor.names[setdiff(1:nfactors,blocks)]
 
 ## up to version 0.96-1 nfactors was ntreat+1
 #         design.info <- list(type="FrF2.blocked", block.name=block.name, 
 #           nruns=nruns, nfactors=ntreat+1, nblocks=nblocks, blocksize=blocksize, 
 #           ntreat=ntreat,factor.names=factor.names,
 #           aliased.with.blocks=aliased.with.blocks, aliased=aliased,
 #           bbreps=bbreps, wbreps=wbreps)
          design.info <- list(type="FrF2.blocked", block.name=block.name, 
            nruns=nruns, nfactors=ntreat, nblocks=nblocks, block.gen=block.gen, blocksize=blocksize, 
            ntreat=ntreat,factor.names=factor.names,
            aliased.with.blocks=aliased.with.blocks, aliased=aliased,
            bbreps=bbreps, wbreps=wbreps, 
            FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version)
        if (!is.null(generators)) {
           if (g>0) 
           design.info <- c(design.info, 
                list(base.design=paste("generator columns:", paste(cand[[1]]$gen, collapse=", "))))
           else 
           design.info <- c(design.info, 
                list(base.design="full factorial"))
                }
           else design.info <- c(design.info, list(catlg.name = catlg.name, base.design=names(cand[1])))
           design.info <- c(design.info, list(map=map))
        
        if (bbreps>1) {
            hilflev <- paste(rep(levels(desdf[,1]), each=bbreps), rep(1:bbreps, nblocks), sep=".")
            desdf[,1] <- factor(paste(desdf[,1], rep(1:bbreps, each=nruns*wbreps),sep="."), levels=hilflev)
        }
            ## make block names reflect the between block replication
        }   ## end of blocked designs

        if (WPs > 1){
            
            ## output for split-plot designs
            if (alias.info==3)
            aliased <- aliases(lm((1:nrow(desmat))~(.)^3,data=data.frame(desmat)))$aliases
            else
            aliased <- aliases(lm((1:nrow(desmat))~(.)^2,data=data.frame(desmat)))$aliases
            aliased <- aliased[which(sapply(aliased,length)>1)]
            ## update: same format with aliased entry of other FrF2 designs
            if (length(aliased) > 0){ 
                if (nfactors<=50) leg <- paste(Letters[1:nfactors],names(factor.names),sep="=")
                else leg <- paste(paste("F",1:nfactors,sep=""),names(factor.names),sep="=")
                aliased <- struc.aliased(recalc.alias.block(aliased, leg), nfactors, alias.info)
                }

            design.info <- list(type="FrF2.splitplot", 
                nruns=nruns, nfactors=nfactors, nfac.WP=nfac.WP, nfac.SP=nfactors-nfac.WP, 
                      factor.names=factor.names,
                nWPs=WPs, plotsize=nruns/WPs, 
                res.WP=res.WP, aliased=aliased, FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version)
            ## below: changed entry in base.design to original generators, combined with map
            ## instead of cand[[1]]$gen from before
            if (!is.null(generators)) 
                design.info <- c(design.info, 
                     list(base.design=paste("generator columns:", 
                     paste(which(names(Yates)[1:(nruns-1)] %in% names(generators)), collapse=", ")), map=map, 
                     orig.fac.order = c(orignew, setdiff(1:nfactors,orignew))))
                else design.info <- c(design.info, list(catlg.name = catlg.name, base.design=names(cand[1]), map=map, 
                     orig.fac.order = c(orignew, setdiff(1:nfactors,orignew))))
                }

    if (is.null(estimable) & is.null(generators) & !(is.list(blocks) | WPs > 1))
        design.info <- list(type="FrF2", nruns=nruns, nfactors=nfactors, factor.names=factor.names, 
            catlg.name = catlg.name,
            catlg.entry=cand[1], aliased = alias3fi(k,cand[1][[1]]$gen,order=alias.info), 
            FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version)
        ## incorporate reasonable further info (blocks, WPs etc.)
    if ((!is.null(generators)) & !is.list(blocks) & !WPs > 1) {
         if (nfactors <= 50)
         names(generators) <- Letters[(k+1):nfactors]
         else names(generators) <- paste("F",(k+1):nfactors, sep="")
         gen.display <- paste(names(generators),sapply(generators,function(obj){
               if (nfactors <= 50)
               paste(if (obj[1]<0) "-" else "", paste(Letters[abs(obj)],collapse=""),sep="")
               else 
               paste(if (obj[1]<0) "-" else "", paste(paste("F",abs(obj),sep=""),collapse=":"),sep="")}), sep="=")
         design.info <- list(type="FrF2.generators", nruns=nruns, nfactors=nfactors, factor.names=factor.names, generators=gen.display, 
               aliased = alias3fi(k,generators,order=alias.info), 
               FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version)
         }
    aus <- desdf
      rownames(aus) <- rownames(desmat) <- 1:nrow(aus)
    class(aus) <- c("design","data.frame")
    attr(aus,"desnum") <- desmat
      ## change 27 Jan 2011: leave orig.no as a factor, but with better-ordered levels
      orig.no <- factor(orig.no, levels=unique(orig.no[orig.no.levord]))
      ## added 14 Feb 2013
      orig.no.rp <- factor(orig.no.rp, levels=unique(orig.no.rp[orig.no.levord]))
    if (!(is.list(blocks) | WPs > 1)) 
       attr(aus,"run.order") <- data.frame("run.no.in.std.order"=orig.no,"run.no"=1:nrow(desmat),"run.no.std.rp"=orig.no.rp)
    else attr(aus,"run.order") <- data.frame("run.no.in.std.order"=orig.no,"run.no"=1:nrow(desmat),"run.no.std.rp"=orig.no.rp)
    ## make repeat.only reflect the calculated status instead of the status requested by the user 
    ## (which can be seen in the creator element)
    if (design.info$type=="FrF2.blocked") {
        if (design.info$wbreps==1) repeat.only <- FALSE
        nfactors <- ntreat
    }
    if (nfactors<=50) design.info$aliased <- c(list(legend=paste(Letters[1:nfactors],names(factor.names),sep="=")),design.info$aliased)
       else design.info$aliased <- c(list(legend=paste(paste("F",1:nfactors,sep=""),names(factor.names),sep="=")),design.info$aliased)
    attr(aus,"design.info") <- c(design.info, list(replications=replications, repeat.only=repeat.only,
      randomize=randomize, seed=seed, creator=creator))
    ## add center points, if requested
    if (ncenter>0) aus <- add.center(aus, ncenter, distribute=center.distribute)
    aus
} 