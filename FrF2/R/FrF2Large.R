nrunsV <- function(nfactors){
   ## make sure to include the latest knowledge
   ## Ryan and Bulutoglu 2010 Technometrics:
   ## Xu 25 to 28 in 1024 runs are MA
   ## Xu 24 to 27 in 2048 runs are MA
   ## Xu 28 to 32 in 2048 runs have been improved by Ryan and Bulutoglu 
       ### change bound in fiveFrF2MA as soon as incorporated in catalogue
   ## Xu 25 factors in 4096 runs is MA
   ## Xu 26 factors has been improved by Ryan and Bulutoglu
       ### change bound in fiveFrF2MA as soon as incorporated in catalogue
   ## 
   fiveFrF2MA <- data.frame(nfactors=c(2,3,5,6,8,11,17,23,28,27,25),nruns=2^(2:12))
   fiveFrF2 <- data.frame(nfactors=c(2,3,5,6,8,11,17,23,33,47,65),nruns=2^(2:12))
   fiveSanchezSanchez <- data.frame(nfactors=c(2,3,5,6,8,11,17,21,29,38,52,69,92,120),nruns=2^(2:15))
   oldwarn <- options("warn")
   options(warn=-1)
   FrF2MA <- min(fiveFrF2MA$nruns[which(fiveFrF2MA$nfactors >= nfactors)])
   FrF2 <- min(fiveFrF2$nruns[which(fiveFrF2$nfactors >= nfactors)])
   SanchezSanchez <- min(fiveSanchezSanchez$nruns[which(fiveSanchezSanchez$nfactors >= nfactors)])
   options(warn=oldwarn$warn)
   if (SanchezSanchez == Inf) message("no regular resolution V design can be created with this software")
   else if (SanchezSanchez<FrF2) message("Only function FrF2Large with nruns=", SanchezSanchez, " creates a regular resolution V design.")
         else if (FrF2MA<=FrF2) message("Function FrF2 with nruns=", FrF2, " creates a minimum aberration regular resolution V design.")
              else message("Function FrF2 with nruns=", FrF2, " creates a regular resolution V design (good, but not necessarily MA).")
   invisible(min(FrF2MA, FrF2, SanchezSanchez))
}

FrF2Large <- function(nruns, nfactors=NULL,
                 factor.names = if(!is.null(nfactors)) {if(nfactors<=50) Letters[1:nfactors]
                                else paste("F",1:nfactors,sep="")} else NULL,
                 default.levels = c(-1,1), ncenter=0, center.distribute=NULL,
                 generators=NULL, replications=1, repeat.only=FALSE,
                 randomize=TRUE, seed=NULL, alias.info=2, ...){
     creator <- sys.call()
     ## provide the SanchezSanchez column numbers
     ## which are NOT a recipe for MA
     ## these can be used with function YatesFly
      SanchezSanchez <- data.frame(nfactors=1:120,
          gen=c(1, 2, 4, 8, 15, 16, 32, 51, 64, 85, 106, 128,
               150, 171, 219, 237, 247, 256, 279, 297, 455, 512, 537,
               557, 594, 643, 803, 863, 998, 1024, 1051, 1070, 1112,
               1169, 1333, 1345, 1620, 1866, 2048, 2076, 2085, 2185,
               2372, 2456, 2618, 2800, 2873, 3127, 3284, 3483, 3557,
               3763, 4096, 4125, 4135, 4174, 4435, 4459, 4469, 4497,
               4752, 5255, 5732, 5804, 5915, 6100, 6369, 6907, 7069,
               8192, 8263, 8351, 8422, 8458, 8571, 8750, 8858, 9124,
               9314, 9500, 10026, 10455, 10556, 11778, 11885, 11984,
               13548, 14007, 14514, 14965, 15125, 15554, 16384, 16457,
               16517, 16609, 16771, 16853, 17022, 17453, 17891, 18073,
               18562, 18980, 19030, 19932, 20075, 20745, 21544, 22633,
               23200, 24167, 25700, 26360, 26591, 26776, 28443, 28905,
               29577, 32705))
    ##check nruns
    k <- round(log2(nruns))
       if (!2^k==nruns) stop("nruns must be a power of 2.")
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
## check that no incompatible options are used together
    if (ncenter>0) if (center.distribute > nruns + 1)
        stop("center.distribute must not be larger than nruns+1")
## simple checks for individual option entries
    if (default.levels[1]==default.levels[2]) stop("Both default levels are identical.")
    if (!(is.logical(repeat.only) & is.logical(randomize)))
         stop("repeat.only and randomize must be logicals (TRUE or FALSE).")
    if (!is.numeric(replications)) stop("replications must be an integer number.")

    if (!alias.info %in% c(2,3))
         stop("alias.info can be 2 or 3 only.")
    if (!(is.numeric(default.levels) | is.character(default.levels)))
         stop("default.levels must be a numeric or character vector of length 2")
    if (!length(default.levels) ==2)
         stop("default.levels must be a numeric or character vector of length 2")

## check factor specifications
    if (is.null(factor.names) & is.null(nfactors) & is.null(generators))
         stop("The number of factors must be specified via nfactors, via factor.names, or via generators.")
    if (!is.null(factor.names) & !(is.character(factor.names) | is.list(factor.names)) )
         stop("factor.names must be a character vector or a list.")
    if (is.null(nfactors)) {if (!is.null(factor.names)) nfactors <- length(factor.names)
                          else if (!is.null(generators)) nfactors <- length(generators)+k
                          }
    ## from here on, nfactors is known
    if (!nfactors==floor(nfactors))
         stop("nfactors must be an integer number.")
    if (!is.null(factor.names) & !length(factor.names)==nfactors)
         stop("There must be nfactors factor names, if any.")
    if (is.null(factor.names))
         if(nfactors<=50) factor.names <- Letters[1:nfactors] else factor.names <- paste("F",1:nfactors,sep="")

## warnings to prevent users from using inferior designs
#    if ((nruns >= 4 & nruns <= 512) | (nruns == 1024 & nfactors <= 24)
#        | (nruns == 2048 & nfactors <= 23)  | (nruns == 4096 & nfactors <= 24))
#        stop("Please use function FrF2, it will guarantee an optimal design.")
#    if (nruns >= 1024 & nruns <= 4096)
#        warning("Usage of function FrF2 might yield a better design \n(up to 33 factors in 1024 runs, \nup to 47 factors in 2048 runs, \nup to 65 factors in 4096 runs)")

    if (nruns <= 4096) stop("Please use function FrF2; it does everything that FrF2Large can do and more.") 

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
      generators <- gen.check(k, generators)   ## check entries for validity
      hilf <- YatesFly(c(2^((1:k)-1), gencalc(generators)), k=k)
      }
    else hilf <- YatesFly(SanchezSanchez$gen[1:nfactors], k=k)

    gen <- hilf$gen    ## vector of column numbers for non-base factors
    if (is.null(generators)) generators <- hilf$Yates.brief[-(1:k)]
                       ## list of column combinations

    ## number of generators
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

               gen <- gencalc(generators)
               gen <- gen*sapply(generators, function(obj) sign(obj[1]))

               cand <- list(custom=list(res=res, nfac=nfactors, nruns=nruns,
                    gen=gen,
                    WLP=WLP, nclear.2fis=nclear.2fis, clear.2fis=clear.2fis, all.2fis.clear=all.2fis.clear))
                   ## needs to be list of list, because access later is always via cand[1]
               class(cand) <- c("catlg","list")
      ## prepare output
      desmat <- hilf$mat
      rownames(desmat) <- 1:nruns
      if (randomize & !is.null(seed)) set.seed(seed)
      rand.ord <- rep(1:nruns,replications)
      if (replications > 1 & repeat.only) rand.ord <- rep(1:nruns,each=replications)
      if (randomize & !repeat.only) for (i in 1:replications)
                  rand.ord[((i-1)*nruns+1):(i*nruns)] <- sample(nruns)
      if (randomize & repeat.only) rand.ord <- rep(sample(1:nruns), each=replications)
        orig.no <- rownames(desmat)
        orig.no <- orig.no[rand.ord]
    ## row added 27 01 2011 (for proper ordering of design)
    orig.no.levord <- sort(as.numeric(orig.no),index=TRUE)$ix
        rownames(desmat) <- NULL
        desmat <- desmat[rand.ord,]
        
      ## order information for output
      orig.no.rp <- orig.no
      if (replications > 1){
                if (repeat.only)
                orig.no.rp <- paste(orig.no.rp, rep(1:replications,nruns),sep=".")
                else
                orig.no.rp <- paste(orig.no.rp, rep(1:replications,each=nruns),sep=".")
           }

      ## prepare generators to be character for output
      if (nfactors <=50)
        generators <- paste(Letters[(k+1):nfactors],colnames(desmat)[(k+1):nfactors],sep="=")
      else
        generators <- paste(paste("F",1:nfactors,sep=""),colnames(desmat)[1:nfactors],sep="=")
      ## factor names for design matrix
      colnames(desmat) <- names(factor.names)

      ## create output data frame
      desdf <- as.data.frame(desmat)
      quant <- rep(FALSE,nfactors)
      for (i in 1:nfactors) {
        ## make all variables into factors
        desdf[,i] <- des.recode(desdf[,i],"-1=factor.names[[i]][1];1=factor.names[[i]][2]")
        quant[i] <- is.numeric(desdf[,i])
        desdf[,i] <- factor(desdf[,i],levels=factor.names[[i]])
        contrasts(desdf[,i]) <- contr.FrF2(2)
        }
      attr(desdf, "desnum") <- desmat
      ## change 27 Jan 2011: leave orig.no as a factor, but with better-ordered levels
      orig.no <- factor(orig.no, levels=unique(orig.no[orig.no.levord]))
      attr(desdf, "run.order") <- data.frame("run.no.in.std.order"=orig.no,
           "run.no"=1:nrow(desmat),"run.no.std.rp"=orig.no.rp)
      design.info <- list(type="FrF2.large", nruns=nruns,
          nfactors=nfactors, factor.names=factor.names, generators=generators,
          aliased = alias3fi(k,cand[1][[1]]$gen,order=alias.info),
          replications=replications,
          repeat.only=repeat.only, randomize=randomize, seed=seed, creator=creator, 
          FrF2.version = sessionInfo(package="FrF2")$otherPkgs$FrF2$Version)
      if (nfactors<=50) design.info$aliased <- c(list(legend=paste(Letters[1:nfactors],names(factor.names),sep="=")),design.info$aliased)
          else design.info$aliased <- c(list(legend=paste(paste("F",1:nfactors,sep=""),names(factor.names),sep="=")),design.info$aliased)
      attr(desdf, "design.info") <- design.info
      class(desdf) <- c("design", "data.frame")
      ## add center points, if requested
      if (ncenter>0) desdf <- add.center(desdf, ncenter, distribute=center.distribute)
      ## hier weiter: sicherstellen, dass dieses alle Voraussetzungen an FrF2-artige Objekte erfüllt
      desdf
}