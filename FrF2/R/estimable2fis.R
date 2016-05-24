#require(gtools)
## replace alias3fi by alias2fi ?
## incorporate map into output object of FrF2
## thus enable alias matrix calculation for all designs, even if columns have been swapped

## functions estimable determine design with some additional info
## function mapcalc determines a map - here, it should be possible to switch between really requesting "clear" or requesting "distinct"

#estimable <- function(estimable, ...){
#     UseMethod("estimable")
#}

mapcalc <- function (estimable, nfac, nruns, res3 = FALSE, select.catlg = catlg, 
    method="VF2", sort = "natural") 
{
## not clear whether map is correct for lad
    if (ncol(estimable) > nruns - nfac - 1) 
        stop("too many interactions requested for this number of runs and factors")
    if (nruns > 32 & res3) {
        warning("res3=TRUE has no effect for designs with more than 32 runs.")
        res3 <- FALSE
    }
    catlg <- select.catlg
    go2 <- graph.empty(n = nfac, directed = FALSE)
    go2 <- add.edges(go2, estimable)
    ## optionally sort vertices by degree, 20 Jul 2012
   deg2 <- degree(go2)
    if (sort %in% c("high", "low")) {
        if (sort == "low") 
            ord2 <- order(deg2)
        else ord2 <- order(deg2, decreasing = TRUE)
        go2 <- permute.vertices(go2, invperm(ord2))
    }
    degree2 <- rev(cumsum(rev(table(deg2))))   ## make it faster to reject non-isomorphic cases
                                      ## 7 Feb 2011
    degs2 <- as.numeric(names(degree2))        ## required minimum degrees
    ## added further pre-filtering criteria 9 July 2012
    indep2 <- independence.number(go2)         ## required maximum independence number
    clique2 <- clique.number(go2)              ## required minimum clique size
    ## reduced attention to dominating designs, if applicable (29 June 2012)
    tobechecked <- catlg[which(nfac(catlg) == nfac & nruns(catlg) == 
        nruns & dominating(catlg))]
    if (length(tobechecked) > 0) 
        if (!res3) 
            tobechecked <- tobechecked[which(res(tobechecked) >= 
                4)]
    if (length(tobechecked) > 0) 
        tobechecked <- tobechecked[which(nclear.2fis(tobechecked) >= 
            ncol(estimable))]
    if (length(tobechecked) == 0) 
        if (res3) 
            stop("The required interactions cannot be accomodated clear of aliasing in ", 
                nruns, " runs.")
        else {
            if (nruns <= 64) 
                stop("The required interactions cannot be accomodated clear of aliasing in ", 
                  nruns, " runs with resolution IV or higher.")
            else stop("No resolution IV or higher design from the current catalogue can accomodate the required interactions clear of aliasing in ", 
                nruns, " runs.")
        }
    map <- NULL
    ## .FrF2.currentlychecked is stored in private environment
    ## for enabling user to access it after aborting without success
    for (i in 1:length(tobechecked)) {
        putFrF2(".FrF2.currentlychecked", names(tobechecked[i]))
        go1 <- graph.empty(n = nfac, directed = FALSE)
          ## previous version subtracted 1 from tobechecked[[i]]$clear.2fis 
          ## for previous igraph node definition; changed 29/06/2012
        go1 <- add.edges(go1, tobechecked[[i]]$clear.2fis)
       ## optionally sort vertices by degree, 20 Jul 2012
        deg1 <- degree(go1)
        if (sort %in% c("high", "low")) {
            if (sort == "low") 
                ord1 <- order(deg1)
            else ord1 <- order(deg1, decreasing = TRUE)
            go1 <- permute.vertices(go1, invperm(ord1))
        }
        degree1 <- rev(cumsum(rev(table(deg1))))
        degs1 <- as.numeric(names(degree1))
        if (max(degs2) <= max(degs1)) {
             ## check for no chance, 7.2.2011
             ## if max(degs2)>max(degs1), subgraph isomorphism is impossible
            comp <- sapply(degs2, function(obj) degree1[min(which(degs1 >= 
                obj))])
            comp[is.na(comp)] <- 0
            if (any(comp < degree2)) 
                next
          ## added further pre-filtering criteria 9 July 2012
            if (independence.number(go1) > indep2) 
                next
            if (clique.number(go1) < clique2) 
                next
            if (method=="LAD") 
            erg <- graph.subisomorphic.lad(go2, go1)
            else erg <- graph.subisomorphic.vf2(go1, go2)
          ## +1 removed from map21 because of adapting to igraph (29 June 2012)
            if (erg$iso) {
                if (method=="LAD") map <- list(erg$map)
                else map <- list(erg$map21)
                if (sort %in% c("high", "low")) 
                  map <- list(ord1[map[[1]]][invperm(ord2)])
                names(map) <- getFrF2(".FrF2.currentlychecked")
                putFrF2(".FrF2.currentlychecked", names(tobechecked[i]))
                break
            }
        }
    }
    if (is.null(map)) 
        if (res3) 
            stop("The required interactions cannot be accomodated clear of aliasing in ", 
                nruns, " runs.")
        else {
            if (nruns <= 64) 
                stop("The required interactions cannot be accomodated clear of aliasing in ", 
                  nruns, " runs with resolution IV or higher.")
            else stop("No resolution IV or higher design from the current catalogue can accomodate the required interactions clear of aliasing in ", 
                nruns, " runs.")
        }
    map
}

##map <-FrF2:::mapcalc(matrix(c(1,2,1,3,1,4,2,3,3,5,4,8,3,7),2,7),9,64)
#print(FrF2:::mapcalc(matrix(c(1,2,1,3,1,4,1,5,3,5,4,5,5,7),2,7),9,32,sort="low"))
#print(FrF2:::mapcalc(matrix(c(1,2,1,3,1,4,1,5,3,5,4,5,5,7),2,7),9,32,sort="high"))
#print(FrF2:::mapcalc(matrix(c(1,2,1,3,1,4,1,5,3,5,4,5,5,7),2,7),9,32))

mapcalc.distinct <- function(estimable, nfac, nruns, res3=FALSE, max.time=60, select.catlg=catlg, perm.start=1:nfac, perms=NULL){
## words of length three not correctly represented yet
      ## called with checked inputs only --> no checks needed
      ## estimable is a matrix with two rows
      if (ncol(estimable)>nruns - nfac -1) stop("too many interactions requested for this number of runs and factors")

      ## implement further structural checks of impossibility (e.g. G1^2+G2^2, if union of G1 and G2 is everything)

      hilf1 <- combn(nfac,2)
      catlg <- select.catlg

      ## work through candidate designs
      if (res3) tobechecked <- catlg[which(nfac.catlg(catlg)==nfac & nruns.catlg(catlg)==nruns)]
      else tobechecked <- catlg[which(nfac.catlg(catlg)==nfac & nruns.catlg(catlg)==nruns & res.catlg(catlg)>=4)]
      if (length(tobechecked)==0){
              if (res3)
              stop("The required interactions cannot be accomodated on distinct columns in ", nruns, " runs.")
              else
              stop("The required interactions cannot be accomodated on distinct columns in ", nruns, " runs with resolution IV or higher.")
              }      
      map <- NULL

   begin_time <- Sys.time()
      ## start time to be handed to check.subisomorphism.special
      ## for stopping calculations in case of long run times
   for (i2 in 1:length(tobechecked)){
      ## determine words of length up to four
      hilf2 <- words.all(round(log2(nruns)),tobechecked[[i2]]$gen,max.length=4)
      n3 <- hilf2$WLP["3"]
      if (!is.na(n3)) hilf3 <- hilf2$words.up.to.length.4[which(sapply(hilf2$words.up.to.length.4,length)==3)]
           else hilf3 <- NULL
          ## list of vectors of length 3
      n4 <- hilf2$WLP["4"]
      if (!is.na(n4)) hilf2 <- hilf2$words.up.to.length.4[which(sapply(hilf2$words.up.to.length.4,length)==4)] 
           else hilf2 <- NULL
          ## list of vectors of length 4
      if (is.na(3) & is.na(n4)) stop("Resolution V+ design, use any allocation of experiment factors to design factors.")
      if (is.null(perms)) erg <- check.subisomorphic.special(estimable,nfac,hilf2,hilf3=hilf3,max.time=max.time,begin_time=begin_time,
               name=names(tobechecked[i2]),res3=res3,perm.start=perm.start)
      else erg <- check.subisomorphic.matrix(estimable,nfac,hilf2,hilf3=hilf3,max.time=max.time,begin_time=begin_time,
               name=names(tobechecked[i2]),res3=res3,perms=perms)
          if (erg$iso) {
                    ## "+1" removed from map21 because of adapting to igraph (29 June 2012)
                    map <- list(erg$map21[1:nfac])
                    names(map) <- names(tobechecked[i2])
                    break}
     }
     if (is.null(map)) {
              if (res3)
              stop("The required interactions cannot be accomodated on distinct columns in ", nruns, " runs.")
              else
              stop("The required interactions cannot be accomodated on distinct columns in ", nruns, " runs with resolution IV or higher.")            
            }
     map
}

check.subisomorphic.special <- function(estimable, nfac, hilf2, hilf3=NULL, res3=FALSE, max.time=60, perm.start=1:nfac, 
    begin_time=Sys.time(), name=NA){
    ## hilf2 contains the list of length 4 words, hilf3 length 3 words (if applicable)
    ## begin_time is current time, unless handed over from previously determined timing
    map <- perm.start

    iso <- TRUE ## default: unpermuted design does it

    if (is.null(hilf3)) res3.cur <- FALSE else res3.cur <- TRUE
    if (!is.null(hilf2)) colpairs <- matrix(NA,4,choose(ncol(estimable),2))

    if (res3.cur) {
       ## if any 2fi from requirement set in length three word --> not adequate
       if (any(apply(estimable,2,function(objcp) any(sapply(hilf3,function(obj) all(objcp %in% obj)))))) 
            iso <- FALSE    
    } 
    ## if not yet discarded, check four letter words
    if (ncol(estimable) > 1) hilf <- combn(ncol(estimable),2)
    if (iso & !is.null(hilf2) & ncol(estimable) > 1){
       ## check estimable for two columns being in same word
       ## words are sorted in ascending order
        for (i in 1:ncol(colpairs)) colpairs[,i] <- sort(c(estimable[,hilf[,i]]))
        if (any(apply(colpairs,2,function(objcp) any(sapply(hilf2,function(obj) all(obj==objcp))))))
            iso <- FALSE
    }
    if (!iso){
       ## if iso FALSE for the unpermuted design, loop through others

    ## loop through all permutations of map
    ## adjust estimable accordingly
    ## repeat above check
    for (i in 2:factorial(nfac)) {
      ## loop makes sure that getNext does not continue
      ## could be omitted if getNext would understand that it is finished
       if (difftime(Sys.time(),begin_time,units="sec") > max.time)               
              stop("A solution could not be found in the allocated max.time of ",max.time," seconds. \n",
                   if (!is.na(name)) paste("Final design:",name), ", Final permutation:", paste(map,collapse=",") )
      map <- getNext(map)
      iso <- TRUE ## start by assuming that current permutation does it

      ### adjust numbers in estimable to permutation
      estcur <- matrix(map[estimable],nrow=2)

    ## check here
    if (res3.cur) {
       ## if any 2fi from requirement set in length three word --> not adequate
       if (any(apply(estcur,2,function(objcp) any(sapply(hilf3,function(obj) all(objcp %in% obj)))))) 
            iso <- FALSE    
    } 
    ## bug fix with version 1.3: restricted to ncol(estimable) > 1
    if (iso & !is.null(hilf2) & ncol(estimable) > 1){
       for (i in 1:ncol(colpairs)) colpairs[,i] <- sort(c(estcur[,hilf[,i]]))
       if (any(apply(colpairs,2,function(objcp) any(sapply(hilf2,function(obj) all(obj==objcp)))))){
            iso <- FALSE
            }
    }
    if (iso) break
    }
    }
    ## previous version subtracted 1 from map for previous igraph node definition; changed 29/06/2012
    if (iso) aus <- list(iso=iso, map21 = map) 
        else aus <- list(iso=iso)
    aus
}

check.subisomorphic.matrix <- function(estimable, nfac, hilf2, hilf3=NULL, res3=FALSE, max.time=60, 
    begin_time=Sys.time(), name=NA, perms=perms){
    ## hilf2 contains the list of length 4 words, hilf3 length 3 words (if applicable)
    ## begin_time is current time, unless handed over from previously determined timing

    ## loop through perms
    ## adjust estimable accordingly
    ## repeat above check
    if (is.null(hilf3)) res3.cur <- FALSE else res3.cur <- TRUE
    if (ncol(estimable)>=2) hilf <- combn(ncol(estimable),2)
    else hilf <- matrix(1,1,1)
    if (!(is.null(hilf2) | ncol(estimable)<2)) {
        colpairs <- matrix(NA,4,choose(ncol(estimable),2))
        }

    for (i in 1:nrow(perms)) {
      ## loop makes sure that getNext does not continue
      ## could be omitted if getNext would understand that it is finished
       if (difftime(Sys.time(),begin_time,units="sec") > max.time)               
              stop("A solution could not be found in the allocated max.time of ",max.time," seconds. \n",
                   if (!is.na(name)) paste("Final design:",name), ", Final permutation:", paste(map,collapse=",") )
      map <- perms[i,]
      iso <- TRUE ## start by assuming that current permutation does it

      ### adjust numbers in estimable to permutation
      estcur <- matrix(map[estimable],nrow=2)

    ## check here
    if (res3.cur) {
       ## if any 2fi from requirement set in length three word --> not adequate
       if (any(apply(estcur,2,function(objcp) any(sapply(hilf3,function(obj) all(objcp %in% obj)))))) 
            iso <- FALSE    
    } 
    if (iso & !(is.null(hilf2) | ncol(estimable)<2)){
       ## create colpairs (=forbidden words of length four)
       for (i in 1:ncol(colpairs)) colpairs[,i] <- sort(c(estcur[,hilf[,i]]))
       if (any(apply(colpairs,2,function(objcp) any(sapply(hilf2,function(obj) all(obj==objcp)))))){
            iso <- FALSE
            }
    }
    if (iso) break
    }
    ## previous version subtracted 1 from map for previous igraph node definition; changed 29/06/2012
    if (iso) aus <- list(iso=iso, map21 = map) 
        else aus <- list(iso=iso)
    aus
}

map2design <- function(map, select.catlg=catlg){
          hilf <- select.catlg[[names(map)]]
          hilf <- FrF2(nruns=hilf$nruns,nfactors=hilf$nfac,generators=hilf$gen,randomize=FALSE)[,map[[1]]]
          colnames(hilf) <- Letters[1:ncol(hilf)]
          hilf
    }
## test <- map2design(map)
## aliases(lm(rnorm(nrow(test))~.^2,test))

formula2matrix <- function(formula) {
          aus <- apply(attr(terms(formula(formula)),"factors")[,
                attr(terms(formula(formula)),"order")==2,drop=FALSE],2,function(obj) which(obj==1))
          if (!is.matrix(aus)) stop("The formula for option estimable must contain main effects for all factors that occur in interactions!")
          aus
    }

char2num <- function(char){
          char <- strsplit(char,"")
          matrix(sapply(char, function(obj) match(obj,Letters)),nrow=2)
    }

estimable.check <- function(estimable,nfac,factor.names){
    ## function to check validity of estimable entries
    ## returns estimable in matrix form
    ## in case of formula, factor.names is used for making sure that correct factors are picked
    if (is.character(factor.names)) fn <- factor.names else fn <- names(factor.names)
    if (!(is.matrix(estimable) | inherits(estimable,"formula") | is.character(estimable)))
         stop("estimable must be a numeric matrix with two rows or a formula or a character vector with length two entries")
    if (is.character(estimable)){
         if (any(!nchar(estimable)==2)) stop("All entries of character vector estimable must be of length 2.")
         if (is.null(nfac)) {if (!length(setdiff(hilf<-unique(unlist(strsplit(estimable,""))),Letters))==0)
             stop("Only factor letters from constant Letters are permitted in estimable.")
             nfac <- max(which(Letters %in% hilf))}
         else { 
         if (!length(setdiff(unique(unlist(strsplit(estimable,""))),Letters[1:nfac]))==0)
             stop("Only factor letters for nfac factors (Letters[1:nfac]) are permitted in estimable.")
             }
         estimable <- list(estimable=char2num(unique(estimable)),nfac=nfac)
    }
    if (inherits(estimable,"formula")){
          if (any(attr(terms(formula(estimable)),"order")>2)) stop("terms of order higher than 2 are not allowed in formula estimable")
          if (!is.null(nfac))
            if (nrow(attr(terms(formula(estimable)),"factors"))>nfac) stop("more than nfac main effects in formula estimable")
          if (all(!attr(terms(formula(estimable)),"order")==2)) {warning("no interaction requested in formula estimable, estimable is ignored.")
               estimable <- NULL}
          ### implement comparison to factor.names ###
          hilf <- row.names(attr(terms(formula(estimable)),"factors"))
          hilf2 <- 1:length(hilf)
          if (!length(setdiff(hilf,fn))==0){
                 if (!length(setdiff(hilf,Letters))==0)
                      stop("formula must refer to elements of Letters or to factor names.")
                 if (!is.null(nfac)) if (!length(setdiff(hilf,Letters[1:nfac]))==0)
                      stop("Factor letters are higher than permitted for number of factors.")
                 for (i in 1:length(hilf))
                     hilf2[i] <- which(Letters == hilf[i])
                 if (is.null(nfac)) nfac <- max(hilf2)
          }
          else { 
                 for (i in 1:length(hilf)){
                     hilf2[i] <- which(fn == hilf[i])
                     ## nfac is known from factor.names
               }
          }
         if (!is.null(estimable)) estimable <- list(estimable = formula2matrix(estimable), nfac=nfac)
         ## adjust factor numbers in case of factor names
         estimable[[1]] <- matrix(hilf2[estimable[[1]]],nrow=2)
    }
    if (is.matrix(estimable)){
         if (!is.numeric(estimable)) stop("matrix estimable must be numeric")
         if (!nrow(estimable)==2) stop("matrix estimable must have two rows")
         if (any(!round(estimable)==estimable)) stop("entries of matrix estimable must be integer numbers") 
         if (max(estimable)>nfac | min(estimable)<1) stop("entries of estimable must be between 1 and nfac")
         if (is.null(nfac)) nfac <- max(estimable)
         estimable <- list(estimable=estimable, nfac=nfac)
    }
    ## output the matrix of estimable interactions
    ## order as in factor.names or factor letters in case of formula
    estimable
}

estimable <- function(estimable, nfac, nruns, 
                      clear=FALSE, res3=FALSE, max.time=60, 
                      select.catlg=catlg, 
                      perm.start=1:nfac, perms=NULL,
                      order = 3, method="VF2", sort="natural"){
        if (clear) {
            ## perhaps relabel estimable entries such that most frequently occurring ones in front
            estimable <- 
            map <- mapcalc(estimable,nfac,nruns,res3=res3, select.catlg=select.catlg, method=method, sort=sort)
        }
           else map <- mapcalc.distinct(estimable,nfac,nruns,res3=res3, max.time=max.time, 
                                select.catlg=select.catlg, perm.start=perm.start, perms=perms)
        test <- map2design(map, select.catlg=select.catlg)
        hilf <- select.catlg[names(map)][[1]]
        ## determine aliased here because of potential remap of factors
        hilf <- alias3fi(hilf$nfac-length(hilf$gen),hilf$gen, order=order)
        if (is.list(hilf)) aus <- list(map=map,design=test,
            aliased=list(main=sort(chartr(paste(Letters[map[[1]]],collapse=""),paste(Letters[1:nfac],collapse=""),hilf$main)),
                         fi2=sort(chartr(paste(Letters[map[[1]]],collapse=""),paste(Letters[1:nfac],collapse=""),hilf$fi2))))
        else aus <- list(map=map,design=test, aliased=hilf)
        aus
    }

#estimable.formula <- function(estimable, nfac, nruns,  clear=FALSE, res3=FALSE, max.time=60, select.catlg=catlg, perm.start=1:nfac, perms=NULL){
#        if (clear) map <- mapcalc(formula2matrix(estimable),nfac,nruns,res3=res3)
#           else map <- mapcalc.distinct(formula2matrix(estimable),nfac,nruns,res3=res3,max.time=max.time, select.catlg=select.catlg, 
#                  perm.start=perm.start, perms=perms)
#        test <- map2design(map)
#        hilf <- catlg[names(map)][[1]]
#        hilf <- alias3fi(hilf$nfac-length(hilf$gen),hilf$gen)      
#        list(map=map,design=test,aliased=list(main=sort(chartr(paste(Letters[map[[1]]],collapse=""),paste(Letters[1:nfac],collapse=""),hilf$main)),
#                                         fi2=sort(chartr(paste(Letters[map[[1]]],collapse=""),paste(Letters[1:nfac],collapse=""),hilf$fi2))))
#    }

#estimable.character <- function(estimable, nfac, nruns,  clear=FALSE, res3=FALSE, max.time=60, select.catlg=catlg, perm.start=1:nfac, perms=NULL){
#        if (clear) map <- mapcalc(char2num(estimable),nfac,nruns,res3=res3)
#           else map <- mapcalc.distinct(char2num(estimable),nfac,nruns,res3=res3,max.time=max.time, select.catlg=select.catlg, 
#                    perm.start=perm.start, perms=perms)
#        test <- map2design(map)
#        hilf <- catlg[names(map)][[1]]
#        hilf <- alias3fi(hilf$nfac-length(hilf$gen),hilf$gen)      
#        list(map=map,design=test,aliased=list(main=sort(chartr(paste(Letters[map[[1]]],collapse=""),paste(Letters[1:nfac],collapse=""),hilf$main)),
#                                         fi2=sort(chartr(paste(Letters[map[[1]]],collapse=""),paste(Letters[1:nfac],collapse=""),hilf$fi2))))
#    }


#estimable(formula2matrix(formula("~(A+B+C+D)^2+(E+F+G+H)+(A+B+C+D):(E+F+G+H)")), 10, 32)
#estimable(formula2matrix(formula("~(A+B+C+D)^2+(E+F+G+H)+(A+B+C+D):(E+F+G+H)")), 10, 64)
#estimable(formula2matrix(formula("~(A+B+C+D+E+F)+A:B+B:C+C:D+C:F+(D+E+F)^2")), 6, 16,res3=TRUE)


