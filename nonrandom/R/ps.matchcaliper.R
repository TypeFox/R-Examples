ps.matchcaliper <- function(vect1,
                            vect2,
                            ratio           = 1,
                            caliper         = "logit",
                            x               = 0.2,
                            givenTmatchingC = TRUE,
                            bestmatch.first = TRUE,
                            setseed         = FALSE)
{ 
  ## #################
  ## Set random number
  if(!identical(setseed, FALSE))  
    set.seed(seed=setseed)  

  
  ## ####################################################

  ## Treated -> untreated (default: givenTmatchingC=TRUE)
  ## Untreated -> treated (givenTmatchingC=FALSE)
  
  if(!givenTmatchingC){ ## swap vect1 (treated) and vect2 (untreated)

    temp  <- vect1  
    vect1 <- vect2  ## in general: vect1 includes treated; vect2
                    ## includes untreated
    vect2 <- temp   ## length(vect2) > length(vect1)

    message(paste("Argument 'givenTmatchingC'=",
                  givenTmatchingC,
                  ": Treated elements were matched to each untreated element.",
                  sep=""))
  }

  ## ##############
  ## Define caliper
  if (!is.numeric(caliper)){
    if (is.character(caliper)){
      if (!is.na(charmatch(caliper, "logit"))){        
        caliper <- x * sd(log(c(vect1,vect2)/(1-c(vect1,vect2))))
      }else{
        stop(paste("Argument 'caliper'= ",caliper," is not permitted.",sep=""))
      }
    }else{
      stop("Argument 'caliper' must be numeric or a string.")
    }
  }

  ## ##############################
  ## Matrix of absolute differences
  D <- abs(outer(vect1, vect2, FUN="-"))
  
  ## Check, whether there are enough controls for each element in vect1=nrow(D)
  list.candidates <- integer(nrow(D))

  ## How many potential matching partner does each treated obs (vect1) have?
  for(i in 1:nrow(D)){
    list.candidates[i] <- length(D[i, D[i,]<=caliper])
  }

  ## check whether obs have not enough potential matching partners?
  if(any((list.candidates-ratio) < 0))
    message(paste("It is not possible to find",
                  ratio , "matching elemtent(s) for",
                  sum((list.candidates-ratio) < 0),
                  "element(s);\n Check the caliper size",
                  round(caliper,3), "or the ratio", ratio))

  
  mv1 <- matrix(nrow=nrow(D),ncol=ratio)  
  lm  <- integer(nrow(D))           ## vector of length nrow(D) = number treated obs.
                                          
  M <- D  
  M[D > caliper] <- NA              ## Set NA to controls which are not matching candidate
                                    ## (see list.candidate)

  names(list.candidates) <- c(1:nrow(M))
  lm <- sort(list.candidates)       ## obs which the smallest number of match.candidates are 
                                    ## the first in matching procedure

  
  ## ##################
  ## Matching algorithm
  for(j in 1:nrow(M)){            
    
    i <- as.integer(names(lm[j])) ## start with the obs with fewest
                                  ## matching partners
    
    l <- rep(NA,times=ncol(M))    ## length(l) = number of untreated
                                  ## obs; maximal number of matching
                                  ## partners

    for(k in 1:ncol(M)){
      if(is.na(M[i,k])==FALSE){     ## obs not potential match partner
                                    ## are set to NA
        
        l[k] <- k                   
      }
    }
    
    l <- l[is.na(l)==FALSE]         ## only potential matching elements included

    if(length(l)>ratio){            ## More than ratio potential matching elements? 

      if(bestmatch.first){          ## TRUE: Take the |ratio| smallest 

        c        <- M[i,]                         ## Difference of treated/control for obs.i
        names(c) <- c(1:ncol(M))
        c        <- sort(c)                       ## Sort differences without NAs
        l        <- as.integer(names(c[1:ratio])) ## Take the first |ratio| of
                                                  ## possible match partners
      }else{
        ## FALSE: Take a sample of size |ratio| 
        l <- sort(sample(c(l),ratio))
      }
      
      mv1[i,]=l  ## Save the matched element(s) l to treated element j 

    }else{ ## length(l) =< ratio

      mv1[i,]=c(l,rep(0,times=ratio-length(l))) ## Set 0 for missing matching
                                                ## elements 
    }

    if(length(l)!= 0){     ## If there are matching elements found, ...
      M[i,-l] <- NA        ## Set the rest of row i to NA except for 
                           ## 'l' row(s) including match partners
      for(k in 1:nrow(M)){        
        if(k!=i){          ## The k treated obs. can not found matching 
          M[k,l] <- NA     ## partners 'l' because there were matched  
        }                  ## already to treated element i ==> set to NA
      }
    } ## length(l)!= 0
    
  } ## end loop over all treated (for (j))


  mv2 <- integer(ncol(M)) ## List of matched elements for each control obs.
  pv1 <- integer(nrow(M)) ## List of matching partners w.r.t treated obs
  pv2 <- integer(ncol(M)) ## List of matching partners w.r.t control obs
  
 
  for(j in 1:nrow(M)){            
    mv2[mv1[j, mv1[j, ] != 0]] <- j  ## non-matched obs are zero
  }

  names(mv2) <- 1:ncol(M)
  
 
  if (any(mv1 == 0))  
    message("Some elements have no matching element.")
  
  if( length(mv2[mv2[] != 0]) != nrow(M) * ratio )
    message(paste("Some elements have no",
                  ratio,
                  "matching element(s) as desired."))

  ## mv1: indices for matched untreated obs 
  ## mv2: index for matched treated obs


  ## define vectors of matching tuples (or 'pairs', if 'ratio' = 1)
  k <- 1               
  for(j in 1:nrow(M)){   

    if(any(mv1[j,]!=0)){
      pv1[j] <- k
      pv2[mv1[j,mv1[j,]!=0]] <- k
      k <- k+1      
    }
    
    if(is.na(vect1[j])){          ## reinsert the 'NA's 
      mv1[j,] <- NA               
      pv1[j]  <- NA
    }
  }

  for(j in c(1:ncol(M)))          ## reinsert the 'NA's 
    if(is.na(vect2[j])){
      mv2[j] <- NA
      pv2[j] <- NA
    }

  
 
  if(!givenTmatchingC){
    temp <- pv1
    pv1  <- pv2
    pv2  <- temp
    temp <- mv1
    mv1  <- mv2
    mv1  <- temp
  }

  matchedvectors <- list(caliper         = caliper,
                         ratio           = ratio,
                         pairvect1       = pv1,
                         pairvect2       = pv2,
                         matchvect1      = mv1,
                         matchvect2      = mv2,
                         givenTmatchingC = givenTmatchingC,
                         bestmatch.first = bestmatch.first)
    
  return(matchedvectors)
  
}
  
