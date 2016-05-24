# Automatically generated from all.nw using noweb


#$Log: pedTrim.q,v $
#Revision 1.4  2009/11/19 15:00:31  sinnwell
#*** empty log message ***
#
#Revision 1.3  2009/11/19 14:57:05  sinnwell
#*** empty log message ***
#
#Revision 1.2  2009/11/17 23:11:09  sinnwell
#change for ped object
#
#Revision 1.1  2008/07/16 20:23:29  sinnwell
#Initial revision
#


findUnavailable <-function(ped, avail) {

  ## find id within pedigree anyone who is not available and
  ## does not have an available descendant
  
  ## avail = TRUE/1 if available, FALSE/0 if not
  
  ## will do this iteratively by successively removing unavailable
  ## terminal nodes
  ## Steve Iturria, PhD, modified by Dan Schaid
  
  cont <- TRUE                  # flag for whether to keep iterating
  
  is.terminal <- (is.parent(ped$id, ped$findex, ped$mindex) == FALSE)
  ## JPS 3/10/14 add strings check in case of char ids
  pedData <- data.frame(id=ped$id, father=ped$findex, mother=ped$mindex,
                        sex=ped$sex, avail, is.terminal, stringsAsFactors=FALSE)  
  iter <- 1

  while(cont)  {
    ##print(paste("Working on iter", iter))
    
    num.found <- 0
    idx.to.remove <- NULL
    
    for(i in 1:nrow(pedData))
      {
        
        if(pedData$is.terminal[i])
          {
            if( pedData$avail[i] == FALSE )   # if not genotyped         
              {
                idx.to.remove <- c(idx.to.remove, i)
                num.found <- num.found + 1
                
                ## print(paste("  removing", num.found, "of", nrow(pedData)))
              }
          }
       
      }

    if(num.found > 0) {

      pedData <- pedData[-idx.to.remove, ]
      ## re-index parents, which varies depending on if the removed indx is
      ## prior to parent index
      for(k in 1:nrow(pedData)){
        if(pedData$father[k] > 0) {
          pedData$father[k] <- pedData$father[k] -
              sum(idx.to.remove < pedData$father[k])
        }
        if(pedData$mother[k]+0) {
          pedData$mother[k] <- pedData$mother[k] -
              sum(idx.to.remove < pedData$mother[k])
        }
      }
      pedData$is.terminal <-
        (is.parent(pedData$id, pedData$father, pedData$mother) == FALSE)
      
    }
    else {
      cont <- FALSE
    }
    iter <- iter + 1   
    
  }
  
  ## A few more clean up steps

  ## remove unavailable founders
  tmpPed <- excludeUnavailFounders(pedData$id, 
                        pedData$father, pedData$mother, pedData$avail)

  ## 
  tmpPed <- excludeStrayMarryin(tmpPed$id, tmpPed$father, tmpPed$mother)

  
  id.remove <- ped$id[is.na(match(ped$id, tmpPed$id))]

  return(id.remove)
  
}


excludeStrayMarryin <- function(id, father, mother){
  # get rid of founders who are not parents (stray available marryins
  # who are isolated after trimming their unavailable offspring)
  ## JPS 3/10/14 add strings check in case of char ids
  trio <- data.frame(id=id, father=father, mother=mother, stringsAsFactors=FALSE)
  parent <- is.parent(id, father, mother)
  founder <- is.founder(father, mother)

  exclude <- !parent & founder
  trio <- trio[!exclude,,drop=FALSE]
  return(trio)

}

excludeUnavailFounders <- function(id, father, mother, avail)
  {
    nOriginal <- length(id)
    idOriginal <- id   
    zed <- father!=0 & mother !=0
    ## concat ids to represent marriages. 
    ## Bug if there is ":" in char subj ids
    marriage <- paste(id[father[zed]], id[mother[zed]], sep=":" )

    sibship <- tapply(marriage, marriage, length)
    nm <- names(sibship)

    splitPos <- regexpr(":",nm)
    dad <- substring(nm, 1, splitPos-1)
    mom <- substring(nm, splitPos+1,  nchar(nm))
    
    ##  Want to look at parents with only one child.
    ##  Look for parents with > 1 marriage.  If any
    ##  marriage has > 1 child then skip this mom/dad pair.
    
    nmarr.dad <- table(dad)
    nmarr.mom <- table(mom)
    skip <- NULL
    
    if(any(nmarr.dad > 1)) {
      ## Dads in >1 marriage
      ckdad <- which(as.logical(match(dad,
                      names(nmarr.dad)[which(nmarr.dad > 1)],nomatch=FALSE)))
      skip <- unique(c(skip, ckdad))
    }
    
    if(any(nmarr.mom > 1)) {
      ## Moms in >1 marriage
      ckmom <- which(as.logical(match(mom,
                      names(nmarr.mom)[which(nmarr.mom > 1)],nomatch=FALSE)))
      skip <- unique(c(skip, ckmom))
    }
      
    if(length(skip) > 0) {
      dad <- dad[-skip]
      mom <- mom[-skip]
      zed <- (sibship[-skip]==1) 
    } else {
      zed <- (sibship==1)
    }

    
    n <- sum(zed)
    idTrimmed <- NULL
    if(n>0)
      {
        
        # dad and mom are the parents of sibships of size 1
        dad <- dad[zed]
        mom <- mom[zed]
        for(i in 1:n){
          ## check if mom and dad are founders (where their parents = 0)
          dad.founder <- (father[id==dad[i]] == 0) & (mother[id==dad[i]] == 0)
          mom.founder <- (father[id==mom[i]] == 0) & (mother[id==mom[i]] == 0)
          both.founder <- dad.founder & mom.founder

          ## check if mom and dad have avail
          dad.avail <- avail[id==dad[i]]
          mom.avail <- avail[id==mom[i]]

          ## define not.avail = T if both mom & dad not avail
          not.avail <- (dad.avail==FALSE & mom.avail==FALSE)
        
          if(both.founder & not.avail)   {
              ## remove mom and dad from ped, and zero-out parent 
              ## ids of their child
                        
            child <- which(father==which(id==dad[i]))          
            father[child] <- 0
            mother[child] <- 0
            
            idTrimmed <- c(idTrimmed, dad[i], mom[i])
            
            excludeParents <- (id!=dad[i]) & (id!=mom[i])
            id <- id[excludeParents]
            father <- father[excludeParents]
            mother <- mother[excludeParents]

            ## re-index father and mother, assume len(excludeParents)==2
            father <- father - 1*(father > which(!excludeParents)[1]) -
              1*(father > which(!excludeParents)[2])
            
            mother <- mother - 1*(mother > which(!excludeParents)[1]) -
              1*(mother > which(!excludeParents)[2])

            avail <- avail[excludeParents]
          } 
        }
      }
    
    nFinal <- length(id)
    nTrimmed = nOriginal - nFinal 
  
    
    return(list(nTrimmed = nTrimmed, idTrimmed=idTrimmed,
                id=id, father=father, mother=mother))
  }



