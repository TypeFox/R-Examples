gd.smouse <- function(population, verbose=TRUE){
  # check to see if the passed data is of the right type
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }
  
  # get initial estimates of number of individuals and the maximnum number of loci
  initmaxind<-(dim(population@tab))[1]
  maxloci<-length(locNames(population))
  
  # create a local matrix with the genomes so we can make changes if necessary and not affect other parts of script
  genomes<-population@tab

  # determine the initial population details
  indpop<-population@pop[1:initmaxind,drop=TRUE] # gets a list of the individual's population
  initnumpops<-length(attributes(indpop)$levels) # number of populations
  poplist<-attributes(indpop)$levels

  # create an empty list which will contain all pops 
  indsrmv<-c()
  popsgone<-c()
  
  # need to find the intervals for each loci...
  cs <- cumsum(population@loc.n.all) # this determines the end of each loci frame for each ind's genotype
  cs.l <- c(1,cs[-length(cs)]+1) # this determine the beginning of each loci frame for each ind's genotype

  # check to see that the populations have more than one individual (probably need to make it more 
  # sophisticated later on e.g. check to see that at least one individual has loci values)
  for (i in 1:initnumpops){
    indcnt<-length(indpop[which(indpop==poplist[i])])
    if (indcnt<2){
      indrmv<-which(indpop==poplist[i])
      indsrmv<-c(indsrmv, indrmv)
      popsgone<-c(popsgone,i)
    }
    if(indcnt>=2){
      popi<-which(indpop==poplist[i])
      genos<-population@tab[popi,]
      cntlow<-0
      for (j in 1:maxloci){
        totalleles<-sum(genos[,cs.l[j]:cs[j]],na.rm=TRUE)
        if(totalleles<1) cntlow<-cntlow+1
      }
      if (cntlow>0) {
        indsrmv<-c(indsrmv, popi)
        popsgone<-c(popsgone, i)
      }
    }
  }

  if(length(popsgone)>0){
    message("WARNING!!! The following populations were dropped due to insufficient numbers of individuals")
    for (i in 1:length(popsgone)){
        message(poplist[popsgone[[i]]])
    }
  }

  # remove individuals if the population only has one individual
  if(length(indsrmv)>0){
    indsrmv<-indsrmv*(-1)
    indpop<-factor(indpop[indsrmv])
    genomes<-genomes[indsrmv,]
    poplist<-attributes(indpop)$levels
  }

  # what is the total number of individuals
  maxind<-(dim(genomes))[1]

  # calculate the number of pops after cleaning list
  numpops<-length(attributes(indpop)$levels) # number of populations

  # matrix to put smoused distances in
  smoused<-matrix(NA,nrow=maxind,ncol=maxind)



  for (i in 1:numpops){   #numpops this is looping over pops
    for (j in i:numpops){ #numpops this is looping over pops
      if(verbose) cat("\r","Comparing population ",levels(population@pop)[i]," with population ",levels(population@pop)[j])
      pop1<-which(indpop==poplist[i]) # gets a list of individuals from pop1
      pop2<-which(indpop==poplist[j]) # gets a list of individuals from pop2
      smoused.loci<-array(NA,c(length(pop2),length(pop1),maxloci))
    
      # this if statement calculates genetic distances between individuals in same population, only half the 
      # calculations need to be done as the other half are repeats (e.g. calculate distance between a and b and
      # then calculate the distance between b and a)
      if (i==j){
        if (length(pop1)>1 && length(pop2)>1){
          for (l in 1:(length(pop1)-1)){    # this is looping over columns, skip last column because no comparison to be made
            for (k in (l+1):length(pop2)){  # this is looping over rows, skip l=k and assign value later to save calcs 
              # extract the genotypes of the two individuals
              i1 <- genomes[pop1[l],]
              i2 <- genomes[pop2[k],]          
              # calculate genetic distances between individuals using Smouse and Peakall 1999
              pairdist<-(i1-i2)^2
              sdists<-sapply(1:maxloci, function(x,cs,cs.l,pairdist) sum(pairdist[cs.l[x]:cs[x]])*0.5  , cs, cs.l,pairdist)
              for (m in 1:maxloci){
                if(is.na(sdists[m])) smoused.loci[k,l,m]<-(-99) else smoused.loci[k,l,m]<-sdists[m]    
              }    
            }
          }
        }
      } 
      # this if statement calculates genetic distances between individuals in different populations. All 
      # calculations have to be done here because there aren't the repeated calculations in the above if
      # statement (e.g. no a = b and then b = a)
      if (i!=j){
        for (k in 1:length(pop1)){
          for (l in 1:length(pop2)){
            # extract the genotypes of the two individuals
            i1 <- genomes[pop1[k],]
            i2 <- genomes[pop2[l],]          
            # calculate genetic distances between individuals using Smouse and Peakall 1999
            pairdist<-(i1-i2)^2
            sdists<-sapply(1:maxloci, function(x,cs,cs.l,pairdist) sum(pairdist[cs.l[x]:cs[x]])*0.5  , cs, cs.l,pairdist)
            for (m in 1:maxloci){
              if(is.na(sdists[m])) smoused.loci[l,k,m]<-(-99) else smoused.loci[l,k,m]<-sdists[m]    
            }    
          }
        }
      }
      # this step ID's all array elements for which we have a distance and allows us to mask other elements 
      # e.g. (the upper triangle elements) so that we can solve missing values
      check<-smoused.loci[,,,drop=FALSE]>=0  
      # now we come up with a value for missing genetic distances between individuals in same population
      if (i==j){
        if (length(pop1)>1 && length(pop2)>1){
          for(k in 1:maxloci){ 
            set<-as.matrix(smoused.loci[,,k]) # get the matrix for a particular loci
            set[upper.tri(set)]<-NA           # block the upper tri angle
            frame<-check[,,k,drop=FALSE]
            locimean<-sum(set[frame],na.rm=TRUE)/sum(frame[lower.tri(frame)],na.rm=TRUE)  
            for(m in 1:(length(pop1)-1)){ # loop over columns
              for (l in (m+1):length(pop2)){ # loop over rows
                if (smoused.loci[l,m,k]<0) smoused.loci[l,m,k]=locimean
              }
            }
          }
        }
        for (k in 1:length(pop1)){
          for (l in k:length(pop2)){
            if(k==l) smoused[pop2[l],pop1[k]]<-0
            if(k!=l){
              smoused[pop2[l],pop1[k]]<-sum(smoused.loci[l,k,])
              smoused[pop1[k],pop2[l]]<-sum(smoused.loci[l,k,])
            }
          }
        }    
      }
      # come up with a value for missing genetic distances between individuals in different populations
      if(i!=j){
        for (k in 1:maxloci){
          set<-as.matrix(smoused.loci[,,k])
          subset<-set[1:length(pop2),1:length(pop1)]
          frame<-as.matrix(check[,,k])
          locimean<-sum(subset[frame[,]],na.rm=TRUE)/sum(frame,na.rm=TRUE)      
          for(l in 1:length(pop2)){
            for (m in 1:length(pop1)){
              if (smoused.loci[l,m,k]<0) smoused.loci[l,m,k]=locimean
            }
          }  
        }
        for (k in 1:length(pop2)){
          for (l in 1:length(pop1)){
            smoused[pop2[k],pop1[l]]<-sum(smoused.loci[k,l,])
            smoused[pop1[l],pop2[k]]<-sum(smoused.loci[k,l,])
          }
        }
      }
    }
  }

  if(length(indsrmv)>0){
    # put names on the rows and columns on d.fast (!!!!make sure to use only non-removed individuals)
    colnames(smoused)<-indNames(population)[indsrmv]
    rownames(smoused)<-indNames(population)[indsrmv]
    # calculate geographical distance
    #geodist<-as.matrix(dist(population@other$utm[indsrmv,]))/1000
  } else if(is.null(indsrmv)){
    colnames(smoused)<-indNames(population)
    rownames(smoused)<-indNames(population)
  }

  # force upper triangle to NA
  smoused[upper.tri(smoused,diag=FALSE)]<-NA
  smoused<-as.dist(smoused)
  return(smoused)
}                                                         