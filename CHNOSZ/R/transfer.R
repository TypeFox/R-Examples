# CHNOSZ/transfer.R
# plot reaction paths on activity diagrams
# 20090325 jmd

transfer <- function(nsteps=500,dmode='coupled',devmax=0.1,
  plot=NULL,ibalance=1,fmode='one',buffers=NULL,alphamax=-2,
  alphastart=-10,T=25,P='Psat',do.title=TRUE,beta=0) {
  # calculate mass transfer among species of interest,
  # conserving the number of moles of balanced thing
  # (i.e., the "conservant"; e.g., Al+3 or PBB)
  # 1 - dissolve all the species (to log0)
  # 2 - precipitate the species that is (meta)stable

  ## dmode (destruction mode)
  # none - don't react seconday species (implies open system)
  # coupled - destruction of old secondary reactant is coupled to formation
  ## fmode (formation mode)
  # one - only form one most stable species
  # many - form various species, relative velocities from reaction affinities

  # keep track of whether each step was successful
  # and mode of destruction and most stable specie
  didwork <- logical()
  dmodes <- character()
  istables <- character()
  myaffs <- list()
  sout <- NULL
  # logarithm of a very small number (approaching zero)
  log0 <- -999
  # logarithm of molality above which something
  # is considered possible for reaction (is present)
  logpresent <- -50
  # the starting (basis0) and current (basis)
  # species and basis conditions
  thermo <- get("thermo")
  basis <- basis0 <- thermo$basis
  species <- species0 <- thermo$species
  dmode0 <- dmode
  buffargs <- buffstuff <- NULL
  # what basis species can tolerate coupling
  #if(missing(icouple)) if(!is.null(plot)) icouple <- plot 
  #if(is.null(icouple)) icouple <- 1:nrow(basis)
  if(!is.null(plot)) icouple <- plot 
  else icouple <- 1:nrow(basis)[-ibalance]

  if(ibalance=='PBB') {
    # balance destruction and formation reactions 
    # on protein backbone group  jmd 20080510
    isprotein <- grep('_',species$name)
    if(length(isprotein) != length(species$name))
      stop('transfer: for balance=PBB, all species must be proteins')
    ibalance <- nrow(basis) + 1
    # add row to basis dataframe
    basis <- rbind(basis,basis[nrow(basis),])
    rownames(basis)[ibalance] <- 'PBB'
    # the stoichiometric matrix becomes not square; so be careful
    basis[ibalance,1:(ibalance-1)] <- 0
    basis$ispecies[ibalance] <- 9999
    basis$logact[ibalance] <- 0
    bs <- as.character(basis$state)
    bs[ibalance] <- 'gas'
    basis$state <- bs
    # add column to species dataframe
    # (insert new column at end of stoic matrix)
    species <- cbind(species[,1:(ibalance-1)],species[,1],
      species[ibalance:ncol(species0)])
    colnames(species)[ibalance] <- 'PBB'
    species[,ibalance] <- protein.length(species$name)
    # done!
    cat('transfer: balancing reactions on PBB\n')
  }

  # howto balance formation with destruction: number of moles of the conservant
  # otherwise we take ibalance supplied by the user
  if(missing(ibalance)) ibalance <- which.balance(species)[1]
  if(is.na(ibalance)) stop('transfer: no conservant specified or found')

  # exponent of destruction
  alpha <- alphastart
  alphas <- numeric()
  # adjustment to exponent of buffer reactions
  # the current basis species and number of
  # candidates for coupled mode 
  ncb <- 0
  icb <- 1
  # we starting with nothing formed
  species$logact <- log0
  buffargs <- buffstuff <- NULL

  # are there buffers? deal with them as follows
  molbasisbuffer <- function(basis,buffers,buffargs=NULL,buffstuff=NULL) {
    # returns the logarithms of activities of basis
    # species, modified by buffer reactions if present
    basis2 <- basis1 <- basis
    if(!is.null(buffers)) {
      for(i in 1:length(buffers$basis)) {
        ibuf <- match(buffers$basis[i],rownames(basis)) 
        basis1$logact[ibuf] <- buffers$buffer[i]
      }
      ibuf <- which(!can.be.numeric(basis1$logact))
      # get the buffered activities the slow way (first)
      # or by a quicker approach
      if(is.null(buffargs)) {
        # the slow way, short version
        # we have to use basis species w/o the PBB
        tempbasis <- basis1[rownames(basis1)!="PBB",]
        thermo$basis <- tempbasis
        assign("thermo", thermo, "CHNOSZ")
        # the slow way, long-winded version
        # (recreating affinity's call to buffer so
        # we can store the intermediate results)
        logact.basis <- as.list(basis$logact)
        # temporarily activate the buffer
        is.buffer <- buffer(logK=NULL)
        # and save corresponding basis and species definitions
        buffstuff <- list(bufbasis=thermo$basis,bufspecies=thermo$species,is.buffer=is.buffer)
        # the logKs of the buffer species
        logK <- affinity(property='logK',T=T,P=P)$values
        # the indices of buffer species
        # would need to add balance argument here for PBB
        bresult <- list()
        myib <- numeric()
        for(i in 1:length(ibuf)) {
          ib <- is.buffer[[i]]
          myib <- c(myib,ib)
          ibasis <- ibuf[i]
          buffargs <- list(logK=logK,ibasis=ibasis,logact.basis=logact.basis,is.buffer=ib)
          bresult <- do.call("buffer",buffargs)$logact.basis[ibuf[i]]
          if(i==1) br <- bresult else br <- c(br,bresult)
        }
        bresult <- br
        buffargs$is.buffer <- is.buffer
        # because that buffer(logK=NULL) command adds to the species, restore them
        species(myib[!myib %in% 1:nrow(species)],delete=TRUE)
      } else {
        # the 'quick' way
        #ibuf <- buffargs$ibasis
        oldbasis <- thermo$basis
        oldspecies <- thermo$species
        thermo$basis <- buffstuff$bufbasis
        thermo$species <- buffstuff$bufspecies
        assign("thermo", thermo, "CHNOSZ")
        is.buffer <- buffargs$is.buffer
        for(i in 1:length(ibuf)) {
          ib <- is.buffer[[i]]
          ibasis <- ibuf[i]
          buffargs <- list(logK=buffargs$logK,ibasis=ibasis,
            logact.basis=buffargs$logact.basis,is.buffer=ib)
          bresult <- do.call("buffer",buffargs)$logact.basis[ibuf[i]]
          if(i==1) br <- bresult else br <- c(br,bresult)
        }
        bresult <- br
        thermo$basis <- oldbasis
        thermo$species <- oldspecies
        assign("thermo", thermo, "CHNOSZ")
      }
      for(i in 1:length(ibuf)) {
        # reference to the moles of species 
        # in the previous step
        molbprev <- 10^basis$logact[ibuf[i]]
        # how much buffered species to transfer at each step
        # refer to the staring point of the simulation?
        #molb0 <- 10^basis0$logact[ibuf[i]]
        # or to the last step?
        molb0 <- molbprev
        # what we get if the buffer goes to completion (beta = 0)
        basis2$logact[ibuf[i]] <- as.numeric(bresult[i])
        b2 <- basis2$logact[ibuf[i]]
        molb2 <- 10^b2
        # the change we're looking at
        mybeta <- alpha + beta
        moldb <- (molb2 - molb0) * 10^(mybeta)
        molb1 <- molbprev + moldb
        # disallow negative amounts of basis species
        if(molb1 < 0) moldb <- molb2 - molbprev
        ii <- match(rownames(basis)[ibuf[i]],buffers$basis)
        cat('transfer: buffer (',mybeta,buffers$buffer[ii],') adds',
          moldb,'moles of',rownames(basis)[ibuf[i]],'\n')
        basis1$logact[ibuf[i]] <- log10(molbprev + moldb)
        basis1$logact[is.infinite(basis1$logact)] <- log0
      }
    }
    # stuff for PBB (not-square stoichiometric matrix)
    if('PBB' %in% rownames(basis1)) dd <- 2 else dd <- 3
    return(list(basis=basis1[1:nrow(basis1),1:(nrow(basis1)+dd)],buffargs=buffargs,buffstuff=buffstuff))
  }

  # major loop
  for(j in 1:nsteps) {

    # 000 - plot title (number of each step) 20090408
    if(do.title) {
      # strikeout the previous number
      title(main=as.character(j-1),col.main="white")
      # plot the current number
      title(main=as.character(j))
    }
     
    # 00 - calculate the exponent of 
    # destruction (alpha) for this cycle
    # basic idea: if last step failed, decrease alpha
    # otherwise increase it (with some exceptions)
    dalpha <- 0
    if(j>1) if(didwork[j-1]) {
      dalpha <- dalpha + 1
    } else dalpha <- dalpha - 1

    # logic for coupled mode
    # (coupled - all)
    if(dmode0=='coupled') {
      if(j>1) {
        # restore all to coupled if last step worked
        if(didwork[length(didwork)]) {
          # cycle through the basis species - 
          # even if they are working
          if(icb < ncb) icb <- icb + 1
          else icb <- 1
          dmode <- dmode0
        }
        else {
          if(dmode=='coupled') {
            # if it didn't work, go to all mode
            dmode <- "all"
            # none of this affects the exponent of destruction
            dalpha <- 0
          } else if(dmode=='all') {
            dmode <- 'coupled'
            dalpha <- -1
            # going back to coupled mode, let us reset
            # the counter for coupled basis species
            icb <- 1
          }
        }
      }  
    }
    # logic for all mode
    # (all - only)
    if(dmode0=='all') {
      if(j>1) {
        if(didwork[length(didwork)]) {
          if(dmode!='only') dmode <- dmode0
          else if(alpha==alphamax) dmode <- dmode0
        }
        else {
          if(alpha < logpresent) dmode <- 'only'
          else if(dmode=='only') if(j>2) 
            if(didwork[j-2] & dmodes[j-2]=='only') dmode <- dmode0
        }
      }
    }
    # update alpha
    alpha <- alpha + dalpha
    # this is kind of arbitrary (so we let the user set it)
    if(alpha > alphamax) alpha <- alphamax
    cat(paste('--- step',j,'---\n'))
    cat(paste('transfer: destruction exponent is',alpha,'\n'))
    # and record the destruction coefficient
    alphas <- c(alphas,alpha)
    dmodes <- c(dmodes,dmode)
    if(dmode0 %in% c('coupled','all')) 
      dmtext <- paste('(',dmode,')') else dmtext <- ''
    cat(paste('transfer: destruction is',dmode0,dmtext,'\n'))
    cat(paste('transfer: formation is ',fmode,' species, balanced on ',rownames(basis)[ibalance],'\n',sep=''))

    # 0 - if this is the first step or the previous step worked
    # calculate chemical activities of basis species
    # and affinities of formation reactions
    # do the buffer calculations starting on the second step,
    mybl <- basis$logact
    if(j>1) {
      # prevent the buffer from being applied if
      # the previous step worked and was coupled (to avoid backtracking
      # into a preceeding stability field, which can lead
      # to impossibility of an accessory conservant at next step)
      if(!(dmodes[length(dmodes)-1]=='coupled' & didwork[j-1])) {
        if(is.null(buffargs)) {
          mbb <- molbasisbuffer(basis=basis,buffers=buffers)
          buffargs <- mbb$buffargs
          buffstuff <- mbb$buffstuff
        } else {
          mbb <- molbasisbuffer(basis=basis,buffers=buffers,buffargs=buffargs,buffstuff=buffstuff)
        }
        # keep track of the changes imposed by the buffer
        # so they can be excluded from devmax checking
        logbuffdev <- as.numeric(mbb$basis$logact) - basis$logact
        mybl <- as.numeric(mbb$basis$logact)
      }
    } 
    # now use this number of moles of basis species
    molbasis0 <- 10^as.numeric(mybl)
    # get the affinities for the first step
    getaff <- function(mybl,sout=NULL) {
      # do it for unit activities of minerals (and proteins?)
      thermo$species$logact <- 0
      # prevent the PBB from getting in here
      thermo$basis$logact <- mybl[1:nrow(basis0)]
      assign("thermo", thermo, "CHNOSZ")
      if(is.null(sout)) {
        # on the first step only, grab the intermediate results
        # they are kept around for reasons of speed
        aff <- affinity(T=T,P=P)
        # store output of subcrt
        sout <- aff$sout
      } else aff <- affinity(sout=sout,T=T,P=P)
      # numerical values of the formation affinities of the species
      myaff <- as.numeric(aff$values)/as.numeric(species[,ibalance])
      return(list(myaff=myaff,sout=sout))
    }
    if(j==1) {
      aff <- getaff(mybl,sout)
      myaff <- aff$myaff
      sout <- aff$sout
    }

    # 1 - calculate the number of moles of primary reacting species
    # we are being fed from a pool of constant composition 
    sl1 <- species0$logact
    # but we only take so many moles of species
    alpha1 <- alpha
    # if only the secondary minerals react
    if(dmode=='only') sl1 <- alpha1 <- log0
    # the moles of primary reacting species
    molspecies1 <- 10^(sl1 + alpha1)
    # the number of moles of basis species from primary reactants
    ipresent <- which(log10(molspecies1) > logpresent)
    molbasis1 <- as.numeric(0 * species[1,1:nrow(basis)])
    if(length(ipresent)>0) for(i in 1:length(ipresent))
      molbasis1 <- molbasis1 + as.numeric(molspecies1[ipresent[i]] *
        species[ipresent[i],1:nrow(basis)])
    # the number of moles of conservant that we form
    molconservant1 <- molbasis1[ibalance]
    # stop when empty but not in only mode
    if(identical(molconservant1,0) & dmode!='only') {
      cat('transfer: stopping: no primary conservant (try more moles of reactant?)\n')
      didwork <- c(didwork,FALSE)
      return(invisible(list(basis=basis,species=species,
        alphas=alphas[didwork],dmodes=dmodes[didwork])))
    }

    # 2 - get the number of moles secondary reacting species
    # (i.e., those previously formed)
    sl2 <- species$logact
    molprevspecies <- 10^sl2
    if(j > 1) {
      molspecies2 <- molprevspecies  # i.e. dmode0='coupled'
      # if(dmode0=="all") molspecies2 <- molprevspecies * 10^alpha
    } else {
      molprevspecies <- molspecies2 <- rep(0,length(sl2))
    }
    # if an open system don't react the previously formed stuff
    if(dmode=='none') {
      molunusedspecies <- molspecies2
      cat(paste('transfer: open system lost moles of species\n'))
      molprevspecies <- molspecies2 <- rep(0,length(sl2))
    }
    # moles of basis species from secondary reactions
    ipresent <- which(log10(molspecies2) > logpresent)
    molbasis2 <- as.numeric(0 * species[1,1:nrow(basis)])
    if(length(ipresent)>0) for(i in 1:length(ipresent))
      molbasis2 <- molbasis2 + as.numeric(molspecies2[ipresent[i]] * 
        species[ipresent[i],1:nrow(basis)])
    # the number of moles of conservant that is formed
    molconservant2 <- molbasis2[ibalance]

    # 3 - calculate the formation of metastable products
    sl3 <- species$logact
    molspecies3 <- sl3 * 0
    # identify the single most stable species
    istable <- which.max(myaff)
    if(j>1) istableprev <- istable else istableprev <- NULL
    if(fmode=='one') {
      # the number of moles of species precipitated
      molspecies3[istable] <- species[istable,ibalance]
    } else { 
      # form various species in some proportion to their
      # chemical affinities (near-equilibrium rates)
      # 20090409 use abundance here -- relative
      # abundances of species in equilibrium
      molspecies3 <- 10^as.numeric(equil.boltzmann(myaff,rep(1,length(myaff)),0))
    }
    # the number of moles of basis species used
    ipresent <- which(log10(molspecies3) > logpresent)
    molbasis3 <- as.numeric(0 * species[1,1:nrow(basis)])
    if(length(ipresent) > 0) for(i in 1:length(ipresent))
      molbasis3 <- molbasis3 + as.numeric(molspecies3[ipresent[i]] * 
        species[ipresent[i],1:nrow(basis)])
    # the number of moles of conservant
    molconservant3 <- molbasis3[ibalance]
    # the conservation ratio
    if(molconservant3 != 0) {
      r <- (molconservant1 + molconservant2) / molconservant3
      # the final formation parameters
      molspecies3 <- r * molspecies3
      molbasis3 <- r * molbasis3
    } else {
      cat('transfer: unbalanced product formation (continuing anyway)\n')
    }
    # a function to aid in an informative exit
    # from the current step
    nextfun <- function() {
      cat(paste('transfer: failed step ( would form',species$name[istable],')\n'))
    }
    # 'only' mode fails if different products were made
    # in the last step (needs work for fmode='many')
    if(j>1) if(dmode=='only') {
      if(didwork[j-1] & dmodes[j-1]=='only') {
        if(!identical(istable,istableprev)) {
          cat('transfer: only mode stops: different products formed\n')
          # it was the step before the previous
          # that caused the change, so we set up to redo it
          # in the original mode.
          dmode <- dmode0
          didwork <- c(didwork,FALSE)
          # use the last known good values of activities
          # (undo the last step!)
          mysl <- species[,ncol(species)-1]
          species$logact <- species[,ncol(species)] <- mysl
          mybl <- basis[,ncol(basis)-1]
          basis$logact <- basis[,ncol(basis)] <- mybl
          # make the last increment very small
          alphas[length(alphas)-1] <- log0
          didwork <- c(didwork,FALSE)
          nextfun()
          next
        }
      }
    }

    # coupled destruction mode
    cmode <- "one"
    if(dmode=='coupled') {
      # A is our first conservant e.g. Al+3 or PBB
      # primary = formation1 (A1,B)
      # secondary = formation2 (A2,B)
      # which are written for two conservants
      # and give the reaction path along field boundaries
      # fraction of conservant A in primary reaction
      r1 <- molconservant1 / (molconservant1 + molconservant2)
      # moles of basis species in primary reaction
      mymolbasis1 <- molbasis1 - molbasis3 * r1
      # moles of basis species in secondary reaction (replacement)
      mymolbasis2 <- molbasis2 - molbasis3 * (1-r1)
      # to compute the coefficients in the coupled reactions
      couplefun <- function(r2,molspecies2,molbasis2,molspecies3,molbasis3,r1) {
        # scale the reaction parmaters by the ratios
        molspecies2 <- -r2 * molspecies2
        molbasis2 <- -r2 * molbasis2
        ms3 <- molspecies3
        mb3 <- molbasis3
        molspecies3 <- ms3 * r1 - sign(r2) * ms3 * (1-r1) * (abs(r2))
        molbasis3 <- mb3 * r1 - sign(r2) * mb3 * (1-r1) * (abs(r2))
        return(list(molspecies2=molspecies2,molbasis2=molbasis2,
          molspecies3=molspecies3,molbasis3=molbasis3))
      }
      # test if we can replace the secondary mineral
      if(!any(mymolbasis2 != 0 & mymolbasis1 != 0)) {
        cat('transfer: coupling: nothing to react\n')
        didwork <- c(didwork,FALSE)
        nextfun()
        next
      } else {
        # to couple the primary and secondary reaction, find a basis species
        # (not primary conservant) with opposite coefficients in the reactions
        iworked <- FALSE
        # calculate the coupling
        # that corresponds to one of the candidates
        # for conservants (usu. the plotting variables)
        # or a mixture of two (diagonal lines)
        # here the conservants don't have to have opposite coefficients
        iworked <- FALSE
        iconservantB <- which(mymolbasis1 != 0 & mymolbasis2 != 0 & 
          1:nrow(basis) != ibalance &
          1:nrow(basis) %in% icouple[1:2])
        onencb <- length(iconservantB)
        if(onencb < 2) {
          #cat('transfer: coupling: cmode two --> one (only one conservant found)\n')
        } else {
          cmode <- "two"
          cat('transfer: found two secondary conservants\n')
          iB1 <- iB2 <- 0
          iB1 <- iconservantB[1]
          # ratio of conservant B1 in primary vs secondary reaction
          rB1 <- mymolbasis1[iB1] / mymolbasis2[iB1]
          myout.1 <- couplefun(rB1,molspecies2,molbasis2,molspecies3,molbasis3,r1)
          iB2 <- iconservantB[2]
          # ratio of conservant B2 in primary vs secondary reaction
          rB2 <- mymolbasis1[iB2] / mymolbasis2[iB2]
          myout.2 <- couplefun(rB2,molspecies2,molbasis2,molspecies3,molbasis3,r1)
          # okay here we are. scale each of the coupled reactions
          # to the coefficient of the conservant appearing in the
          # replacement reaction (the one whose boundary we are at)
          rrB1 <- mymolbasis2[iB1] / myout.1$molbasis2[iB1]
          rrB2 <- mymolbasis2[iB2] / myout.2$molbasis2[iB2]
          if(is.infinite(rrB1) | is.infinite(rrB2) | identical(rrB1,0) | identical(rrB2,0)) {
            cat('transfer: fall back to one secondary conservant (some conservant does not appear)\n')
          } else if(abs(log10(abs(rrB1/rrB2))) > 10) {
            cat('transfer: fall back to one secondary conservant (conservants differ too much)\n')
          } else {
            cat('transfer: coupling on',mymolbasis2[iB1],rownames(basis)[iB1],
              mymolbasis2[iB2],rownames(basis)[iB2],'\n')
            for(i in 1:length(myout.1)) myout.1[[i]] <- rrB1 * myout.1[[i]]
            for(i in 1:length(myout.2)) myout.2[[i]] <- rrB2 * myout.2[[i]]
            # now we sum them up
            myout.out <- myout.1
            for(i in 1:length(myout.out)) myout.out[[i]] <- myout.1[[i]] + myout.2[[i]]
            # and scale again, this time to moles of conservant
            mc3 <- myout.out$molbasis3[ibalance]
            mc2 <- myout.out$molbasis2[ibalance]
            mc23 <- mc2 - mc3
            # this value should be opposite and equal to molconservant1
            rr <- - molconservant1 / mc23
            for(i in 1:length(myout.out)) myout.out[[i]] <- rr * myout.out[[i]]
            myout <- myout.out
            # check if that worked
            if(mc23==0)
              cat('transfer: fall back to one secondary conservant (unconstrained secondary formation)\n')
            else if(any(myout$molspecies2 < 0)) 
              cat('transfer: fall back to one secondary conservant (negative species used)\n')
            else if(any(myout$molspecies2 > molspecies2))
              cat('transfer: fall back to one secondary conservant (negative species left)\n')
            else if(any(myout$molspecies3 < 0))
              cat('transfer: fall back to one secondary conservant (negative species formed)\n')
            else {
              iworked <- TRUE
            }
          }
        }

        if(!iworked) {
          cmode <- "one"
          # coming from two, the conservant is one of the endmembers
           iconB <- which(sign(mymolbasis1) != sign(mymolbasis2) & 
             mymolbasis1 != 0 & mymolbasis2 != 0 & 
             1:nrow(basis) != ibalance &
             1:nrow(basis) %in% icouple &
             abs(log10(abs(mymolbasis1/mymolbasis2))) < abs(logpresent/2))
          ncb <- length(iconB)
          # if the last step worked and the possible
          # conservants have changed, reset the count
          if(!icb %in% 1:ncb) icb <- 1
          iconservantB <- iconB
            if(ncb==0) {
              cat('transfer: coupling: no conservant\n')
              didwork <- c(didwork,FALSE)
              nextfun()
              next
            } else {
              iconservantB <- iconservantB[icb]
              # ratio of conservant B in primary to secondary reactions
              r2 <- mymolbasis1[iconservantB] / mymolbasis2[iconservantB]
              ctext <- rownames(basis)[iconservantB]
            }
            if(r2 > 1) {
              cat(paste('transfer: coupling: insufficient ',
                rownames(basis)[iconservantB],' in reactant\n',sep=''))
              didwork <- c(didwork,FALSE)
              nextfun()
              next
            } else {
              cat(paste('transfer: coupling replacement on',
                ctext,'\n')) 
              myout <- couplefun(r2,molspecies2,molbasis2,molspecies3,molbasis3,r1)
              if(any(myout$molspecies2 < 0)) 
                cat('transfer: coupling failed (negative species used)\n')
              else if(any(myout$molspecies2 > molspecies2))
                cat('transfer: coupling failed (negative species left)\n')
              else if(any(myout$molspecies3 < 0))
                cat('transfer: coupling failed (negative species formed)\n')
              else {
                iworked <- TRUE
              }
            }
          }
      }
      if(iworked) {
        molspecies2 <- myout$molspecies2
        molbasis2 <- myout$molbasis2
        molconservant2 <- molbasis2[ibalance]
        molspecies3 <- myout$molspecies3
        molbasis3 <- myout$molbasis3
      } else {
        didwork <- c(didwork,FALSE)
        nextfun()
        next
      }
    }  # done with coupled mode

    # 4 - calculate the change in moles of basis species
    # and loop if the change is negative or too large
    # summing up the primary and secondary reactants
    molspecies <- molspecies1 + molspecies2
    molconservant <- molconservant1 + molconservant2
    # these species are present
    ipresent <- which(log10(molspecies) > logpresent)
    if(length(ipresent)==0) {
      cat('transfer: nothing to destroy\n')
      # for only mode, increase alpha for next step
      # (actual effect will be alpha + 2 + dalpha)
      if(dmode=='only') alpha <- alpha + 2
      didwork <- c(didwork,FALSE)
      nextfun()
      next
    }
    # the final number of moles of basis species
    molbasis <- molbasis0 + molbasis1 + molbasis2 - molbasis3
    # (negative value here is not good ... try smaller step size)
    wm <- which(molbasis < 0)
    if(any(molbasis < 0)) {
      for(i in 1:length(wm))
        cat(paste('transfer: negative moles (',molbasis[wm[i]],') of basis species',
          c2s(rownames(basis)[wm[i]]),sep=' '),'\n')
      didwork <- c(didwork,FALSE)
      nextfun()
      next
    }

    # slow things down if our deviations are becoming huge
    if(!is.null(devmax)) {
      # j - the current step, j1 - the previous step
      logmolbasis.j <- log10(molbasis)
      if(!any(didwork)) logmolbasis.j1 <- basis$logact else logmolbasis.j1 <- basis[,ncol(basis)]
      if(j > 1) {
        if(any(logbuffdev > devmax)) {
          idev <- which(logbuffdev > devmax)
          cat(paste('transfer: change from buffer of',rownames(basis)[idev],'exceeded deviation limits\n'))
          print(logbuffdev)
          didwork <- c(didwork,FALSE)
          nextfun()
          next
        }
        dev <- abs(logmolbasis.j - logmolbasis.j1 - logbuffdev)
        dev.orig <- logmolbasis.j - logmolbasis.j1 - logbuffdev
      } else {
        dev <- abs(logmolbasis.j - logmolbasis.j1)
        dev.orig <- logmolbasis.j - logmolbasis.j1
      }
      # here we find that setting 999 as a dummy logact
      # for the conservant will become Inf and flag
      # a deviation limit...
      if(any(dev > devmax)) {
        idev <- dev > devmax
        cat(paste('transfer: change of',rownames(basis)[idev],
          'by',dev.orig[idev],'exceeds deviation limit\n'))
        didwork <- c(didwork,FALSE)
        nextfun()
        next
      }
    }
    # the secondary destruction reactions
    ipresent <- which(log10(molspecies2) > logpresent)
    if(length(ipresent) > 0) {
      mysl <- 10^species$logact[ipresent] - molspecies2[ipresent]
      if(any(mysl < 0)) {
        cat('transfer: too much secondary reactant; looping\n')
        didwork <- c(didwork,FALSE)
        nextfun()
        next
      }
    }

    # calculating new affinities and checking that
    # we didn't cross a boundary in coupled mode
    aff <- getaff(log10(molbasis),sout)
    mynewaff <- aff$myaff
    if(dmode=='coupled' & cmode=='one') {
      inewstable <- which.max(mynewaff)
      if(inewstable!=istable) {
         cat('transfer: coupling: crossed boundary\n')
         didwork <- c(didwork,FALSE)
         nextfun()
         next
      }
    }

    # from now on the step is ensured
    myaff <- mynewaff
    didwork <- c(didwork,TRUE)
    # the primary destruction reactions
    ipresent <- which(log10(molspecies1) > logpresent)
    if(length(ipresent) > 0) for(i in 1:length(ipresent)) {
      cat(paste('transfer: reacting ',molspecies1[ipresent[i]],' moles of ',
        species$name[ipresent[i]],' (primary) \n',sep=''))
      # we don't actually do anything here
    }
    # report the secondary reactions
    # and update species activities
    ipresent <- which(log10(molspecies2) > logpresent)
    if(length(ipresent) > 0) {
      mysl <- 10^species$logact[ipresent] - molspecies2[ipresent]
      cat(paste('transfer: reacting ',molspecies2[ipresent],' moles of ',
        species$name[ipresent],'\n',sep=''))
      species$logact[ipresent] <- log10(mysl)
    }
    # formation reaction report and update
    ipresent <- which(log10(molspecies3) > logpresent)
    if(length(ipresent) > 0) {
      ipresent <- ipresent[sort(molspecies3[ipresent],index.return=TRUE,decreasing=TRUE)$ix]
      for(i in 1:length(ipresent)) {
        cat(paste('transfer: forming',molspecies3[ipresent[i]],'moles of',
          species$name[ipresent[i]],'\n'))
        species$logact[ipresent[i]] <- log10(10^species$logact[ipresent[i]] + 
          molspecies3[ipresent[i]])
      }
    }
    # clean up logarithms
    species$logact[is.infinite(species$logact)] <- log0
    # finally the new values for moles of basis species
    basis$logact <- log10(molbasis)

    # copy those columns to the end
    species <- cbind(species,species$logact)
    basis <- cbind(basis,basis$logact)
    colnames(species)[ncol(species)] <- colnames(basis)[ncol(basis)] <-
      paste('X',j,sep='')
    istables <- c(istables,istable)
    myaffs <- c(myaff,myaffs)

    # plot the activities of the basis species
    # if identified (by e.g. plot=c(2,3))
    if(!is.null(plot)) {
      myact <- as.numeric(basis[,ncol(basis)])
      basisnames <- rownames(basis)
      elementnames <- colnames(basis)[1:nrow(basis)]
      iZelement <- match('Z',elementnames)
      nH1 <- nH2 <- Hact <- 0
      # calculate activity ratios for charged basis species e.g. aK+/aH+ aMg+2/a2H+
      if(!is.na(iZelement)) {
        for(i in 1:length(plot)) {
          hasZ <- basis[plot[i],iZelement]!=0
          if(hasZ & !basisnames[plot[i]] %in% c('H+','e-')) {
            myact[plot[i]] <- myact[plot[i]] - myact[basisnames=='H+'] * basis$Z[plot[i]]
          }
        }
      }
      # 20090328 take care of pH and pe
      iHe <- which(basisnames %in% c("H+","e-"))
      if(length(iHe) > 0) myact[iHe] <- -myact[iHe]
      points(myact[plot[1]],myact[plot[2]])
    }

  } # end major loop

  # we've finished all that; restore basis and species definitions
  species$logact <- species0$logact
  # this is the second place to be careful of PBB
  if('PBB' %in% rownames(basis)) basis$logact <- c(basis0$logact,0)
  else basis$logact <- basis0$logact
  thermo$basis <- basis0
  thermo$species <- species0
  assign("thermo", thermo, "CHNOSZ")

  # report the success rate and total progress
  aaa <- alphas
  alphas <- alphas[didwork]
  dmodes <- dmodes[didwork]
  progress <- 100*sum(10^alphas)
  cat(paste('transfer: ',length(which(didwork)),'of',nsteps,'steps succeeded (',progress,'% overall progress)\n'))

  # return the results
  return(invisible(list(basis=basis,species=species,alphas=alphas,
    dmodes=dmodes,istables=istables,myaffs=myaffs,didwork=didwork)))

}

draw.transfer <- function(t,ylim=c(-10,1),ylimbasis=c(-12,-2),logprogress=FALSE) {
  # plot the logarithms of activities of basis species
  # and of number of moles of species as a function
  # of reaction progress (in percent completion)
  progress <- 100*sum(10^t$alphas)
  par(mfrow=c(2,1))
  # cumulative sum of progress
  myprogress <- numeric()
  for(i in 1:length(t$alphas)) {
    if(i==1) myoldprogress <- 0 else
      myoldprogress <- myprogress[length(myprogress)]
    myprogress <- c(myprogress,myoldprogress+100*10^t$alphas[i])
  }
  if(logprogress) {
    xlab <- "log percent reaction progress"
    myprogress <- log10(myprogress)
    xlim <- range(myprogress)
  } else {
    xlab <- "percent reaction progress"
    xlim <- c(0,progress)
  }
  # the basis species
  thermo.plot.new(xlim=xlim,ylim=ylimbasis,xlab=xlab,ylab='logact')
  istatecol <- match('state',colnames(t$basis))
  for(i in 1:nrow(t$basis)) {
    lines(myprogress,t$basis[i,(istatecol+1):ncol(t$basis)],lty=i)
  }
  legend(legend=rownames(t$basis),lty=1:nrow(t$basis),x='bottomright')
  # the species of interest
  thermo.plot.new(xlim=xlim,ylim=ylim,xlab=xlab,ylab='logmol')
  inamecol <- match('name',colnames(t$species))
  for(i in 1:nrow(t$species))
    lines(myprogress,t$species[i,(inamecol+1):ncol(t$species)],lty=i)
  legend(legend=t$species$name,lty=1:nrow(t$species),x='bottomright')
}

feldspar <- function(which="closed",plot.it=FALSE) {
  # open- and closed-system reaction paths
  # for weathering of k-feldspar
  # after Steinmann, Lichtner, and Shotyk, 1994
  # call feldspar("closed")
  # or feldspar("open")
  # setup conditions for feldspar reaction
  #basis(c('Al+3','SiO2','K+','H2O','H+','O2'))
  basis(delete=TRUE)
  # SLS89 use H4SiO4 instead of SiO2 - use the secondary database
  add.obigt()
  basis(c('Al+3','H4SiO4','K+','H2O','H+','O2'))
  # some of SLS89's initial conditions
  basis(c('K+','H4SiO4'),c(-6,-6))
  # the candidate species
  species(delete=TRUE)
  species(c('k-feldspar','muscovite','pyrophyllite','kaolinite','gibbsite'))
  # setup a diagram on which to plot a reaction path
  basis('pH',0)
  a <- affinity(H4SiO4=c(-6,-2),'K+'=c(-3,8))
  diagram(a)
  # identify the basis species whose activities will
  # be plotted by transfer()
  plot <- c(2,3) 
  # return to SLS89 initial conditions
  basis('pH',4)
  # start with miniscule amounts of all species
  species(1:5,-999)
  # except feldspar (the primary reactant)
  species('k-feldspar',-4)

  if(which=='closed') {
    # closed system diagram (SLS94 Fig. 2)
    tr <- transfer(550,dmode='coupled',plot=plot,devmax=0.2)
    if(plot.it) draw.transfer(tr)
  } else if(which=='open') {
    # open system (SLS94 Fig. 3)
    # A* - B* - C* - D*
    tr <- transfer(450,dmode='none',plot=plot)     
    # E* - F*
    species(c('k-felsdspar','kaolinite'),c(-999,-4))
    tr <- transfer(150,dmode='none',plot=plot)
    # F* - H*
    species(c('k-feldspar','kaolinite'),c(-4,-999))
    basis('H4SiO4',tr$basis[rownames(tr$basis)=='H4SiO4',ncol(tr$basis)])
    tr <- transfer(420,dmode='none',plot=plot)
    if(plot.it) draw.transfer(tr)
  }
  return(invisible(tr))
}

apc <- function(which="open",basis="CO2",plot.it=FALSE) {
  # react APC2 to other proteins in the anaphase-promoting complex, e.g.
  # apc("open")
  # apc("closed")
  # apc("many")
  # apc("buffer")
  # assign basis species
  basis(delete=TRUE)
  if(basis=="CO2") basis(c("CO2","H2O","NH3","H2","H2S"),c(-10,0,-4,-10,-7))
  else if(basis=="acetic") basis(c("acetic acid","H2O","NH3","H2","H2S"),c(-5.5,0,-4,-10,-7))
  basis("H2","aq")
  # load the proteins
  species(c("APC1","APC2","APC5","CDC16","APC10","SWM1"),"YEAST")
  # create a diagram
  if(basis=="CO2") a <- affinity(CO2=c(-10,0),H2=c(-10,0))
  else if(basis=="acetic") a <- affinity(C2H4O2=c(-10,-2),H2=c(-10,-4))
  diagram(a, normalize=TRUE)
  # set APC2 to react
  species(1:nrow(species()),-999)
  species("APC2_YEAST",0)
  if(which=="open") {
    tr <- transfer(220,ibalance="PBB",plot=c(1,4),dmode="none",devmax=0.2)
    if(plot.it) draw.transfer(tr,ylim=c(-22,-6),logprogress=TRUE)
  } else if(which=="closed") {
    tr <- transfer(510,ibalance="PBB",plot=c(1,4),dmode="coupled",devmax=0.15)
    if(plot.it) draw.transfer(tr,ylim=c(-22,-6),logprogress=TRUE)
  } else if(which=="many") {
    tr <- transfer(250,ibalance="PBB",plot=c(1,4),dmode="none",devmax=0.15,fmode="many")
    if(plot.it) draw.transfer(tr,ylim=c(-15,-8),logprogress=TRUE)
  } else if(which=="buffer") {
    mod.buffer("H2S","H2S","aq",0)
    species("APC2_YEAST",0)
    tr <- transfer(700,ibalance="PBB",plot=c(1,4),dmode="coupled",devmax=0.15,
      buffers=list(basis="H2S",buffer="H2S"),beta=4)
    if(plot.it) draw.transfer(tr,ylim=c(-22,-6),logprogress=TRUE)
  }
  return(invisible(tr))
}

