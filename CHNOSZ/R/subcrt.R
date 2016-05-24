# CHNOSZ/subcrt.R
# calculate standard molal thermodynamic propertes
# 20060817 jmd

subcrt <- function(species, coeff=1, state=NULL, property=c('logK','G','H','S','V','Cp'),
  T=seq(273.15,623.15,25), P='Psat', grid=NULL, convert=TRUE, exceed.Ttr=FALSE,
  logact=NULL, action.unbalanced='warn', IS=0) {

  # revise the call if the states have 
  # come as the second argument 
  if(!is.null(coeff[1])) {
    if(is.numeric(state[1])) newcoeff <- state else newcoeff <- 1
    if(is.character(coeff[1])) newstate <- coeff else newstate <- NULL
    if(is.character(coeff[1])) {
      if(missing(T)) {
        if(identical(newcoeff,1) & !(identical(newcoeff,state))) 
          return(subcrt(species,state=coeff,property=property,P=P,grid=grid,
            convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
          else return(subcrt(species,coeff=newcoeff,state=coeff,property=property,
            P=P,grid=grid,convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
      } else {
        if(identical(newcoeff,1) & !(identical(newcoeff,state))) 
          return(subcrt(species,state=coeff,property=property,T=T,P=P,grid=grid,
            convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
          else return(subcrt(species,coeff=newcoeff,state=coeff,property=property,
            T=T,P=P,grid=grid,convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
      }
    }
  }

  do.reaction <- FALSE
  #if(!missing(coeff) & coeff!=1) do.reaction <- TRUE
  if(!missing(coeff)) do.reaction <- TRUE

  # species and states are made the same length
  if(!is.null(state[1])) {
    if(length(state) > length(species)) species <- rep(species,length.out=length(state))
    if(length(species) > length(state)) state <- rep(state,length.out=length(species))
    state <- state.args(state)
  }

  # allowed properties
  properties <- c('rho','logK','G','H','S','Cp','V','kT','E')
  # property checking
  prop <- tolower(property)
  notproperty <- property[!prop %in% tolower(properties)]
  if(length(notproperty) > 0) stop(paste(notproperty,
    'are not valid properties\ntry rho, logK, G, H, S, V, Cp, kT, or E (or their lowercase equivalents)'))
  # length checking
  if(do.reaction & length(species)!=length(coeff)) 
    stop('coeff must be same length as the number of species.')
  if(length(IS)>1) if(!identical(grid,'IS')) {
    if(is.null(grid)) grid <- 'IS'
    else stop('if you want length(IS) > 1, set grid=\'IS\'')
  }
  if(!is.null(logact)) logact <- rep(logact,length.out=length(coeff))
  # normalize temperature units
  if(!missing(T)) {
    if(convert) T <- envert(T,'K')
    else if(!missing(convert) & convert) T <- envert(T,'K')
  }
  if(is.numeric(P[1])) {
    if(convert) P <- envert(P,'bar')
  }

  # gridding?
  do.grid <- FALSE
  if(!is.null(grid)) if(!is.logical(grid)) do.grid <- TRUE
  newIS <- IS
  if(do.grid) {
    if(grid=='T') {
      newT <- numeric()
      for(i in 1:length(T)) newT <- c(newT,rep(T[i],length(P)))
      newP <- rep(P,length(T))
      T <- newT; P <- newP
    }
    if(grid=='P') {
      newP <- numeric()
      for(i in 1:length(P)) newP <- c(newP,rep(P[i],length(T)))
      newT <- rep(T,length(P))
      T <- newT; P <- newP
    }
    if(grid=='IS') {
      ll <- length(T)
      if(length(P) > 1) ll <- length(P)
      newIS <- numeric()
      for(i in 1:length(IS)) newIS <- c(newIS,rep(IS[i],ll))
      tpargs <- TP.args(T=T,P=P)
      T <- rep(tpargs$T,length.out=length(newIS))
      P <- rep(tpargs$P,length.out=length(newIS))
    }
  } else {
    # expansion of Psat and equivalence of argument lengths
    tpargs <- TP.args(T=T,P=P)
    T <- tpargs$T; P <- tpargs$P
  }

  # get species information
  thermo <- get("thermo")
  # pre-20110808, we sent numeric species argument through info() to
  # get species name and state(s)
  # but why slow things down if we already have a species index?
  # so now phase species stuff will only work for character species names
  if(is.numeric(species[1])) {
    ispecies <- species
    species <- as.character(thermo$obigt$name[ispecies])
    state <- as.character(thermo$obigt$state[ispecies])
    newstate <- as.character(thermo$obigt$state[ispecies])
    sinfo <- ispecies
  } else {
    # from names, get species indices and states and possibly
    # keep track of phase species (cr1 cr2 ...)
    sinfo <- numeric()
    newstate <- character()
    for(i in 1:length(species)) {
      mysearch <- species[i]
      if(can.be.numeric(mysearch)) mysearch <- thermo$obigt$name[as.numeric(mysearch)]
      si <- info.character(mysearch, state[i])
      # that could have the side-effect of adding a protein; re-read thermo
      thermo <- get("thermo", "CHNOSZ")
      if(is.na(si[1])) stop('no info found for ',species[i],' ',state[i])
      if(!is.null(state[i])) is.cr <- state[i]=='cr' else is.cr <- FALSE
      if(thermo$obigt$state[si[1]]=='cr1' & (is.null(state[i]) | is.cr)) {
        newstate <- c(newstate,'cr')
        sinfo <- c(sinfo,si[1])
      } else {
        newstate <- c(newstate,as.character(thermo$obigt$state[si[1]]))
        sinfo <- c(sinfo,si[1])
      }
    }
  }

  # to make the following more readable and maybe save
  # run time, keep some parts of thermo$obigt handy
  ton <- thermo$obigt$name
  tos <- thermo$obigt$state

  # stop if species not found
  noname <- is.na(sinfo)
  if(TRUE %in% noname)
    stop(paste('species',species[noname],'not found.\n'))

  # take care of mineral phases
  state <- as.character(tos[sinfo])
  name <- as.character(ton[sinfo])
  # a counter of all species considered
  # inpho is longer than sinfo if cr1 cr2 ... phases are present
  # sinph shows which of sinfo correspond to inpho
  # pre-20091114: the success of this depends on there not being duplicated aqueous or other
  # non-mineral-phase species (i.e., two entries in obigt for Cu+ screw this up
  # when running the skarn example).
  # after 20091114: we can deal with duplicated species (aqueous at least)
  inpho <- sinph <- coeff.new <- numeric()
  for(i in 1:length(sinfo)) {
     if(newstate[i]=='cr') {
       searchstates <- c('cr','cr1','cr2','cr3','cr4','cr5','cr6','cr7','cr8','cr9') 
       tghs <- thermo$obigt[(ton %in% name[i]) & tos %in% searchstates,]
       # we only take one if they are in fact duplicated species and not phase species
       if(all(tghs$state==tghs$state[1])) tghs <- thermo$obigt[sinfo[i],]
     } else tghs <- thermo$obigt[sinfo[i],]
     inpho <- c(inpho,as.numeric(rownames(tghs))) 
     sinph <- c(sinph,rep(sinfo[i],nrow(tghs)))
     coeff.new <- c(coeff.new,rep(coeff[i],nrow(tghs)))
  }

  # where we keep info about the species involved
  reaction <- data.frame( coeff=coeff.new,name=ton[inpho],
    formula = thermo$obigt$formula[inpho],state=tos[inpho],
    ispecies=inpho, stringsAsFactors=FALSE)
  # make the rownames readable ... but they have to be unique
  if(length(unique(inpho))==length(inpho)) rownames(reaction) <- as.character(inpho)

  # wetness etc.
  isH2O <- reaction$name=='water' & reaction$state=='liq'
  isaq <- reaction$state=='aq'

  #if(length(T)==1) T.text <- paste(T,units('T')) else T.text <- paste(length(T),'values of T')
  #if(length(P)==1) P.text <- paste(P,units('P')) else P.text <- paste(length(P),'values of P')
  ut <- T
  if(identical(grid,'IS')) ut <- unique(ut)
  if(length(ut)==1) T.text <- paste(ut,'K') else {
    T.text <- paste(length(ut),'values of T')
  }
  if(length(P)==1) {
    if(can.be.numeric(P)) P.text <- paste(round(as.numeric(P),2),'bar')
    else P.text <- "P"
  } else P.text <- 'P'
  #} else P.text <- paste(length(P),'values of P')
  if(identical(P[[1]],'Psat')) P.text <- P
  if(any(c(isH2O,isaq))) P.text <- paste(P.text,' (wet)',sep='')
  if(length(species)==1 & convert==FALSE) {
    # we don't think we want messages here
  } else {
    msgout(paste('subcrt:',length(species),'species at',T.text,'and',P.text,'\n'))
  }

  # inform about unbalanced reaction
  if(do.reaction) {
    # the mass balance ... is zero for a balanced reaction
    mss <- makeup(sinfo, coeff, sum=TRUE)
    # take out very small numbers
    mss[abs(mss) < 1e-7] <- 0
    # report and try to fix any non-zero mass balance
    if(any(mss!=0) & !is.null(action.unbalanced)) {
      # the missing composition: the negative of the mass balance
      miss <- -mss
      # drop elements that are zero
      miss <- miss[miss!=0]
      msgout("subcrt: reaction is not balanced; it is missing this composition:\n")
      # we have to do this awkward dance to send a formatted matrix to msgout
      msgout(paste(capture.output(print(miss)), collapse="\n"), appendLF=TRUE)
      # look for basis species that have our compositoin
      tb <- thermo$basis
      if(!is.null(tb)) {
        if(all(names(miss) %in% colnames(tb)[1:nrow(tb)])) {
          # the missing composition as formula
          ft <- as.chemical.formula(miss)
          # the basis species needed to supply it
          bc <- species.basis(ft)
          # drop zeroes
          bc.new <- bc[,(bc[1,]!=0),drop=FALSE]
          # and get the states
          b.state <- as.character(thermo$basis$state)[bc[1,]!=0]
          bc <- bc.new
          # special thing for Psat
          if(P.text=='Psat') P <- P.text
          else P <- outvert(P,"bar")
          # add to logact values if present
          if(!is.null(logact)) {
            ila <- match(colnames(bc),rownames(thermo$basis))
            nla <- !(can.be.numeric(thermo$basis$logact[ila]))
            if(any(nla)) warning('subcrt: logact values of basis species',
              c2s(rownames(thermo$basis)[ila]),'are NA.')
            logact <- c(logact,thermo$basis$logact[ila])
          }
          # warn user and do it!
          ispecies.new <- tb$ispecies[match(colnames(bc),rownames(tb))]
          b.species <- thermo$obigt$formula[ispecies.new]
          if(identical(species,b.species) & identical(state,b.state))
            msgout("subcrt: balanced reaction, but it is a non-reaction; restarting...\n")
          else msgout('subcrt: adding missing composition from basis definition and restarting...\n')
          newspecies <- c(species, tb$ispecies[match(colnames(bc), rownames(tb))])
          newcoeff <- c(coeff, as.numeric(bc[1, ]))
          newstate <- c(state, b.state)
          return(subcrt(species=newspecies, coeff=newcoeff, state=newstate,
            property=property, T=outvert(T, "K"), P=P, grid=grid, convert=convert, logact=logact, exceed.Ttr=FALSE))
        } else if(identical(action.unbalanced,'warn')) 
            warning(paste('reaction was unbalanced, missing', as.chemical.formula(miss)),call.=FALSE)
      } else {
        if(identical(action.unbalanced,'warn')) 
          warning(paste('reaction was unbalanced, missing', as.chemical.formula(miss)),call.=FALSE)
      }
    }
  }

  # calculate the properties
  # if we want affinities we must have logK
  if(!is.null(logact)) if(!'logk' %in% prop) prop <- c('logk',prop)
  # if logK but not g was requested, get g ...
  if('logk' %in% prop & ! 'g' %in% prop) eprop <- c(prop,'g') else eprop <- prop
  # don't request logk from the eos ...
  eosprop <- eprop[!eprop %in% c('logk','rho')]
  # also get g if we are dealing with mineral phases
  if(!'g' %in% eprop & length(inpho) > length(sinfo)) eosprop <- c(eosprop,'g')
  # the reaction result is in out
  out <- list()
  # aqueous species
  if(TRUE %in% isaq | 'rho' %in% eprop) {
    # load the water properties (better here, once,
    # than possible many times in hkf()).
    wprop.PT <- character()
    wprop.PrTr <- 'rho'
    dosupcrt <- thermo$opt$water != "IAPWS95"
    if(TRUE %in% (prop %in% c('logk','g','h','s'))) wprop.PrTr <- c(wprop.PrTr,'YBorn')
    if(dosupcrt | TRUE %in% (prop %in% c('logk','g','h'))) wprop.PrTr <- c(wprop.PrTr,'diel')
    H2O.PrTr <- water(wprop.PrTr,T=thermo$opt$Tr,P=thermo$opt$Pr)
    if(TRUE %in% (prop %in% c('cp'))) {wprop.PT <- c(wprop.PT,'XBorn','YBorn')}
    if(TRUE %in% (prop %in% c('v'))) {wprop.PT <- c(wprop.PT,'QBorn')}
    if(TRUE %in% (prop %in% c('kt'))) {wprop.PT <- c(wprop.PT,'NBorn')}
    if(TRUE %in% (prop %in% c('e'))) {wprop.PT <- c(wprop.PT,'UBorn')}
    # get additional properties required for omega derivatives
    if(dosupcrt) wprop.PT <- c(wprop.PT,'alpha','daldT','beta','diel')
    H2O.PT <- water(c(wprop.PrTr,wprop.PT),T=T,P=P)
    if(TRUE %in% isaq) {
      # now the species stuff
      # 20110808 if inpho are the species indices let's avoid
      # the overhead of info() and use new obigt2eos() instead
      #si <- info(inpho[isaq],quiet=TRUE)
      si <- obigt2eos(thermo$obigt[inpho[isaq],], "aq", fixGHS = TRUE)
      domega <- thermo$obigt$name[inpho[isaq]] != 'H+'
      p.aq <- hkf(eosprop,T=T,P=P,ghs=si,eos=si,H2O.PT=H2O.PT,H2O.PrTr=H2O.PrTr,domega=domega)
      if(any(IS!=0)) p.aq <- nonideal(inpho[isaq],p.aq,newIS,T)
      out <- c(out,p.aq)
    }
  }
  # crystalline, gas, liquid (except water) species
  iscgl <- reaction$state %in% c('liq','cr','gas','cr1','cr2','cr3',
    'cr4','cr5','cr6','cr7','cr8','cr9') & reaction$name != 'water'

  if(TRUE %in% iscgl) {
    #si <- info(inpho[iscgl],quiet=TRUE)
    si <- obigt2eos(thermo$obigt[inpho[iscgl],], "cgl", fixGHS = TRUE)
    p.cgl <- cgl(eosprop,T=T,P=P,ghs=si,eos=si)
    # replace Gibbs energies with NA where the
    # phases are beyond their temperature range
    if('g' %in% eosprop) {
      # 20080304 this code is weird and hard to read - needs a lot of cleanup!
      # 20120219 cleaned up somewhat; using exceed.Ttr and NA instead of do.phases and 999999
      # the numbers of the cgl species (becomes 0 for any that aren't cgl)
      ncgl <- iscgl
      ncgl[iscgl] <- 1:nrow(si)
      for(i in 1:length(iscgl)) {
        # not if we're not cgl
        if(!iscgl[i]) next
        # not if dPdTtr is calling us
        if(identical(caller.name(3),"dPdTtr")) next
        # name and state
        myname <- reaction$name[i]
        mystate <- reaction$state[i]
        # check if we're below the transition temperature
        if(!(reaction$state[i] %in% c('cr1','liq','cr','gas'))) {
          Ttr <- Ttr(inpho[i]-1,P=P,dPdT=dPdTtr(inpho[i]-1))
          if(any(T < Ttr)) {
            status.Ttr <- "(extrapolating G)"
            if(!exceed.Ttr) {
              # put NA into the value of G
              p.cgl[[ncgl[i]]]$G[T<Ttr] <- NA
              status.Ttr <- "(using NA for G)"
            } 
            msgout(paste('subcrt: some points below transition temperature for',myname, mystate, status.Ttr, '\n'))
          }
        }
        # check if we're above the transition temperature
        if(!(reaction$state[i] %in% c('cr','liq','gas')))
          Ttr <- Ttr(inpho[i],P=P,dPdT=dPdTtr(inpho[i]))
        else {
          Ttr <- thermo$obigt$z.T[inpho[i]]
          if(is.na(Ttr)) next
        }
        if(all(Ttr==0)) next
        if(any(T >= Ttr)) {
          status.Ttr <- "(extrapolating G)"
          if(!exceed.Ttr) {
            p.cgl[[ncgl[i]]]$G[T>=Ttr] <- NA
            status.Ttr <- "(using NA for G)"
          }
          msgout(paste('subcrt: some points above transition temperature for',myname, mystate, status.Ttr, '\n'))
        }
      }
    }
    out <- c(out,p.cgl)
  }

  # water
  if(TRUE %in% isH2O) {
    if(!exists('H2O.PT',inherits=FALSE)) H2O.PT <- water('rho',T=T,P=P)
    if(length(eosprop)==0) eosprop <- 'rho'
    #msgout(paste('subcrt: water equation of state:',c2s(eosprop),'\n'))
    p.H2O <- list(tmp=water(eosprop,T=T,P=P))
    out <- c(out,rep(p.H2O,length(which(isH2O==TRUE))))
  }

  # use variable-pressure standard Gibbs energy for gases
  isgas <- reaction$state %in% "gas" 
  if(TRUE %in% isgas & "g" %in% eprop & thermo$opt$varP) {
    for(i in which(isgas)) out[[i]]$G <- out[[i]]$G - convert(log10(P), "G", T=T)
  }

  # logK
  if('logk' %in% prop) {
    for(i in 1:length(out)) {
      # NOTE: the following depends on the water function renaming g to G
      out[[i]] <- cbind(out[[i]],data.frame(logK=convert(out[[i]]$G,'logK',T=T)))
      colnames(out[[i]][ncol(out[[i]])]) <- 'logK'
    }
  }

  # ordering the output
  # the indices of the species in out thus far
  ns <- 1:nrow(reaction)
  is <- c(ns[isaq],ns[iscgl],ns[isH2O])
  if(length(ns)!=length(is)) stop('subcrt: not all species are accounted for.')
  v <- list()
  for(i in 1:length(is))  v[[i]] <- out[[match(ns[i],is)]]
  out <- v

  # deal with phases (cr1 cr2) here
  # we have to eliminate rows from out, 
  # reaction and values from isaq, iscgl, isH2O
  out.new <- list()
  reaction.new <- reaction
  isaq.new <- logical()
  iscgl.new <- logical()
  isH2O.new <- logical()
  #print(sinfo)
  #print(sinph)
  #print(reaction)
  for(i in 1:length(sinfo)) {
    iphases <- which(sinfo[i]==sinph)
    # deal with repeated species here ... divide iphases 
    # by the number of duplicates
    #print(iphases)
    if(TRUE %in% duplicated(inpho[iphases])) {
      iphases <- iphases[length(which(sinfo==sinfo[i]))]
    }
    if(length(iphases)>1) {
      msgout(paste('subcrt:',length(iphases),'phases for',thermo$obigt$name[sinfo[i]],'... '))
      # assemble the Gibbs energies for each species
      for(j in 1:length(iphases)) {
        G.this <- out[[iphases[j]]]$G
        if(length(which(is.na(G.this))) > 0 & exceed.Ttr) warning(paste('subcrt: NAs found for G of ',
          reaction$name[iphases[j]],' ',reaction$state[iphases[j]],' at T-P point(s) ',
          c2s(which(is.na(G.this)),sep=' '),sep=''),call.=FALSE)
        if(j==1) G <- as.data.frame(G.this)
        else G <- cbind(G,as.data.frame(G.this))
      }
      # find the minimum-energy phase at each T-P point
      phasestate <- numeric()
      out.new.entry <- out[[1]]
      for(j in 1:nrow(G)) {
        ps <- which.min(as.numeric(G[j,]))
        if(length(ps)==0) {
          # minimum not found: NAs have crept in (like something wrong with Psat?)
          # (or no non-NA value of G to begin with, e.g. aegerine)
          ps <- 1
          if(exceed.Ttr) warning('subcrt: stable phase for ',reaction$name[iphases[ps]],' at T-P point ',j,
          ' undetermined (using ',reaction$state[iphases[ps]],')',call.=FALSE)
        } 
        phasestate <- c(phasestate,ps)
        out.new.entry[j,] <- out[[ iphases[ps] ]][j,]
      }

      # update our objects
      out.new[[i]] <- cbind(out.new.entry,data.frame(state=phasestate))
      reaction.new[i,] <- reaction[iphases[phasestate[1]],]
      # mark the minerals with multiple phases
      rs <- as.character(reaction.new$state)
      rs[i] <- 'cr*'
      reaction.new$state <- rs
      isaq.new <- c(isaq.new,isaq[iphases[phasestate[1]]])
      iscgl.new <- c(iscgl.new,iscgl[iphases[phasestate[1]]])
      isH2O.new <- c(isH2O.new,isH2O[iphases[phasestate[1]]])
      # info for the user
      up <- unique(phasestate)
      if(length(up)>1) { word <- 'are'; p.word <- 'phases' }
      else { word <- 'is'; p.word <- 'phase' }
      msgout(paste(p.word,c2s(unique(phasestate)),word,'stable\n'))
    } else {
      # multiple phases aren't involved ... things stay the same
      out.new[[i]] <- out[[iphases]]
      # hmm.. this could mess up our coefficients 20091103
      #reaction.new[i,] <- reaction[iphases,]
      coeff.orig <- reaction$coeff
      reaction.new[i,] <- reaction[iphases,]
      reaction.new$coeff <- coeff.orig
      rs <- as.character(reaction.new$state)
      rs[i] <- as.character(reaction$state[iphases])
      reaction.new$state <- rs
      isaq.new <- c(isaq.new,isaq[iphases])
      iscgl.new <- c(iscgl.new,iscgl[iphases])
      isH2O.new <- c(isH2O.new,isH2O[iphases])
    }
  }
  out <- out.new
  reaction <- reaction.new[1:length(sinfo),]
  isaq <- isaq.new
  iscgl <- iscgl.new
  isH2O <- isH2O.new

  #print(out)

  newprop <- eprop[eprop!='rho']
  # the order of the properties
  #if(ncol(out[[1]])>1) for(i in 1:length(out)) {
  if(length(newprop)>1) for(i in 1:length(out)) {
    # keep state/loggam columns if they exists
    ipp <- match(newprop,tolower(colnames(out[[i]])))
    if('state' %in% colnames(out[[i]])) ipp <- c(ipp,match('state',colnames(out[[i]]))) 
    if('loggam' %in% colnames(out[[i]])) ipp <- c(ipp,match('loggam',colnames(out[[i]]))) 
    out[[i]] <- out[[i]][,ipp,drop=FALSE]
  }

  # add up reaction properties
  if(do.reaction) {
    o <- 0
    statecols <- NULL
    # do our affinity calculations here
    if(!is.null(logact)) {
      logQ <- logK <- rep(0,length(T))
      for(i in 1:length(coeff)) {
        logK <- logK + out[[i]]$logK * coeff[i]
        logQ <- logQ + logact[i] * coeff[i]
      }
      reaction <- cbind(reaction,logact)
      A <- logK - logQ
      # convert A/2.303RT (no dims) to cal mol-1
      # then to the user's units (outvert) from cal
      A <- outvert(convert(-A,'G',T=T),'cal')
    }
    # the addition of properties
    for(i in 1:length(coeff)) {
      # assemble state columns if they exist
      if('state' %in% colnames(out[[i]])) {
         sc <- as.data.frame(out[[i]]$state)
         out[[i]] <- out[[i]][,-match('state',colnames(out[[i]]))]
         colnames(sc) <- as.character(reaction$name[i])
         if(is.null(statecols)) statecols <- sc
         else statecols <- cbind(statecols,sc)
      }
      # include a zero loggam column if we need it
      # for those species that are ideal
      o.i <- out[[i]]
      if('loggam' %in% colnames(o.i)) if(!'loggam' %in% colnames(o))
        o <- cbind(o,loggam=0)
      if('loggam' %in% colnames(o)) if(!'loggam' %in% colnames(o.i))
        o.i <- cbind(o.i,loggam=0)
      o <- o + o.i * coeff[i]
    }
    # output for reaction (stack on state columns if exist)
    if(!is.null(statecols)) out <- list(reaction=reaction,out=o,state=statecols)
    else out <- list(reaction=reaction,out=o)
  } else {
    # output for species: strip the coeff column from reaction
    reaction <- reaction[,-match('coeff',colnames(reaction))]
    out <- c(list(species=reaction),out)
  }
  # append T,P,rho, A, logQ columns and convert units
  for(i in 2:length(out)) {
    # affinity and logQ
    if(i==2) if(do.reaction & !is.null(logact)) {
      out[[i]] <- cbind(out[[i]],data.frame(logQ=logQ,A=A))
    }
    # 20120114 only prepend T, P, rho columns if we have more than one T
    if(length(T) > 1) {
      # 20090329 added checks for converting T, P units
      if(convert) T.out <- outvert(T,"K") else T.out <- T
      if(convert) P.out <- outvert(P,"bar") else P.out <- P
      # try to stuff in a column of rho if we have aqueous species
      # watch out! supcrt-ish densities are in g/cc not kg/m3
      if('rho' %in% prop | (missing(property) & any(c(isaq,isH2O))) & (names(out)[i])!='state') 
        out[[i]] <- cbind(data.frame(T=T.out,P=P.out,rho=H2O.PT$rho/1000),out[[i]])
      else
        out[[i]] <- cbind(data.frame(T=T.out,P=P.out,out[[i]]))
    }
    if(convert) {
      for(j in 1:ncol(out[[i]])) {
        if(colnames(out[[i]])[j] %in% c('G','H','S','Cp')) out[[i]][,j] <- outvert(out[[i]][,j],'cal')
      }
    }
  }
  # put ionic strength next to any loggam columns
  for(i in 2:length(out)) if('loggam' %in% colnames(out[[i]])) out[[i]] <- cbind(out[[i]],IS=newIS)
  # more fanagling for species
  if(!do.reaction) {
    out <- list(species=out$species,out=out[2:length(out)])
    # add names to the output
    names(out$out) <- as.character(reaction$name)
  }
  return(out)
}

