# CHNOSZ/info.R
# search database for species names or formulas
# and retrieve thermodynamic properties of species
# 20061024 extraced from species.R jmd
# 20120507 code rewrite and split into info.[character,approx,numeric];
#   these functions expect arguments of length 1; 
#   info() handles longer arguments

info.text <- function(ispecies) {
  # a textual description of species name, formula, source, e.g.
  # CO2 [CO2(aq)] (SSW01, SHS89, 11.Oct.07)
  this <- get("thermo")$obigt[ispecies, ]
  sourcetext <- this$ref1
  ref2 <- this$ref2
  if(!is.na(ref2)) sourcetext <- paste(sourcetext, ref2, sep=", ")
  date <- this$date
  if(!is.na(date)) sourcetext <- paste(sourcetext, date, sep=", ")
  out <- paste(this$name, " [", this$formula, "(", this$state, ")", "] (", sourcetext, ")", sep="")
  return(out)
}

info.character <- function(species, state=NULL, check.protein=TRUE) {
  # returns the rownumbers of thermo$obigt having an exact match of 'species' to
  # thermo$obigt$[species|abbrv|formula] or NA otherwise
  # a match to thermo$obigt$state is also required if 'state' is not NULL
  # (first occurence of a match to species is returned otherwise)
  thermo <- get("thermo")
  # all matches for the species
  matches.species <- thermo$obigt$name==species | 
    thermo$obigt$abbrv==species | thermo$obigt$formula==species
  # since thermo$obigt$abbrv contains NAs, convert NA results to FALSE
  matches.species[is.na(matches.species)] <- FALSE
  # turn it in to no match if it's a protein in the wrong state
  ip <- suppressMessages(iprotein(species))
  if(any(matches.species) & !is.na(ip) & !is.null(state)) {
    matches.state <- matches.species & grepl(state, thermo$obigt$state)
    if(!any(matches.state)) matches.species <- FALSE
  }
  # no match, not available
  if(!any(matches.species)) {
    # unless it's a protein
    if(check.protein) {
      # did we find a protein? add its properties to obigt
      if(!is.na(ip)) {
        # here we use a default state from thermo$opt$state
        if(is.null(state)) state <- thermo$opt$state
        # retrieve the amino acid composition
        aa <- ip2aa(ip)
        # add up protein properties
        eos <- aa2eos(aa, state)
        # the real assignment work 
        nrows <- suppressMessages(mod.obigt(eos))
        thermo <- get("thermo", "CHNOSZ")
        matches.species <- rep(FALSE, nrows)
        matches.species[nrows] <- TRUE
      } else return(NA)
    } else return(NA)
  }
  # do we demand a particular state
  if(!is.null(state)) {
    # special treatment for H2O: aq retrieves the liq
    if(species %in% c("H2O", "water") & state=="aq") state <- "liq"
    # the matches for both species and state
    matches.state <- matches.species & grepl(state, get("thermo")$obigt$state)
    if(!any(matches.state)) {
      # the requested state is not available for this species
      available.states <- thermo$obigt$state[matches.species]
      if(length(available.states)==1) a.s.verb <- "is" else a.s.verb <- "are"
      a.s.text <- paste("'", available.states, "'", sep="", collapse=" ")
      msgout("info.character: requested state '", state, "' for ", species, 
        " but only ", a.s.text, " ", a.s.verb, " available\n")
      return(NA)
    }
    matches.species <- matches.state
  }
  # all of the species that match
  ispecies <- which(matches.species)
  # we return only the first species that matches
  # unless they are 'cr1' 'cr2' etc. and we requested state 'cr'
  if(identical(state, "cr")) ispecies.out <- ispecies
  else ispecies.out <- ispecies[1]
  # let user know if there is more than one state for this species
  if(length(ispecies) > length(ispecies.out)) {
    ispecies.other <- ispecies[!ispecies %in% ispecies.out]
    othertext <- paste(thermo$obigt$state[ispecies.other], collapse=", ")
    msgout("info.character: found ", species, "(", thermo$obigt$state[ispecies.out], 
      "), also available in ", othertext, "\n")
  }
  return(ispecies.out)
}

info.numeric <- function(ispecies, check.it=TRUE) {
  # from a numeric species index in 'ispecies' return the 
  # thermodynamic properties and equations-of-state parameters
  thermo <- get("thermo")
  # if we're called with NA, return an empty row
  if(is.na(ispecies)) {
    this <- thermo$obigt[1,]
    this[] <- NA
    return(this)
  }
  this <- thermo$obigt[ispecies,]
  # species indices must be in range
  ispeciesmax <- nrow(thermo$obigt)
  if(ispecies > ispeciesmax | ispecies < 1) 
    stop(paste("species index", ispecies, "not found in thermo$obigt\n"))
  # remove scaling factors on EOS parameters depending on state
  # use new obigt2eos function here
  this <- obigt2eos(this, this$state)
  # identify any missing GHS values
  naGHS <- which(is.na(this[8:10]))
  # a missing one of G, H or S can cause problems for subcrt calculations at high T
  if(length(naGHS)==1) {
    # calculate a single missing one of G, H, or S from the others
    GHS <- as.numeric(GHS(as.character(this$formula), G=this[,8], H=this[,9], S=this[,10]))
    msgout("info.numeric: ", colnames(this)[8:10][naGHS], " of ",
      this$name, "(", this$state, ") is NA; set to ", round(GHS[naGHS],2), "\n")
    this[, naGHS+7] <- GHS[naGHS]
  } 
  # now perform consistency checks for GHS and EOS parameters if check.it=TRUE
  if(check.it) {
    # check GHS if they were all present
    if(length(naGHS)==0) calcG <- checkGHS(this)
    # check tabulated heat capacities against EOS parameters
    calcCp <- checkEOS(this, this$state, "Cp")
    # fill in NA heat capacity
    if(!is.na(calcCp) & is.na(this$Cp)) {
      msgout("info.numeric: Cp of ", this$name, "(", this$state, ") is NA; set by EOS parameters to ", round(calcCp, 2), "\n")
      this$Cp <- as.numeric(calcCp)
    }
    # check tabulated volumes - only for aq (HKF equation)
    if(identical(this$state, "aq")) {
      calcV <- checkEOS(this, this$state, "V")
      # fill in NA volume
      if(!is.na(calcV) & is.na(this$V)) {
        msgout("info.numeric: V of ", this$name, "(", this$state, ") is NA; set by EOS parameters to ", round(calcV, 2), "\n")
        this$V <- as.numeric(calcV)
      }
    }
  } # done checking
  # all done!
  return(this)
}

info.approx <- function(species, state=NULL) {
  # returns species indices that have an approximate match of 'species'
  # to thermo$obigt$[name|abbrv|formula], 
  # possibly restricted to a given state
  thermo <- get("thermo")
  if(!is.null(state)) this <- thermo$obigt[thermo$obigt$state==state, ]
  else this <- thermo$obigt
  # only look for fairly close matches
  max.distance <- 0.1
  approx.name <- agrep(species, this$name, max.distance)
  approx.abbrv <- agrep(species, this$abbrv, max.distance)
  approx.formula <- agrep(species, this$formula, max.distance)
  approx.species <- unique(c(approx.name, approx.abbrv, approx.formula))
  if(!is.na(approx.species[1])) {
    # show the names of the species
    if(length(approx.species)==1) {
      msgout("info.approx: '", species, "' is similar to ", info.text(approx.species), "\n")
    } else {
      napprox.max <- 25
      exttext <- ":"
      if(length(approx.species) > napprox.max) exttext <- paste(" (showing first ", napprox.max, ")", sep="")
      msgout("info.approx: '", species, "' is ambiguous; has approximate matches to ", 
        length(approx.species), " species", exttext, "\n")
      printout <- capture.output(print(thermo$obigt$name[head(approx.species, napprox.max)]))
      msgout(paste(printout, collapse="\n"), "\n")
    }
    return(approx.species)
  }
  # if we got here there were no approximate matches
  msgout("info.approx: '", species, "' has no approximate matches\n")
  return(NA)
}

info <- function(species=NULL, state=NULL, check.it=TRUE) {
  ## return information for one or more species in thermo$obigt
  ## if no species are requested, summarize the available data  20101129
  thermo <- get("thermo")
  if(is.null(species)) {
    msgout("info: 'species' is NULL; summarizing information about thermodynamic data...\n")
    msgout(paste("thermo$obigt has", nrow(thermo$obigt[thermo$obigt$state=="aq", ]), "aqueous,",
      nrow(thermo$obigt), "total species\n"))
    msgout(paste("number of literature sources: ", nrow(thermo$refs), ", elements: ",
      nrow(thermo$element), ", buffers: ", length(unique(thermo$buffers$name)), "\n", sep=""))
    msgout(paste("number of proteins in thermo$protein is", nrow(thermo$protein), "from",
      length(unique(thermo$protein$organism)), "organisms\n"))
    # print information about SGD.csv, ECO.csv, HUM.csv
    more.aa(organism="Sce")
    more.aa(organism="Eco")
    #pdata.aa(organism="HUM")
    # print information about yeastgfp.csv
    yeastgfp()
    return()
  }
  ## run info.numeric or info.character depending on the input type
  if(is.numeric(species)) {
    out <- lapply(species, info.numeric, check.it)
    # if we different states the column names could be different
    if(length(unique(unlist(lapply(out, names)))) > ncol(thermo$obigt)) {
      # make them the same as thermo$obigt
      out <- lapply(out, function(row) {
        colnames(row) <- colnames(thermo$obigt); return(row)
      }) 
    }
    # turn the list into a data frame
    out <- do.call(rbind, out)
  } else {
    # state and species should be same length
    if(!is.null(state)) {
      lmax <- max(length(species), length(state))
      state <- rep(state, length.out=lmax)
      species <- rep(species, length.out=lmax)
    }
    # loop over the species
    out <- sapply(seq_along(species), function(i) {
      # first look for exact match
      ispecies <- info.character(species[i], state[i])
      # if no exact match and it's not a protein, show approximate matches (side effect of info.approx)
      if(identical(ispecies, NA) & !grepl("_", species[i])) ispecies.notused <- info.approx(species[i], state[i])
      # do not accept multiple matches
      if(length(ispecies) > 1) ispecies <- NA
      return(ispecies)
    })
  }
  ## all done!
  return(out)
}

