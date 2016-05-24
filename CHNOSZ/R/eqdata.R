# eqdata.R
# extract data from EQ6 output files
# 20091028 just for aqueous species
# 20110410 also get results for solid phases
# 20110516 use grand summary of solid phases, handle data blocks
#   with blank rows, add progress messages, get activity of water
# 20110805 add mineral saturation states and speciation summaries

eqdata <- function(file, species, property="log act", outfile=TRUE) {

  # the available properties for different types of data:
  # solid phases, aqueous species, mineral saturation states, speciation summary
  props <- list(
    solid = c("product", "log moles", "moles", "grams", "volume, cc"),
    aqueous = c("species", "moles", "grams", "conc", "log conc", "log g", "log act"),
    mineral = c("affinity, kcal"),
    speciation = c("molal conc", "per cent")
  )

  ## process the 'property' argument to
  ## figure out the data type for the requested property
  # match property to each list element of props
  iprops <- sapply(props, function(x,y) match(y,x), y=property)
  # is this property available?
  iok <- which(!is.na(iprops))
  if(length(iok) == 0) stop(paste('property "', property, 
    '" is not available for any data type',sep=''))
  # only take the first match
  # (moles and grams will be matched to solid but not aqueous)
  itype <- iok[1]
  # the position (column number) of this property in the data block
  iprop <- iprops[itype]
  # the name of our data type (solid, aqueous, mineral, or speciation)
  mytype <- names(iprop)
  cat(paste("eqdata: data type for", property, "is", mytype, "\n"))
  # if we're extracting speciation data, can only do it for one basis species
  if(mytype == "speciation") {
    if(length(species) > 1) 
      stop("speciation data can only be had for a single basis species")
  }
  # header lines that begin the result block for this data type
  headers <- list(
    # summary of solid product phases
    #solid = " product                    log moles        moles",
    # grand summary of solid phases
    solid = "      phase/end-member       log moles        moles",
    aqueous = "   species                moles        grams",
    mineral = "   mineral        affinity, kcal    state",
    speciation = paste("aqueous species accounting for 99% or more of",species)
  )
  header <- headers[[itype]]
  # lines that let us know the data block has ended
  enders <- list(
    solid = "--- mineral saturation state summary ---",
    aqueous = "--- major aqueous species contributing to mass balances ---",
    mineral = "--- summary of gas species ---",
    speciation = "- - - - - - - - - - - - - - - - - - - - - - -"
  )
  ender <- enders[[itype]]
  ## done processing 'property' argument

  # string constants that show up in the file
  zistring <- " stepping to zi="
  Tstring <- "                     temperature    ="
  H2Ostring <- "                              activity of water ="

  # to find the lines identifying each step of xi
  zilines <- function(lines) grep(zistring, lines)
  # to find the lines indicating the temperature
  Tlines <- function(lines) grep(Tstring, lines)
  # to find the lines with activity of water
  H2Olines <- function(lines) grep(H2Ostring, lines)

  # to read the moles etc of aqueous species or solid phases
  getdata <- function(lines, ihead, property, species) {
    # the property of a species is NA unless we find it
    # so we make a list with an entry for each species
    # the entries are NAs repeated to length of ihead
    prop.out <- lapply(1:length(species),function(x) rep(NA,length(ihead)))
    names(prop.out) <- species
    # loop over each header line
    for(i in 1:length(ihead)) {
      # the line number of the header
      n <- ihead[i]
      # start reading two lines below the header
      dn <- 2
      # or starting 4 lines below the header for 'speciation'
      if(mytype == "speciation") dn <- 4
      myline <- lines[n+dn]
      # read until we reach the end of the block
      while( length(grep(ender, myline)) == 0 ) {
        # read all the values that occur on this line
        myvalues <- strsplit(myline, " ")[[1]]
        myvalues <- myvalues[!myvalues==""]
        # logic below depdends on the data type
        if(mytype == "mineral") {
          # look for names of minerals whose saturation states are desired
          jspecies <- which(myvalues %in% species)
          if(length(jspecies) > 0) {
            # a loop, in case two desired minerals are on the same line
            for(j in 1:length(jspecies)) {
              ispecies <- match(myvalues[jspecies[j]],species)
              prop.out[[ispecies]][i] <- as.numeric(myvalues[jspecies[j]+1])
            }
          }
        } else if(mytype == "speciation") {
          # speciation summary
          # what's the name of this species and value of the property?
          myspec <- myvalues[[1]]
          myval <- myvalues[[iprop + 1]]
          # if this is the first line of the first data block,
          # then this is the first species in prop.out
          if(i==1 & dn==4) names(prop.out)[1] <- myspec
          else if(!myspec %in% names(prop.out)) {
            # if this species isn't already in prop.out, add it
            prop.new <- list(rep(NA, length(ihead)))
            names(prop.new) <- myspec
            prop.out <- c(prop.out, prop.new)
          }
          # save the value for this species
          ispec <- match(myspec, names(prop.out))
          prop.out[[ispec]][i] <- as.numeric(myval)
        } else {
          # properties of aqueous species or solid phases
          ispecies <- match(myvalues[1], species)
          # did we hit one of the desired species or phases?
          if(!is.na(ispecies)) prop.out[[ispecies]][i] <- as.numeric(myvalues[iprop])
        }
        # go to the next line
        dn <- dn + 1
        myline <- lines[n+dn]
      }
    }
    prop.out <- data.frame(prop.out, check.names=FALSE)
    return(prop.out)
  }

  # to get the values of zi or T for each of the data blocks identified by ihead
  getziT <- function(lines, iziT, ihead, Astring, Bstring=NULL) {
    # the linenumbers of zi/T before each of the data blocks
    myziT <- sapply(ihead,function(x, y) {max(iziT[x > y])}, iziT)
    # to get values of zi/T:
    # first split at Astring (beginning of line)
    # then at Bstring (after value)
    # then convert to numeric
    if(is.null(Bstring)) ziT <- sapply(myziT,function(x) {
      as.numeric(strsplit(lines[x], Astring)[[1]][2]) } )
    else ziT <- sapply(myziT,function(x) {
      as.numeric(strsplit(strsplit(lines[x], Astring)[[1]][2], Bstring)[[1]][1]) } )
    return(ziT)
  }

  # put it all together
  # first read the entire file
  lines <- readLines(file)
  cat(paste("eqdata: read", length(lines), "lines from", file, "\n"))
  # get the line numbers where the data blocks start
  # without fixed=TRUE this fails for e.g. zn+2 !!
  ihead <- grep(header, lines, fixed=TRUE)
  if(length(ihead)==0) stop(paste("no data blocks found for", mytype))
  cat(paste("eqdata: found", length(ihead), "data blocks starting at line", ihead[1], "\n"))
  # get values of zi for these data blocks
  izi <- zilines(lines)
  zi <- getziT(lines, izi, ihead, zistring, ",")
  # get values of T for these data blocks
  iT <- Tlines(lines)
  T <- getziT(lines, iT, ihead, Tstring, " degrees")
  # get activity of water
  iH2O <- H2Olines(lines)
  aH2O <- getziT(lines, iH2O, ihead, H2Ostring)
  # get properties of aqueous species or solid phases
  mydata <- getdata(lines, ihead, property, species)
  # make a single table
  out <- cbind(zi, T, aH2O, mydata)
  # save it to a file
  if(isTRUE(outfile) | is.character(outfile)) {
    if(!is.character(outfile)) outfile <- paste(file, property, "csv", sep=".")
    write.csv(out, outfile, row.names=FALSE)
    cat(paste("eqdata: saved results to", outfile, "\n"))
  } 
  # done!
  return(invisible(out))

}
