# CHNOSZ/util.data.R
# add or change entries in the thermodynamic database

today <- function() {
  # write today's date in the format used in SUPCRT data files
  # e.g. 13.May.12 for 2012-05-13
  t <- date()
  tt <- unlist(strsplit(t, " "))
  # for single-digit days there is an extra space
  tt <- tt[!tt==""]
  tday <- tt[3]
  tmonth <- tt[2]
  tyear <- substr(tt[5], start=3, stop=4)
  return(paste(tday, tmonth, tyear, sep="."))
}

mod.obigt <- function(...) {
  # add or modify species in thermo$obigt
  thermo <- get("thermo")
  # the names and values are in the arguments
  # this works for providing arguments via do.call
  args <- list(...)
  # this is needed if we are called with a list as the actual argument
  if(is.list(args[[1]])) args <- args[[1]]
  if(length(args) < 2) stop("please supply at least a species name and a property to update")
  if(is.null(names(args))) stop("please supply named arguments")
  # if the first argument is numeric, it's the species index
  if(is.numeric(args[[1]][1])) {
    ispecies <- args[[1]]
  } else {
    # if the name of the first argument is missing, assume it's the species name
    if(names(args)[1]=="") names(args)[1] <- "name"
    # search for this species, use check.protein=FALSE to avoid infinite loop when adding proteins
    # and suppressMessages to not show messages about matches of this name to other states
    if("state" %in% names(args)) ispecies <- suppressMessages(mapply(info.character, 
      species=args$name, state=args$state, check.protein=FALSE, SIMPLIFY=TRUE, USE.NAMES=FALSE))
    else ispecies <- suppressMessages(mapply(info.character, 
      species=args$name, check.protein=FALSE, SIMPLIFY=TRUE, USE.NAMES=FALSE))
  }
  # the column names of thermo$obigt, split at the "."
  cnames <- c(do.call(rbind, strsplit(colnames(thermo$obigt), ".", fixed=TRUE)), colnames(thermo$obigt))
  # the columns we are updating
  icol <- match(names(args), cnames)
  if(any(is.na(icol))) stop(paste("properties not in thermo$obigt:", paste(names(args)[is.na(icol)], collapse=" ")) )
  # the column numbers for properties that matched after the split
  icol[icol > 40] <- icol[icol > 40] - 40
  icol[icol > 20] <- icol[icol > 20] - 20
  # which species are new and which are old
  inew <- which(is.na(ispecies))
  iold <- which(!is.na(ispecies))
  # the arguments as data frame
  args <- data.frame(args, stringsAsFactors=FALSE)
  if(length(inew) > 0) {
    # the right number of blank rows of thermo$obigt
    newrows <- thermo$obigt[1:length(inew), ]
    # if we don't know something it's NA
    newrows[] <- NA
    # put in a default state
    newrows$state <- thermo$opt$state
    # the formula defaults to the name
    newrows$formula <- args$name[inew]
    # fill in the columns
    newrows[, icol] <- args[inew, ]
    # now check the formulas
    e <- tryCatch(makeup(newrows$formula), error=function(e) e)
    if(inherits(e, "error")) {
      warning("please supply a valid chemical formula as the species name or in the 'formula' argument")
      # transmit the error from makeup
      stop(e)
    }
    # assign to thermo$obigt
    thermo$obigt <- rbind(thermo$obigt, newrows)
    rownames(thermo$obigt) <- NULL
    assign("thermo", thermo, "CHNOSZ")
    # update ispecies
    ntotal <- nrow(thermo$obigt)
    ispecies[inew] <- (ntotal-length(inew)+1):ntotal
    # inform user
    msgout(paste("mod.obigt: added ", newrows$name, "(", newrows$state, ")", sep="", collapse="\n"), "\n")
  }
  if(length(iold) > 0) {
    # loop over species
    for(i in 1:length(iold)) {
      # the old values and the state
      oldprop <- thermo$obigt[ispecies[iold[i]], icol]
      state <- thermo$obigt$state[ispecies[iold[i]]]
      # tell user if they're the same, otherwise update the data entry
      if(isTRUE(all.equal(oldprop, args[iold[i], ], check.attributes=FALSE))) 
        msgout("mod.obigt: no change for ", args$name[iold[i]], "(", state, ")\n")
      else {
        thermo$obigt[ispecies[iold[i]], icol] <- args[iold[i], ]
        assign("thermo", thermo, "CHNOSZ")
        msgout("mod.obigt: updated ", args$name[iold[i]], "(", state, ")\n")
      }
    }
  }
  return(ispecies)
}

add.obigt <- function(file=system.file("extdata/thermo/OBIGT-2.csv",package="CHNOSZ"),
  force=FALSE,E.units="cal") {
  # add/replace entries in thermo$obigt from values saved in a file
  # only replace if force==TRUE
  if(missing(file)) {
    # we use force=TRUE for the default data file
    if(missing(force)) force <- TRUE
  }
  thermo <- get("thermo")
  to1 <- thermo$obigt
  id1 <- paste(to1$name,to1$state)
  to2 <- read.csv(file,as.is=TRUE)
  id2 <- paste(to2$name,to2$state)
  # check if the file is compatible with thermo$obigt
  tr <- try(rbind(to1,to2),silent=TRUE)
  if(identical(class(tr),'try-error')) stop(paste(file,"is not compatible with thermo$obigt data table."))
  # match the new species to existing ones
  does.exist <- id2 %in% id1
  ispecies.exist <- na.omit(match(id2, id1))
  nexist <- sum(does.exist)
  # convert from J if necessary
  if(tolower(E.units)=="j") {
    # loop over each row
    for(i in 1:nrow(to2)) {
      # GHS and EOS parameters
      icol <- (8:18)
      # if it's aqueous, also include omega
      if(to2$state[i]=="aq") icol <-(8:19)
      # don't touch volume (in column 12)
      icol <- icol[icol!=12]
      # convert to calories
      to2[i,icol] <- convert(to2[i,icol],"cal")
    }
  }
  # keep track of the species we've added
  inew <- numeric()
  if(force) {
    # replace existing entries
    if(nexist > 0) {
      to1[ispecies.exist, ] <- to2[does.exist, ]
      to2 <- to2[!does.exist, ]
      inew <- c(inew, ispecies.exist)
    }
  } else {
    # ignore any new entries that already exist
    to2 <- to2[!does.exist, ]
    nexist <- 0
  }
  # add new entries
  if(nrow(to2) > 0) {
    to1 <- rbind(to1, to2)
    inew <- c(inew, (length(id1)+1):nrow(to1))
  }
  # commit the change
  thermo$obigt <- to1
  rownames(thermo$obigt) <- 1:nrow(thermo$obigt)
  assign("thermo", thermo, "CHNOSZ")
  # message about file, if file argument is missing (default)
  if(missing(file)) {
    msgout("add.obigt: using default file:\n") 
    msgout(file, "\n")
  }
  msgout("add.obigt: read ", length(does.exist), " rows; made ", 
    nexist, " replacements, ", nrow(to2), " additions, units = ", E.units, "\n")
  msgout("add.obigt: use data(thermo) to restore default database\n")
  return(invisible(inew))
}

browse.refs <- function(key=NULL) {
  ## browse to web page associated with a given source
  ## of thermodynamic data. first version: 20110615
  # 'key' can be
  # NULL: show a table of all sources in a browser
  # character: open a web page for each listed source
  # numeric: open one or two web pages for each listed species
  # list: the output of subcrt()
  ## first retrieve the sources table
  thermo <- get("thermo")
  x <- thermo$refs
  ## show a table in the browser if 'key' is NULL 
  if(is.null(key)) {
    # create the html links
    cite <- x$citation
    x$citation <- sprintf("<a href='%s' target='_blank'>%s</a>", x$URL, cite)
    notlinked <- x$URL=="" | is.na(x$URL)
    x$citation[notlinked] <- cite[notlinked]
    # remove the last (URL) component
    #x$URL <- NULL
    x <- x[1:4]
    # count the times each source is listed in OBIGT.csv
    ns1 <- sapply(x$key, function(x) length(which(thermo$obigt$ref1==x)) )
    ns1.2 <- sapply(x$key, function(x) length(which(thermo$obigt$ref2==x)) )
    ns1 <- ns1 + ns1.2
    ns1[ns1==0] <- ""
    # count the times each source is listed in OBIGT-2.csv
    o2 <- read.csv(system.file("extdata/thermo/OBIGT-2.csv", package = "CHNOSZ"))
    ns2 <- sapply(x$key, function(x) length(which(o2$ref1==x)) )
    ns2.2 <- sapply(x$key, function(x) length(which(o2$ref2==x)) )
    ns2 <- ns2 + ns2.2
    ns2[ns2==0] <- ""
    # count the times each source is listed in protein.csv
    npr <- sapply(x$key, function(x) length(which(thermo$protein$ref==x)) )
    npr[npr==0] <- ""
    # count the times each source is listed in stress.csv
    stressfile <- system.file("extdata/abundance/stress.csv", package="CHNOSZ")
    stressdat <- read.csv(stressfile, check.names=FALSE, as.is=TRUE)
    nst <- sapply(x$key, function(x) length(which(stressdat[2,]==x)) )
    nst[nst==0] <- ""
    # append the counts to the table to be shown
    x <- c(x,list(ns1=ns1,ns2=ns2,npr=npr,nst=nst))
    # title to display for web page
    title <- "Sources of Thermodynamic Data in CHNOSZ"
    ### the following is adapted from print.findFn in package 'sos'
    f0 <- tempfile()
    File <- paste(f0, ".html", sep="")
    Dir <- dirname(File)
    js <- system.file("extdata/js", "sorttable.js", package = "CHNOSZ")
    file.copy(js, Dir)
    ## Sundar's original construction:
    con <- file(File, "wt")
    on.exit(close(con))
    .cat <- function(...)
      cat(..., "\n", sep = "", file = con, append = TRUE)
    ## start
    cat("<html>", file = con)
    .cat("<head>")
    .cat("<title>", title, "</title>")
    .cat("<script src=sorttable.js type='text/javascript'></script>")
    .cat("</head>")
    ### boilerplate text
    .cat("<h1>Listing of all entries in thermo$refs</h1>")
    .cat("<h3>Click on hyperlinked references to open URL in new window</h3>")
    .cat("<h3>Click on column headers to sort</h3>")
    .cat("<h3>Columns 'n..' give number of times each reference appears in data tables:</h3>")
    .cat("ns1: 'ref1' and 'ref2' in data/OBIGT.csv<br>")
    .cat("ns2: 'ref1' and 'ref2' in extdata/thermo/OBIGT-2.csv<br>")
    .cat("npr: 'ref' in data/protein.csv<br>")
    .cat("nst: second row in data/stress.csv<br><p>")
    ### start table and headers
    .cat("<table class='sortable' border='1'>\n<thead>")
    .cat("<tr>")
    .cat(sprintf("  <th>%s</th>\n</tr>",
                 paste(names(x), collapse = "</th>\n  <th>")))
    .cat("</thead>\n<tbody>")
    ### now preparing the body of the table
    paste.list <- c(lapply(x, as.character), sep = "</td>\n  <td>")
    tbody.list <- do.call("paste", paste.list)
    tbody <- sprintf("<tr>\n  <td>%s</td>\n</tr>", tbody.list)
    tbody <- sub("<td><a", "<td class=link><a", tbody, useBytes = TRUE)
    .cat(tbody)
    ### finish it!
    .cat("</tbody></table></body></html>")
    ### end adaptation from print.findFn
    # show table in browser
    browseURL(File)
    cat("browse.refs: table of references is shown in browser\n")
  } else if(is.character(key)) {
    # open the URL(s) of the given source(s)
    for(i in seq_along(key)) {
      ix <- match(key[i],x$key)
      if(is.na(ix)) {
        cat(paste("browse.refs: reference key",key[i],"not found\n"))
        next
      } 
      URL <- x$URL[ix]
      if(URL=="" | is.na(URL)) {
        cat(paste("browse.refs: no URL available for reference key",key[i],"\n"))
        next
      }
      cat(paste("browse.refs: opening URL for ",key[i]," (",x$author[ix],", ",x$year[ix],")\n",sep=""))
      browseURL(x$URL[ix])
    }
    return(invisible(URL))
  } else if(is.numeric(key)) {
    # open the URL(s) of sources associated with the indicated species
    sinfo <- suppressMessages(info(key))
    mysources <- unique(c(sinfo$ref1,sinfo$ref2))
    mysources <- mysources[!is.na(mysources)]
    return(browse.refs(mysources))
  } else if(is.list(key)) {
    if("species" %in% names(key)) ispecies <- key$species$ispecies
    else if("reaction" %in% names(key)) ispecies <- key$reaction$ispecies
    else stop("list does not appear to be a result from subcrt()")
    if(is.null(ispecies)) stop("list does not appear to be a result from subcrt()")
    return(browse.refs(ispecies))
  }
}

obigt2eos <- function(obigt,state,fixGHS=FALSE) {
  # remove scaling factors from EOS parameters
  # and apply column names depending on the EOS
  if(identical(state, "aq")) {
    obigt[,13:20] <- t(t(obigt[,13:20]) * 10^c(-1,2,0,4,0,4,5,0))
    colnames(obigt)[13:20] <- c('a1','a2','a3','a4','c1','c2','omega','Z') 
  } else {
    obigt[,13:20] <- t(t(obigt[,13:20]) * 10^c(0,-3,5,0,-5,0,0,0))
    colnames(obigt)[13:20] <- c('a','b','c','d','e','f','lambda','T')
  }
  if(fixGHS) {
    # fill in one of missing G, H, S
    # for use esp. by subcrt because NA for one of G, H or S 
    # will hamper calculations at high T
    # which entries are missing just one
    imiss <- which(rowSums(is.na(obigt[,8:10]))==1)
    if(length(imiss) > 0) {
      for(i in 1:length(imiss)) {
        # calculate the missing value from the others
        ii <- imiss[i]
        GHS <- as.numeric(GHS(as.character(obigt$formula[ii]),G=obigt[ii,8],H=obigt[ii,9],S=obigt[ii,10]))
        icol <- which(is.na(obigt[ii,8:10]))
        obigt[ii,icol+7] <- GHS[icol]
      }
    }
  }
  return(obigt)
}

checkEOS <- function(eos, state, prop, ret.diff=FALSE) {
  # compare calculated properties from equation-of-state
  # parameters with reference (tabulated) values
  # print message and return the calculated value
  # if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  thermo <- get("thermo")
  # get calculated value based on EOS
  if(identical(state, "aq")) {
    if(prop=="Cp") {
      # value of X consistent with IAPWS95
      X <- -2.773788E-7
      # we use the value of X consistent with SUPCRT
      X <- -3.055586E-7
      refval <- eos$Cp
      calcval <- eos$c1 + eos$c2/(298.15-thermo$opt$Theta)^2 + eos$omega*298.15*X
      tol <- thermo$opt$Cp.tol
      units <- "cal K-1 mol-1"
    } else if(prop=="V") {
      # value of Q consistent with IAPWS95
      Q <- 0.00002483137
      # value of Q consistent with SUPCRT92
      Q <- 0.00002775729
      refval <- eos$V
      calcval <- 41.84*eos$a1 + 41.84*eos$a2/2601 + 
        (41.84*eos$a3 + 41.84*eos$a4/2601) / (298.15-thermo$opt$Theta) - Q * eos$omega
      tol <- thermo$opt$V.tol
      units <- "cm3 mol-1"
    }
  } else {
    # all other states
    if(prop=="Cp") {
      refval <- eos$Cp
      Tr <- thermo$opt$Tr
      calcval <- eos$a + eos$b*Tr + eos$c*Tr^-2 + eos$d*Tr^-0.5 + eos$e*Tr^2 + eos$f*Tr^eos$lambda
      tol <- thermo$opt$Cp.tol
      units <- "cal K-1 mol-1"
    }
  }
  # calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else {
    # return the calculated value
    # if the difference is higher than tol
    if(!is.na(calcval)) {
      if(!is.na(refval)) {
        if(abs(diff) > tol) {
          msgout(paste("checkEOS: ", prop, " of ", eos$name, " ", eos$state, " (", rownames(eos),
            ") differs by ", round(diff,2), " ", units, " from tabulated value\n", sep=""))
          return(calcval)
        }
      } else return(calcval)
    }
  }
  # return NA in most cases
  return(NA)
}

checkGHS <- function(ghs, ret.diff=FALSE) {
  # compare calculated G from H and S with reference (tabulated) values
  # print message and return the calculated value if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  thermo <- get("thermo")
  # get calculated value based on H and S
  ina <- is.na(ghs$formula)
  if(any(ina)) {
    msgout("checkGHS: formula of ", ghs$name[ina], "(", ghs$state[ina], ") is NA\n")
    Se <- NA
  } else Se <- entropy(as.character(ghs$formula))
  refval <- ghs[,8]
  DH <- ghs[,9]
  S <- ghs[,10]
  calcval <- DH - thermo$opt$Tr * (S - Se)
  # now on to the comparison
  # calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else if(!is.na(calcval)) {
    if(!is.na(refval)) {
      diff <- calcval - refval
      if(abs(diff) > thermo$opt$G.tol) {
        msgout(paste("checkGHS: G of ", ghs$name, " ", ghs$state, " (", rownames(ghs),
          ") differs by ", round(diff), " cal mol-1 from tabulated value\n", sep=""))
        return(calcval)
      }
    } else return(calcval)
  } else {
    # calculating a value of G failed, perhaps b/c of missing elements
    return(NULL)
  }
  # return NA in most cases
  return(NA)
}

check.obigt <- function() {
  # function to check self-consistency between
  # values of Cp and V vs. EOS parameters
  # and among G, H, S values
  # 20110808 jmd replaces 'check=TRUE' argument of info()
  checkfun <- function(what) {
    # looking at thermo$obigt or OBIGT-2.csv
    if(what=="OBIGT") to <- get("thermo")$obigt
    else if(what=="OBIGT-2") {
      file <- system.file("extdata/thermo/OBIGT-2.csv",package="CHNOSZ")
      to <- read.csv(file,as.is=1:7)
    }
    ntot <- nrow(to)
    # where to keep the results
    DCp <- DV <- DG <- rep(NA,ntot)
    # first get the aqueous species
    isaq <- to$state=="aq"
    eos.aq <- obigt2eos(to[isaq,],"aq")
    DCp.aq <- checkEOS(eos.aq,"aq","Cp",ret.diff=TRUE)
    DV.aq <- checkEOS(eos.aq,"aq","V",ret.diff=TRUE)
    cat(paste("check.obigt: GHS for",sum(isaq),"aq species in",what,"\n"))
    DG.aq <- checkGHS(eos.aq,ret.diff=TRUE)
    # store the results
    DCp[isaq] <- DCp.aq
    DV[isaq] <- DV.aq
    DG[isaq] <- DG.aq
    # then other species, if they are present
    if(sum(!isaq) > 0) {
      eos.cgl <- obigt2eos(to[!isaq,],"cgl")
      DCp.cgl <- checkEOS(eos.cgl,"cgl","Cp",ret.diff=TRUE)
      cat(paste("check.obigt: GHS for",sum(!isaq),"c,g,l species in",what,"\n"))
      DG.cgl <- checkGHS(eos.cgl,ret.diff=TRUE)
      DCp[!isaq] <- DCp.cgl
      DG[!isaq] <- DG.cgl
    }
    # put it all together
    out <- data.frame(table=what,ispecies=1:ntot,name=to$name,state=to$state,DCp=DCp,DV=DV,DG=DG)
    return(out)
  }
  # check both databases in CHNOSZ
  out1 <- checkfun("OBIGT")
  out2 <- checkfun("OBIGT-2")
  out <- rbind(out1,out2)
  # set differences within a tolerance to NA
  out$DCp[abs(out$DCp) < 1] <- NA
  out$DV[abs(out$DV) < 1] <- NA
  out$DG[abs(out$DG) < 500] <- NA
  # take out species where all reported differences are NA
  ina <- is.na(out$DCp) & is.na(out$DV) & is.na(out$DG)
  out <- out[!ina,]
  # round the values
  out$DCp <- round(out$DCp,2)
  out$DV <- round(out$DV,2)
  out$DG <- round(out$DG)
  # how to make the file at extdata/thermo/obigt_check.csv
  # write.csv(out,"obigt_check.csv",na="",row.names=FALSE)
  # return the results
  return(out)
}

RH2obigt <- function(compound=NULL, state="cr", file=system.file("extdata/thermo/RH98_Table15.csv", package="CHNOSZ")) {
  # get thermodynamic properties and equations of state parameters using 
  # group contributions from Richard and Helgeson, 1998   20120609 jmd
  # read the compound names, physical states, chemical formulas and group stoichiometry from the file
  # we use check.names=FALSE because the column names are the names of the groups,
  # and are not syntactically valid R names, and stringsAsFactors=FALSE
  # so that formulas are read as characters (for checking with as.chemical.formula)
  dat <- read.csv(file, check.names=FALSE, stringsAsFactors=FALSE)
  # "compound" the compound names and states from the file
  comate.arg <- comate.dat <- paste(dat$compound, "(", dat$state, ")", sep="")
  # "compound" the compound names and states from the arguments
  if(!is.null(compound)) comate.arg <- paste(compound, "(", state, ")", sep="")
  # identify the compounds
  icomp <- match(comate.arg, comate.dat)
  # check if all compounds were found
  ina <- is.na(icomp)
  if(any(ina)) stop(paste("compound(s)", paste(comate.arg[ina], collapse=" "), "not found in", file))
  # initialize output data frame
  out <- get("thermo")$obigt[0, ]
  # loop over the compounds
  for(i in icomp) {
    # the group stoichiometry for this compound
    thisdat <- dat[i, ]
    # take out groups that are NA or 0
    thisdat <- thisdat[, !is.na(thisdat)]
    thisdat <- thisdat[, thisdat!=0]
    # identify the groups in this compound
    igroup <- 4:ncol(thisdat)
    ispecies <- info(colnames(thisdat)[igroup], state=thisdat$state)
    # check if all groups were found
    ina <- is.na(ispecies)
    if(any(ina)) stop(paste("group(s)", paste(colnames(thisdat)[igroup][ina], collapse=" "), "not found in", thisdat$state, "state"))
    # group additivity of properties and parameters: add contributions from all groups
    thiseos <- t(colSums(get("thermo")$obigt[ispecies, 8:20] * as.numeric(thisdat[, igroup])))
    # group additivity of chemical formula
    formula <- as.chemical.formula(colSums(i2A(ispecies) * as.numeric(thisdat[, igroup])))
    # check if the formula is the same as in the file
    if(!identical(formula, thisdat$formula)) 
      stop(paste("formula", formula, "of", comate.dat[i], "(from groups) is not identical to", thisdat$formula, "(listed in file)" ))
    # build the front part of obigt data frame
    thishead <- data.frame(name=thisdat$compound, abbrv=NA, formula=formula, state=thisdat$state, 
      ref1=NA, ref2=NA, date=today(), stringsAsFactors=FALSE)
    # insert the result into the output
    out <- rbind(out, cbind(thishead, thiseos))
  }
  return(out)
}
