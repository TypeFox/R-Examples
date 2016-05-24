##==============================================================================
##                                                                            ##
## LINEAR INVERSE MODELLING INPUT   -   LIM                                   ##
##    Karline Soetaert                                                        ##
##                                                                            ##
## -----------------------------                                              ##
## part 1: LIM problem, input as a list                                       ##
## -----------------------------                                              ##
##                                                                            ##
## read.limfile :Reads inverse input text file; returns a list                ##
## setup.limfile:Creates inverse model using list returned by read.limfile    ##
##                                                                            ##
##==============================================================================


##==============================================================================
## Reads an inverse input file, and puts information in a list.
## The file contains declarations denoted inbetween two headers
## Header is a a line starting with "## "
## and then the declarations
## Each header is terminated with ### END
## the "!" sign is a comment
##==============================================================================

Setup <- function(...) UseMethod("Setup")
Setup.character <- function(...) Setup.limfile(...)


Read <-  function (file, verbose = FALSE, checkLinear = TRUE, remtabs = TRUE) {

  if (is.null(file))
    return(NULL)

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## internal function:
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
## find the start and end of a section - only 1st 4 char
## of section name are considered
##------------------------------------------------------------------------------

  findsection <- function(Names) {

    pos <- NULL
    for (i in seq(1, numhead)) {
      sec1 <- toupper(gsub("[ #]", "", Lines[headline[i]]))
      if (substr(sec1, 1, 4) %in% Names)
        pos <- c(pos, i)
    }
    if (length(pos) > 1)
      stop("error: section is present more than once",  Names)
    if (is.null(pos))
      return(NA)
    sec2 <- toupper(gsub("[ #]", "", Lines[headline[pos + 1]]))
    if (length(grep("END", sec2)) == 0)
       warning("heading not properly ended")
    sec2 <- gsub("END", "", sec2)
    if (!substr(sec2, 1, 4) %in% Names)
      warning("heading not properly ended")
    if (headline[pos] == headline[pos + 1] - 1)
      return(NA)
    findsection <- c(headline[pos] + 1, headline[pos + 1] - 1)
  }

##------------------------------------------------------------------------------
## replace everything between pos[1] and pos[2]
##------------------------------------------------------------------------------

  Replace <- function(Strings, replacement, pos, All = TRUE) {

    for (i in 1:length(Strings)) {
      if (class(pos) == "character") {
        i1 <- regexpr(pos[1], Strings[i], fixed = TRUE)[1]
        i2 <- regexpr(pos[2], Strings[i], fixed = TRUE)[1]
      }
      else {
        i1 <- pos[1]
        i2 <- pos[2]
        All <- FALSE
      }

      if (i1 > 0 & i2 > 0) {
        pattern <- substr(Strings[i], i1, i2)
        Strings[i] <- gsub(pattern, replacement, Strings[i], fixed = TRUE)
      } else if (i1>0 & i2 <0) {
        Strings[i] <- gsub(pos[1], replacement, Strings[i], fixed = TRUE)
      } else if (i1<0 & i2 >0)
        Strings[i] <- gsub(pos[2], replacement, Strings[i], fixed = TRUE)
    }
    if (All & length(grep(pos[1], Strings, fixed = TRUE)) > 0)
      Strings <- Replace(Strings, replacement, pos)
    if (All & length(grep(pos[2], Strings, fixed = TRUE)) > 0)
      Strings <- Replace(Strings, replacement, pos)
    return(Strings)
  }

##------------------------------------------------------------------------------
## remove everything between pos[1] and pos[2]
##------------------------------------------------------------------------------
  Remove <- function(Strings, pos, All = TRUE)
    Replace(Strings, "", pos, All = All)

##------------------------------------------------------------------------------
## remove everything after first occurrence of pos
##------------------------------------------------------------------------------
  RemoveToEnd <- function(Strings, pos) {
    len <- length(Strings)
    while (length(grep(pos, Strings, fixed = TRUE)) > 0) {
      for (i in 1:len) {
        cpos <- regexpr(pos, Strings[i], fixed = TRUE)[1]
        if (cpos > 0)
          Strings[i] <- substr(Strings[i], 1, cpos - 1)
      }
    }
    return(Strings)
  }

##------------------------------------------------------------------------------
## split a sentence in a name and a value
##------------------------------------------------------------------------------
  Splitvaluename <- function(Strings) {
    nlen <- length(Strings)
    if (nlen == 0)
      return(NULL)
    pp <- strsplit(Strings, "=")
    res <- NULL
    for (i in 1:nlen) {
      pa <- 0
      vv <- NA
      if (!is.na(pp[[i]][2])) {
         vv <- suppressWarnings(as.double(pp[[i]][2]))
         if (is.na(vv))
            vv <- 1
         pa <- which(toupper(pars$name) == toupper(pp[[i]][2]))
      }
      res <- rbind(res, data.frame(name = pp[[i]][1], val = vv, parnr = pa))
    }
    return(res)
  }
##------------------------------------------------------------------------------
## finding parameters, components, externals etc..
##------------------------------------------------------------------------------
  Findpar <- function(string) {
        ss <- which(toupper(string) == toupper(parnames))
        if (length(ss) == 0)
            ss <- NA
        return(ss)
  }
  Findvar <- function(string) {
        ss <- which(toupper(string) == toupper(varnames))
        if (length(ss) == 0)
            ss <- NA
        return(ss)
  }
  Findcomp <- function(string) {
        ss <- which(toupper(string) == toupper(compnames))
        if (length(ss) == 0)
            ss <- NA
        return(ss)
  }
  Findexternal <- function(string) {
        ss <- which(toupper(string) == toupper(externnames))
        if (length(ss) == 0)
            ss <- NA
        return(ss)
  }
  Findflow <- function(string) {
        if (is.null(nrow(flows)))
            return(NA)
        ss <- toupper(string)
        if (length(grep("FLOWTO(", ss, fixed = TRUE)) > 0) {
            ss <- gsub(")", "", ss)
            ss <- gsub("FLOWTO(", "", ss, fixed = TRUE)
            cc <- Findcomp(ss)
            if (is.na(cc))
                cc <- -Findexternal(ss)    # externals are neg. numbers
            flownr <- which(flows$to == cc)
            # check if not bidirectional flows
            if (any(!posreac[flownr]))
               stop(paste("Cannot expand", string,"some flows are bidirectional"))
            return(flownr)
        }
        if (length(grep("FLOWFROM(", ss, fixed = TRUE)) > 0) {
            ss <- gsub(")", "", ss)
            ss <- gsub("FLOWFROM(", "", ss, fixed = TRUE)
            cc <- Findcomp(ss)
            if (is.na(cc))
                cc <- -Findexternal(ss)
            flownr <- which(flows$from == cc)
            if (any(!posreac[flownr]))
               stop(paste("Cannot expand", string,"some flows are bidirectional"))
            return(flownr)
        } 
        ss <- gsub(")", "", ss)
        ss <- gsub("FLOW(", "", ss, fixed = TRUE)
        ss <- gsub("@@", ",", ss)
        ss <- gsub("@", ",", ss)        
        flownr <- which(toupper(flows$name) == ss)[1]
        if (is.na(flownr)) flownr <- which(toupper(flows$fname) == ss)[1]
        return(flownr)
  }
  FindMass <- function(string) {
    if (is.null(nrow(marker)))
      return(NA)
    ss <- toupper(string)
    if (length(grep("MASSBALANCE(", ss, fixed = TRUE)) >  0) {
      ss <- gsub(")", "", ss)
      ss <- gsub("MASSBALANCE(", "", ss, fixed = TRUE)
      cc <- Findcomp(ss)
      vcc <- marker$val[cc]
      if (is.na(cc))
        return(NA)
      flowto <- which(flows$to == cc)
      cfrom <- flows[flowto, ]$from
      cfrom[cfrom < 0] <- -cfrom[cfrom < 0] + length(compnames)
      vals <- marker$val[cfrom] - vcc
      return(list(ff = flowto, vv = vals))
    }
  }
  Findreact <- function(string) {
    string <- gsub("@","->",string)
    ss <- which(toupper(string) == toupper(reacnames[,1]))
    if (length(ss) == 0)
      ss <- which(toupper(string) == toupper(reacnames[,2]))
    if (length(ss) == 0)
      ss <- NA
    return(ss)
  }
##------------------------------------------------------------------------------
## splitting flows into components from and components to
##------------------------------------------------------------------------------
  Splitflows <- function(Flowlines,createcomp=FALSE) {
      nr <- 1
      Flowlines <- RemoveToEnd(Flowlines, "=")
      pos <- rep(TRUE,length(Flowlines))
      pos[grep("@@", Flowlines, fixed = TRUE) ]<-FALSE

      ss <- gsub(")", "", toupper(Flowlines))
      ss <- gsub("FLOW(", "", ss, fixed = TRUE)
      ss <- gsub("@@", ",", ss)
      SS <- gsub("@", ",", ss)
      comps <- strsplit(SS, ",")
      Flows <- NULL
      for (i in 1:length(Flowlines)) {
          ss <- comps[[i]]
          from <- Findcomp(ss[1])
          if (is.na(from))
              from <- -Findexternal(ss[1])
          if (is.na(from) && createcomp )
             {ES <- emptystruct(ss[1])
              ES$nr <- nr
              nr <- nr+1
              comp <<- rbind(comp,ES)
              compnames <<- c(compnames,toupper(ss[1]))
              from       <- nrow(comp)}

          to <- Findcomp(ss[2])
          if (is.na(to))
              to <- -Findexternal(ss[2])
          if (is.na(to) && createcomp )
             {ES <- emptystruct(ss[2])
              ES$nr <- nr
              nr <- nr+1
              comp <<- rbind(comp,ES)
              compnames <<- c(compnames,toupper(ss[2]))
              to <- nrow(comp)}

          Flows <- rbind(Flows, data.frame(from = from, to = to))
      }
      Flows <- cbind(Flows, fname = SS)
      return(list(Flows=Flows,pos=pos))
  }
##------------------------------------------------------------------------------
## parsing items  (e.g. -1*A*rB)
##------------------------------------------------------------------------------
  finditem <- function(ss, ct = 1, res = NULL) {
    if (ss == "")
      return(res)
    ss <- unlist(strsplit(ss, "\\-"))
      if (length(ss) > 1) {
        if (ss[1] != "")
          res <- finditem(ss[1], ct = ct, res = res)
        for (j in 2:length(ss)) res <- finditem(ss[j], ct = -1 *
              ct, res = res)
          return(res)
      }
      string <- unlist(strsplit(ss, "\\*"))  # split in multiplicative terms
      vv <- ff <- cc <- ee <- rr <- NA       #vars,flows,components,externals,rates
      pp <- c(NA, NA, NA, NA)                #parameters, up to 4 can be multiplied
      np <- 0                                #number of parameters
      Val <- ct                              #constant value
      GG <- NA                               #mass balances
      for (SS in string) {                   # multiplication with a NA -> NA
        if (SS == "XNA") {
           val <- NA
           break()
        }
        val <- Val * suppressWarnings(as.double(SS))
        if (is.na(val)){
          Val <- val <- Val * 1
          PP <- Findpar(SS)
            if (is.na(PP)) {
              FF <- Findflow(SS)
              if (any(is.na(FF))) {
                VV <- Findvar(SS)
                if (is.na(VV)) {
                  CC <- Findcomp(SS)
                  if (is.na(CC)) {
                    EE <- Findexternal(SS)
                    if (is.na(EE)) {
                      RR <- Findreact(SS)
                      if (is.na(RR)) {
                        GG <- FindMass(SS)
                        if (any(is.na(GG)) && !createcomp) {
                          stop(paste("in string", ss, "cannot find item",SS))
                        } else
                        if (any(is.na(GG))) {   # create new component
                          ES <- emptystruct(SS)
                          ES$nr <- length(compnames) +1
                          comp <<- rbind(comp,ES)
                          compnames <<- c(compnames,toupper(SS))
                          cc <- CC <- ES$nr
                        }
                      }
                      else {
                        if (!is.na(rr))
                        stop(paste("in string", ss, "there is more than one reaction"))
                        rr<-RR
                      }
                    }
                    else {
                      if (!is.na(ee))
                      stop(paste("in string", ss, "there is more than one external"))
                      ee <- EE
                    }
                  }
                  else {
                    if (!is.na(cc))
                      stop(paste("in string", ss, "there is more than one component"))
                    cc <- CC
                  }
                }
                else {
                  if (!is.na(vv))
                    stop(paste("in string", ss, "there is more than one variable"))
                  vv <- VV
                }
              }
              else {
                if (!any(is.na(ff)))
                  stop(paste("in string", ss, "there is more than one flow"))
                ff <- FF
              }
            }
            else {
              np <- np + 1
              if (np > 4)
                stop(paste("in string", ss, "there are more than 4 parameters"))
              pp[np] <- PP
            }
         }
         else Val <- val
      }
      if (!any(is.na(GG))) {
          ff <- GG$ff
          val <- val * GG$vv
      }
      if (checkLinear) {
          if (!is.na(vv)) {
              if (varunknown[vv]) {
                if (colunknown == 9 & !is.na(ff))   # unknowns are flows
                  stop(paste("system non linear: in string",
                    ss, ": a flow is multiplied with a variable containing a flow"))
                if (colunknown == 10 & !is.na(cc))  # unknowns are components
                  stop(paste("system non linear: in string",
                    ss, ": a component is multiplied with a variable containing a component"))
                if (colunknown == 12 & !is.na(rr))  # unknowns are reactions
                  stop(paste("system non linear: in string",
                    ss, ": a reaction is multiplied with a variable containing a reaction"))
              }
          }
      }
      return(rbind(res, cbind(val = val, par1 = pp[1], par2 = pp[2],
          par3 = pp[3], par4 = pp[4], var = vv, flow = ff,
          comp = cc, external = ee, reaction = rr)))
  }
##------------------------------------------------------------------------------
## splitting strings in left and right section
##------------------------------------------------------------------------------
  Splitleftright <- function(string, del = "=", check = FALSE) {
    hasdel <- grep(del, string)
    if (check & length(hasdel) != length(string))
      stop(paste("cannot split in left and right: delimiter not found", del))
    leftright <- matrix(data = string, nrow = length(string), ncol = 2)

    for (i in hasdel)
      leftright[i, ] <- unlist(strsplit(string[i], del, fixed = TRUE))[1:2]

    colnames(leftright) <- c("left", "right")
    return(leftright)
  }
##------------------------------------------------------------------------------
## merging strings on several lines
##------------------------------------------------------------------------------
  Mergestrings <- function(strings, sep = "&") {
    if (length(strings) == 0)
       return(NULL)
    ii <- grep(sep, strings)
    if (length(ii) <= 0)
       return(strings)
    for (i in rev(ii)) {
      strings[i] <- paste(strings[i], strings[i + 1], sep = "")
      strings <- strings[-(i + 1)]
    }
    strings <- gsub("&", "", strings)
    return(strings)
  }
##------------------------------------------------------------------------------
## if no values declared: creates an empty structure
##------------------------------------------------------------------------------
  emptystruct <- function(strings) {
    if (length(strings) == 0)
      return(NULL)

    return(data.frame(name = strings, nr = 1:length(strings),
           val = NA, par1 = NA, par2 = NA, par3 = NA, par4 = NA,
           var = NA, flow = NA, comp = NA, external = NA, reaction =NA))
  }
##------------------------------------------------------------------------------
## splitting an equation in a name and equation sector
##------------------------------------------------------------------------------
  SplitNameEquation <- function(strings) {
    if (length(strings) == 0)
      return(NULL)
    if (length(grep("=", strings)) <= 0)
      return(emptystruct(strings))
    varLR <- Splitleftright(strings, "=")
    varnames <- varLR[, 1]
    vars <- ParseLine(varLR[, 2], names=varnames)
    return(vars)
  }
##------------------------------------------------------------------------------
## splitting an equation in left and right section for = , >, and <
##------------------------------------------------------------------------------
  SplitEquation <- function(string) {
    equal <- grep("=", string)
    large <- grep(">", string)
    small <- grep("<", string)
    leftright <- matrix(nrow = length(string), ncol = 2)
    if (length(equal) > 0)
      leftright[equal, ] <- matrix(ncol = 2, byrow = TRUE,
            data = unlist(strsplit(string[equal], "=", fixed = TRUE)))
    if (length(large) > 0)
      leftright[large, ] <- matrix(ncol = 2, byrow = TRUE,
            data = unlist(strsplit(string[large], ">", fixed = TRUE)))
    if (length(small) > 0)
      leftright[small, ] <- matrix(ncol = 2, byrow = TRUE,
            data = unlist(strsplit(string[small], "<", fixed = TRUE)))
    type <- vector(length = length(string))
    type[large] <- ">"
    type[small] <- "<"
    type[equal] <- "="
    ii <- which (type=="FALSE")
    if (length(ii)>0)
      stop(paste("cannot proceed: following is not an equation: ",string[ii]))
    leftright <- cbind(leftright, type)
    colnames(leftright) <- c("left", "right", "type")
    return(leftright)
    }
##------------------------------------------------------------------------------
## splitting inequalities containing [] in > and <
##------------------------------------------------------------------------------
    CleanInequality <- function(Strings) {
      Strings <- gsub(">=", ">", Strings)
      Strings <- gsub("<=", "<", Strings)
      Strings <- gsub("==", "=", Strings)
      IneqString <- NULL
      for (ss in Strings) {
        eqsign <- regexpr("=", ss, fixed = TRUE)[1]
        if (eqsign < 0) {
          IneqString <- c(IneqString, ss)
          next()
        }
        s2 <- s1 <- ss
        while (regexpr("[", s1, fixed = TRUE)[1] > 0) {
          opened <- regexpr("[", s1, fixed = TRUE)[1]
          closed <- regexpr("]", s1, fixed = TRUE)[1]
          comma <- opened - 1 + regexpr(",", substr(s1,
                opened, closed), fixed = TRUE)[1]
          if (eqsign < opened) {
             s1 <- Remove(s1, c(opened, comma))
             s1 <- sub("]", "", s1, fixed = TRUE)
             opened <- regexpr("[", s2, fixed = TRUE)[1]
             closed <- regexpr("]", s2, fixed = TRUE)[1]
             comma <- opened - 1 + regexpr(",", substr(s2,
               opened, closed), fixed = TRUE)[1]
             s2 <- Remove(s2, c(comma, closed))
             s2 <- sub("[", "", s2, fixed = TRUE)
           }
           else {
             s1 <- Remove(s1, c(comma, closed))
             s1 <- sub("[", "", s2, fixed = TRUE)
             opened <- regexpr("[", s2, fixed = TRUE)[1]
             closed <- regexpr("]", s2, fixed = TRUE)[1]
             comma <- opened - 1 + regexpr(",", substr(s2,
                  opened, closed), fixed = TRUE)[1]
             s2 <- Remove(s2, c(opened, comma))
             s2 <- sub("]", "", s2, fixed = TRUE)
         }
       }
       s1 <- sub("=", "<", s1, fixed = TRUE)
       s2 <- sub("=", ">", s2, fixed = TRUE)
       IneqString <- c(IneqString, s1, s2)
       IneqString <- gsub("]", "", IneqString, fixed = TRUE)
       IneqString <- gsub("[", "", IneqString, fixed = TRUE)
     }
     return(IneqString)
  }
##------------------------------------------------------------------------------
## splitting reaction containing <-> into one -> and one <-
##------------------------------------------------------------------------------
  CleanReaction <- function(Strings) {
    Strings <- gsub("@@", "<->", Strings)
    Strings <- gsub("@", "->", Strings)
    Names    <- Splitleftright(Strings, ":")
    nn <- which (grep(":",Strings)<0)
    Names [nn,1] <- paste("f",1:length(nn),sep="")
    Strings<- Names[,2]
    pos <- regexpr("<->", Strings, fixed = TRUE) <=0
    Strings <- gsub("<->", "->", Strings)
    Strings <- gsub("-", "", Strings)
    return(list(S=Strings,pos=pos,names=Names))
  }

##------------------------------------------------------------------------------
## parsing an equation in single items
##------------------------------------------------------------------------------
  ParseEquation <- function(LR, pref = "EQ",
      names = paste(pref, 1:length(LR), sep = "")) {

    parsed <- NULL
    for (i in 1:nrow(LR)) {
      ifelse(LR[i, 3] %in% c("=", "<"), ct <- c(-1, 1), ct <- c(1, -1))
      eq <- gsub("@@", ",", LR[i, ])
      eq <- gsub("@", ",", eq)
      for (ii in 1:2) {
        ss <- unlist(strsplit(eq[ii], "\\+"))
        for (j in 1:length(ss)) {
          if (ss[j]!="") {
            new <- finditem(ss[j], ct = ct[ii])
            if (is.null(new) && createcomp ) {
              ES <- emptystruct(ss[j])
              ES$nr <- nr
              nr <- nr+1
              comp <<- rbind(comp,ES)
              compnames <<- c(compnames,toupper(ss[j]))
              new   <- nrow(comp)
            }
            if (!is.null(new))
              parsed <- rbind(parsed, data.frame(name = names[i],
                      nr = i, new))
          }
        }
      }
    }
    return(parsed)
  }
##------------------------------------------------------------------------------
## parsing a line equation in single items
##------------------------------------------------------------------------------
  ParseLine <- function(String, pref="var",
      names=paste(pref,1:length(String), sep = "")) {

    cost <- NULL
    for (i in 1:length(String)) {
      eq <- String[i]
      if (is.na(eq))
        eq <- "XNA"
      ss <- unlist(strsplit(eq, "\\+"))
      for (j in 1:length(ss))
        cost <- rbind(cost, data.frame(name = names[i],  nr = i,
                                       finditem(ss[j], ct = 1)))
    }
    return(cost)
  }
##------------------------------------------------------------------------------
## checking variable declarations: should be a function of unknowns
##------------------------------------------------------------------------------
  checkvar <- function(col) {
    nvars <- length(varnames)
    for (i in 1:nvars) {
      term <- vars[vars$nr == i, ]
      for (ii in 1:nrow(term)) {
        if (!is.na(term[ii, col])) {
          varunknown[i] <- TRUE
          next()
        }
        if (!is.na(term$var[ii])) {
          if (varunknown[term$var[ii]]) {
             varunknown[i] <- TRUE
             next()
          }
        }
      } # ii

    }  # i
    ii <- which(!varunknown)
    if (length(ii) > 0) {
      print(varnames[ii])
      stop("the above variable(s) are  not a true variable, but a parameter")
    }
    return(varunknown)
  }

#==============================================================================#
#==============================================================================#
#=                          MAIN FUNCTION STATEMENTS                          =#
#==============================================================================#
#==============================================================================#

##-------------------------------------------------------------------------------
##    CLEANING UP TEXT LINES
##-------------------------------------------------------------------------------
  if (verbose)
    print("reading file lines")
  Lines <- readLines(file)
  if (verbose)
    print("cleaning up file lines")
  Lines <- RemoveToEnd(Lines, "!")
  len <- length(Lines)
  Lines <- Remove(Lines, c("{", "}"))
  Lines <- Mergestrings(Lines)
  if (remtabs)
    Lines <- gsub("\t", "", Lines)
  Lines <- gsub(" ", "", Lines)
  Lines <- gsub("<->", "@@", Lines)
  Lines <- gsub("->", "@", Lines)
  Lines <- gsub("\"", "'", Lines)
  Lines <- Lines[nchar(Lines) != 0]
##------------------------------------------------------------------------------
##  FIND HEADING LINES AND READ ELEMENTS
##------------------------------------------------------------------------------
  headline <- grep("#", Lines)
  numhead <- length(headline)
  pars <- NULL
  comp <- NULL
  extern <- NULL
  rate <- NULL
  marker <- NULL
  flows <- NULL
  cost <- NULL
  profit <- NULL
  vars <- NULL
  equis <- NULL
  inequis <- NULL
  reac <- NULL
  posreac <- NULL
  Type <- "simple"
  createcomp <- FALSE
  parnames <- varnames <- compnames <- externnames <- reacnames <- NULL
  if (length(numhead) == 0)
    stop("there are no headings in file")
## --------------------------   Parameters    ----------------------------------
  if (verbose)
    print("reading parameters")
  paramsec <- findsection("PARA")
  if (length(paramsec) > 1) {
    parnames <- Splitleftright(Lines[paramsec[1]:paramsec[2]], "=")[, 1]
    pars <- SplitNameEquation(Lines[paramsec[1]:paramsec[2]])
    errors <- c(which(is.na(pars$val)), which(pars$nr ==
          pars$par1), which(pars$nr == pars$par3), which(pars$nr ==
          pars$par4))
   if (length(errors) > 0)
     stop(paste("parameter values not defined", parnames[errors]))

  }
## --------------------------   Components    ----------------------------------
  createcomp <- TRUE
  nr <- 0
  if (verbose)
    print("reading components")
  compsec <- findsection(c("STOC", "COMP", "DECI", "STAT", "UNKN"))
  if (length(compsec) > 1) {
    createcomp <- FALSE
    compnames <- Splitleftright(Lines[compsec[1]:compsec[2]], "=")[, 1]
    comp <- SplitNameEquation(Lines[compsec[1]:compsec[2]])
  }
  compnames <- toupper(compnames)
## --------------------------    Externals    ----------------------------------
  if (verbose)
    print("reading externals")
  externsec <- findsection("EXTE")
  if (length(externsec) > 1) {
    externnames <- Splitleftright(Lines[externsec[1]:externsec[2]], "=")[, 1]
    extern <- SplitNameEquation(Lines[externsec[1]:externsec[2]])
  }
  Externnames <- toupper(externnames)
## --------------------------      Rates      ----------------------------------
  if (verbose)
    print("reading rates")
  ratesec <- findsection("RATE")
  if (length(ratesec) > 1) {
    rate <- SplitNameEquation(Lines[ratesec[1]:ratesec[2]])
    Ratenames <- toupper(rate$name)
    if (sum(!Ratenames %in% compnames) > 0)
      stop(paste("error: a rate is defined, but no corresponding component"),
              rate$name[which(!Ratenames %in% compnames)])
    ii <- which(!compnames %in% Ratenames)
    if (length(ii) > 0) {
      rate <- rbind(rate, emptystruct(comp$name[ii]))
      Ratenames <- c(Ratenames, compnames[ii])
    }
    rate <- rate[match(compnames, Ratenames), ]
    rate$nr <- 1:nrow(rate)
  }
## ----------------------  Markers, e.g isotopes   -----------------------------
  marksec <- findsection("MARK")
  if (length(marksec) > 1) {
    if (verbose)
      print("reading markers")
    marker <- SplitNameEquation(Lines[marksec[1]:marksec[2]])
    Markernames <- toupper(marker$name)
    compext <- c(compnames, Externnames)
    if (sum(!Markernames %in% compext) > 0)
      stop(paste("error: a marker is defined, but no corresponding component or external"),
           marker$name[which(!Markernames %in% compext)])
    ii <- which(!compext %in% Markernames)
    if (length(ii) > 0) {
      marker <- rbind(marker, emptystruct(compext[ii]))
      Markernames <- c(Markernames, compext[ii])
    }
    marker <- marker[match(compext, Markernames), ]
    marker$nr <- 1:nrow(marker)
  }
## --------------------------     Flows       ----------------------------------
  if (verbose)
    print("reading flows")
  flowsec <- findsection("FLOW")
  if (length(flowsec) > 1)  {
    ff <- Splitleftright(Lines[flowsec[1]:flowsec[2]],":")
    flownames<- ff[,1]
    FF <- Splitflows(ff[,2],createcomp)
    flows <- FF$Flows
    posreac <- FF$pos   # flows that should be positive...
    flows$name <- flownames
    Type <- "flow"
  }
  ifelse(length(flowsec) > 1, colunknown <- 9, colunknown <- 10)
## ---------------------------   Reactions    ----------------------------------
  if (verbose)
    print("reading reactions")
  reacsec <- findsection(c("REAC"))
  if (length(reacsec) > 1) {
    if(Type == "flow") stop("cannot specify flows AND reactions")
    Type <- "reaction"
    colunknown <- 12
    reacstrings <- Lines[reacsec[1]:reacsec[2]]
    reacst <- CleanReaction(reacstrings)
    reacnames <- reacst$names
    equiLR <- SplitEquation(reacst$S)
    reac <- ParseEquation(equiLR, names = reacnames[,1])

    reac$val <- -1*reac$val
    posreac <-reacst$pos  # reactions that should be positive...
  }
## Check for duplicated names
  allnames <- c(parnames, varnames, externnames, compnames,reacnames[,1])
  ii <- duplicated(allnames)
  if (sum(ii) > 0) {
    print(allnames[ii])
    stop("cannot proceed: the above name(s) are declared twice:")
  }
## --------------------------    Variables    ----------------------------------
  if (verbose)
    print("reading variables")
  varsec <- findsection("VARI")
  if (length(varsec) > 1) {
    varnames <- Splitleftright(Lines[varsec[1]:varsec[2]], "=")[, 1]
    varunknown <- rep(FALSE, length(varnames))
    vars <- SplitNameEquation(Lines[varsec[1]:varsec[2]])
    if (checkLinear)
      varunknown <- checkvar(col = colunknown)
  }
## -----------------------    Cost to minimise    ------------------------------
  if (verbose)
    print("reading cost")
  costsec <- findsection(c("COST", "MINI"))
  if (length(costsec) > 1) {
    coststr <- Lines[costsec[1]:costsec[2]]
    if (length(grep(":", coststr, fixed = TRUE)) > 0) {
      eqnames <- Splitleftright(coststr, ":")
      cost    <- ParseLine(eqnames[,2], names = eqnames[,1])
    } else  cost <- ParseLine(coststr, "cost")
  }
## -----------------------    Cost to minimise    ------------------------------
  if (verbose)
    print("reading profit")
  profitsec <- findsection(c("PROF", "MAXI"))
  if (length(profitsec) > 1)  {
    profstr <- Lines[profitsec[1]:profitsec[2]]
    if (length(grep(":", profstr, fixed = TRUE)) > 0) {
      eqnames <- Splitleftright(profstr, ":")
      profit  <- ParseLine(eqnames[,2], names = eqnames[,1])
    } else profit <- ParseLine(profstr, "cost")
  }
  createcomp <- FALSE
## --------------------------  Equalities (=)  ---------------------------------
  if (verbose)
    print("reading equalities")
  equisec <- findsection("EQUA")
  if (length(equisec) > 1) {
    eqstrings <- Lines[equisec[1]:equisec[2]]
    if (length(grep(":", eqstrings, fixed = TRUE)) > 0) {
      eqnames <- Splitleftright(eqstrings, ":")
      equiLR <- SplitEquation(eqnames[, 2])
      equis <- ParseEquation(equiLR, names = eqnames[, 1])
    }
    else equis <- ParseEquation(SplitEquation(eqstrings), "eq")
  }
## -------------------------- Inequalities (>) ---------------------------------
  if (verbose)
    print("reading inequalities")
  inequisec <- findsection(c("CONS", "INEQ"))
  if (length(inequisec) > 1) {
    ineq <- CleanInequality(Lines[inequisec[1]:inequisec[2]])
    if (length(grep(":", ineq, fixed = TRUE)) > 0) {
      eqnames <- Splitleftright(ineq, ":")
      equiLR <- SplitEquation(eqnames[, 2])
      inequis <- ParseEquation(equiLR, names = eqnames[, 1])
    }
    else inequis <- ParseEquation(SplitEquation(ineq), "ineq")
  }
  if (length(flowsec) > 1) {
    flows$fname <- gsub(",", "->", flows$fname)
    flows$name <- gsub("@@", "<->", flows$name)
    flows$name <- gsub("@", "->", flows$name)
    Type <- "flow"
  }
  res <- list(file = file, pars = pars, comp = comp, rate = rate,
      extern = extern, flows = flows, vars = vars, cost = cost,
      profit = profit, equations = equis, constraints = inequis,
      reactions = reac, posreac=posreac, marker = marker,
      parnames = parnames, varnames = varnames,
      compnames = compnames, externnames = externnames, Type = Type)
  class(res) <- "liminput"
  return(res)
}

################################################################################
###########            CREATING LIM INPUT MATRICES                 ############# 
################################################################################

Setup.limfile <- function (file, verbose=TRUE, ...) {

  Inverse.input <- Read(file, verbose=verbose,...)
  if (is.null(Inverse.input))
    stop("error: cannot setup lim: file does not exist")
  Setup.liminput(Inverse.input)
}
################################################################################

Setup.liminput <- function (liminput, ...)  {

  calcterm <- function(term) {
    calc <- term$val
    if (!is.na(term$par1))
      calc <- calc * parval[term$par1]
    if (!is.na(term$par2))
      calc <- calc * parval[term$par2]
    if (!is.na(term$par3))
      calc <- calc * parval[term$par3]
    if (!is.na(term$par4))
      calc <- calc * parval[term$par4]
    if (!is.na(term$comp))
      calc <- calc * compval[term$comp]
    if (!is.na(term$external))
      calc <- calc * externval[term$extern]
    return(calc)
  }
  calcvalues <- function(pars, npar) {
    if (npar == 0)
      return(NULL)
    parval <- rep(0, length = npar)
    for (i in 1:npar) {
      term <- pars[pars$nr == i, ]
      for (ii in 1:nrow(term)) {
        parval[i] <- parval[i] + calcterm(term[ii, ])
      }
    }
    return(parval)
  }
  calcvariable <- function(vars, ii, vec = rep(0, nunkn + 1)) {
    term <- vars[vars$nr == ii, ]
    for (i in 1:nrow(term)) {
      ct <- calcterm(term[i, ])
      iv <- term$var[i]
      ifl <- term[i, unkncol]
      if (!is.na(iv))
        vec[1:nunkn] <- vec[1:nunkn] + varmat[iv, 1:nunkn] * ct
        if (!is.na(ifl)) {
           vec[ifl] <- vec[ifl] + ct
        }
        if (is.na(iv) & is.na(ifl))
          vec[nunkn + 1] <- vec[nunkn + 1] - ct
    }
    return(vec)
  }

## begin main body...
  ncomp <- length(liminput$compnames)
  nextern <- length(liminput$externnames)
  ifelse (is.null(liminput$rate), nrate <- 0, nrate <- ncomp)
  ifelse (is.null(liminput$marker), nmarker <- 0, nmarker <- nrow(liminput$marker))
  npar <- length(liminput$parnames)
  nvars <- length(liminput$varnames)
  flows <- liminput$flows
  reactions <- liminput$reactions
  ispos <- TRUE
  ifelse (is.null(reactions),nreac <-0, nreac <- max(reactions$nr))
  nflow <- nrow(flows)
   neqs <- nineqs <- 0
  if (is.null(nflow))
    nflow <- 0
  if(liminput$Type == "flow")     Unknown <- "F" else
  if(liminput$Type == "reaction") Unknown <- "R" else Unknown <- "S"
  if (Unknown == "F") {
    nunkn <- nflow
    if (nunkn == 0)
      stop("cannot create inverse matrices: number of flows=0")
    unkncol <- 9
    Unknownnames <- liminput$flows$name
    if (any(!liminput$posreac))
      ispos<-FALSE
  } else if (Unknown == "R") {
    nunkn <- max(reactions$nr)
    if (nunkn == 0)
      stop("cannot create inverse matrices: number of unknowns =0")
    unkncol <- 12
    Unknownnames <- unique(liminput$reactions$name)
    if (any(!liminput$posreac))
      ispos<-FALSE
    } else {
        nunkn <- ncomp
        unkncol <- 10
        if (nunkn == 0)
            stop("cannot create inverse matrices: number of components=0")
        Unknownnames <- liminput$compnames
    }
    parval <- rep(NA, npar)
    compval <- rep(NA, ncomp)
    if (!is.null(liminput$pars)) {
        while (any(is.na(parval))) parval <- calcvalues(liminput$pars,
            npar)
    }
    ifelse(Unknown%in% c("F","R"), compval <- calcvalues(liminput$comp,
        ncomp), compval <- rep(1, ncomp))
    externval <- calcvalues(liminput$extern, nextern)
    rateval   <- calcvalues(liminput$rate, nrate)
    markerval <- calcvalues(liminput$marker, nmarker)
    if (length(rateval) == 0)
        rateval <- rep(0, ncomp)
    rateval[is.na(rateval)] <- 0
    VarA <- VarB <- NULL
    if (!is.null(liminput$vars)) {
        varmat <- matrix(data = 0, nrow = nvars, ncol = nunkn + 1)
        for (i in 1:nvars) varmat[i, ] <- calcvariable(liminput$vars,
            i)
        VarA <- varmat[, 1:nunkn]
        VarB <- varmat[, nunkn + 1]
        if (nvars > 1)
            ii <- which(rowSums(abs(VarA)) == 0)
        else ii <- NULL
        if (length(ii) > 0) {
            print(liminput$varnames[ii])
            stop(paste("cannot proceed: the above variable(s) are  not a true variable, but a parameter:"))
        }
    }
    eqmat <- NULL
    if (!is.null(liminput$equations)) {
        neqs <- max(liminput$equations$nr)
        eqmat <- matrix(data = 0, nrow = neqs, ncol = nunkn + 1)
        for (i in 1:neqs) eqmat[i, ] <- calcvariable(liminput$equations,
            i)
    }       
    ineqmat <- NULL
    if (!is.null(liminput$constraints)) {
        nineqs <- max(liminput$constraints$nr)
        ineqmat <- matrix(data = 0, nrow = nineqs, ncol = nunkn + 1)
        for (i in 1:nineqs) ineqmat[i, ] <- calcvariable(liminput$constraints,
            i)
    }
    cost <- NULL
    if (!is.null(liminput$cost)) {
        ncost <- max(liminput$cost$nr)
        costmat <- matrix(data = 0, nrow = ncost, ncol = nunkn + 1)
        for (i in 1:ncost) costmat[i, ] <- calcvariable(liminput$cost,
            i)
        cost <- costmat[, 1:nunkn]
    }
    profit <- NULL
    if (!is.null(liminput$profit)) {
        nprofit <- max(liminput$profit$nr)
        profitmat <- matrix(data = 0, nrow = nprofit, ncol = nunkn + 1)
        for (i in 1:nprofit) profitmat[i, ] <- calcvariable(liminput$profit,
            i)
        profit <- profitmat[, 1:nunkn]
    }
    A <- B <- G <- H <- NULL
    if (nflow > 0) {
        A <- matrix(data = 0, nrow = ncomp, ncol = nunkn)
        for (i in 1:nflow) {
            if (flows$from[i] > 0)
                A[flows$from[i], i] <- A[flows$from[i], i] -
                  1
            if (flows$to[i] > 0)
                A[flows$to[i], i] <- A[flows$to[i], i] + 1
        }
        B <- rateval
    } else 
    if (nreac>0)
    {
        A <- matrix(data = 0, nrow = ncomp, ncol = nunkn)
# KARLINE: HERE IS THE NEW SECTION
        for (i in 1:nrow(reactions)) {
            st <- reactions[i,"comp"]
            rr <- reactions[i,"nr"]
            A[st,rr] <- A[st,rr] + reactions[i,"val"]
        }
        B <- rateval
    }
    if (!is.null(eqmat)) {
        A <- rbind(A, eqmat[, 1:nunkn])
        B <- c(B, eqmat[, nunkn + 1])
    }
    if (nflow > 0 || nreac>0 ) {
        # CHECK A-matrix for possible EXTERNALS
        Aext <- cbind(A[1:ncomp,],rateval)
        Aext[Aext!=0] <- 1
        RS <- rowSums(Aext)
        if (any (RS<=1))
         {warning("following component(s) is possibly an external:")
          warning(paste(" ",liminput$compnames[which(RS<=1)]))}
        # not all reaction rates are positive (i.e. some can go both ways <->)
        G <- diag(nunkn)
        G<- G[liminput$posreac,]
        H <- rep(0,sum(liminput$posreac))
    }
    if (!is.null(ineqmat)) {
        G <- rbind(ineqmat[, 1:nunkn], G)
        H <- c(ineqmat[, nunkn + 1], H)
    }
    if (nflow > 0) {
        nel <- ncomp + nextern
        elnames <- c(liminput$compnames, liminput$externnames)
        Flowmatrix <- matrix(data = 0, nrow = nel, ncol = nel, 
            dimnames = list(elnames, elnames))
        flowsft <- as.matrix(flows[, 1:2])
        flowsft[flowsft < 0] <- -flowsft[flowsft < 0] + ncomp
        Flowmatrix[flowsft] <- 1:nflow
    }
    else Flowmatrix <- NULL
    if (nmarker == 0) {
        markers = NULL
    }
    else markers = data.frame(name = c(liminput$compnames,
        liminput$externnames), val = markerval)
    res <- list(file = liminput$file, NUnknowns = nunkn,
        NEquations = neqs, NConstraints = nineqs, NComponents = ncomp,
        NExternal = nrow(liminput$extern), NVariables = nvars,
        A = A, B = B, G = G, H = H, Cost = cost, Profit = profit,
        Flowmatrix = Flowmatrix, VarA = VarA, VarB = VarB, 
        Parameters = data.frame(name = liminput$parnames, val = parval), 
        Components = data.frame(name = liminput$compnames, val = compval), 
        Externals = data.frame(name = liminput$externnames, val = externval), 
        rates = data.frame(name = liminput$compnames, val = rateval), 
        markers = markers, Variables = liminput$varnames,
        costnames = unique(liminput$cost$name),
        profitnames = unique(liminput$profit$name),
        eqnames = unique(liminput$equations$name), 
        ineqnames = unique(liminput$constraints$name),
        Unknowns = Unknownnames, ispos = ispos)
    class(res) <- "lim"
    return(res)
}
