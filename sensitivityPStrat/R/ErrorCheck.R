.CheckEmptyPrincipalStratum <- function(empty.principal.stratum, ...) {
  ## empty.principal.stratum must be
  ## 1. not missing
  ## 2. not NULL
  ## 3. a 2 element vector
  ## 4. elements must be diffrent values
  ## 5. cannot be NA

  if(missing(empty.principal.stratum))
    return("'empty.principal.stratum' argument must be supplied")

  ErrMsg <- NULL
  if(is.null(empty.principal.stratum) ||
     length(empty.principal.stratum) != 2L)
    ErrMsg <- c(ErrMsg,
                "'empty.principal.stratum' argument must be a two element vector")
  else if(empty.principal.stratum[1] == empty.principal.stratum[2])
    ErrMsg <- c(ErrMsg,
                "'empty.principal.stratum' argument cannot have the same values for both elements")
  

  if(length(empty.principal.stratum) &&
     any(is.na(empty.principal.stratum)))
    ErrMsg <- c(ErrMsg,
                "'empty.principal.stratum' may not contain a NA")

  return(ErrMsg)
}


.CheckSelection <- function(selection, s, empty.principal.stratum, ...) {
  ## selection must be
  ## 1. not missing
  ## 2. a single element vector
  ## 3. not NA
  ## 4. contained in empty.principal.stratum if it exists or contained in
  ##     s if it exists
  
  if(missing(selection))
    return("'selection' argument must be supplied")

  ErrMsg <- NULL
  if(length(selection) > 0L && any(is.na(selection)))
    ErrMsg <- c(ErrMsg,
                "'selection' may not be NA")

  if(length(selection) != 1L)
    ErrMsg <- c(ErrMsg,
                "'selection' argument must be a single element vector")

  if(length(selection) > 0L) {
    if(!missing(empty.principal.stratum)) {
      if(!all((selection %in% empty.principal.stratum) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principal.stratum'")

    } else if(!missing(s)) {
      if(!all((selection %in% s) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principal.stratum'")
      
    }
  }
  
  return(ErrMsg)
}

.CheckGroupings <- function(groupings, ...) {
  ## groupings must be
  ## 1. not missing
  ## 2. cannot contain NAs
  
  if(missing(groupings))
    return("'groupings' argument must be supplied")

  ErrMsg <- NULL
  if(is.null(groupings) || length(groupings) != 2L)
    ErrMsg <- c(ErrMsg,
                "'groupings' argument must be a two element vector")

  if(length(groupings) && any(is.na(groupings)))
    ErrMsg <- c(ErrMsg,
                "'groupings' may not contain a NA")
  
  return(ErrMsg)
}

.CheckTrigger <- function(trigger, d) {
  ## trigger must be
  ## 1. not missing
  ## 2. cannot contain NAs

  if(missing(trigger))
    return("'trigger' argument must be supplied")

  ErrMsg <- NULL
  if(is.null(trigger) || length(trigger) < 1L)
    ErrMsg <- c(ErrMsg,
                "'trigger' argument must be a single element vector")

  if(length(trigger) && any(is.na(trigger)))
    ErrMsg <- c(ErrMsg,
                "'trigger' may be NA")

  if(length(trigger) > 0L && !missing(d) &&
     !all((trigger %in% d) | is.na(trigger)))
        ErrMsg <- c(ErrMsg,
                    "'trigger' value does not appear in given values of 'd'")
  return(ErrMsg)
}

.CheckTau <- function(tau) {
  ## for argument tau
  ## 1. it must not be missing
  ## 2. it must be numeric
  ## 3. it must have at least 1 element
  ## 4. it can have no more then 2 elements
  ## 5. no elements may be NA
  
  if(missing(tau))
    return("'tau' argument must be supplied")

  ErrMsg <- NULL
  if(!is.numeric(tau))
    ErrMsg <- c(ErrMsg, "'tau' must be a numeric")

  if(length(tau) < 1 || length(tau) > 2)
    ErrMsg <- c(ErrMsg, "'tau' must be a vector of 1 or 2 elements")

  if(any(is.na(tau)))
    ErrMsg <- c(ErrMsg, "'tau' may not contain any NA values")

  return(ErrMsg)
}

.CheckPhiPiPsi <- function(phi, Pi, psi, p0, p1) {
  ## for arguments phi, Pi, and psi
  ## 1. Only one of these argument may be specified
  ## 2. the specified argument may not be missing
  ## 3. the specified argument must have 1 or more values
  ## 4. Pi must be in range of 0 to 1
  ## 5. phi must be in range of 0 to 1
  ## 6. phi, Pi, psi must be numeric

  if(sum(c(missing(phi), missing(Pi), missing(psi))) == 0)
    return("one of arguments 'phi', 'Pi', or 'psi' must be specified")

  ErrMsg <- NULL
  if(sum(c(missing(phi), missing(Pi), missing(psi)) > 1))
    ErrMsg <- c(ErrMsg,
                "only one of the arguments 'phi', 'Pi', or 'psi' may be specified")

  if(!missing(phi) && !is.null(phi)) {
    if(length(phi) < 1)
      ErrMsg <- c(ErrMsg, "'phi' must have 1 or more elements")
                  
    if(any(is.na(phi)))
      ErrMsg <- c(ErrMsg,
                  "'phi' may not contain NA values")

    if(!is.numeric(phi))
      ErrMsg <- c(ErrMsg, "'phi' must be a numeric vector")
    
    if(!missing(p0) && !missing(p1) && any(phi > min(p0, p1)/p1 | phi < max(0, p0 + p1 - 1)/p1))
        ErrMsg <- c(ErrMsg, "values of 'phi' must be constitant with data provided 'phi > min(p0, p1)/p1, phi < max(0, p0 + p1 - 1)/p1'")
    else if(any(phi > 1 | phi < 0))
        ErrMsg <- c(ErrMsg, "values of 'phi' must be between 0 and 1")

  }
    
  if(!missing(Pi) && !is.null(Pi)) {
    if(length(Pi) < 1)
      ErrMsg <- c(ErrMsg, "'Pi' must have 1 or more elements")
                  
    if(any(is.na(Pi)))
      ErrMsg <- c(ErrMsg,
                  "'Pi' may not contain NA values")

    if(!is.numeric(Pi))
      ErrMsg <- c(ErrMsg, "'Pi' must be a numeric vector")

    if(!missing(p0) && !missing(p1) && any(Pi > min(p0,p1) | Pi < max(0,p0 + p1 - 1)))
      ErrMsg <- c(ErrMsg, "values of 'Pi' must be constitant with data provided 'Pi > min(p0,p1), Pi < max(0,p0 + p1 - 1)'")
    else if(any(Pi > 1 | Pi < 0))
      ErrMsg <- c(ErrMsg, "values of 'Pi' between 0 and 1")
  }

  if(!missing(psi) && !is.null(psi)) {
    if(length(psi) < 1)
      ErrMsg <- c(ErrMsg, "'psi' must have 1 or more elements")
                  
    if(any(is.na(psi)))
      ErrMsg <- c(ErrMsg,
                  "'psi' may not contain NA values")

    if(!is.numeric(psi))
      ErrMsg <- c(ErrMsg, "'psi' must be a numeric vector")
  }

  return(ErrMsg)
}
  
.CheckLength <- function(z, s, d, y, v) {
  ## for vectors z, s, d, and y
  ## 1. must be same length

  vectorLengths <- c(if(!missing(z)) length(z),
                     if(!missing(s)) length(s),
                     if(!missing(d)) length(d),
                     if(!missing(y)) length(y),
                     if(!missing(v)) length(v))
  
  vectorMissing <- c(missing(z), missing(s), missing(d), missing(y), missing(v))
  
  if(any(!vectorMissing) && length(unique(vectorLengths)) > 1)
    return(sprintf("Vectors %s are of unequal lengths", paste(c('z','s','d','y')[!vectorMissing], collapse=',')))

  return(NULL)
}

.CheckZ <- function(z, groupings, na.rm, ...) {
  ## z must be
  ## 1. not missing
  ## 2. cannot contain NAs unless na.rm == TRUE
  ## 3. all z values must match the values of groupings unless groupings is
  ##      missing

  if(missing(z))
    return("'z' argument must be supplied")

  ErrMsg <- NULL
  if(any(is.na(z)) && !na.rm)
    ErrMsg <- c(ErrMsg,
                "'z' cannot contain any NA values")

  if(!missing(groupings) && !all(z %in% groupings | is.na(z)))
    ErrMsg <- c(ErrMsg,
                "All values of 'z' must match one of the two values in 'groupings'")
  return(ErrMsg)
}

.CheckS <- function(s, empty.principal.stratum, na.rm, ...) {
  ## s must be
  ## 1. not missing
  ## 2. cannot contain NAs unless na.rm == TRUE
  ## 3. all s values must match the values of empty.principal.stratum if exists

  if(missing(s))
    return("'s' argument must be supplied")

  ErrMsg <- NULL
  if(any(is.na(s)) && !na.rm)
    ErrMsg <- c(ErrMsg,
                "argument 's' cannot contain any NA values")

  if(!missing(empty.principal.stratum) &&
     !all(s %in% empty.principal.stratum | is.na(s)))
    ErrMsg <- c(ErrMsg,
                "All values of 's' must match one of the two values in 'empty.principal.stratum'")

  return(ErrMsg)
}



.CheckY <- function(y, s, selection, na.rm, ...) {
  ## y must be
  ## 1. not missing
  ## 2. cannot be NA is s is selected.

  if(missing(y))
    return("'y' argument must be supplied")

  ErrMsg <- NULL
  if(!missing(selection) && !missing(s) && !na.rm &&
     length(s) == length(y) && any(selection %in% s) &&
     any(s %in% selection & !is.na(s) & is.na(y)))
    ErrMsg <- c(ErrMsg,
                sprintf("argument 'y' cannont contain a NA value if the corresponding 's' is %s",
                        paste(selection, collapse=",")))
  return(ErrMsg)
}

.CheckD <- function(d, s, selection, na.rm) {
  ## d must be
  ## 1. not missing
  ## 2. cannot be NA is s is selected.

  if(missing(d))
    return("'d' argument must be supplied")

  ErrMsg <- NULL
  if(!missing(selection) && !missing(s) && !na.rm &&
     length(s) == length(d) && any(selection %in% s) &&
     any(s %in% selection & !is.na(s) & is.na(d)))
    ErrMsg <- c(ErrMsg,
                sprintf("argument 'd' cannot contain a NA value if the corresponding 's' is %s", paste(selection, collapse=",")))

  return(ErrMsg)
}

.CheckV <- function(v, followup.time, na.rm) {
  ## v must be
  ## 1. present if followup.time is non-NULL or non-missing
  ## 2. cannot be pressent if followup.time if NULL
  ## 2. a numeric vector
  ## 3. cannot contain NA unless na.rm == TRUE
  if(missing(followup.time) || is.null(followup.time)) {
    if(!missing(v))
      return("'v' argument cannot be supplied without a 'followup.time' argument")
  } else {
    if(missing(v))
      return("'followup.time' argument cannot be supplied without a 'v' argument")

  }
    
}

.CheckCi <- function(ci, ci.type) {
  ErrMsg <- NULL
  if(!missing(ci.type) && ((missing(ci) && length(ci.type == 1L)) ||
                            (length(ci.type) != length(ci))))
    ErrMsg <- c(ErrMsg, "if 'ci.type' is specified it must be the same length as 'ci'")
  if(!missing(ci.type) && any(fail <- pmatch(ci.type, c('upper','lower','twoSided'), nomatch=0, duplicates.ok=TRUE) == 0))
    ErrMsg <- c(ErrMsg,
                paste("'ci.type' values ",
                      paste(dQuote(ci.type[which(fail)]), collapse=', '),
                      " are not valid 'ci.type' values"))
  return(ErrMsg)
}
