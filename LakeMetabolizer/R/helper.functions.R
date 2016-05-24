#'@title tests data.frame for column names
#'@description 
#'tests \code{data} for data column names
#'
#'@usage
#'has.vars(data, var.names)
#'@param data Object of class data.frame
#'@param var.names A character vector of names to test against \code{data}
#'
#'@return a boolean vector of same length as \code{var.names} 
#'
#'@keywords methods
#'@import plyr
#'@importFrom methods is
#'@author
#'Luke A. Winslow
#'@seealso 
#'\link{get.vars}
#'\link{rmv.vars}
#'@export
has.vars <- function(data, var.names){
  
  if(!is(data, 'data.frame')){
    stop('Data must be of class data.frame')
  }
  
  header = names(data)
  
  has.var = rep(FALSE, times=length(var.names))
  
  for(i in 1:length(var.names)){
    has.var[i] = any(grepl(paste(var.names[i], '[\\w_]?', sep=''), header, ignore.case=TRUE))
  }
  
  return(has.var)  
}

#'@title subsets data.frame according to header names
#'@description 
#'subsets \code{data} according to header names
#'
#'@usage
#'get.vars(data, var.names)
#'@param data Object of class data.frame
#'@param var.names A character vector of names to get from \code{data}
#'
#'@return An object of class data.frame
#'
#'@keywords methods

#'@author
#'Luke A. Winslow
#'@seealso 
#'\link{has.vars}
#'\link{rmv.vars}
#'@export
get.vars <- function(data, var.names){
  
  if(!is(data, 'data.frame')){
    stop('Data must be of class data.frame')
  }
  
  header = names(data)
  
  datetimeI = grepl('datetime', header, ignore.case=TRUE)
  
  if(!any(datetimeI)){
    stop("Can't find datetime column in supplied data.")
  }
  
  varI = rep(FALSE, times=length(var.names))
  
  for(i in 1:length(var.names)){
    varI = varI | grepl(paste(var.names[i], '[\\w_]?', sep=''), header, ignore.case=TRUE)
  }
  
  if(!any(varI)){
    stop("No variable pattern matches:", paste(var.names, collapse=' '))
  }
  
  return(data[, varI | datetimeI])
}

#'@title gets surface water temperatures
#'@description 
#'grabs best available data for surface water temperature
#'
#'@param data Object of class data.frame
#'@param s.range a numeric vector of length=2 with the range for depth measurements to still be considered 'surface'
#'@return An object of class data.frame
#'
#'@keywords methods

#'@author
#'Jordan S. Read
#'@seealso 
#'\link{has.vars}
#'\link{get.vars}
#'\link{rmv.vars}
#'@export
get.Ts <- function(data,s.range=c(0,1)){
  wtr <- get.vars(data,'wtr') 
  datetimeI <- var.indx(wtr,'datetime')
  depths <- c(NA, get.offsets(wtr)) # append NA for datetime col
  depths[depths > s.range[2] | depths < s.range[1]] <- NA
  varI <- which.min(depths)
  
  if (length(varI)==0){
    stop(paste0('no matches for water temperatures within depth range ',s.range[1],' to ',s.range[2]))
  }
  
  # want to use this here, but it returns wtr_0 and wtr_0.5 etc.
  #Ts <- get.vars(wtr,head.nm[varI])
  
  Ts <- wtr[, c(datetimeI,varI)]
  return(Ts)
}

#'@title subsets data.frame according to header names
#'@description 
#'subsets \code{data} according to header names. Excludes all matches to \code{var.name}
#'
#'@usage
#'rmv.vars(data, var.name, ignore.missing=TRUE, ignore.offset=FALSE)
#'@param data Object of class data.frame
#'@param var.name A character vector of names to remove from \code{data}
#'@param ignore.missing Boolean, should an error be thrown if no matching data found
#'@param ignore.offset Should the numerical offset be ignored in the match, (e.g. all \code{wtr} columns removed, or \code{wtr_0} specifically)
#'
#'@return An object of class data.frame
#'
#'@keywords methods

#'@author
#'Luke A. Winslow
#'@seealso 
#'\link{has.vars}
#'\link{get.vars}
#'@export
rmv.vars <- function(data, var.name, ignore.missing=TRUE, ignore.offset=FALSE){
  if(ignore.offset){
    varI = var.indx(data, var.name)
  }else{
    varI = grep(var.name, names(data), ignore.case=TRUE)
  }
  
  if(length(varI) > 0){
    varI = varI * -1
    return(data[, varI])
  }else{
    if(!ignore.missing){
      stop('No variable by that name found')
    }
  }
  
}

#'@title finds matching column names in data.frame
#'@description 
#'returns index of column matches for \code{data} according to header names matches with \code{var.names}.
#'
#'@usage
#'var.indx(data, var.name)
#'@param data Object of class data.frame
#'@param var.name A character vector of names to find matches with \code{data}
#'
#'@return a boolean vector with same length as \code{var.names}
#'
#'@keywords methods

#'@author
#'Luke A. Winslow
#'@seealso 
#'\link{has.vars}
#'\link{get.vars}
#'\link{rmv.vars}
#'@export
var.indx <- function(data, var.name){
  if(length(var.name) != 1){
    stop('var.indx only operates on one variable at a time')
  }
  
  header = names(data)  
  
  indx = grep(paste(var.name, '[\\w_]?', sep=''), header, ignore.case=TRUE)
  
  return(indx)
}



# test edit
# =============================================
# = Function to predict dimensions of merge() =
# =============================================
# RDB 07-May-2014
# Predict merge() dimensions
# Both desired (no duplicates) and expected (how merge will behave) dimensions
# 'all' argument corresponds to the argument of the same name in merge()
pred.merge <- function(x1, x2, all=FALSE){
  common.names <- intersect(names(x1), names(x2))
  
  if(length(common.names)>1){
    fact1 <- do.call(paste, as.list(x1[,common.names])) #'factors' from x1
    fact2 <- do.call(paste, as.list(x2[,common.names])) # factors from x2	
  }else{
    fact1 <- x1[,common.names]
    fact2 <- x2[,common.names]
  }
  
  
  fInt <- intersect(fact1, fact2) # common elements of fact1 and fact2, same as desired.aF (see below)
  
  o1 <- table(fact1[fact1%in%fInt])
  o2 <- table(fact2[fact2%in%fInt])
  out.aF <- sum(o1*o2)
  
  if(all){ # if you used all=TRUE in merge()		
    just1 <- sum(fact1%in%setdiff(fact1, fact2))
    just2 <- sum(fact2%in%setdiff(fact2, fact1))
    
    out.aT <- just1 + just2 + out.aF # This is what merge will give when all=TRUE
    desired.aT <- sum(!duplicated(c(fact1, fact2))) # non-duplicated output if all=TRUE
    
    list(desired=desired.aT, merge=out.aT)
    
  }else{
    out.aF # This is what merge *will* output
    desired.oF <- length(fInt) # Length of desired output, assuming you don't want duplicated join rows
    
    list(desired=desired.oF, merge=out.aF)
  }
}



# =======================
# = round.time Function =
# =======================
# RDB 16May2014
# example: round.time(x, "5 minutes")
# example: round.time(x, "90 min")
# example: round.time(x, "0.5 hours")
# x is a time format, preferably POSIXct, or can be coerced w/ as.POSIXct
# if x needs to be converted to POSIX, define input.format if x currently isn't in a 'standard unambiguous format'
# default output.format=NULL leads to output of class POSIXct, character otherwise
round.time <- function(x, units, input.format=NULL, output.format=NULL){
  # x = head(t.sonde0.na2[,"date"], 20) + 120
  # units = "df.345 min"
  # Check for invalid input classes
  stopifnot(
    is.character(units) & 
      (is.null(output.format) | is.character(output.format)) &
      (is.null(input.format) | is.character(input.format))
  )
  
  # Determine time unit
  unit.choices <- c("sec", "min", "hour", "day")
  choices.or <- paste(unit.choices, collapse="|")
  unit.pattern <- paste(".*(", choices.or, ").*", sep="")
  unit <- gsub(unit.pattern, "\\1", units)
  if(is.na(unit)){stop("not a valid unit, use sec, min, hour, or day")}
  which.choice <- which(unit==unit.choices)
  
  # Determine time interval
  u.time.pattern <- "(?:[0-9]+\\.[0-9]+)|(?:[0-9]+\\.)|(?:\\.[0-9]+)|(?:[0-9]+)"
  u.time.char <- regmatches(units, regexpr(u.time.pattern, units, perl=TRUE))
  u.time <- as.numeric(u.time.char)
  u.time <- ifelse(is.na(u.time), 1, u.time)
  
  unit.cutoff <- switch(unit, sec=60, min=60, hour=24, day=1)
  
  # =========================================================================
  # = Check for invalid input (before slow [attempted] conversion to POSIX) =
  # =========================================================================
  if(sign(u.time)==-1L){
    stop("time interval must be positive")
  }
  # Deal with case where units are 1 second (or less)
  if(unit=="sec" & u.time<=1L){
    return(format.Date(x, format=output.format))
  } else
    
    # Fractional time intervals  convert to smaller unit
    if((trunc(u.time)-u.time)!=0){
      if(sign(u.time)==1L){
        while((trunc(u.time)-u.time)!=0){
          if(unit=="sec"){stop("time interval must be an integer when converted to units of seconds")}
          unit <- unit.choices[which.choice-1]
          which.choice <- which(unit==unit.choices)
          unit.cutoff <- switch(unit, sec=60, min=60, hour=24)
          u.time <- unit.cutoff*u.time
        }
      }else{
        stop("time interval must be positive")
      }
    } else 
      
      # Deal with case where units are days
      if(unit=="day"){
        if(u.time==1){
          return(format.Date(trunc.POSIXt(x + 43200, units = units), format=output.format))
        }else{
          stop("units must be <= 1 day")
        }
      } else 
        
        # Deal w/ cases where time interval is 1 unit
        if(u.time==1){
          unit <- unit.choices[which.choice-1]
          which.choice <- which(unit==unit.choices)
          unit.cutoff <- switch(unit, sec=60, min=60, hour=24)
          u.time <- unit.cutoff
        } 
  
  # Deal with cases where time interval is > 1 of a larger unit
  # Note that this follows up on case where u.time is > 1 and is not an integer
  if(u.time>unit.cutoff){
    u.time <- u.time%%unit.cutoff
    mod.mess <- paste("Rounding to units =", u.time, unit) # may or may not want to make this a warning ...
    warning(mod.mess)
  }
  
  # =============================================
  # = Convert to POSIX, or if can't, give error =
  # =============================================
  if(!"POSIXct"%in%class(x)){
    if(is.null(input.format)){
      x <- as.POSIXct(x)
    }else{
      x <- as.POSIXct(x, format=input.format)
    }
  }
  
  # ===========================================================
  # = Matching units (e.g., min) and unit multiples (e.g., 5) =
  # ===========================================================
  
  which.choice <- which(unit==unit.choices)
  form.unit <- c("%S", "%M", "%H", "%d")[which.choice]
  mult <- as.integer(format.Date(x, format=form.unit))/u.time
  after <- round(mult, 0)*u.time
  # direction <- sign(after-before)
  
  # trunc.unit <- unit.choices[min(which.choice+1, length(unit.choices))]
  trunc.unit <- unit.choices[min(which.choice+1, length(unit.choices))]
  rounded <- trunc.POSIXt(x, trunc.unit) + switch(unit, sec = 1, min = 60, hour = 3600, day = 86400)*after
  if(!is.null(output.format)){
    return(format.Date(rounded, format=output.format))
  }else{
    return(rounded)
  }	
}



# ==================
# = Conquer a List =
# ==================
#RDB
conquerList <- function(x, naming=NULL){
  # If x is not a list, don't bother
  if(!is.list(x) | is.data.frame(x)){return(x)}
  
  s1 <- length(x)
  s2 <- length(x[[1]])
  u1 <- unlist(x, recursive=FALSE)
  stopifnot(length(u1)==s1*s2)
  
  
  if(is.data.frame(x[[1]])){
    single.row <- nrow(x[[1]]) == 1L
  }else{
    single.row <- FALSE
  }
  
  # return value from ldply() if it will work (e.g., if each element of list x contains a row of a data frame)
  if(single.row & is.list(x)){ # the checking for is.list() is a bit redundant with earlier check
    return(cbind(naming, ldply(x)))
  }
  
  #
  s2C <- unlist(lapply(x[[1]], class))
  cqd <- vector("list", s2)
  for(i in 1:s2){
    ti <- seq(i, s1*s2, s2)
    tl <- vector("list", s1)
    for(j in 1:s1){
      tl[[j]] <- u1[[ti[j]]]
    }
    if(is.data.frame(tl[[1]])|!is.list(tl[[1]])){
      if(!is.null(naming)){
        cqd[[i]] <- cbind(naming,ldply(tl))
      }else{
        cqd[[i]] <- ldply(tl)
      }
    }else{
      cqd[[i]] <- llply(tl)
    }
  }
  return(cqd)
}

# =======================================================
# = Simple way of estimating watts entering water layer =
# =======================================================
#RDB
#'@export
watts.in <- function(top, bot, irr, z1perc){
  # top = the top of the layer in meters (e.g., 2)
  # bottom = the bottom of the layer in meters (e.g., 4)
  # irr = PAR, measured in uE
  # z1perc = depth of 1 percent surface light, measured in meters (e.g., 4)
  
  watts <- 0.2174*irr # convert PAR to watts/m^2 (0.2174)
  kd <- log(0.01)/-z1perc # calculate an average kd for the photic zone
  watts*exp(-kd*top) - watts*exp(-kd*bot) # Estimate watts gained as the difference between watts entering at the top and exiting at the bottom
}

# ==================
# = Calculate Mode =
# ==================
#RDB
Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ================================
# = Calculate sampling frequency =
# ================================
calc.freq <- function(datetime){
  freq <- round(Mode(1/diff(date2doy(datetime))))
}

# =================================
# = Checks all inputs for NA vals =
# =================================
# if error=TRUE, will throw an error when 
# NA in input is found
complete.inputs <- function(..., error=FALSE){
  
  inputs = list(...)
  for(i in 1:length(inputs)){
    if(!is.null(inputs[[i]]) && any(is.na(inputs[[i]]))){
      if(error){
        stop('Input ', names(inputs[i]), ' contains NA value')
      }else{
        return(FALSE)
    
      }
    }
  }
  return(TRUE)
}
