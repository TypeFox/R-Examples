# functions for variable filtering

#
# define generic functions
#
Vinclude <- function( this, ... ) {
  UseMethod("Vinclude", this)
}
Vexclude <- function( this, ... ) {
  UseMethod("Vexclude", this)
}
Vtime <- function( this, ... ) {
  UseMethod("Vtime", this)
}

#
# define functions
#
Vinclude.Connect <- function(this, name.patterns=NULL, label.patterns=NULL, and=TRUE, ...){
  # start with a set of variables provided in the variables argument, or all variables in the query
  # variables that match the conditions provided in name, type will be returned as a vector
  if (is.null(name.patterns) & is.null(label.patterns)){
    # no arguments given, return all non-archived variable system names
    return(GetNames(this, kind="system"))
  } else {
    # one or more arguments given, apply filtering
    # find the matching variables from the arguments
    return(VariableMatcher(this, name.patterns, label.patterns, and))
  }
}


Vexclude.Connect <- function(this, variable.names=NULL, name.patterns=NULL, label.patterns=NULL, and=TRUE, ...){
  if (is.null(variable.names)){
    # no variable names were provided, so assume that we start with all variables as eligible
    variable.names <- GetNames(this, kind="system")
  } 
  
  # find variables matching patterns, then omit the matches
  excluded.variables <- VariableMatcher(this, name.patterns, label.patterns, and)
  return(variable.names[ !(variable.names %in% excluded.variables) ])
}


Vtime.Connect <- function(this, variable.names, domains, ...){
  # Given an array of variable names, this appends periods to them.
  # Mixes of time range and time point variables are handled automatically
  # For example, a time range variable will not have "Focal Point" appended to it.
  
  variable.names.all <- GetNames(this, kind="system")
  idx.valid <- variable.names %in% variable.names.all
  if (any(!idx.valid)){
    stop("Invalid variable(s): ", paste(variable.names[!idx.valid], collapse=", "))
  }
  
  separator = "$" # this will be inserted between the variable name and the time domain
  
  # set time domains labels
  # note that time independent variables will not have a domain appended
  time.range.label <- "Time Range"
  time.point.label <- "Time Point"
  time.indep.label <- "Time Independent"
  
  valid.time.ranges.abbr  <- this$timeRanges$ABBREVIATION
  valid.time.ranges.short <- this$timeRanges$SHORT_NAME
  valid.time.points.abbr  <- this$timePoints$ABBREVIATION
  valid.time.points.short <- this$timePoints$SHORT_NAME
  idx.valid <- domains %in% c(valid.time.ranges.abbr, valid.time.ranges.short, valid.time.points.abbr, valid.time.points.short)
  if (any(!idx.valid)){
    stop("Invalid time range(s) or time point(s): ", paste(domains[!idx.valid], collapse=", "))
  }
  
  # break variables into three groups: time range, time point, and time independent
  # build a list mapping variable names to periods
  names.to.domain.type.list <- as.list(this$variables$TIME_DOMAIN_TYPE)
  names(names.to.domain.type.list) <- this$variables$NAME
  
  # extract the domain type for each variable
  domains.types.for.variables <- names.to.domain.type.list[variable.names]
  
  # build indices of variables from each type of time period
  idx.timerange <- domains.types.for.variables == time.range.label
  idx.timepoint <- domains.types.for.variables == time.point.label
  idx.timeindep <- domains.types.for.variables == time.indep.label
  
  # build indices of domain types for each time domain
  idx.domains.timerange <- domains %in% c(valid.time.ranges.abbr, valid.time.ranges.short)
  idx.domains.timepoint <- !idx.domains.timerange # timerange is not timepoint
  
  # test for missing domains
  if (any(idx.timerange) & !any(idx.domains.timerange)){
    warning(sum(idx.timerange), " time range variables do not have time range domains and will not be returned.")
  }
  if (any(idx.timepoint) & !any(idx.domains.timepoint)){
    warning(sum(idx.timepoint), " time point variables do not have time point domains and will not be returned.")
  }
  
  # initialize list of output variables, start with time independent
  if (any(idx.timeindep)){
    variables.out <- variable.names[idx.timeindep]  
  } else {
    variables.out <- c() # initialize to empty
  }
  
  # loop for each time domain
  for (domain in domains[idx.domains.timerange]){
    newvars <- paste(variable.names[idx.timerange], rep(domain, sum(idx.timerange)), sep=separator)
    variables.out <- c(variables.out, newvars)
  }
  for (domain in domains[idx.domains.timepoint]){
    newvars <- paste(variable.names[idx.timepoint], rep(domain, sum(idx.timepoint)), sep=separator)
    variables.out <- c(variables.out, newvars)
  }
  return(unique(variables.out))
}


VariableMatcher <- function(conn, variable.patterns, labels, and){
  # this function is called by Including and Excluding
  # Given a vector of variables in the variables argument, this
  # returns a subset vector that match the given criteria
  
  matchlist <- list()
  variable.names.system  <- GetNames(conn, kind="system")
  #variable.names.display <- GetNames(conn, kind="display")

  # match variables argument against all variables
  if (is.null(variable.patterns)){
    # no name filter argument, so select all variables
    matchlist[[1]] <- rep(TRUE, length(variable.names.system))
  } else {
    # a name filter argument was given, search for variable pattern matches
    matchlist[[1]] <- GrepLoop(variable.patterns, variable.names.system, boolean=TRUE)
  }
  
  # match labels against all variables in data frame
  if (!is.null(labels)){
    matchlist[[ length(matchlist)+1 ]] <- VariableLabelMatcher(conn, variable.names.system, labels)
  }
  
  # return results
  if (length(matchlist)==1){
    # the simple case is that there is one matching criteria
    matches <- variable.names.system[ matchlist[[1]] ]
  } else {
    # we have more than one matching criteria, so we need to AND or OR together the criteria
    # convert criteria into array of booleans, with a row for each criteria and column for each variable
    match.matrix <- matrix(unlist(matchlist), nrow=length(matchlist), ncol=length(variable.names.system), byrow=TRUE)
    
    # apply matching criteria
    if (and) {
      # all conditions must be met
      idx <- apply(match.matrix, 2, all)
    } else {
      # any condition must be met
      idx <- apply(match.matrix, 2, any)
    }
    matches <- variable.names.system[idx]
  }
  return(matches)
}


VariableLabelMatcher <- function(conn, system.names, labels){
  # return an index of TRUE/FALSE values indicating if the variables in system.names match any
  # labels provided
  all.labels <- conn$labels$LABEL_NAME
  all.ids    <- conn$labels$VARIABLE_ID
  # build a map of IDs to variable system names
  id.to.name.list <- as.list(conn$variables$NAME)
  names(id.to.name.list) <- as.character(conn$variables$ID)
  
  # build a vector of unique variable IDs that match labels
  unique.matching.ids <- as.character(unique(all.ids[GrepLoop(labels, all.labels, boolean=TRUE)]))
  
  # Return an index of variable names
  return(system.names %in% id.to.name.list[unique.matching.ids])
}