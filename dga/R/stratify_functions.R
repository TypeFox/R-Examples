require(chron)
make.strata <-  function(overlaps, dates = NULL, locations = NULL, 
                                demographics = NULL, date.defs = "monthly", loc.defs = NULL, 
                                demog.defs = NULL, start.date = NULL, end.date = NULL){
  #This function creates list intersection counts that go into MSE and counts by data source, aggregated
  # (optionally) by date, location, and demographics. The only required argument is
  # "overlaps", which is individual-level records where a 1 in the ij'th spot in the
  # matrix indicates that individual i was recorded on list j. Inputting additional 
  # information, such as the date, will aggregate by those variables. 
  
  #overlaps: a matrix of list appearance indicatators. Each row is a person, 
  # each column is a list.
  #dates: a chron object of the dates of death. 
  #date.defs: "daily" aggregates only by day, "weekly" aggregates weekly
  # beginning with the first Sunday in the dataset, "monthly", "yearly". 
  # "dayofweek" gives counts by day of the week.
  # If an integer is given, this will aggregate in k-day intervals beginning
  # with the first day that appears. Defaults to monthly.
  #locations: factors indicating geographic location
  #loc.defs: a list partitioning the locations into groups. Defaults to all 
  # unique values treated separately. 
  #demographics: factor indicating demographic groups
  #demog.defs: a list partitioning the demographics into. Defaults to all unique
  # values treated separately.
  
  #Examples:
  # Aggregate only by week: aggregate_overlaps(overlaps, dates, date.def = "weekly)
  # Aggregate by year and location, where locations are not grouped: 
  #       aggregate_overlaps(overlaps, dates, date.def = "yearly", locations)
  # Aggregate by 2 day increments and location, where there are unique location levels
  #       1:10 and locations 1:5 are in group a and locations 6:10 are in group b.
  #       loc.defs <- list("a" = 1:5, "b" = 6:10)
  #       aggregate_overlaps(overlaps, dates, date.def = 2, locations, loc.defs = loc.defs)
  # Aggregate by demographic (sex) only, where sex takes values M, F, A, NA, and U
  #       and we would like to group these as M, F, and other.
  #       demog.defs <- list("M" = "M", "F" = "F", "Other" = c("A", NA, "U"))
  #       aggregate_overlaps(overlaps, demographics, demog.defs = demog.defs)
  
  
  date.cat <- NULL
  loc.cat <- NULL
  demog.cat <- NULL
  
  
  # get date categories
  if(!is.null(dates)){
    if(!is.null(start.date)){
      min.date <- start.date
    }else{
      min.date <- min(dates)
    }
    if(!is.null(end.date)){
      max.date <- end.date
    }else{
      max.date <- max(dates)
    }
    if(date.defs == "yearly"){
      date.cat <- years(dates)
    }
    if(date.defs == "monthly"){
      date.cat <- paste(months(dates), years(dates))
    }
    if(date.defs == "weekly"){
      wd <- weekdays(dates)
      min.sun <- (min.date - 0:6)
      min.sun <- min.sun[weekdays(min.sun) == "Sun"]
      max.sun <- max.date + 0:6
      max.sun <- max.sun[weekdays(max.sun) == "Sun"]
      intervals <- c( seq(0, max.sun - min.sun, 7))
      cuts <- c(min.sun - 7, min.sun + intervals, max.sun + 7)-1
      cut.labels <- cuts + 1
      date.cat <- cut(dates, cuts, labels = cut.labels[1:(length(cut.labels)-1)])
    }
    if(date.defs == "daily"){
      date.cat <- dates
    }
    if(date.defs == "dayofweek"){
      date.cat <- weekdays(dates)
    }
    if(is.numeric(date.defs)){
      intervals <- seq(0, max.date-min.date, date.defs)- 1
      cuts <- c(min.date + intervals, max.date + date.defs)
      cut.labels <- cuts + 1
      date.cat <- cut(dates, cuts, labels = cut.labels[1:(length(cuts)-1)])
    }
  }
  
  if(!is.null(locations)){
    if(is.null(loc.defs)){
      loc.cat <- locations
    }else{
      loc.cat <- nrow(overlaps)
      for(loc in names(loc.defs)){
        loc.cat[is.element(locations, loc.defs[[loc]])] <- loc
      }
    }
  }
  
  if(!is.null(demographics)){
    if(is.null(demog.defs)){
      demog.cat <- demographics
    }else{
      demog.cat <- nrow(overlaps)
      for(demog in names(demog.defs)){
        demog.cat[is.element(demographics, demog.defs[[demog]])] <- demog
      }
    }
  }
  
  strata <- paste(date.cat, loc.cat, demog.cat)
  nlist <- ncol(overlaps)
  if(nlist > 1){
    capt.history <- apply(t(overlaps)*2^((nlist-1):0), 2, sum)
    strata.history <- sapply(split(capt.history, strata), cfunction,nlist = nlist)
    source.counts <- t(sapply(split(overlaps, strata), sfunction))
  }else{
    strata.history <- t(table(strata))
    source.counts <- strata.history
  }
  
  # get by-list region counts
  out <- list()
  out$overlap.counts <- t(strata.history)
  out$source.counts <- source.counts
  return(out)
}

cfunction <- function(x, nlist){
  out <- table(c(x, 0:(2^nlist-1)))-1
}
sfunction <- function(x){
  out <- apply(x, 2, sum)
}

check.strata <- function(strata){  
  flag <- TRUE
  #first check that this has the right stuff in it
  stopifnot(names(strata)==c('overlap.counts', 'source.counts'))
  
  #check to make sure we have a number of lists we can handle
  num.lists <- ncol(strata$source.counts)
  if(num.lists <3){
    print('You have fewer than 3 lists! Are you sure you want to do this?')
    flag <- FALSE
  }
  if(num.lists > 5){
    print('Sorry! We have only pre-computed all of the graphs for three, four, and five lists! 
          Come back later! Or sub-select three, four, or five of your lists to use this package.')
    flag <- FALSE
  }
  
  #check to make sure that none of the lists are empty
  zeroes <- apply(strata$source.counts==0, 1, sum)>0
  if(sum(zeroes)>0){
    print(paste("these strata have empty lists:", rownames(strata$source.counts)[zeroes]))
    flag <- FALSE
  }
  
  small <- apply(strata$source.counts<10, 1, sum)>0
  if(sum(small)>0){
    print(paste("Proceed with caution. These strata have lists with very few (less than 10) records:", 
                rownames(strata$source.counts)[small]))
    flag <- FALSE  
  }
  
  #check to make sure that every list overlaps with at least one other list
  X <- integer.base.b(0:(2^num.lists-1))
  bad.strata <- NULL
  for(i in 1:num.lists){
    inds <- X[,i]==1 & apply(X[,-i], 1, sum)>0
    two.way.overlaps <- apply(strata$overlap.counts[,inds], 1, sum)
    bad.strata <- c(bad.strata, rownames(strata$overlap.counts)[two.way.overlaps==0])
  }
  
  if(length(bad.strata)>0){
    print(paste("These strata have lists that do not intersect with any other lists: ", bad.strata))
    flag <- FALSE
  }
  
  #Eventually, check to see if there are any "islands"
  
  if(flag){
    print("I didn't find any obvious problems. Onward!")
  }
  return(flag)
}