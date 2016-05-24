combineCORE <-
function (z) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

# check number of members
# identities may not match even if numbers of members match

 if (length(mem <- unique(sapply(z, function(x) dim(x$weights)[1]))) != 1) 
      stop("merge components have differing numbers of ensemble members")

# check member names

 nam <- data.frame(lapply(z, function(x) sort(dimnames(x$weights)[[1]])))
 if (length(sum(!duplicated(nam))) != 1)
   stop("member names differ among merge components")
 nam <- sort(dimnames(z[[1]]$weights)[[1]])

 if (length(forecastHour <- unique(sapply(z, function(x) attr(x, "forecastHour")))) != 1) 
      stop("merge components differ in forecastHour")
 if (length(initializationTime <- 
    unique(sapply(z, function(x) attr(x, "initializationTime")))) != 1) 
      stop("merge components differ in initialization time")
 if (length(lag <- unique(sapply(z, function(x) x$training$lag))) != 1) 
      stop("merge components differ training lag")

# handle exchangeble designation
# complicated by the fact that the representation is not unique
 
 exch <- lapply(z, function(x) x$exhangeable)
 exch <- lapply( exch, function( x) 
          if (!is.null(x)) sort(as.factor(as.character(x))) else NULL)
 exch <- lapply( exch, function( x,m) {if (is.null(x) || all(x == 1:m)) NULL else x},
                 m = mem)
 if (length(unique(sapply(exch,length))) != 1) 
   stop("merge components differ in exchangeable members")
 if (!all(sapply(exch,is.null))) {
   if (sum(!duplicated(as.data.frame(exch))) != 1)
     stop("merge components differ in exchangeable members")
     exch <- exch[[which(!duplicated(as.data.frame(exch)))]]
 }
 else exch <- NULL

 l <- sapply(z, function(x) length(x$training$days) == 1 && !is.list(x$training$days))
 days <- if (all(l))  unique(sapply( z, function(x) x$training$days)) else NULL

 ttableLIST <- lapply(z, function(x) as.list(x$training$table))
 ttableDATES <- unlist(lapply(ttableLIST, names))
 ttableLENGTH <- unlist(lapply(ttableLIST, function(x) sapply(x, length)))

 dates <- rep(ttableDATES, ttableLENGTH)

 table <- lapply(split(cbind.data.frame( length = unlist(ttableLIST), dates = dates), as.factor(dates)),
                 function(x) unique(x$length))
 if (all(sapply(table, length)) == 1) table <- unlist(table)
   
 nIterLIST <- lapply( z, function(x) x$nIter)
 dateLIST <- lapply( nIterLIST, names)

 dateALL <- unlist(dateLIST)
 notdup <- !duplicated(dateALL)
 nIterALL <- unlist(nIterLIST)

if (any(sapply(split( cbind.data.frame( nIterALL, dateALL), dateALL), function(x) 
               length(unique(x)))) != 1)
  stop("cannot merge - different specs for the same date")

 DATES <- sort(dateALL[notdup])

 nIter <- nIterALL[notdup][DATES]
 
 weights <- matrix( NA,  nrow = mem, ncol = length(nIterALL))

 j <- 0
 for (i in 1:length(z)) {
    l <- dim(z[[i]]$weights)[2]
    weights[,j+(1:l)] <- z[[i]]$weights[nam,]
    j <- j+l
 }

 weights <- weights[,notdup]

 dimnames(weights) <- list( nam, names(nIter))

# use member ordering as in x 
 nam <- dimnames(z[[1]]$weights)[[1]]
 weights <- weights[nam,DATES]

 structure(list( training = list(days = days, lag = lag, table = table), 
                 weights = weights, nIter = nIter, exchangeable = exch, 
                 sortedMemberNames = nam,
                 forecastHour = forecastHour,
                 initializationTime = initializationTime, notdup = notdup)) 
}

