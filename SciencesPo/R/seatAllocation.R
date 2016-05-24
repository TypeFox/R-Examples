#' @encoding UTF-8
#' @title The Hamilton Method of Allocating Seats Proportionally
#'
#' @description Computes the Alexander Hamilton's apportionment method (1792), also known as Hare-Niemeyer method or as Vinton's method. The Hamilton method is a largest-remainder method which uses the Hare Quota.
#'
#' @param parties A vector containig parties labels or candidates in the same order of \code{votes}.
#' @param votes A vector with the formal votes received by the parties/candidates.
#' @param seats An integer for the number of seats to be returned.
#' @param \dots Additional arguements (currently ignored)
#' @return A \code{data.frame} of length \code{parties} containing apportioned integers (seats) summing to \code{seats}.
#' @details The Hamilton/Vinton Method sets the divisor as the
#' proportion of the total population per house seat.
#' After each state's population is divided by the divisor,
#' the whole number of the quotient is kept and the fraction
#' dropped resulting in surplus house seats. Then, the first
#' surplus seat is assigned to the state with the largest
#' fraction after the original division. The next is assigned to
#' the state with the second-largest fraction and so on.
#' @references
#'  Lijphart, Arend (1994). \emph{Electoral Systems and Party Systems: A Study of Twenty-Seven Democracies, 1945-1990}. Oxford University Press.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#'
#' @seealso \code{\link{dHondt}}, \code{\link{highestAverages}}, \code{\link{largestRemainders}}, \code{\link{politicalDiversity}}.
#'
#' @importFrom utils head
#' @examples
#' votes <- sample(1:10000, 5)
#' parties <- sample(LETTERS, 5)
#' hamilton(parties, votes, seats = 4)
#'
#' @export
#' @rdname hamilton
`hamilton` <-function(parties=NULL, votes=NULL, seats=NULL,...) UseMethod("hamilton")


#' @export
#' @rdname hamilton
`hamilton` <-function(parties=NULL, votes=NULL, seats=NULL,...){
  # Modified :
  # v0.0 2011-10-25
  # v0.1 2012-07-10
  # v0.2 2016-01-05
  .temp <- data.frame(
    parties = parties,
    scores = votes / sum(votes) * seats,
    perc = round(votes / sum(votes),3));
  integer <- with(.temp, floor(scores));
  fraction <- with(.temp, scores - integer);
  remainder <- seats - sum(integer);
  .temp[,2] <- integer;
  extra <- utils::head(order(fraction, decreasing=TRUE), remainder);
  .temp$scores[extra] <- (.temp$scores[extra] + 1);
  if(sum(.temp$scores) != seats) stop("Allocation error.");
  names(.temp) <-c("Party", "Seats", "\u0025Seats");
  print(.temp, digits = max(3, getOption("digits") - 3))
}
NULL



#' @encoding UTF-8
#' @title The D'Hondt Method of Allocating Seats Proportionally
#'
#' @description The function calculate the seats allotment in legislative house, given the total number of seats and the votes for each party based on the Victor D'Hondt's method (1878), which is mathematically equivalent to the method proposed by Thomas Jefferson few years before (1792).
#'
#' @param parties A vector containig parties labels or candidates accordingly to the \code{votes} vector order.
#' @param votes A vector containing the total number of formal votes received by the parties/candidates.
#' @param seats An integer for the number of seats to be filled (the district magnitude).
#' @param \dots Additional arguements (currently ignored)
#'
#' @return A \code{data.frame} of length \code{parties} containing apportioned integers (seats) summing to \code{seats}.
#'
#' @keywords Electoral
#' @references
#'  Lijphart, Arend (1994). \emph{Electoral Systems and Party Systems: A Study of Twenty-Seven Democracies, 1945-1990}. Oxford University Press.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @seealso \code{\link{highestAverages}}, \code{\link{largestRemainders}},  \code{\link{hamilton}}, \code{\link{politicalDiversity}}.
#'
#' @note Adapted from Carlos Bellosta's replies in the R-list.
#'
#' @examples
#' # Example: 2014 Brazilian election for the lower house in
#' # the state of Ceara. Coalitions were leading by the
#' # following parties:
#'
#' results <- c(DEM=490205, PMDB=1151547, PRB=2449440,
#' PSB=48274, PSTU=54403, PTC=173151)
#'
#' dHondt(parties=names(results), votes=results, seats=19)
#'
#' # The next example is for the state legislative house of Ceara (2014):
#'
#' votes <- c(187906, 326841, 132531, 981096, 2043217,15061,103679,109830, 213988, 67145, 278267)
#'
#' parties <- c("PCdoB", "PDT","PEN", "PMDB", "PRB","PSB","PSC", "PSTU", "PTdoB", "PTC", "PTN")
#'
#' dHondt(parties, votes , seats=42)
#'
#' @importFrom utils head
#' @rdname dHondt
#' @export
`dHondt` <- function(parties=NULL, votes=NULL, seats=NULL, ...) UseMethod("dHondt")

#' @rdname dHondt
#' @export
`dHondt` <-function(parties=NULL, votes=NULL, seats=NULL, ...){
  # Modified :
  # v0.0 2011-10-25
  # v0.1 2012-07-10
  # v0.2 2016-01-05
  # creates a party score object
  .temp <- data.frame(
    parties = rep(parties, each = seats ),
    scores = as.vector(sapply( votes, function(x) x /
                                 1:seats ))
  );
  out <- with(.temp, (parties[order(-scores)][1:seats]))
  out <- freq(out, digits = 3);
  names(out) <-c("Party", "Seats", "\u0025Seats");
  # out <- out[ order(out[,2], decreasing = TRUE),]
  return(out)
}
NULL




#' @encoding latin1
#' @title Highest Averages Methods of Allocating Seats Proportionally
#'
#' @description Computes the highest averages method for a variety of formulas of allocating seats proportionally.
#' @param parties A character vector for parties labels or candidates in the same order as \code{votes}. If \code{NULL}, alphabet will be assigned.
#' @param votes A numeric vector for the number of formal votes received by each party or candidate.
#' @param seats The number of seats to be filled (scalar or vector).
#' @param method A character name for the method to be used. See details.
#' @param threshold A numeric value between (0~1). Default is set to 0.
#' @param \dots Additional arguements (currently ignored)
#'
#' @return A \code{data.frame} of length \code{parties} containing apportioned integers (seats) summing to \code{seats}.
#' @keywords Electoral
#'
#' @details The following methods are available:
#' \itemize{
#' \item {"dh"}{d'Hondt method}
#' \item {"sl"}{Sainte-Lague method}
#' \item {"msl"}{Modified Sainte-Lague method}
#' \item {"danish"}{Danish modified Sainte-Lague method}
#' \item {"hsl"}{Hungarian modified Sainte-Lague method}
#' \item {"imperiali"}{The Italian Imperiali (not to be confused with the Imperiali quota which is a Largest remainder method)}
#' \item {"hh"}{Huntington-Hill method}
#' \item {"wb"}{Webster's method}
#' \item {"jef"}{Jefferson's method}
#' \item {"ad"}{Adams's method}
#' \item {"hb"}{Hagenbach-Bischoff method}
#' }
#'
#' @references
#' Gallagher, Michael (1992). "Comparing Proportional Representation
#' Electoral Systems: Quotas, Thresholds, Paradoxes and Majorities".
#' \emph{British Journal of Political Science}, 22, 4, 469-496.
#'
#'  Lijphart, Arend (1994). \emph{Electoral Systems and Party Systems: A Study of Twenty-Seven Democracies, 1945-1990}. Oxford University Press.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @seealso \code{\link{largestRemainders}}, \code{\link{dHondt}}, \code{\link{hamilton}}, \code{\link{politicalDiversity}}. For more details see the \emph{Indices} vignette: \code{vignette('Indices', package = 'SciencesPo')}.
#'
#' @examples
#' # Results for the state legislative house of Ceara (2014):
#' votes <- c(187906, 326841, 132531, 981096, 2043217, 15061, 103679,109830, 213988, 67145, 278267)
#'
#' parties <- c("PCdoB", "PDT", "PEN", "PMDB", "PRB", "PSB", "PSC", "PSTU", "PTdoB", "PTC", "PTN")
#'
#' highestAverages(parties, votes, seats = 42, method = "dh")
#'
#' # Let's create a data.frame with typical election results
#' # with the following parties and votes to return 10 seats:
#'
#' my_election <- data.frame(
#' party=c("Yellow", "White", "Red", "Green", "Blue", "Pink"),
#' votes=c(47000, 16000,	15900,	12000,	6000,	3100))
#'
#' highestAverages(my_election$party,
#' my_election$votes,
#' seats = 10,
#' method="dh")
#'
#' # How this compares to the Sainte-Lague Method
#'
#'(dat= highestAverages(my_election$party,
#' my_election$votes,
#' seats = 10,
#' method="sl"))
#'
#' # Plot it
#' bar.plot(data=dat, "Party", "Seats") +
#' theme_fte()
#'
#' @rdname highestAverages
#' @export
`highestAverages` <- function(parties=NULL, votes=NULL, seats=NULL, method=c("dh", "sl", "msl", "danish", "hsl", "hh", "imperiali", "wb", "jef", "ad", "hb"), threshold=0, ...) UseMethod("highestAverages")



#' @export
#' @rdname highestAverages
`highestAverages.default` <- function(parties=NULL, votes=NULL, seats=NULL, method=c("dh", "sl", "msl", "danish", "hsl", "hh", "imperiali", "wb", "jef", "ad", "hb"), threshold=0, ...){
  # Modified :
  # v0.0 2013-11-21
  # v0.1 2014-10-02
  # v0.2 2016-01-13
  # local vars for using later
  .ratio <- votes/sum(votes)
  .votes <- ifelse(.ratio < threshold, 0, votes)

  # To deal with  NULL party labels
  if (is.null(parties)){
    parties <- replicate(length(votes),
                         paste(sample(LETTERS, 3,
                                      replace=TRUE), collapse=""))
  }

  # Define Quotient
  switch(method,
         dh = { #d'Hondt
           divisor.vec <- seq(from = 1, by = 1, length.out = seats)
           method.name <- c("d'Hondt")
         },
         sl = { #Sainte-Lague
           divisor.vec <- seq(from = 1, by = 2, length.out = seats)
           method.name <- c("Sainte-Lagu\u00EB")
         },
         msl = { #Modified Sainte-Lague
           divisor.vec <- c(1.4, seq(from = 3, by = 2, length.out = seats-1))
           method.name <- c("Modified Sainte-Lagu\u00EB")
         },
         danish = { #Danish
           divisor.vec <- c(1, seq(from = 4, by = 3, length.out = seats-1))
           method.name <- c("Danish Sainte-Lagu\u00EB")
         },
         hsl = { #Hungarian
           divisor.vec <- c(1.5, seq(from = 3, by = 2, length.out = seats-1))
           method.name <- c("Hungarian Sainte-Lagu\u00EB")
         },
         imperiali = { #Imperiali
           divisor.vec <- c(1, seq(from = 1.5, by = .5, length.out = seats-1))
           method.name <- c("Imperiali")
         },
         hh = { #Huntington-Hill Equal Proportions Method
           divisor.vec0 <- seq(from = 1, by = 1, length.out = seats)
           divisor.vec <- sqrt(divisor.vec0 * (divisor.vec0 - 1))
           method.name <- c("Hungtinton-Hill")
         },
         wb = { #Webster Major Fractions Method
           divisor.vec0 <- seq(from = 1, by = 2, length.out = seats)
           divisor.vec <- (divisor.vec0+(divisor.vec0 - 1))/2
           method.name <- c("Webster")
         },
         jef = { #Jefferson Greatest Divisors or Hagenbach-Bischoff Method
           divisor.vec <- seq(from = 1, by = 1, length.out = seats)
           method.name <- c("Jefferson")
         },
         ad = { #Adam's Method Smallest Devisors
           divisor.vec <- c(0, seq(from = 1, by = 1, length.out = seats-1))
           method.name <- c("Adam's Method")
         },
         hb = { #Hagenbach-Bischoff Method
           divisor.vec <- seq(from = 1, by = 1, length.out = seats)
           method.name <- c("Hagenbach-Bischoff")
         }
  )

  # ratio = as.vector(sapply(votes, function(x) x /
  # sum(votes)))
  .temp <- data.frame(
    parties = rep(parties, each = seats ),
    scores = as.vector(sapply(.votes, function(x) x /
                                divisor.vec ))
  );

  out <- with(.temp, (parties[order(-scores)][1:seats]))

  out <- freq(out, digits = 3);
  names(out) <-c("Party", "Seats", "\u0025Seats");
  # Political diversity indices
  ENP.votes <- 1/sum(.ratio^2)
  ENP.seats <- 1/sum((out$Seats/sum(out$Seats))^2)
  LSq.index <- sqrt(0.5*sum((((votes/sum(votes))*100) - ((out$Seats/sum(out$Seats))*100))^2))

  cat("Method:", method.name, "\n")
  shorten(round(divisor.vec, 2), 4)
  cat(paste("ENP:",round(ENP.votes,2),"(After):",round(ENP.seats,2)),"\n")
  cat(paste("Gallagher Index: ", round(LSq.index, 2)), "\n \n")
  return(out)
}
NULL






#' @encoding latin1
#' @title Largest Remainders Methods of Allocating Seats Proportionally
#'
#' @description Computes the largest remainders method for a variety of formulas of allocating seats proportionally.
#' @param parties A character vector for parties labels or candidates in the order as \code{votes}. If \code{NULL}, a random combination of letters will be assigned.
#' @param votes A numeric vector for the number of formal votes received by each party or candidate.
#' @param seats The number of seats to be filled (scalar or vector).
#' @param method A character name for the method to be used. See details.
#' @param threshold A numeric value between (0~1). Default is set to 0.
#' @param \dots Additional arguements (currently ignored)
#'
#' @return A \code{data.frame} of length \code{parties} containing apportioned integers (seats) summing to \code{seats}.
#' @keywords Electoral
#'
#' @details The following methods are available:
#' \itemize{
#' \item {"dh"}{d'Hondt method}
#' \item {"sl"}{Sainte-Lague method}
#' }
#'
#' @references
#' Gallagher, Michael (1992). "Comparing Proportional Representation
#' Electoral Systems: Quotas, Thresholds, Paradoxes and Majorities".
#' \emph{British Journal of Political Science}, 22, 4, 469-496.
#'
#'  Lijphart, Arend (1994). \emph{Electoral Systems and Party Systems: A Study of Twenty-Seven Democracies, 1945-1990}. Oxford University Press.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @seealso  \code{\link{highestAverages}}, \code{\link{dHondt}}, \code{\link{hamilton}}, \code{\link{politicalDiversity}}. For more details see the \emph{Indices} vignette: \code{vignette('Indices', package = 'SciencesPo')}.
#'
#' @examples
#' # Let's create a data.frame with typical election results
#' # with the following parties and votes to return 10 seats:
#'
#' my_election <- data.frame(
#' party=c("Yellow", "White", "Red", "Green", "Blue", "Pink"),
#' votes=c(47000, 16000,	15900,	12000,	6000,	3100))
#'
#' largestRemainders(my_election$party,
#' my_election$votes, seats = 10,  method="droop")
#'
#' @rdname largestRemainders
#' @export
`largestRemainders` <- function(parties=NULL, votes=NULL, seats=NULL, method=c("dh", "sl", "msl", "danish", "hsl", "hh", "imperiali", "wb", "jef", "ad", "hb"), threshold=0, ...) UseMethod("largestRemainders")



#' @export
#' @rdname largestRemainders
`largestRemainders.default` <- function(parties=NULL, votes=NULL, seats=NULL, method=c("dh", "sl", "msl", "danish", "hsl", "hh", "imperiali", "wb", "jef", "ad", "hb"), threshold=0, ...){
  # Modified :
  # v0.0 2013-11-21
  # v0.1 2014-10-02
  # v0.2 2016-01-13
  # local vars for using later
  .ratio <- votes/sum(votes)
  .votes <- ifelse(.ratio < threshold, 0, votes)

  # To deal with  NULL party labels
  if (is.null(parties)){
    parties <- replicate(length(votes),
                         paste(sample(LETTERS, 3,
                                      replace=TRUE), collapse=""))
  }

  # Define Quotient

}

