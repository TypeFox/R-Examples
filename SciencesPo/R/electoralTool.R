#' @encoding UTF-8
#' @title Political Diversity Indices
#'
#' @description Analyzes political diversity in an electoral unity or across unities. It provides methods for estimating the effective number of parties and other fragmentation/concetration measures. The intuition of these coefficients is to counting parties while weighting them by their relative political--or electoral strength.
#'
#' @param x A data.frame, a matrix-like, or a vector containing values for the number of votes or seats each party received.
#' @param index The type of index desired, one of "laakso/taagepera",  "golosov", "herfindahl", "gini", "shannon", "simpson", "invsimpson".
#' @param margin The margin for which the index is computed.
#' @param base The logarithm base used in some indices, such as the "shannon" index.
#'
#' @details Very often, political analysts say things like \sQuote{two-party system} and \sQuote{multi-party system} to refer to a particular kind of political party system. However, these terms alone does not tell exactly how fragmented--or concentrated a party system actually is. For instance, after the 2010 general election, 22 parties obtained representation in the Lower Chamber in Brazil. Nonetheless, among these 22 parties, nine parties together returned only 28 MPs. Thus, an index to assess the weigh or the \bold{Effective Number of Parties} is important and helps to go beyond the simple count of parties in a legislative branch.
#'
#' A widely accepted algorithm was proposed by M. Laakso and R. Taagepera: \deqn{N = \frac{1}{\sum p_i^2}}{N = 1/ \sum p_i^2}, where \bold{N} denotes the effective number of parties and \bold{p_i} denotes the \eqn{it^h} party's fraction of the seats.
#'
#' In fact, this formula may be used to compute the vote share for each party. This formula is the reciprocal of a well-known concentration index (\bold{the Herfindahl-Hirschman index}) used in economics to study the degree to which ownership of firms in an industry is concentrated. Laakso and Taagepera correctly saw that the effective number of parties is simply an instance of the inverse measurement problem to that one. This index makes rough but fairly reliable international comparisons of party systems possible.
#' \bold{The Inverse Simpson index},
#' \deqn{ 1/ \lambda = {1 \over\sum_{i=1}^R p_i^2} = {}^2D}
#' Where \eqn{\lambda} equals the probability that two types taken at random from the dataset (with replacement) represent the same type. This simply equals true fragmentation of order 2, i.e. the effective number of parties that is obtained when the weighted arithmetic mean is used to quantify average proportional diversity of political parties in the election of interest.
#'
#' Another measure is the \bold{Least squares index (lsq)}, which measures the disproportionality produced by the election. Specifically, by the disparity between the distribution of votes and seats allocation.
#'
#' Recently, Grigorii Golosov proposed a new method for computing the effective number of parties  in which both larger and smaller parties are not attributed unrealistic scores as those resulted by using the Laakso/Taagepera index.I will call this as (\bold{Golosov}) and is given by the following formula: \deqn{N = \sum_{i=1}^{n}\frac{p_{i}}{p_{i}+p_{i}^{2}-p_{i}^{2}}}
#'

#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#'
#' @seealso \code{\link{cox.shugart}}, \code{\link{inv.cox.shugart}}, \code{\link{farina}}, \code{\link{grofman}}, \code{\link{gallagher}}, \code{\link{lijphart}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}
#'
#' @references Gallagher, Michael and Paul Mitchell (2005) \emph{The Politics of Electoral Systems.} Oxford University Press.
#'
#' Golosov, Grigorii (2010) The Effective Number of Parties: A New Approach, \emph{Party Politics,} \bold{16:} 171-192.
#'
#' Laakso, Markku and Rein Taagepera (1979) Effective Number of Parties: A Measure with Application to West Europe, \emph{Comparative Political Studies,} \bold{12:} 3-27.
#'
#' Nicolau, Jairo (2008) \emph{Sistemas Eleitorais.} Rio de Janeiro, FGV.
#'
#' Taagepera, Rein and Matthew S. Shugart (1989) \emph{Seats and Votes: The Effects and Determinants of Electoral Systems.} New Haven: Yale University Press.
#'
#' @keywords Diversity, Basics, Elections
#' @examples
#' # Here are some examples, help yourself:
#' # The wikipedia examples
#'
#' A <- c(.75,.25);
#' B <- c(.75,.10,rep(0.01,15))
#' C <- c(.55,.45);
#'
#' # The index by "laakso/taagepera" is the default
#' politicalDiversity(A)
#' politicalDiversity(B)
#'
#' # Using Grigorii Golosov method gives:
#' politicalDiversity(B, index="golosov")
#' politicalDiversity(C, index="golosov")
#'
#' # The 1980 presidential election in the US (vote share):
#' US1980 <- c("Democratic"=0.410, "Republican"=0.507,
#' "Independent"=0.066, "Libertarian"=0.011, "Citizens"=0.003,
#' "Others"=0.003)
#'
#' politicalDiversity(US1980)
#'
#' politicalDiversity(US1980, index= "herfindahl")
#'
#' politicalDiversity(US1980, index = "H") # will match Herfindahl
#'
#' # The 1999 Finland election:
#' votes_1999 <- c(612963, 600592, 563835,
#' 291675, 194846, 137330, 111835, 28084, 26440, 28549, 20442,
#' 10378, 10104, 5451, 5194, 4481, 3903, 3455, 21734)
#'
#' seats_1999 <- c(51, 48, 46, 20, 11, 11, 10, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#'
#' # 2010 Brazilian legislative election
#'
#' votes_2010 = c("PT"=13813587, "PMDB"=11692384, "PSDB"=9421347,
#' "DEM"=6932420, "PR"=7050274, "PP"=5987670, "PSB"=6553345,
#' "PDT"=4478736, "PTB"=3808646, "PSC"=2981714, "PV"=2886633,
#' "PC do B"=2545279, "PPS"=2376475, "PRB"=1659973, "PMN"=1026220,
#' "PT do B"=605768, "PSOL"=968475, "PHS"=719611, "PRTB"=283047,
#' "PRP"=232530, "PSL"=457490,"PTC"=563145)
#'
#' seats_2010 = c("PT"=88, "PMDB"=79, "PSDB"=53, "DEM"=43,
#' "PR"=41, "PP"=41, "PSB"=34, "PDT"=28, "PTB"=21, "PSC"=17,
#' "PV"=15, "PC do B"=15, "PPS"=12, "PRB"=8, "PMN"=4, "PT do B"=3,
#'  "PSOL"=3, "PHS"=2, "PRTB"=2, "PRP"=2, "PSL"=1,"PTC"=1)
#'
#' politicalDiversity(seats_2010)
#'
#' politicalDiversity(seats_2010, index= "golosov")
#'
#' @export politicalDiversity
#' @docType methods
#' @rdname politicalDiversity-methods
#' @aliases politicalDiversity,numeric,character,integer,numeric,ANY-method
`politicalDiversity`<- setClass("politicalDiversity", representation(x = "numeric",index="character", margin="integer", base="numeric"))
setGeneric("politicalDiversity", def=function(x, index = "laakso/taagepera", margin=1, base = exp(1)){
standardGeneric("politicalDiversity")})

#' @rdname politicalDiversity-methods
setMethod(f="politicalDiversity", definition=function(x, index = "laakso/taagepera", margin = 1, base = exp(1)){
      x <- drop(as.matrix(x))
      index <- .Match(arg = index, choices = c("laakso/taagepera", "golosov", "lsq", "enc",  "enp", "herfindahl", "gini", "simpson", "invsimpson","shannon") )
      if (length(dim(x)) > 1) {
        total <- apply(x, margin, sum)
        x <- sweep(x, margin, total, "/")
      }
      else {
        x <- x/sum(x)
      }
      if (index == "shannon")
        x <- -x * log(x, base)
      else if (index=="golosov")
        x <- sum((x)/((x)+((x[1])^2)-((x)^2)))
      else x <- x * x
      if (length(dim(x)) > 1)
        idx <- apply(x, margin, sum, na.rm = TRUE)
      else idx <- sum(x, na.rm = TRUE)
      if (index == "simpson"||index == "herfindahl")
        idx <- 1 - idx
      else if (index == "laakso/taagepera" || index == "invsimpson"||  index == "enc" || index == "enp")
        idx <- 1/idx
      return(round(idx, 3))
})### end -- politicalDiversity function
NULL






#' @title Gallagher Index
#'
#' @description Calculates the Gallagher index of LSq index.
#'
#' @param v A numeric vector of data values for votes each political party obtained.
#' @param s A numeric vector of data values for seats each political party obtained, the election outcome as seats.
#' @param \dots Additional arguements (currently ignored)
#'
#' @details The representation score is calculated as: sqrt(sum((Z-R)^2)/2).
#' @return A single score (The Gallagher's Representation Score.) given the votes each party received and seats obtained.
#'
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#'
#' @seealso \code{\link{cox.shugart}}, \code{\link{inv.cox.shugart}}, \code{\link{politicalDiversity}}, \code{\link{grofman}}, \code{\link{farina}},  \code{\link{lijphart}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}
#'
#' @references
#'  Gallagher, M. (1991) Proportionality, disproportionality and electoral systems. Electoral Studies 10(1):33-51.
#'
#'  @examples
#' # 2012 Queensland state elecion
#' pvotes= c(49.65, 26.66, 11.5, 7.53, 3.16, 1.47)
#' pseats = c(87.64, 7.87, 2.25, 0.00, 2.25, 0.00)
#'
#' gallagher(pvotes, pseats)

#' @export
#' @rdname gallagher
`gallagher`<- function(v, s, ...) UseMethod("gallagher")

#' @export
#' @rdname gallagher
`gallagher` <-function(v, s, ...){
  idx=sqrt(sum((v-s)^2)/2)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- gallagher function
NULL





#' @title Lijphart Index of Proportionality
#'
#' @description Calculates the Lijphart index of proportionality based on a vector of votes and a vector for the electoral outcome.
#'
#' @param v A numeric vector of data values for votes each political party obtained.
#' @param s A numeric vector of data values for seats each political party obtained, the election outcome as seats.
#' @param \dots Additional arguements (currently ignored)
#'
#' @return A single score given the votes each party received and seats obtained.
#'
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#'
#' @seealso \code{\link{cox.shugart}}, \code{\link{inv.cox.shugart}}, \code{\link{politicalDiversity}}, \code{\link{grofman}}, \code{\link{gallagher}},  \code{\link{farina}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}
#'
#' @examples

#' # 2012 Quebec provincial election:
#' pvotes = c(PQ=31.95, Lib=31.20, CAQ=27.05,QS=6.03,Option=1.89, Other=1.88)
#' pseats = c(PQ=54, Lib=50, CAQ=19, QS=2, Option=0, Other=0)
#'
#' lijphart(pvotes, pseats)
#'
#' @rdname lijphart
#' @export
`lijphart`<- function(v, s, ...) UseMethod("lijphart")


#' @export
#' @rdname lijphart
`lijphart` <-function(v, s, ...){
  idx=max(s-v)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- lijphart function
NULL




#' @title Grofman Index
#'
#' @description Calculates the Grofman index of proportionality based on a vector of votes and a vector for the electoral outcome.
#'
#' @param v A numeric vector of data values for votes each political party obtained.
#' @param s A numeric vector of data values for seats each political party obtained, the election outcome as seats.
#' @param \dots Additional arguements (currently ignored)
#'
#' @return A single score given the votes each party received and seats obtained.
#'
#'  @references
#' Taagepera, R., and B. Grofman. Mapping the indices of seats-votes disproportionality and inter-election volatility. Party Politics 9, no. 6 (2003): 659-77.
#'
#' @seealso \code{\link{cox.shugart}}, \code{\link{inv.cox.shugart}}, \code{\link{politicalDiversity}}, \code{\link{farina}}, \code{\link{gallagher}},  \code{\link{lijphart}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}
#'
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#' @examples
#'
#' # 2012 Quebec provincial election:
#' pvotes = c(PQ=31.95, Lib=31.20, CAQ=27.05,QS=6.03,Option=1.89, Other=1.88)
#' pseats = c(PQ=54, Lib=50, CAQ=19, QS=2, Option=0, Other=0)
#'
#' grofman(pvotes, pseats)
#'
#' @export
#' @rdname grofman
`grofman`<- function(v, s, ...) UseMethod("grofman")

#' @rdname grofman
#' @export
`grofman` <- function(v, s, ...){
  N <- politicalDiversity(s, index = "laakso/taagepera")
  idx=(1/N) * sum(abs(v-s))/2
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- grofman function
NULL



#' @title Farina Index
#'
#' @description Calculates the Farina index also referred to as the cosine proportionality score based on a vector of votes and a vector for the electoral outcome.
#' @param v A numeric vector of data values for votes each political party obtained.
#' @param s A numeric vector of data values for seats each political party obtained, the election outcome as seats.
#' @param \dots Additional arguements (currently ignored)
#'
#' @return A single score given the votes each party received and seats obtained.
#'
#' @seealso \code{\link{cox.shugart}}, \code{\link{inv.cox.shugart}}, \code{\link{politicalDiversity}}, \code{\link{grofman}}, \code{\link{gallagher}},  \code{\link{lijphart}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}.
#'
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#' @references
#' Koppel, M., and A. Diskin. (2009) Measuring disproportionality, volatility and malapportionment: axiomatization and solutions. Social Choice and Welfare 33, no. 2: 281-286.
#'
#' @examples
#'
#' # 2012 Queensland state elecion
#' pvotes= c(49.65, 26.66, 11.5, 7.53, 3.16, 1.47)
#' pseats = c(87.64, 7.87, 2.25, 0.00, 2.25, 0.00)
#'
#' farina(pvotes, pseats)
#'
#' @export farina
#' @rdname farina
`farina`<- function(v, s, ...) UseMethod("farina")

#' @export
#' @rdname farina
`farina` <- function(v, s, ...){
  idx= acos(sum(v*s)/(sum(v^2)*sum(s^2))^.5)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- farina function
NULL





#' @title Cox-Shugart Measure of Proportionality
#'
#' @description Calculate the Cox and Shugart measure of
#'  proportionalitybased on a vector of votes and a vector for
#'  the electoral outcome. This measure is also referred to as the regression index.
#'
#' @param v A numeric vector of data values for votes each political party obtained.
#' @param s A numeric vector of data values for seats each political party obtained, the election outcome as seats.
#' @param \dots Additional arguements (currently ignored)

#'
#' @return A single score given the votes each party received and seats obtained.
#'
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#'
#' @seealso \code{\link{inv.cox.shugart}}, \code{\link{farina}}, \code{\link{politicalDiversity}}, \code{\link{grofman}}, \code{\link{gallagher}},  \code{\link{lijphart}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}.
#'
#' @examples
#' if (interactive()) {
#' # 2012 Queensland state elecion:
#' pvotes= c(49.65, 26.66, 11.5, 7.53, 3.16, 1.47)
#' pseats = c(87.64, 7.87, 2.25, 0.00, 2.25, 0.00)
#'
#' cox.shugart(pvotes, pseats)
#'
#' # 2012 Quebec provincial election:
#' pvotes = c(PQ=31.95, Lib=31.20, CAQ=27.05, QS=6.03, Option=1.89, Other=1.88)
#' pseats = c(PQ=54, Lib=50, CAQ=19, QS=2, Option=0, Other=0)
#'
#' cox.shugart(pvotes, pseats)
#' }
#'
#' @export
#' @rdname cox.shugart
`cox.shugart` <- function(v, s, ...) UseMethod("cox.shugart")


#' @rdname cox.shugart
`cox.shugart` <- function(v, s, ...){
  S <- mean(s)
  V <- mean(v)
  idx <- sum((s-S) * (v-V))/sum((v-V)^2)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- cox.shugart function
NULL



#' @title Inverse Cox-Shugart Measure of Proportionality
#'
#' @description Calculate the inverse Cox and Shugart measure of
#'  proportionality based on votes and seats,
#'  the electoral outcome.
#'
#' @param v A numeric vector of data values for votes each political party obtained.
#' @param s A numeric vector of data values for seats each political party obtained, the election outcome as seats.
#' @param \dots Additional arguements (currently ignored)

#'
#' @return A single score given the votes each party received and seats obtained.
#'
#' @author Daniel Marcelino \email{dmarcelino@@live.com}
#'
#' @examples
#'
#' # 2012 Queensland state elecion:
#' pvotes= c(49.65, 26.66, 11.5, 7.53, 3.16, 1.47)
#' pseats = c(87.64, 7.87, 2.25, 0.00, 2.25, 0.00)
#'
#' inv.cox.shugart(pvotes, pseats)
#'
#' # 2012 Quebec provincial election:
#' pvotes = c(PQ=31.95, Lib=31.20, CAQ=27.05, QS=6.03, Option=1.89, Other=1.88)
#' pseats = c(PQ=54, Lib=50, CAQ=19, QS=2, Option=0, Other=0)
#'
#' inv.cox.shugart(pvotes, pseats)
#'
#' @seealso \code{\link{cox.shugart}}, \code{\link{farina}}, \code{\link{politicalDiversity}}, \code{\link{grofman}}, \code{\link{gallagher}},  \code{\link{lijphart}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}.
#'
#' @export
#' @rdname inv.cox.shugart
`inv.cox.shugart` <-function(v, s, ...) UseMethod("inv.cox.shugart")


#' @export
#' @rdname inv.cox.shugart
`inv.cox.shugart` <-function(v, s, ...){
  V <- mean(v)
  S <- mean(s)
  idx <- sum((v-V) * (s-S))/sum((s-S)^2)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- inv.cox.shugart function
NULL




#' @title Atkinson Index of Inequality
#'
#' @description Calculates the Atkinson Index. This inequality measure is espcially good at determining which end of the distribution is contributing most to the observed inequality.
#'
#' @param x A vector of data values of non-negative elements.
#' @param n A vector of frequencies of the same length as \code{x}.
#' @param parameter A parameter of the inequality measure (if set to \code{NULL} the default parameter of the respective measure is used).
#' @param na.rm A logical. Should missing values be removed? The Default is set to \code{na.rm=FALSE}.
#' @param \dots Additional arguements (currently ignored)
#'
#' @references
#' Cowell, F. A. (2000) Measurement of Inequality in Atkinson, A. B. / Bourguignon, F. (Eds): \emph{Handbook of Income Distribution}. Amsterdam.
#'
#' Cowell, F. A. (1995) \emph{Measuring Inequality} Harvester Wheatshef: Prentice Hall.
#'
#' @seealso \code{\link{herfindahl}}, \code{\link{rosenbluth}},  \code{\link{gini}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}.
#' @examples
#' if (interactive()) {
#' # generate a vector (of incomes)
#' x <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012)
#'
#' # compute Atkinson coefficient with parameter=0.5
#' atkinson(x, parameter=0.5)
#'}
#' @export
#' @rdname atkinson
`atkinson` <-function(x, n = rep(1, length(x)), parameter=0.5, na.rm=FALSE, ...) UseMethod("atkinson")
NULL

#' @export
#' @rdname atkinson
`atkinson` <- function(x, n = rep(1, length(x)), parameter = 0.5, na.rm = FALSE, ...){
  x <- rep(x, n)    # same handling as Lc and Gini
  if(na.rm) x <- na.omit(x)
  if (any(is.na(x)) || any(x < 0)) return(NA_real_)

  if(is.null(parameter)) parameter <- 0.5
  if(parameter==1)
    idx <- 1 - (exp(mean(log(x)))/mean(x))
  else
  {
    x <- (x/mean(x))^(1-parameter)
    idx <- 1 - mean(x)^(1/(1-parameter))
  }
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- atkinson function
NULL






#' @title Rosenbluth Index of Concentration
#'
#' @description Calculates the Rosenbluth Index of concentration, also known as Hall or Tiedemann Indices.
#'
#' @param x A vector of data values of non-negative elements.
#' @param n A vector of frequencies of the same length as \code{x}.
#' @param na.rm A logical. Should missing values be removed? The Default is set to \code{na.rm=FALSE}.
#' @param \dots Additional arguements (currently ignored)
#'
#' @references
#' Cowell, F. A. (2000) Measurement of Inequality in Atkinson, A. B. / Bourguignon, F. (Eds): \emph{Handbook of Income Distribution}. Amsterdam.
#'
#' Cowell, F. A. (1995) \emph{Measuring Inequality} Harvester Wheatshef: Prentice Hall.
#'
#' @seealso \code{\link{atkinson}}, \code{\link{herfindahl}},  \code{\link{gini}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}.
#'
#' @examples
#' # generate a vector (of incomes)
#' x <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012)
#'
#' # compute Rosenbluth coefficient
#' rosenbluth(x)
#'
#' @export
#' @rdname rosenbluth
`rosenbluth` <-function(x, n = rep(1, length(x)), na.rm=FALSE, ...)  UseMethod("rosenbluth")
NULL

#' @export
#' @rdname rosenbluth
`rosenbluth` <-function(x, n = rep(1, length(x)), na.rm = FALSE, ...){
  x <- rep(x, n)
  if(na.rm) x <- na.omit(x)
  if (any(is.na(x)) || any(x < 0)) return(NA_real_)
  n <- length(x)
  x <- sort(x)
  idx <- (n:1)*x
  idx <- 2*sum(idx/sum(x))
  idx <- 1/(idx-1)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- rosenbluth function
NULL









#' @encoding UTF-8
#' @title Herfindahl Index of Concentration
#'
#' @description Calculates the Herfindahl Index of concentration.
#'
#' @param x A vector of data values of non-negative elements.
#' @param n A vector of frequencies of the same length as \code{x}.
#' @param parameter A parameter of the concentration measure (if set to \code{NULL} the default parameter of the respective measure is used).
#' @param na.rm A logical. Should missing values be removed? The Default is set to \code{na.rm=FALSE}.
#' @param \dots Additional arguements (currently ignored)
#'
#' @details This index is also known as the \emph{Simpson Index} in ecology, the \emph{Herfindahl-Hirschman Index (HHI)} in economics, and as the \emph{Effective Number of Parties (ENP)} in political science.
#'
#' @references
#' Cowell, F. A. (2000) Measurement of Inequality in Atkinson, A. B. / Bourguignon, F. (Eds): \emph{Handbook of Income Distribution}. Amsterdam.
#'
#' Cowell, F. A. (1995) \emph{Measuring Inequality} Harvester Wheatshef: Prentice Hall.
#'
#' @seealso \code{\link{atkinson}}, \code{\link{rosenbluth}}, \code{\link{politicalDiversity}}, \code{\link{gini}}. For more details see the Indices vignette: \code{vignette("Indices", package = "SciencesPo")}.
#' @examples
#' # generate a vector (of incomes)
#' x <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012)
#'
#' # compute the Herfindahl coefficient with parameter=1
#' herfindahl(x, parameter=1)
#'
#'
#' @export
#' @rdname herfindahl
`herfindahl` <- function(x, n = rep(1, length(x)), parameter = 1, na.rm = FALSE, ...) UseMethod("herfindahl")


#' @export
#' @rdname herfindahl
`herfindahl` <- function(x, n = rep(1, length(x)), parameter = 1, na.rm = FALSE, ...){
  x <- rep(x, n)
  if(na.rm) x <- na.omit(x)
  if (any(is.na(x)) || any(x < 0)) return(NA_real_)

  if(is.null(parameter))
    m <- 1
  else
    m <- parameter
  idx <- x/sum(x)
  idx <- idx^(m+1)
  idx <- sum(idx)^(1/m)
  print(idx, digits = max(3, getOption("digits") - 3))
}### end -- herfindahl function
NULL




#' # hareQuota <- function(votes, seats, try.quota, droop.quota){}
#'
#'# electoralTool <- function(party,
#'#  district, seats, polar, time, blocks,
#'#   verbose=TRUE){
#'#   if(verbose)}

#'# votes86 = c(54109, 365376, 14325, 241402, 478405, 11788, 3055097, 2668055, 7973080, 50590, 13840, 1236493, 81424, 84156, 19520016, 41280, 6601, 28194, 27399, 35768, 1679, 27202, 50669, 387733, 166179, 2537691, 1722801 ,5844, 30907)


#'# 2010: ALL
#'# votes10=c(6932420,2580925, 22936, 1425, 4579006, 720154, 11867988, 1052218, 6862846, 2376475, 7097712, 1659973, 231282, 291341, 6582622, 2965055, 9500548, 177773, 464180, 969954, 54252, 14251798, 605768, 3912093, 565409, 155434, 2885915)

#'# only Elected:
#'# votes10=c(4550184,1823077,2596556, 96060, 8492579, 446162, 5197306, 1356469, 5889206, 871461, 47488, 50488, 4737508, 1859443, 6362262, 40093, 442756, 10472717, 178734, 2349527, 104015, 1005770)


#'#ro <- c(50740, 80954)

#'# results <- c(DEM=490205, PMDB=1151547, PRB=2449440, PSB=48274, PSTU=54403, PTC=173151)

#'# votes14=c(DEM=3868200, PCdoB=1799619, PCB=37253, PCO=8267, PDT=3141818, PEN=634682, PHS=903968, PMDB=10053108, PMN=432807, PP=6158835, PPL=103606, PPS=1875826, PR=5448721, PRB=4297373, PROS=1879940, PRP=655107, PRTB=430995, PSB=5574401, PSC=2420581, PSD=5637961, PSDB=9145950, PSDC=491280, PSL=772628, PSOL=1486393, PSTU=151353, PT=11803985, PTdoB=807509, PTB=3703639, PTC=312548, PTN=682854, PV=1808991, SD=2621639)
