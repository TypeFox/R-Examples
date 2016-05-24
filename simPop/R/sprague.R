## create Sprague multiplier data





#' Sprague index (multipliers)
#' 
#' Using the Sprague multipliers, the age counts are estimated for each year
#' having 5-years interval data as input.
#' 
#' The input is population counts of age classes 0-4, 5-9, 10-14, ... , 77-74,
#' 75-79, 80+.
#' 
#' @name sprague
#' @param x numeric vector of age counts in five-year intervals
#' @return Population counts for age 0, 1, 2, 3, 4, ..., 78, 79, 80+.
#' @author Matthias Templ
#' @seealso \code{\link{whipple}}
#' @references G. Calot and J.-P. Sardon.  Methodology for the calculation of
#' Eurostat's demographic indicators.  Detailed report by the European
#' Demographic Observatory
#' 
#' \url{http://www.cedefop.europa.eu/EN/Files/Methodology_for_the_calculation_of_Eurostats_demographic_indicators.pdf}
#' @keywords arith
#' @export
#' @examples
#' 
#' ## example from the world bank
#' x <- data.frame(age=as.factor(c(
#'   "0-4",
#'   "5-9","10-14","15-19", "20-24",
#'   "25-29","30-34","35-39","40-44","45-49",
#'   "50-54","55-59","60-64","65-69","77-74","75-79","80+"
#'     )),
#'   pop=c(1971990, 2095820,2157190, 2094110,2116580,   2003840, 1785690,
#'         1502990, 1214170, 796934,  627551,  530305, 488014,
#'         364498, 259029,158047,  125941)
#' )
#' 
#' s  <- sprague(x[,2])
#' s
#'   
#' all.equal(sum(s), sum(x[,2]))
#' 
sprague <- function(x){
  if(length(x) != 17) stop("input writen for original sprague with age groups 0-4,
  5-9, 10-14, ... ,77-74,75-79,80+")
  breaks=c(seq(0,80,5),150)
#  if(tabular==FALSE){
#    x <- cut(x, breaks, right=FALSE)
#    Ns <- table(x)
#  } else {
    Ns <- x
#    if( length(Ns) != 18 ) {
#      warning("input table with five-yers age groups expected but not provided.")
#    }
#  }

  plus80 <- Ns[length(Ns)]

  multipliers <- data.frame(G1 = c(0.3616, 0.264, 0.184, 0.12, 0.0704,
            0.0336, 0.008, -0.008, -0.016, -0.0176,
            -0.0128,  -0.0016, 0.0064, 0.0064, 0.0016,
            rep(0,10)),
    G2 = c(-0.2768,  -0.096, 0.04, 0.136, 0.1968,
            0.2272, 0.232, 0.216, 0.184, 0.1408,
            0.0848, 0.0144, -0.0336, -0.0416, -0.024,
            -0.0144,  -0.008, 0, 0.008, 0.0144,
            0.0176, 0.016, 0.008, -0.008, -0.0336),
    G3 = c(0.1488, 0.04, -0.032, -0.072, -0.0848,
           -0.0752, -0.048, -0.008, 0.04, 0.0912,
            0.1504, 0.2224, 0.2544, 0.2224, 0.1504,
            0.0912, 0.04, -0.008, -0.048,  -0.0752,
            -0.0848, -0.072,  -0.032, 0.04, 0.1488),
    G4 = c(-0.0336, -0.008, 0.008,  0.016, 0.0176,
            0.0144,  0.008, 0, -0.008, -0.0144,
            -0.024,  -0.0416, -0.0336,  0.0144, 0.0848,
            0.1408,  0.184,  0.216, 0.232, 0.2272,
            0.1968, 0.136,  0.04, -0.096, -0.2768),
    G5 = c(0,  0, 0,  0, 0,
            0,  0, 0, 0, 0,
            0.0016,  0.0064, 0.0064, -0.0016, -0.0128,
            -0.0176, -0.016, -0.008, 0.008, 0.0336,
            0.0704, 0.12,  0.184, 0.264, 0.3616)
  )
  multipliers <- cbind(groups=rep(c("lowest","low","normal","high", "highest"),
                                  each=5), multipliers)

  # 27 years:
  ## correspond to n3 since it is
  ## in group 25,26,27,28,29
  infoGroup <- function(n, mult=multipliers, mybreaks=breaks, popN=Ns){
    ## from a five years group, which one of the five years:
    groups <- (n) %% 5
    ## extreme group or normal:
    if(n < 5){
      tab <- subset(mult, subset=groups=="lowest")
    } else if(n >= 5 & n < 10){
      tab <- subset(mult, subset=groups=="low")
    } else if(n >= 75 & n < 80){
      tab <- subset(mult, subset=groups=="high")
    } else if(n >= 79){
      tab <- subset(mult, subset=groups=="highest")
    } else{
      tab <- subset(mult, subset=groups=="normal")
    }
    ## my groups:
    ng <- cut(n, mybreaks, right=FALSE)
    mygroup <- which(levels(ng) %in% ng)
    ## cohort
    rowsm <- tab[groups+1,2:ncol(tab)]
    if(mygroup==1){
      s <- seq(mygroup,mygroup+4,1)
    } else if(mygroup==2){
      s <- seq(mygroup-1, mygroup+3,1)
    } else if(mygroup==17){
      s <- seq(mygroup-4,mygroup,1)
    } else if(mygroup==16){
      s <- seq(mygroup-3,mygroup+1,1)
    } else{
      s <- seq(mygroup-2,mygroup+2, 1)
    }


    cohort <- sum(rowsm * popN[s])
    return(cohort)
  }

   cohorts <- sapply(0:79, infoGroup)
  ## group 80+
   cohorts <- c(cohorts, plus80)
   names(cohorts) <- c(0:79,"80+")
   return(cohorts)
}
