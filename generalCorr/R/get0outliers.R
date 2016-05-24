#' Function to compute outliers and their count using Tukey method
#' using 1.5 times interquartile range (IQR) to define boundarirs.
#' 
#' 
#' @param x {vector of data.}
#' @param verbo {set to TRUE(default) assuming printed details are desired.}
#' @param mult {=1.5(default), the number of times IQR is used 
#'  in defining outlier boundaries.}
#' @return 
#' \item{below}{which items are lower than the lower limit} 
#' \item{above}{which items are larger than the upper limit} 
#' \item{low.lim}{the lower boundary for outlier detection} 
#' \item{up.lim}{the upper boundary for outlier detection} 
#' \item{nUP}{count of number of data points above upper boundary} 
#' \item{nLO}{count of number of data points below lower boundary} 
#' @note The function removes the missing data before checking for outliers.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @examples
#' 
#' set.seed(101);x=sample(1:100)[1:15];x[16]=150;x[17]=NA
#' get0outliers(x)#correctly identifies outlier=150
#' 
#' 
#' @export

get0outliers <-
function(x,verbo=TRUE,mult=1.5) {
    # function to compute the number of outliers automatically author
    # H. D. Vinod, Fordham university, New York, 24 March, 2006
    x2 = x[!is.na(x)]
    su = summary(x2)
    if (ncol(as.matrix(x)) > 1) {
        print("Error: input to get.outliers function has 2 or more columns")
        return(0)
    }
    iqr = su[5] - su[2]
    dn = su[2] - mult * iqr
    up = su[5] + mult * iqr
    LO = which(x2 < dn)  #vector of values below the lower limit
    nLO = length(LO)
    UP = which(x2 > up)
    nUP = length(UP)
if (verbo) print(c(" Q1-1.5*(inter quartile range)=", as.vector(dn), 
"number of outliers below it are=", 
        as.vector(nLO)), quote = F)
    
    if (nLO > 0) {
if(verbo) print(c("Actual values below the lower limit are:", x2[LO]), 
            quote = F)
    }
    
    # print(up, nUP)
if(verbo) print(c(" Q3+1.5*(inter quartile range)=", as.vector(up), 
" number of outliers above it are=", 
        as.vector(nUP)), quote = F)
    
    if (nUP > 0) {
if(verbo)   print(c("Actual values above the upper limit are:", x2[UP]), 
            quote = F)
    }
    
    list(below = LO, nLO = nLO, above = UP, nUP = nUP, low.lim = dn, 
        up.lim = up)
}
