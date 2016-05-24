#' @title Expected frequencies using margins 
#' @export expected_margin_cfa
#' @description Calculates the expected frequencies of counts based on the margins of the k-dimensional contingency table. 
#' @details only main effects are considered.
#' 
#' @param Pfreq Object of class "Pfreq" (see. function \code{\link{dat2fre}}).
#' @param blank Either (1) character vector defining the pattern (with spaces between variable categories), which will be ignored for calculation of expected frequencies; or (2) a numeric vector defining the position(s) of the pattern in object of class \code{"Pfreq"}, which will be ignored for calculation of expected frequencies. At default (\code{blank=NULL}) all possible pattern, as listed in object of class \code{"Pfreq"}, are included for calculation of expected frequencies.
#' 
#' @return An vector object giving the expected counts.
#' @references No references in the moment
#' @examples #######################################
#' # expected counts for LienertLSD data example.
#' data(LienertLSD) # load example data
#' expected_margin_cfa(Pfreq = LienertLSD) # calculation of expected counts (only main effects).
#' ####################################### 
expected_margin_cfa <- function(Pfreq, blank=NULL){
  # Pfreq object of class "Pfreq"
  # blank either (1) character vector defining the pattern  
  if(class(blank)!="NULL"){
    if(class(blank)!="character"){
      Pfreq$Freq[blank] <- NA  
    }
    if(class(blank)=="character"){
      posp <- base::apply(Pfreq[,1:(ncol(Pfreq)-1)],1,paste, collapse=" ") 
      Freq_ind <- sapply(blank, function(x){which(x==posp)})
      Pfreq$Freq[Freq_ind] <- NA
    }  
  }
  
  # helper functions
  ntable <- function(Pfreq){
    ndim <- base::apply(Pfreq[,1:(ncol(Pfreq)-1)],2,function(x){length(unique(x))})
    ndimnames <- lapply(as.list(Pfreq[,1:(ncol(Pfreq)-1)]),FUN = function(x){as.character(unique(x))})
    # names(ndimnames) <- names(ndim)
    res <- array(data = Pfreq$Freq ,dim = rev(ndim), dimnames = rev(ndimnames)) #hier war noch ein fehler -es mus rev() verwendet werden
    #      res <- Pfreq$Freq
    #      dim(res)  <- rev(ndim) 
    #res <- table(fre2dat(Pfreq,fact = T))
    res
  }
  margin.table.na <- function (x, margin = NULL,na.rm=T) {
    if (!is.array(x)) 
      stop("'x' is not an array")
    if (length(margin)) {
      z <- apply(x, margin, sum, na.rm=na.rm)
      dim(z) <- dim(x)[margin]
      dimnames(z) <- dimnames(x)[margin]
    }
    else return(sum(x,na.rm=na.rm))
    class(z) <- oldClass(x)
    z
  }  
  #### start function
  
  res_tab <- ntable(Pfreq)
  res_marg <- lapply((1:(ncol(Pfreq)-1)),FUN = function(x){margin.table.na(res_tab,x)})
  res_pm <- as.matrix(expand.grid(... = res_marg))
  result <- apply(res_pm,1,prod) / (margin.table.na(res_tab))^(length(dim(res_tab))-1)### das sind die expected
  # result <- rev(result)
  return(result)
}
