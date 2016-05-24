#' scifix --- corrects missing e or E in 
#' scientific notation
#' INPUT = scalar either in or not in 
#' scientific notation
#' OUTPUT = corrected numeric value
#' M. Lau | July 2012
#' ------------------------------------

scifix <- function(x){
  x <- as.character(x)
                                        #e/E check
  e.check <- grepl('e',x)|grepl('E',x)
                                        #+/- check
  pm.check <- grepl('\\+',x)|grepl('\\-',x)
  if (pm.check&e.check==FALSE){
    if (grepl('\\+',x)){
      x <- sub('\\+','E+',x)
    }else{
      x <- sub('\\-','E-',x)
    }
  }else{}
  return(as.numeric(x))
}
