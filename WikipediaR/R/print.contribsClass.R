#' S3method print contribs
#' 
#' @description S3method to print contribsClass objects, results of the contribs function
#' 
#' @param x the contribsClass object, result of the contribs function
#' @param maxprint the maximal number of pages to print
#' @param ... others parameters
#' 
#' @method print contribsClass 
#' @seealso contribs
#'  
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @export
#' 
print.contribsClass <-function(x, maxprint = 10,...)
{
  if(class(x)!="contribsClass") stop("The argument Object must be a class contribsClass")
  
  cl <- x$call
  if (!is.null(cl)){
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  if(any(names(x)=="page"))
  {
    cat("Page title:", x$page[1],"\n")
    cat("Page id:", x$page[2],"\n")
    cat("Domain:", x$page[3],"\n")
    cat("\n")  
    
    if(any(names(x)=="contribs"))
    {
      cat("Total number of contributions:",nrow(x$contribs))
      cat("\n")   
      cat("List of contributions :","\n")
      nbprint <- min(maxprint,nrow(x$contribs))
      print(x$contribs[1:nbprint,])
      if(nbprint <nrow(x$contribs) ) {cat("... last",nrow(x$contribs) -maxprint,"pages not printed")}
    } 
    else{cat("0 contribution")}
    
    cat("\n")
    
  } # end if the page exists
  
  
}