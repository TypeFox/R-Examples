#' S3method print backLinks
#' 
#' @description S3method to print backLinksClass objects, results of the backLinks function
#' 
#' @param x the backLinksClass object, result of the backLinks function
#' @param maxprint the maximal number of pages to print
#' @param ... others parameters
#' 
#' @method print backLinksClass 
#' @seealso backLinks
#'  
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 

#' 
#' @export
#' 
print.backLinksClass <-function(x, maxprint = 10,...)
{
  if(class(x)!="backLinksClass") stop("The argument Object must be a class backLinksClass")
  
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
    
    if(any(names(x)=="backLinks"))
    {
      cat("Total number of back links:",nrow(x$backLinks))
      cat("\n")   
      cat("List of back links :","\n")
      nbprint <- min(maxprint,nrow(x$backLinks))
      print(x$backLinks[1:nbprint,])
      if(nbprint <nrow(x$backLinks) ) {cat("... last",nrow(x$backLinks) -maxprint,"pages not printed")}
    } 
    else{cat("0 back link")}
    cat("\n")  
  } # end if the page exists
  
  
}