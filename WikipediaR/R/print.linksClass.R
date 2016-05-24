#' S3method print links
#' 
#' @description S3method to print linksClass objects, results of the links function
#' 
#' @param x the linksClass object, result of the links function
#' @param maxprint the maximal number of pages to print
#' @param ... others parameters
#' 
#' @method print linksClass 
#' @seealso links
#'  
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou
#' 
#' @export
#' 
print.linksClass <-function(x, maxprint = 10,...)
{
  if(class(x)!="linksClass") stop("The argument Object must be a class linksClass")
  
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
    
    if(any(names(x)=="links"))
    {
      cat("Total number of links:",nrow(x$links))
      cat("\n")   
      cat("List of links :","\n")
      nbprint <- min(maxprint,nrow(x$links))
      print(x$links[1:nbprint,])
      if(nbprint <nrow(x$links) ) {cat("... last",nrow(x$links) -maxprint,"pages not printed \n")}
    } 
    else{cat("0 link \n")}
    
    if(any(names(x)=="extlinks"))
    {
      cat("\n")
      cat("Total number of external links:",length(x$extlinks))
      cat("\n")   
      cat("List of external links :","\n")
      nbprint <- min(maxprint,length(x$extlinks))
      print(x$extlinks[1:nbprint])
      if(nbprint <length(x$extlinks) ) {cat("... last",length(x$extlinks) -maxprint,"external links not printed")}
    } 
    else{cat("0 external link")}
    cat("\n")  
  } # end if the page exists
  
  
}