#' S3method print userContribs
#' 
#' @description S3method to print userContribsClass objects, results of the userContribs function
#' 
#' @param x the userContribsClass object, result of the userContribs function
#' @param maxprint the maximal number of contributions to print
#' @param ... others parameters
#' 
#' @method print userContribsClass
#' @seealso userContribs
#' 
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou

#' 
#' @export
#' 
print.userContribsClass <-function(x, maxprint = 10,...)
{
  if(class(x)!="userContribsClass") stop("The argument Object must be a class userContribsClass")
  
  cl <- x$call
  if (!is.null(cl)){
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }  
  if(any(names(x)=="user"))
  {
    cat("User name:", x$user[1],"\n")
    cat("User id:", x$user[2],"\n")
    cat("Domain:", x$user[3],"\n")
    cat("\n")  
    
    if(any(names(x)=="contribs"))
    {
      cat("Total number of contributions:",nrow(x$contribs))
      cat("\n")   
      cat("List of contributions :","\n")
      nbprint <- min(maxprint,nrow(x$contribs))
      print(x$contribs[1:nbprint,])
      if(nbprint <nrow(x$contribs) ) {cat("... last",nrow(x$contribs) -maxprint,"contributions not printed")}
    }   
    else{cat("0 contribution")}
    cat("\n")  
  } # end user exists
}