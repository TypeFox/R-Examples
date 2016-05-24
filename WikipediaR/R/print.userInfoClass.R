#' S3method print userInfo
#' 
#' @description S3method to print userInfoClass objects, results of the userInfo function
#' 
#' @param x the userInfoClass object, result of the userInfo function
#' @param ... others parameters
#' 
#' @method print userInfoClass
#' @seealso userInfo
#' 
#' @author Avner Bar-Hen, Louise Baschet, Francois-Xavier Jollois, Jeremie Riou

#' 
#' @export
#' 
print.userInfoClass <-function(x,...)
{
  if(class(x)!="userInfoClass") stop("The argument Object must be a class userInfoClass")
  
  cl <- x$call
  if (!is.null(cl)){
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }  
  if(any(names(x)=="info"))
  {
    cat("User name:", as.vector(x$info$name),"\n")
    cat("User id:", as.vector(x$info$userid),"\n")
    cat("Domain:", as.vector(x$info$domain),"\n")
    cat("\n")  
    if(length(x$info) > 3 ) {cat("others informations: \n") ; print(x$info[4:length(x$info)])}
  }
  if(any(names(x)=="rights"))
  {
    cat("\n")   
    cat("Total number of rights:",length(x$rights))
    cat("\n")   
    cat("List of rights:","\n")
    print(x$rights)
  }   
  
  if(any(names(x)=="groups"))
  {
    cat("\n")   
    cat("Total number of groups:",length(x$groups))
    cat("\n")   
    cat("List of groups:","\n")
    print(x$groups)
  }  
  if(any(names(x)=="implicitgroups"))
  {
    cat("\n")   
    cat("Total number of implicit groups:",length(x$implicitgroups))
    cat("\n")   
    cat("List of implicit groups:","\n")
    print(x$implicitgroups)
  }    
  cat("\n")  
}