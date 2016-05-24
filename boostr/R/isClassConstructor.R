#' @title check if a function is a(n S3) class constructor
#' @description takes a function and returns a boolean indicating whether its
#' output gets assigned a class.
#' 
#' @export
#' 
#' @param func any function
#' 
#' @details The body of \code{func} is search for one of three idioms:
#' \enumerate{
#'  \item \code{UseMethod("className")}
#'  \item \code{class(output) <- classes}
#'  \item \code{attr(output, "class") <- classes}
#' }
#' If either is found, the assigned class (or classes) are returned as the 
#' \code{classes} attribute of the output. If none are found, a value of 
#' \code{FALSE} is returned (with no attributes). 
#' 
#' @return a boolean. If the return value is \code{TRUE} the boolean has
#' attribute \code{classes} which returns the (potential) classes for the
#' output of \code{func}
#' 
#' @examples
#' 
#' isClassConstructor(mean) # FALSE
#'  
#' # simple output
#' library(randomForest)
#' isClassConstructor(randomForest) # TRUE
#' 
#' # complicated output (multiple values in "classes")
#' isClassConstructor(glm) # TRUE 
#' isClassConstructor(lm) # TRUE
#' 

isClassConstructor <- function(func) {
# if I was really slick, I'd look at the last line (or return value)
# and then backwards track to try and infer class.
  
  # easy route: find UseMethod
  funcBody <- as.character(body(func))
  if (sum(grepl(pattern="UseMethod", x=funcBody)) != 0) {
    out <- TRUE
    # pattern is ...UseMethod..."class"... 
    # where ... are non-alpha-numeric characters
    className <- gsub(pattern='[^[:alnum:]]*(UseMethod)[^[:alnum:]]*([[:alnum:]]*)[^[:alnum:]]*',
                      replacement="\\2",
                      perl=TRUE,
                      x=paste(body(func), collapse="#"))
    attr(out, "classes") <- className
    return(out)
    
    # harder route: class(something) <- ...
  } else if (sum(str_count(funcBody, "class[(][[:alnum:]]*[)] <-")) != 0) {
    out <- TRUE
    # or class(something) <- "XXX"
    lineNum <- as.logical(str_count(funcBody, "class[(][[:alnum:]]*[)] <-"))
    relevantLine <- funcBody[lineNum]
    
    if (sum(lineNum) > 1) {
      warning('class() gets called multiple times inside constructor')
    }
    
    classes <- unlist(str_extract_all(relevantLine, '["][[:alnum:]]*["]'))
    classes <- sapply(str_split(classes, '\"'), function(x){x[2]})
    
    attr(out, "classes") <- classes
    return(out)
    # harest route: out(something ... "class") <- ...
  } else if (sum(str_count(funcBody, 'attr[(].*["]class["][)] <-') != 0)) {
    # need to look for out(something, "class") <- "XXX"
    out <- TRUE
    
    lineNum <- as.logical(str_count(funcBody, 'attr[(].*["]class["][)] <-'))
    relevantLine <- funcBody[lineNum]
    
    if (sum(lineNum) > 1) {
      warning('attr(..., "class") gets called multiple times inside constructor')
    }
    strings <- unlist(str_extract_all(relevantLine, '["][[:alnum:]]*["]'))
    classes <- sapply(str_split(strings, '\"'), function(x){x[2]})
    classes <- classes[sapply(classes, function(class) class != "class")]
    
    attr(out, "classes") <- classes
    return(out)
  } else {
    return(FALSE)
  }
}
