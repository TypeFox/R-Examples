# created on Saturday, July 11, 2009 at 12:57
# by MS
# now exist in R
#family <- function(object, ...)
#UseMethod("family")
#----------------------------------------------------------------------------------------
#  family and gamlss.family functions
#----------------------------------------------------------------------------------------
as.family <- function(object)
if(inherits(object, "family")) object else family(object)
#family.default<-function(object)
#{
#    cl <- data.class(object)[1]
#    return(switch(cl,
#        family = object,
#        "function" = family(object()),
#        character = family(get(object)),
#        name = family(eval(object)),
#        call = family(eval(object)),
#        NULL = gaussian(),
#        stop("The object argument is invalid")))
#}
#---------------------------------------------------------------------------------------
gamlss.family <-function (object, ...) 
UseMethod("gamlss.family")
#---------------------------------------------------------------------------------------
gamlss.family.default<-function(object,...)
{
    cl <- data.class(object)[1]
    return(switch(cl,
        gamlss.family = object,
        "function" = gamlss.family(object()),
        character = gamlss.family(get(object)),
        name = gamlss.family(eval(object)),
        call = gamlss.family(eval(object)),
        NULL = NO(),
        stop("The object argument is invalid")))
}

as.gamlss.family<-function(object)
if(inherits(object, "gamlss.family")) object else gamlss.family(object)
#----------------------------------------------------------------------------------
print.gamlss.family<-function (x, ...) 
{
    cat("\nGAMLSS Family:", x$family, "\n")
    if ("mu"%in%names(x$parameters))    
           cat("Link function for mu   :", x$mu.link, "\n")
    if ("sigma"%in%names(x$parameters)) 
           cat("Link function for sigma:", x$sigma.link, "\n")
    if ("nu"%in%names(x$parameters))    
           cat("Link function for nu   :", x$nu.link, "\n")
    if ("tau"%in%names(x$parameters))   
           cat("Link function for tau  :", x$tau.link, "\n")
}
#------------------------------------------------------------------------------------
#######################################################################################
#                 checklink
#######################################################################################
#---------------------------------------------------------------------------------------- 
#     Kalliope Akantziliotou
checklink <- function (which.link = NULL, 
                       which.dist = NULL, 
                             link = NULL, 
                        link.List = NULL
                        )  
{
    if (is.null(which.link)) 
       stop(paste("The parameter link name has not been defined."))
    if (is.null(which.dist)) 
       stop(paste("The distribution has not been defined."))
    if (is.null(link)) 
       stop(paste("The link has not been defined."))
    if (is.null(link.List)) 
       stop(paste("The list of links has not been defined."))
   #-----
   linktemp <- link
    if (!is.character(linktemp)) 
    {     # if not a character make it one
        linktemp <- deparse(linktemp)
    } 
   if (linktemp %in% link.List ) # if permisible link get it 
        stats <-  make.link.gamlss(linktemp)
    else if (is.character(link)) #
    {
        stats <- make.link.gamlss(link)
        linktemp <- link
    }
    else 
    {
        if (inherits(eval(link), "link-gamlss")) 
        {
            stats <- eval(link)
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else 
        {
            stop(gettextf("\"%s\ link \"%s\" not available, for \"%s\ family; available links are %s",  
            which.link, linktemp, which.dist, paste(sQuote(link.List), collapse = ", ")),
                domain = NA)
        }
    }
   link.result <- c(linktemp,stats)
   link.result
}
#----------------------------------------------------------------------------------------
#checklink <- function (which.link = NULL, 
#                       which.dist = NULL, 
#                             link = NULL, 
#                        link.List = NULL
#                        )  
#{
#    if (is.null(which.link)) 
#       stop(paste("The parameter link name has not been defined."))
#    if (is.null(which.dist)) 
#       stop(paste("The distribution has not been defined."))
#    if (is.null(link)) 
#       stop(paste("The link has not been defined."))
#    if (is.null(link.List)) 
#       stop(paste("The list of links has not been defined."))
#    linktemp <- link
#    if (!is.character(linktemp)) 
#    {
#        linktemp <- deparse(linktemp)
#        if (linktemp == which.link) 
#            linktemp <- eval(link)
#    }
#    if (any(linktemp == link.List)) 
#        stats <- make.link.gamlss(linktemp)#ms Sunday, February 20, 2005 
#    else 
#    {
#        availlinks <- ""
#        availlinksArr <- array(0:0, c(1,length(link.List)))
#        for  (i in 1:length(link.List))
#        {
#            if (length(link.List)==1 || i ==length(link.List))
#                availlinksArr[,i] <- paste("`",link.List[i], "'.", sep="")
#            else
#            {
#               comma.and <-ifelse((i==length(link.List)-1), " and", ",")
#               availlinksArr[,i] <- paste("`",link.List[i], "'", comma.and, sep="")
#            }
#            availlinks <-paste(availlinks, availlinksArr[,i])
#        }
#         if (length(link.List) != 1 )
#         {
#         stop(paste("`",linktemp, "' as ", which.link, " not available for ", which.dist, 
#           " family, available links are", availlinks, sep=""))  
#         }
#         else
#         {   
#         stop(paste("`",linktemp, "' as ", which.link, " not available for ", which.dist, 
#           " family, available links is only", availlinks, sep=""))
#         }
#    }
#    link.result <- c(linktemp,stats)
#    link.result
#}
#----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
# numerical derivatives
# taken from the Writing R Extensions p 48
# modified to have a varied delta 
numeric.deriv <- function(expr, theta, delta = NULL,  rho=sys.frame(sys.parent()))
{
  eps <- sqrt(.Machine$double.eps)  
  ans <- eval(substitute(expr), rho)
 grad <- matrix(,length(ans), length(theta), dimnames=list(NULL, theta))
for (i in seq(along=theta)) 
  {
  old <- get(theta[i], envir=rho)
delta <-  if (is.null(delta)) eps * min(1, abs(old)) else delta
 assign(theta[i], old+delta, envir=rho)
 ans1 <- eval(substitute(expr), rho)
 assign(theta[i], old, envir=rho)
 grad[,i] <- (ans1 - ans)/delta
  }
attr(ans, "gradient") <- grad
ans
}
#----------------------------------------------------------------------------------------
