#' Convert between numerical units
#'
#' This function converts between numerical units.
#' It works similar to the \code{iconv} function:
#' You provide vector \code{x} and a \code{from} and a \code{to} unit name
#' and the function converts.
#'
#' Additionally, you may provide a unitset name. Here, the
#' analogy to \code{iconv} ceases. Think of unitset as a
#' dimension of units, or a context for units. Predefined
#' unitsets are "Length", "Mass", "Energy", and "Temperature". 
#' It is recommened to provide the unitset name.
#' A list of available unitsets and the units defined by them
#' can be obtained with \code{uconvlist()}. 
#'
#' This is a proof of concept for the \code{datamart} function.
#'
#' @param x     numerical vector
#' @param from  character, unit to convert from. 
#' @param to    character, unit to convert to
#' @param uset  optional, character, unit set to use.
#'
#' @seealso \code{\link{datamart}}
#'
#' @examples
#' uconv(1, "horse length", "m", "Length")
#' uconv(1:10, "TWh", "PJ", "Energy")
#'
#' @export
uconv <- function(x, from, to, uset=NULL) {
      # Singleton pattern: global variable .UnitSetManager
      if(length(.UnitSetManager)==0) assign("x", unitsetmanager(), envir=.UnitSetManager)
      
      # no uset provided? Iterate
      if(is.null(uset)) {
        for (uset in queries(.UnitSetManager$x)) {
          # cat(uset, "\n")
          res <- try(uconv(x=x, from=from, to=to, uset=uset), silent=TRUE)
          if(!inherits(res, "try-error")) {
            warning("No unitset specified, using '", uset, "'", call.=FALSE)
            return(res)
          }
        }
        stop("uconv doesn't know how to convert from '", from, "' to '", to, "'.")
      }
      
      # get conversion resource
      obj <- try(query(.UnitSetManager$x, uset, verbose=FALSE), silent=TRUE)
      if(!is.numeric(obj)) stop("uconv doesn't know the unitset: '", uset, "'")
      
      # find the units
      cf <- obj[from] / obj[to]
      names(cf) <- NULL
      if(is.na(cf)) stop("Don't know how to convert from '", from, "' to '", to, "'.")
      
      # do conversion
      return(x*cf)
}

#' List unitsets and their units
#'
#' The function lists the currently available
#' unitsets and the units supported by them.
#'
#' @return named list, names=Unitsets, values=Units in these Unitsets
#' @export
uconvlist <- function() {
  if(length(.UnitSetManager)==0) assign("x", unitsetmanager(), envir=.UnitSetManager)
  usets <- queries(.UnitSetManager$x)
  res <- list()
  for (uset in usets) {
    obj <- query(.UnitSetManager$x, uset)
    if(is.numeric(obj))
        res[[uset]] <- names(obj)
    else if(is.function(obj)) {
        # from call to vector
        tmp <- as.list(formals(obj)$from)
        res[[uset]] <- unlist(tmp[2:length(tmp)])
    }
  }
  return(res)
}


#' UnitSetManager -- A class for unit conversion data
#'
#' This is an internal class. It administers the unitsets
#' used by the \code{uconv} method. One instance, usually
#' the only one, is created at startup.
#'
#' @examples
#' getSlots("UnitSetManager2")
#'
#' @name UnitSetManager-class
#' @rdname uconv
setClass(Class="UnitSetManager2", contains="Mashup")

uset.volume <- function(...) { 
    ret <- c(
    siunit("l", value=0.001, extended=TRUE, dimension=1), 
    siunit("m\u00b3", value=1, extended=TRUE, dimension=3),
    quart=1.1365225*0.001, 
    pint=0.56826125*0.001,
    "US gal"=3.785411784*0.001,
    "cb ft"=28.316846592*0.001,
    "cb in"=0.016387064*0.001,
    hl=100 # hektoliter
    )
    comment(ret) <- "m\u00b3"
    ret
}

uset.area <- function(...) { 
    ret <- c(
    siunit("m\u00b2", extended=TRUE, dimension=2), 
    ha=10^4, # Hektar
    a=100, #Ar
    b=10^(-28), # Barn
    sqin=(2.54/100)^2, # square inch
    "sq ft"= 0.3048^2, # square feet
    sqmi=1609.344^2, # square mile
    "mi\u00b2"=1609.344^2, # square mile
    ac=4046.8564224, # acre
    Rai=1600, # thailand
    Pyeong=400/121, # korea
    Tsubo=400/121 # japan
    )
    comment(ret) <- "m\u00b2"
    ret
}


uset.length <- function(...) { 
    ret <- c(
    siunit("m", extended=TRUE), 
    inch=2.54/100,
    "in"=2.54/100,
    bigpts=0.0003527778, # big points = 1/72 in
    thou=2.54/100000,
    ft=0.3048,
    yard=0.9144,
    mile=1609.344,
    fathom= 1.8288,
    "nautical mile"=1852,
    furlong = 201,
    "horse length"=2.4,
    ly=9.46e+15,
    AU=1.5e+11,
    RE=6370000,
    pc=30.8e+15,
    li=500,
    Ken=20/11, # japan
    "chi (PRC)"=1/3.0,
    "chi (Taiwan)"=10/33.,
    chek=0.371475,
    cun=33.3*10^(-3),
    tsun=37.148*10^(-3),
    sun=30.3*10^(-3)
    )
    comment(ret) <- "m"
    ret
}

uset.mass <- function(...) { 
    ret <- c(
    siunit("g", 10^(-3)),
    siunit("t", 10^3),
    siunit("eV", 1.783*10^(-36)), # 1 GeV/c2 = 1.783*10^(-27) kg
    u=1.66*10^(-27),
    sl=14.593903,
    slug=14.593903,
    lb=0.45359237,
    mP=2.1765113*10^(-8)
  )
  comment(ret) <- "kg"
  ret
}

uset.energy <- function(...) { 
    ret <- c(
    siunit("J", 10^(-6)),
    siunit("Wh", 3.6/10^(3)),
    siunit("toe", 41880),
    siunit("cal", 4.184*10^(-6)),
    siunit("tce", 29290),
    "BTU"=1055*10^(-6),
    Therm=105.5,
    
    # natural gas
    "cubic feet NG"=1.076,
    "m\uB3 NG"=38,
    
    # crude oil
    "bbl CL"=5713,
    "l CL"=35948498/10^6
    )
    comment(ret) <- "MJ"
    ret
}

# package-wide variable, lazy initialized, internal
.UnitSetManager <- new.env()

# The function is internal.
#
# param uname       the name of the unit to define, e.g. "g"
# param value       scaling factor, see details.
# param extended    logical (default=FALSE), add not so common prefixes?
# param dimension   allow squared and cubed unit for area and space (default=1)
siunit <- function(uname, value=1.0, extended=FALSE, dimension=1) {
    res1 <- value * 10^(seq(from=0, to=24, by=3))
    names(res1) <- paste(c('', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'), uname, sep="")
    
    res2 <- value *  10^(seq(from=-3, to=-24, by=-3))
    mu <- "\u03BC"
    Encoding(mu) <- "UTF-8"
    names(res2) <- paste(c('m', mu, 'n', 'p', 'f', 'a', 'z', 'y'), uname, sep="")
    
    res <- c(res1, res2)

    if(extended) {
      res3 <- value*c(100, 10, 0.1, 0.01)
      names(res3) <- paste(c("h", "da", "d", "c"), uname, sep="")
      res <- c(res, res3)
    }
    return(res^dimension)
}
  

#' Constructor for UnitSetManager objects
#'
#' Internal function to create an UnitSetManager object.
#'
#' @rdname uconv
unitsetmanager <- function() datamart(
    # Temperature=uset.temperature,
    resfunc(resource="Length", fun=uset.length),
    resfunc(resource="Energy", fun=uset.energy),
    resfunc(resource="Mass", fun=uset.mass),
    resfunc(resource="Area", fun=uset.area),
    resfunc(resource="Volume", fun=uset.volume), 
    clss="UnitSetManager2"
)
