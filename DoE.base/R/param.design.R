param.design <- function(inner, outer, direction="long", responses=NULL, ...){
    if (!"design" %in% class(inner))
       stop ("inner must be a design")
    if (!direction %in% c("long","wide")) stop("direction must be one of long or wide")
    if (!is.null(responses)) {if (!is.character(responses)) 
       stop("if given, responses must be a character vector of long format response names")
       responses <- make.names(responses, unique=TRUE)
       }
    di <- design.info(inner)
    if (!di$randomize) warning("inner array should be randomized")
    if (di$replications>1 | di$repeat.only) 
       stop("inner / outer arrays cannot have replications")
    if (!"design" %in% class(outer))
        if (is.data.frame(outer) | is.matrix(outer) | is.list(outer) | is.array(outer))
             stop("outer must be a vector or a data frame of class design")
    aus <- cross.design(inner, outer, randomize=FALSE)
    ## modify design.info attribute
    design.info(aus)$type <- "param"
    ctypes <- design.info(aus)$cross.types
    if (all(substr(ctypes,1,4)=="FrF2") | (substr(ctypes[1],1,4)=="FrF2" & substr(ctypes[2],1,6)=="vector" ))
        design.info(aus)$type <- "FrF2.param"
    if (!is.null(design.info(aus)$aliased)) names(design.info(aus)$aliased) <- c("inner","outer")
    design.info(aus)$inner <- factor.names(inner)
    design.info(aus)$outer <- factor.names(outer)
    if (!is.null(responses)) response.names(aus) <- responses 
       else response.names(aus) <- NULL
    if (nrow(outer)>8)
       warning("Are you sure ?\nA Taguchi inner/outer array with more than 8 runs in the outer array is very unusual")
    if (direction=="wide") aus <- paramtowide(aus)
    aus
}