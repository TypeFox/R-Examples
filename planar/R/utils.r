

test_complex <- function (x) any(Im(x) != 0)

is.metal <- function(x){
  if(is.character(x)) return(TRUE)
  Im(x) != 0
}

order_metal <- function(x){
  order(unlist(lapply(x, function(m) ifelse(is.character(m), m, mean(m)))))
}

##' epsilon_label
##'
##' characterise the layers of a structure with unique labels for metals and dielectrics
##' @title epsilon_label
##' @param epsilon list of real or complex values
##' @param names optional unique character names in order of appearance
##' @return factor
##' @export
##' @family user_level conversion utility
##' @author baptiste Auguie
epsilon_label <- function(epsilon = list(3.5, 1, 3, 1, "epsAu", 3, 3.5),
                          names=NULL){
  test_metal <- sapply(epsilon, is.metal)
  metals <- unique(epsilon[test_metal])
  dielectrics <- unique(epsilon[!test_metal])
  
  if(length(dielectrics)){
    sort_d <- order(unlist(lapply(dielectrics, mean)))
    dielectrics <- dielectrics[sort_d]
    
    dnames <- paste0("n", seq_along(dielectrics))[match(epsilon, 
                                                        dielectrics)]
  } else dnames <- rep(NA, length(epsilon))
  
  if(length(metals)){
    sort_m <- order_metal(metals)
    metals <- metals[sort_m]
    
    mnames <- paste0("metal", seq_along(metals))[match(epsilon, 
                                                       metals)]
  } else mnames <- c()
  
  f <- ifelse(is.na(dnames), mnames, dnames)
  uf <- sort(unique(f))
  if(!is.null(names) && length(names) == length(uf))
    return(factor(f, levels=uf, labels=names))
  
  factor(f, levels=uf)
}


##' epsilon_dispersion
##'
##' apply a function to a range of wavelength and return dielectric function
##' @title epsilon_dispersion
##' @param epsilon list of real or complex values
##' @param wavelength numeric vector
##' @param envir environment to look for functions
##' @return list
##' @export
##' @family utility
##' @author baptiste Auguie
epsilon_dispersion <- function(epsilon, wavelength=seq(400, 1000),
                               envir = parent.frame()){
  dispersive <- sapply(epsilon, is.character)
  if(!any(dispersive)) return(as.list(epsilon))
  replacement <- lapply(epsilon[dispersive], 
                        do.call, list(wavelength), envir = envir)
  epsilon[dispersive] <- lapply(replacement, "[[", "epsilon")
  epsilon
}

##' relabel factors
##'
##' Wide to long format data.frame with new factor variable(s) describing the original columns
##' @title classify
##' @param d data.frame
##' @param id column id
##' @param vars variables
##' @param ... passed on to melt
##' @return data.frame
##' @export
##' @family helping_functions
##' @author Baptiste Auguie
classify <- function(d, id=NULL, vars=NULL, ...){

  m <- melt(d, id.vars=id, ...)

  id.variables <- list()
  for (ii in seq_along(vars)){
    id.variables[[ii]] <- rep(vars[[ii]], each=nrow(d))
  }
  names(id.variables) <- names(vars)

  data.frame(m, id.variables)
}

##' relabel factors
##'
##' @title modify_levels
##' @param f factor
##' @param modify named list
##' @return factor
##' @export
##' @family helping_functions
##' @author Baptiste Auguie
modify_levels <- function(f, modify=list()){
  f <- factor(f)
  levs = levels(f)
  m = match(modify,levs)
  levs[m] = names(modify)
  factor(f,labels=levs)
}

Curry <- function (FUN, ...) 
{
    .orig = list(...)
    function(...) do.call(FUN, c(.orig, list(...)))
}


##' raman_shift
##'
##' converts Raman shift to wavelength
##' @title sort_factor
##' @param laser vector of laser wavelengths in nm
##' @param shift vector of Raman shifts in cm-1
##' @return matrix of shifted wavelengths (all combinations)
##' @export
##' @author Baptiste Auguie
raman_shift <- function(laser=c(514, 632.8), shift = c(520, 610)){
  
  res <- as.matrix(
    vapply(laser, function(lambda) 1 / (1 / lambda - shift *1e-7 ), shift))
  rownames(res) <- as.character(shift)
  colnames(res) <- as.character(laser)
  res
  
}