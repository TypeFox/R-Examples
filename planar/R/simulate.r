##' simultate the far-field response of a multilayer stack
##'
##' wrapper around recursive_fresnelcpp for a stack structure
##' @title simulate_ff
##' @param fun function returning a stack
##' @param ... further arguments passed to fun
##' @param s stack (optional)
##' @param wavelength numeric vector
##' @param angle incident angle in radians
##' @param polarisation p or s
##' @return data.frame
##' @export
##' @family helping_functions user_level stack
##' @author Baptiste Auguie
simulate_ff <- function(..., s=NULL, fun = tamm_stack, 
                        wavelength = seq(400, 1000), 
                        angle=0, polarisation=c("p","s")){
  
  polarisation <- match.arg(polarisation)
  
  if(is.null(s)) s <- fun(...)
  stopifnot(check_stack(s))
  
  thickness <- s[["thickness"]]
  Nlay <- length(thickness)
  
  ## check that the stack is embedded with specific substrate and superstrate
  if(!(thickness[1] == 0L && thickness[Nlay] == 0L))
    s <- embed_stack(s)
  
  epsilon <- epsilon_dispersion(s[["epsilon"]], wavelength)
  
  results <- recursive_fresnelcpp(wavelength=wavelength,
                                  angle=angle,
                                  epsilon = epsilon, 
                                  thickness = s[["thickness"]], 
                                  polarisation = polarisation)
  
  d <- data.frame(results[c("wavelength", "angle", "R", "T", "A")])
  attr(d, "stack") <- s
  d
}

##' simultate the internal field of a multilayer stack
##'
##' wrapper around multilayer_field for a stack structure
##' @title simulate_nf
##' @param fun function returning a stack
##' @param ... further arguments passed to fun
##' @param s stack (optional)
##' @param wavelength numeric vector
##' @param angle incident angle in radians
##' @param polarisation p or s
##' @param res number of points
##' @param dmax maximum distance from stack boundary
##' @param field logical, return the real electric field
##' @return data.frame
##' @export
##' @family helping_functions user_level stack
##' @author Baptiste Auguie
simulate_nf <- function(..., s=NULL, fun = tamm_stack, 
                        wavelength = 630, 
                        angle=0, polarisation=c("p","s"),
                        dmax=0, res=1e4, field=FALSE){
  
  polarisation <- match.arg(polarisation)
  psi <- if(polarisation == "p") 0 else pi/2
  
  if(is.null(s)) s <- fun(...)
  stopifnot(check_stack(s))
  
  thickness <- s[["thickness"]]
  Nlay <- length(thickness)
  
  ## check that the stack is embedded with specific substrate and superstrate
  if(!(thickness[1] == 0L && thickness[Nlay] == 0L))
    s <- embed_stack(s)
  epsilon <- epsilon_dispersion(s[["epsilon"]], wavelength)
  
  d <- internal_field(wavelength=wavelength, angle=angle, psi=psi,
                      thickness = s[["thickness"]], 
                      dmax=dmax,  res=res, 
                      epsilon=unlist(epsilon), 
                      field = field)

  attr(d, "stack") <- s
  d
}

