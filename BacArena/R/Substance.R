# Substance inherits from Grid and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

#' Structure of the S4 class "Substance"
#' 
#' Structure of the S4 class \code{Substance} representing substances in the environment which can be produced or consumed.
#' @import ReacTran deSolve
#' @export Substance
#' @exportClass Substance
#' @rdname Substance
#'
#' @slot smax A number representing the start concentration of the substance for each grid cell in the environment. 
#' @slot diffmat A sparse matrix containing all concentrations of the substance in the environment.
#' @slot name A character vector representing the name of the substance.
#' @slot difunc A character vector ("pde","cpp" or "r") describing the function for diffusion.
#' @slot difspeed A number indicating the diffusion speed (given by cm^2/s).
#' @slot diffgeometry Diffusion coefficient defined on all grid cells (initially set by constructor).
#' @slot pde R-function that computes the values of the derivatives in the diffusion system
#' @slot boundS A number defining the attached amount of substance at the boundary (Warning: boundary-function must be set in pde!)
setClass("Substance",
         representation(
           smax = "numeric",
           diffmat = "Matrix",
           name = "character",
           difunc = "character",
           difspeed = "numeric",
           diffgeometry = "list",
           pde = "character",
           boundS = "numeric"
         ),
         prototype(
           difunc = "pde",
           pde="Diff2d",
           boundS = 0
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{Substance}
#' 
#' The constructor to get a new object of class \code{Substance}
#' @export
#' @name Substance-constructor
#' 
#' @param smax A number representing the start concentration of the substance for each grid cell in the environment. 
#' @param difspeed A number indicating the diffusion speed (given by cm^2/s).
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param gridgeometry A list containing grid geometry parameter 
#' @param ... Arguments of \code{\link{Substance-class}}
#' @return Object of class \code{Substance}
Substance <- function(n, m, smax, gridgeometry, difspeed=1, ...){
  diffmat <- Matrix::Matrix(smax, nrow=n, ncol=m, sparse=TRUE)
  
  Dgrid <- ReacTran::setup.prop.2D(value = difspeed, grid = gridgeometry$grid2D)
  diffgeometry <- list(Dgrid=Dgrid)
  new("Substance", smax=smax, diffmat=diffmat, difspeed=difspeed, diffgeometry=diffgeometry, ...)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("smax", function(object){standardGeneric("smax")})
setMethod("smax", "Substance", function(object){return(object@smax)})
setGeneric("diffmat", function(object){standardGeneric("diffmat")})
setMethod("diffmat", "Substance", function(object){return(object@diffmat)})
setGeneric("name", function(object){standardGeneric("name")})
setMethod("name", "Substance", function(object){return(object@name)})
setGeneric("difunc", function(object){standardGeneric("difunc")})
setMethod("difunc", "Substance", function(object){return(object@difunc)})
setGeneric("difspeed", function(object){standardGeneric("difspeed")})
setMethod("difspeed", "Substance", function(object){return(object@difspeed)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for naive diffusion (neighbourhood) of the Substance matrix
#'
#' @description The generic function \code{diffuseR} implements the diffusion in the Moore neighbourhood in \code{R}.
#' @export
#' @rdname diffuseR
#'
#' @param object An object of class Substance.
#' @details The diffusion is implemented by iterating through each cell in the grid and taking the cell with the lowest concentration in the Moore neighbourhood to update the concentration of both by their mean.
#' @seealso \code{\link{Substance-class}} and \code{\link{diffusePDE}}
#' @examples
#' arena <- Arena(n=100, m=100, stir=FALSE, Lx=0.025, Ly=0.025)
#' sub <- Substance(n=20,m=20,smax=40,name='test',difunc='r', 
#'                  gridgeometry=arena@gridgeometry) #initialize test substance
#' diffuseR(sub)
setGeneric("diffuseR", function(object){standardGeneric("diffuseR")})
#' @export
#' @rdname diffuseR
setMethod("diffuseR", "Substance", function(object){
  smat <- as(object@diffmat, "matrix")
  smatn <- matrix(NA, nrow=dim(smat)[1]+2, ncol=dim(smat)[2]+2) #define environment with boundary conditions
  smatn[2:(dim(smat)[1]+1), 2:(dim(smat)[2]+1)] <- smat #put the values into the environment
  i <- sample(1:dim(smat)[1], dim(smat)[1])
  j <- sample(1:dim(smat)[2], dim(smat)[2])
  for(ic in seq_along(i)){
    for(jc in seq_along(j)){
      neighbours <- c(smatn[ic,jc], 
                      smatn[ic+1,jc], 
                      smatn[ic+2,jc], 
                      smatn[ic+2,jc+1],
                      smatn[ic+2,jc+2], 
                      smatn[ic+1,jc+2],
                      smatn[ic,jc+2],
                      smatn[ic,jc+1])
      minc <- min(neighbours, na.rm = T)
      if(smat[ic,jc] > minc){
        nmin <- which(neighbours == minc)
        if(length(nmin)!=1){
          nmin <- sample(nmin, 1)
        }
        switch(nmin,
              {mn <- mean(c(smat[ic,jc], smat[ic-1,jc-1])); smat[ic,jc] <- mn; smat[ic-1,jc-1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic,jc-1])); smat[ic,jc] <- mn; smat[ic,jc-1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic+1,jc-1])); smat[ic,jc] <- mn; smat[ic+1,jc-1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic+1,jc])); smat[ic,jc] <- mn; smat[ic+1,jc] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic+1,jc+1])); smat[ic,jc] <- mn; smat[ic+1,jc+1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic,jc+1])); smat[ic,jc] <- mn; smat[ic,jc+1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic-1,jc+1])); smat[ic,jc] <- mn; smat[ic-1,jc+1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic-1,jc])); smat[ic,jc] <- mn; smat[ic-1,jc] <- mn})
      }
    }
  }
  eval.parent(substitute(object@diffmat <- as(smat, "sparseMatrix")))
})

#' @title Function for diffusion of the Substance matrix
#'
#' @description The generic function \code{diffusePDE} implements the diffusion by the solving diffusion equation.
#' @export
#' @rdname diffusePDE
#'
#' @param object An object of class Substance.
#' @param init_mat A matrix with values to be used by the diffusion.
#' @param gridgeometry A list specifying the geometry of the Arena
#' @param tstep A numeric value giving the time step of integration
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in simEnv() by the function estimate_lrw())
#' @details Partial differential equation is solved to model 2d diffusion process in the arena.
#' @seealso \code{\link{Substance-class}} and \code{\link{diffuseR}}
#' @examples
#' arena <- Arena(n=100, m=100, stir=FALSE, Lx=0.025, Ly=0.025)
#' sub <- Substance(n=100,m=100,smax=0,name='test', difspeed=0.1, 
#'                  gridgeometry=arena@gridgeometry) #initialize test substance
#' sub@diffmat[ceiling(100/2),ceiling(100/2)] <- 40
#' diffusePDE(sub, init_mat=as.matrix(sub@diffmat),
#'            gridgeometry=arena@gridgeometry, tstep=arena@tstep)
setGeneric("diffusePDE", function(object, init_mat, gridgeometry, lrw=NULL, tstep){standardGeneric("diffusePDE")})
#' @export
#' @rdname diffusePDE
setMethod("diffusePDE", "Substance", function(object, init_mat, gridgeometry, lrw=NULL, tstep){
  #init_mat <- as.matrix(object@diffmat)
  if(is.null(lrw)){
    lrw=estimate_lrw(gridgeometry$grid2D$x.N, gridgeometry$grid2D$y.N)}
  D <- object@difspeed*3600 # change unit of diff const to cm^2/h
  solution <- deSolve::ode.2D(y = init_mat, func = get(object@pde), times=c(0,0+tstep), parms = c(gridgeometry=gridgeometry, diffgeometry=object@diffgeometry, boundS=object@boundS),
                     dimens = c(gridgeometry$grid2D$x.N, gridgeometry$grid2D$y.N), method="lsodes", lrw=lrw)#160000
  diff_mat <- matrix(data=solution[2,][-1], ncol=ncol(init_mat), nrow=nrow(init_mat))
  return(diff_mat)
})



#show function for class Substance

setMethod(show, signature(object="Substance"), function(object){
  print(paste('Compound ',object@name,' of class Substance with a total concentration of ',
              sum(object@diffmat)/length(c(object@diffmat)),' mmol per gridcell.',sep=''))
})
