#' Estimates the intensity of the point pattern by using the SPDE method from r-INLA.
#' 
#' Estimates the intensity of the point pattern by using the SPDE method from r-INLA.
#'
#' @param coords ppp object or matrix with x and y coordinates of the observed bombs  
#' @param win observation window, either of class owin or a matrix with the x and y coordinates of the boundary,
#' not neccessary if coords is a ppp object
#' @param mesh (optional) a predefined mesh for the spde model
#' @param fine_mesh logical, if FALSE a coarse mesh will be created, if TRUE a fine mesh will be created, 
#' only used if argument mesh is NULL
#' @param weights (optional) integration weights for the spde model, only used if argument mesh is NULL
#' @param alpha (optional) alpha value for the spde model, only used if argument spde is NULL
#' @param ... additional arguments for the construction of the spde model (see \code{\link[INLA]{inla.spde2.matern}})
#' @param npixel number of pixel per dimension (see \code{\link[spatstat]{spatstat.options}})
#' @importFrom methods as
#' @export
#' @return A list of
#'    \item{ intensest }{ Pixel image with the estimated intensities of the random field. }
#'    \item{ mesh }{ The mesh. }
#' @examples
#' \dontrun{
#' data(craterA)
#' est_spde <- est_intens_spde(coords=craterA)
#' image.plot(list(x=est_spde$intensest$xcol, y=est_spde$intensest$yrow, 
#'                 z=log(t(est_spde$intensest$v))), main="estimated logarithmic intensity")
#' points(craterA)
#' }

est_intens_spde <- function(coords, win=NULL, npixel=50, fine_mesh=FALSE, mesh=NULL, weights=NULL, alpha=2, ...){
  if(!requireNamespace("INLA")){
    warning("This function requires R-INLA, but the R-INLA package is not available on CRAN.\n
         Trying to install it from source(\"https://www.math.ntnu.no/inla/givemeINLA.R\") \n 
See www.r-inla.org for more information!
")
    source("https://www.math.ntnu.no/inla/givemeINLA.R")
  }
  if(is.null(coords)){stop("argument coords is missing, with no default")}
  if(!is.ppp(coords)){
    if(!is.matrix(coords)){stop("coords has to be of type ppp or matrix")}
    if(is.null(win)){stop("argument win is missing, with no default")}
    if(!is.owin(win) & !is.matrix(win)){stop("argument win has to be of class owin or a matrix")}
    n <- nrow(coords)
    xy_coords <- coords
  }
  else{
    xy_coords <- cbind(coords$x,coords$y)
    win <- coords$win
    n <- coords$n
  }
  
  if(is.owin(win)){
    if(win$type == "polygonal"){boundary <- cbind(win$bdry[[1]]$x, win$bdry[[1]]$y)}
    else{boundary <- cbind(win$xrange[c(1,2,2,1)], win$yrange[c(1,1,2,2)])}
  }
  
  spatstat::spatstat.options(npixel=npixel)
  
  if(!is.null(mesh) && attributes(mesh)$class != "inla.mesh"){
    warning("argument mesh is not of class inla.mesh, est_intens_spde creates a new mesh")
  }
  if(is.null(mesh)){
    # create mesh
    xrange <- win$xrange[2]- win$xrange[1]
    yrange <- win$yrange[2]- win$yrange[1]
    mesh_boundary <- INLA::inla.mesh.segment(boundary)
    if(fine_mesh){
      mesh <- INLA::inla.mesh.2d(boundary=mesh_boundary, max.edge=min(xrange,yrange)*c(0.025,0.05), cutoff=min(xrange,yrange)*0.0125)
    }
    else{
      mesh <- INLA::inla.mesh.2d(boundary=mesh_boundary, max.edge=min(xrange,yrange)*c(0.1,0.2), cutoff=min(xrange,yrange)*0.05)
    }
  }
  
  spde <- INLA::inla.spde2.matern(mesh=mesh, alpha=alpha,...)
  
  if(is.null(weights)){
    crater_region <- as(boundary, 'gpc.poly')
    tiles <- tile.list(deldir(mesh$loc[,1], mesh$loc[,2]))
    weights <- sapply(tiles, function(p) area.poly(intersect(as(cbind(p$x, p$y), 'gpc.poly'), crater_region)))
  }
  
  mesh_n <- mesh$n
  y <- rep(0:1, c(mesh_n, n))
  e <- c(weights, rep(0, n))
  A <- rBind(Diagonal(mesh_n, rep(1, mesh_n)), INLA::inla.spde.make.A(mesh, xy_coords))
  stack <- INLA::inla.stack(data=list(y=y, e=e), A=list(1,A), tag='crater', effects=list(list(beta0=rep(1,mesh_n+n)), list(i=1:mesh_n)))
  
  
  inla_result <- INLA::inla(y ~ 0 + beta0 + f(i, model=spde), family="poisson", data=INLA::inla.stack.data(stack), 
                      control.predictor=list(A=INLA::inla.stack.A(stack)), E=INLA::inla.stack.data(stack)$e)
  random_field <- INLA::inla.spde2.result(inla_result, 'i', spde) 
  projector <- INLA::inla.mesh.projector(mesh, dims=c(npixel,npixel), xlim=c(win$xrange[1],win$xrange[2]), ylim=c(win$yrange[1],win$yrange[2]))
  intensity_spde <- INLA::inla.mesh.project(projector, random_field$summary.value$mean) + mean(inla_result$marginals.fix[[1]][,1])
  
  pixel_in_poly <- matrix(in.poly(projector$lattice$loc, rbind(boundary,boundary[1,])), ncol=npixel)
  pixel_in_poly[pixel_in_poly==0] <- NA
  intensity_in_poly <- t(pixel_in_poly*exp(intensity_spde))
  
  output <- list(intensest=as.im(X=intensity_in_poly, W=win), mesh=mesh)
  
  return(output)
}
