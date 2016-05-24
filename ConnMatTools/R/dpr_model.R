# Implements methods from the following references:
# 
# @references Grüss, A., Kaplan, D. M., and Lett, C. 2012. Estimating local 
#   settler-recruit relationship parameters for complex spatially explicit 
#   models. Fisheries Research, 127-128: 34-39.
# @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal 
#   per recruit: An efficient method for assessing sustainability in marine 
#   reserve networks. Ecological Applications, 16: 2248-2263.
# @references White, J. W. 2010. Adapting the steepness parameter from 
#   stock-recruit curves for use in spatially explicit models. Fisheries 
#   Research, 102: 330-334.
# @references Grüss A, Kaplan DM, Hart DR (2011) Relative Impacts of Adult
#   Movement, Larval Dispersal and Harvester Movement on the Effectiveness of
#   Reserve Networks. PLoS ONE 6:e19960
# @references Beverton RJH, Holt SJ (1957) On the dynamics of exploited fish
#   populations. H.M.S.O., London. 533 pp.
NULL

#' Correction for slope of settler-recruit relationship
#' 
#' This function corrects the slope of the settler-recruit curve so that the 
#' collapse point of the spatially-explicit population model corresponding to 
#' the connectivity matrix agrees with that of the global non-spatially-explicit
#' model.  Uses the method in White (2010).
#' 
#' @param conn.mat a square connectivity matrix.
#' @param slope slope at the origin of the settler-recruit relationship.  Only 
#'   interesting to fix this argument if it is a vector of length = 
#'   \code{dim(conn.mat)[2]} (i.e., if the slope varies among sites and one
#'   wants to globally scale all slopes so that the collapse point matches the
#'   global collapse point).
#' @param natural.LEP value of lifetime-egg-production (LEP), also known as 
#'   eggs-per-recruit, in the absence of fishing.  Can be a vector of length = 
#'   \code{dim(conn.mat)[2]}.  Defaults to 1.
#' @param critical.FLEP Fraction of natural.LEP at which collapse occurs. 
#'   Defaults to 0.35.
#' @param use.arpack Boolean determining if calculation is to be done with 
#'   \code{\link[igraph]{arpack}} function from the \link[igraph]{igraph} package. This is much 
#'   quicker for large matrices, but requires \link[igraph]{igraph}. Defaults to TRUE, 
#'   but will use eigen instead if \link[igraph]{igraph} is not found.
#'   
#' @return The slope argument corrected so that collapse happens when LEP is 
#'   critical.FLEP * natural.LEP.
#'   
#' @references White, J. W. 2010. Adapting the steepness parameter from 
#'   stock-recruit curves for use in spatially explicit models. Fisheries 
#'   Research, 102: 330-334.
#'   
#' @seealso See also \code{\link{eigs}}, \code{\link[igraph]{arpack}}
#' 
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
#' @include eigs.R
settlerRecruitSlopeCorrection <- function(conn.mat,slope=1,natural.LEP=1,
                                          critical.FLEP=0.35,use.arpack=TRUE) {
  if (class(conn.mat) != "matrix")
    stop("Input conn.mat must be a matrix.")
  
  if (length(natural.LEP)>1)
    C = conn.mat %*% (natural.LEP * critical.FLEP)
  else
    C = conn.mat * (natural.LEP * critical.FLEP)
  
  if (length(slope)>1)
    C = diag(slope) %*% C
  else
    C = slope * C
  
  v = eigs(M=C,nev=1,which="LM",
           use.arpack=use.arpack)$values

  return(slope/Mod(v))
}

#' Beverton-Holt settler-recruit relationship
#' 
#' Calculates recruitment based on the settler-recruit relationship from 
#' Beverton & Holt (1957): \code{ slope * settlers / (1+slope*settlers/Rmax) }
#' 
#' \code{slope} and \code{Rmax} can both either be scalars or vectors of the
#' same length as \code{S}.
#' 
#' @param S a vector of settlement values, 1 for each site.
#' @param slope slope at the origin of the settler-recruit relationship.  Can be
#'   a vector of same length as \code{S}.
#' @param Rmax maximum recruitment value.
#'   
#' @return A vector of recruitment values.
#'   
#' @references Beverton RJH, Holt SJ (1957) On the dynamics of exploited fish 
#'   populations. H.M.S.O., London. 533 pp.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
BevertonHolt <- function(S,slope=1/0.35,Rmax=1) {
  return( (slope * S) / (1 + slope * S / Rmax ) )
}

#' Hockey-stick settler-recruit relationship
#' 
#' Calculates recruitment based on a settler-recruit relationship that increases
#' linearly until it reaches a maximum values.
#' 
#' \code{slope} and \code{Rmax} can both either be scalars or vectors of the 
#' same length as \code{S}.
#' 
#' @param S a vector of settlement values, 1 for each site.
#' @param slope slope at the origin of the settler-recruit relationship.  Can be
#'   a vector of same length as \code{S}.
#' @param Rmax maximum recruitment value.
#'   
#' @return A vector of recruitment values.
#'   
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248-2263.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
hockeyStick <- function(S,slope=1/0.35,Rmax=1) {
  R = slope * S
  R[ R>Rmax] = Rmax
  return( R )
}

#' Population dynamics model based on lifetime-egg-production
#' 
#' This function implements the marine population dynamics model described in 
#' Kaplan et al. (2006).  This model is most appropriate for examining 
#' equilibrium dynamics of age-structured populations or temporal dynamics of 
#' semelparous populations.
#' 
#' @param LEP a vector of lifetime-egg-production (LEP; also known as 
#'   eggs-per-recruit (EPR)) for each site.
#' @param conn.mat a square connectivity matrix.  \code{dim(conn.mat) = 
#'   rep(length(LEP),2)}
#' @param recruits0 a vector of initial recruitment values for each site.
#' @param timesteps a vector of timesteps at which to record egg production, 
#'   settlement and recruitment.
#' @param settler.recruit.func a function to calculate recruitment from the 
#'   number of settlers at each site.  Defaults to \code{\link{hockeyStick}}.
#' @param \dots additional arguments to settler.recruit.func.  Typically
#'   \code{Rmax} and \code{slope}.
#'   
#' @return A list with the following elements: \item{eggs}{egg production for
#'   the timesteps in \code{timesteps}}
#'   
#'   \item{settlers}{Similar for settlement}
#'   
#'   \item{recruits}{Similar for recruitment}
#'   
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248-2263.
#'   
#' @seealso See also \code{\link{BevertonHolt}}, \code{\link{hockeyStick}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.dpr_model.R
#' @export
DispersalPerRecruitModel <- 
  function(LEP,conn.mat,recruits0,
           timesteps=10, settler.recruit.func=hockeyStick,...) {
    
    if (class(conn.mat) != "matrix")
      stop("Input conn.mat must be a matrix.")
    
    r = recruits0
    
    eggs = matrix(NA,nrow=dim(conn.mat)[2],ncol=length(timesteps))
    settlers = eggs
    recruits = settlers
    
    for (t in 1:max(timesteps)) {
      e = r * LEP
      s = conn.mat %*% e
      r = settler.recruit.func(s,...)
      
      I = which( t == timesteps )
      if (length(I)>0) {
        eggs[,I] = e
        settlers[,I] = s
        recruits[,I] = r        
      }
    }
    
    return(list(eggs=eggs,settlers=settlers,recruits=recruits))
  }

#' Extended DPR population dynamics model to include homerange movement
#' 
#' This function implements the marine population dynamics model described in 
#' Gruss et al. (2011).  The model is an extension of the dispersal-per-recruit 
#' model in Kaplan et al. (2006) to include movement in a homerange and a 
#' gravity model for fishing effort redistribution.
#' 
#' @param larval.mat a square larval connectivity matrix.  \code{dim(larval.mat)
#'   = rep(length(recruits0),2)}
#' @param adult.mat a square adult homerange movement matrix. 
#'   \code{dim(adult.mat) = rep(length(recruits0),2)}. adult.mat must
#'   be properly normalized so that each column sums to 1.
#' @param recruits0 a vector of initial recruitment values for each site.
#' @param f0 a vector of initial real fishing mortalities for each site.
#' @param timesteps a vector of timesteps at which to record egg production, 
#'   settlement and recruitment.
#' @param settler.recruit.func a function to calculate recruitment from the 
#'   number of settlers at each site.  Defaults to \code{\link{hockeyStick}}.
#' @param LEP.of.f a function that returns lifetime-egg-productions given a 
#'   vector of fishing rates.
#' @param YPR.of.f a function that returns yields-per-recruit given a vector of 
#'   fishing rates.
#' @param gamma exponent for the gravity model.  Defaults to 0, i.e., no gravity
#'   model.
#' @param gravity.ts.interval number of timesteps between updates of gravity 
#'   model.  Defaults to 1, i.e., every timestep.
#' @param \dots additional arguments to settler.recruit.func.
#'   
#' @return A list with the following elements:
#'   
#'   \item{eggs}{egg production for the timesteps in \code{timesteps}}
#'   
#'   \item{settlers}{Similar for settlement}
#'   
#'   \item{recruits}{Similar for recruitment}
#'   
#'   \item{fishing.mortality}{Real spatial distribution of fishing mortality}
#'   
#'   \item{effective.fishing.mortality}{Effective fishing mortality taking into 
#'   account adult movement}
#'   
#'   \item{yield}{Real spatial distribution of yield}
#'   
#'   \item{effective.yield}{Effective yield indicating where fish biomass caught
#'   originates from}
#'   
#' @references Grüss A, Kaplan DM, Hart DR (2011) Relative Impacts of Adult 
#'   Movement, Larval Dispersal and Harvester Movement on the Effectiveness of 
#'   Reserve Networks. PLoS ONE 6:e19960
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248-2263.
#'   
#' @seealso See also \code{\link{BevertonHolt}}, \code{\link{hockeyStick}}, 
#'   \code{\link{DispersalPerRecruitModel}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
DPRHomerangeGravity <- 
  function(larval.mat,adult.mat,recruits0,f0,
           timesteps=10, settler.recruit.func=hockeyStick,
           LEP.of.f=function(f) 1-f, YPR.of.f=function(f) f, 
           gamma=0, gravity.ts.interval=1, ...) {
    
    if (class(larval.mat) != "matrix")
      stop("Input larval.mat must be a matrix.")
    if (class(adult.mat) != "matrix")
      stop("Input adult.mat must be a matrix.")
    
    r = recruits0
    f = f0
    ftot = sum(f)
    
    eggs = matrix(NA,nrow=dim(larval.mat)[2],ncol=length(timesteps))
    settlers = eggs
    recruits = settlers
    fishing.mortality=eggs
    effective.fishing.mortality=fishing.mortality
    yield=eggs
    effective.yield=yield
    
    for (t in 1:max(timesteps)) {
      feff = t(t(f) %*% adult.mat)
      
      LEP = LEP.of.f(feff)
      
      e = r * LEP
      s = larval.mat %*% e
      r = settler.recruit.func(s,...)
      
      YPR = YPR.of.f(feff)
      Yeff = YPR * r
      
      if (prod(dim(adult.mat))>1)
        Y = (adult.mat %*% (Yeff/feff)) * f
      else
        Y = Yeff
      
      I = which( t == timesteps )
      if (length(I)>0) {
        eggs[,I] = e
        settlers[,I] = s
        recruits[,I] = r        
        fishing.mortality[,I] = f        
        effective.fishing.mortality[,I] = feff        
        yield[,I] =         
        effective.yield[,I] = feff        
      }
      
      if ((gamma>0) && (t%%gravity.ts.interval == 0)) {
        YY = Y^(1/gamma)
        f = ftot * YY / sum(YY)
      }
    }
    
    return(list(eggs=eggs,settlers=settlers,recruits=recruits,
                fishing.mortality=fishing.mortality,
                effective.fishing.mortality=effective.fishing.mortality,
                yield=yield,effective.yield=effective.yield))
  }

#' Uniform Laplacian connectivity matrix
#' 
#' This function generates a connectivity matrix that is governed by a Laplacian
#' distribution: \code{D[i,j]=exp(abs(x[i]-y[i]-shift)/disp.dist)/2/disp.dist}
#' 
#' The \code{boundary} argument can have the following different values: 
#' "nothing" meaning do nothing special with boundaries; "conservative" meaning force
#' columns of matrix to sum to 1; and "circular" meaning wrap edges.
#' 
#' @param num.sites number of sites.  Sites are assumed to be aligned on a 
#'   linear coastline.
#' @param disp.dist dispersal distance in "site" units (i.e., 1 site = 1 unit of
#'   distance)
#' @param shift advection distance in "site" units.  Defaults to 0.
#' @param boundaries string indicating what to do at boundaries.  Defaults to 
#'   "nothing".  Possible values include: "nothing", "conservative" and
#'   "circular"
#'   
#' @return A square connectivity matrix
#'   
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248-2263.
#'   
#' @seealso See also \code{\link{DispersalPerRecruitModel}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.laplacianConnMat.R
#' @export
laplacianConnMat <- function(num.sites,disp.dist,shift=0,boundaries="nothing") {
  # Function integrates exp(-exponent*abs(x))*exponent/2 from start to end.  start <= end
  expint <- function(exponent,start,end) {
    d <- matrix(0,nrow=dim(start)[1],ncol=dim(start)[2])
    
    if (is.infinite(exponent)) {
      if (exponent>0) 
        exponent = 1e20
      else
        exponent = -1e20
    }
    
    posS = start>=0
    posE = end>= 0
    
    pos = posS & posE
    neg = !posS & !posE
    mixed = !pos & !neg
    
    d[pos] = 0.5 * (exp(-exponent*start[pos]) - exp(-exponent*end[pos]))
    d[neg] = 0.5 * (exp(exponent*end[neg]) - exp(exponent*start[neg]))
    d[mixed] = 1-0.5 * (exp(exponent*start[mixed])+exp(-exponent*end[mixed]))
    
    return(d)
  }
  
  x = matrix(1:num.sites,nrow=num.sites,ncol=num.sites)
  x = x - t(x) - shift
  
  dm = expint(1/disp.dist,x-0.5,x+0.5)

  if (boundaries=="conservative") {
    dm = dm %*% diag(apply(dm,2,sum))
  } else if(boundaries=="circular") {
    for (k in c(-10:-1,1:10))
      dm = dm + expint(1/disp.dist,x+k*num.sites-0.5,x+k*num.sites+0.5)
    
    # If circular, force conservative
    dm = dm %*% diag(apply(dm,2,sum))    
  } else if (boundaries=="nothing") {}
  else stop("Bad boundaries argument")
  
  return(dm)
}
