#' California Water Olympics
#' 
#' Digitized results from the California Water Olympics. The variables are as follows:
#' \itemize{
#'   \item t The time (in seconds) since the start of the model run.
#'   \item Q The flow rate [\eqn{ft^3 s^{-1}}].
#'   \item x The distance downstream [\eqn{ft}] at which the hydrograph was recorded.
#' }
#' The data can be used to validate numerical solutions to flood wave routing for a channel 
#' under the following conditions:
#' \itemize{
#'   \item Channel width is 100 feet.
#'   \item Channel slope is 0.001.
#'   \item Channel extent is 150,000 feet.
#'   \item Channel roughness (Manning's n) is 0.045.
#'   \item Channel sideslope is 0 (rectangular channel).
#'   \item Initial flow rate is 250 cfs.
#'   \item Upstream boundary condition is defined as 
#'     \deqn{Q(t < 9000) = 250 + \frac{750}{\pi}(1 - \cos{\frac{\pi t}{4500}})} 
#'     \deqn{Q(t >= 9000) = 250}
#' }
#' @docType data
#' @keywords datasets
#' @name waterolympics
#' @usage data(waterolympics)
#' @format A data frame with 40 rows and 3 variables
#' @references Sobey, Rodney. "H11: Hydrograph Routing." 
#'   Review of One-Dimensional Hydrodynamic and Transport Models. Bay-Delta 
#'   Modeling Forum, 15 June 2001. Web. 13 Mar. 2015. 
#'   <http://www.cwemf.org/1-DReview/>.
NULL


#' @title Flood wave routing
#' @description Route a flood wave down a prismatic channel.
#' @param So Channel slope [\eqn{L L^{-1}}].
#' @param n Manning's roughness coefficient.
#' @param Cm Unit conversion coefficient for Manning's equation. For SI units, Cm = 1.
#' @param g Gravitational acceleration [\eqn{L T^{-2}}].
#' @param B Channel bottom width [\eqn{L}].
#' @param SS Channel sideslope [\eqn{L L^{-1}}].
#' @param initial.condition The initial flow rate [\eqn{L^3 T^{-1}}], assumed constant 
#'   throughout the channel.
#' @param boundary.condition Vector specifying the upstream boundary condition 
#'   for the full duration of the model. If \code{engine = "Kinematic"}, values are 
#'   assumed to be flow [\eqn{L^3 T^{-1}}]. If \code{engine = "Dynamic"}, the form of the 
#'   boundary condition is determined by the argument \code{boundary.type}.
#' @param downstream.condition Only used if \code{engine = "Dynamic"}. Vector specifying 
#'   the upstream boundary condition for the full duration of the model. Must be the same
#'   length as \code{boundary.condition}.
#' @param timestep Temporal resolution of the model. Also the assumed time interval [\eqn{T}] 
#'   between elements of \code{boundary.condition} and \code{downstream.condition}. 
#'   The user is responsible for ensuring numerical stability.
#' @param spacestep the spatial resolution of the model, interpreted as the distance [\eqn{L}]  
#'   between nodes in the model domain. The user is responsible for ensuring numerical stability.
#' @param numnodes The number of nodes used to discretize the channel. The total channel extent is
#'   computed as \code{spacestep*(numnodes - 1)}.
#' @param monitor.nodes the nodes to be monitored every time step. Specified as a vector of node 
#'   indices, with 1 being the upstream boundary and \code{numnodes} being the downstream boundary.
#' @param monitor.times the time steps at which to monitor every node. Specified as a vector of 
#'   indices of \code{boundary.condition}. Defaults to five equally-spaced time steps including 
#'   the first and last time steps.
#' @param engine The engine to be used for routing the flood wave. 
#'   May be either "Kinematic" or "Dynamic".
#' @param scheme Only used if \code{engine = "Dynamic"}. Specifies whether to use the 
#'   Lax Diffusive scheme or the MacCormack predictor-corrector scheme.
#' @param boundary.type Only used if \code{engine = "Dynamic"}. Specifies what boundary data
#'   is supplied. Possible characters are If \code{boundary.type = "QQ"}, both \code{boundary.condition}
#'   and \code{downstream.condition} are assumed to be flows [\eqn{L^3 T^{-1}}]. If 
#'   \code{boundary.type = "Qy"} the upstream boundary is assumed to be flow
#'   while the downstream boundary is assumed to be depth [\eqn{L}]. Other possibilities
#'   are \code{"yQ"} and \code{"yy"}. 
#' @return data.frame with columns:
#'   \item{step}{Time step.}
#'   \item{node}{Node index.}
#'   \item{time}{Time since start.}
#'   \item{distance}{Downstream distance.}
#'   \item{flow}{Flow rate.}
#'   \item{depth}{Flow depth.}
#'   \item{velocity}{Flow velocity.}
#'   \item{area}{Flow area.}
#'   \item{monitor.type}{Row refers to a monitored node ("node") or timestep ("timestep").}
#' @details Provides implementations of a Kinematic Wave Model (KWM) and 
#'   a Dynamic Wave Model (DWM) with the choice of two numerical schemes. The MacCormack 
#'   scheme is a second-order accurate predictor-corrector scheme that provides efficient 
#'   flood wave routing. The Lax diffusive scheme can be used to obtain smooth solutions for
#'   problems with discontinuities in the boundary conditions, e.g. sudden gate closures.
#'   The DWM implementation uses the Method of Characteristics (MOC) to compute the flow 
#'   regime at the model boundaries, and allows the user to specify boundaries in terms of 
#'   depths and/or flows. the KWM implementation assumes the normal depth at the upstream
#'   boundary and is only first-order accurate.
#' @examples
#' \dontrun{ 
#' # kinematic wave routing
#' times = seq(0, 30000, by = 25)
#' floodwave = ifelse(times >= 9000, 250,
#'   250 + (750/pi)*(1 - cos(pi*times/(60*75))))
#' route_wave(0.001, 0.045, 1.486, 32.2, 100, 0, initial.condition = 250, 
#'   boundary.condition = floodwave, timestep = 25, spacestep = 50, 
#'   numnodes=301, monitor.nodes = c(1, 101, 201, 301), 
#'   monitor.times = seq(1, length(times), by = 10), engine = "Kinematic")
#' # dynamic wave routing with zero-gradient downstream condition using MacCormack scheme
#' route_wave(0.001, 0.045, 1.486, 32.2, 100, 0, initial.condition = 250, 
#'   boundary.condition = floodwave, downstream.condition = rep(-1, length(times)), 
#'   timestep = 25, spacestep = 500, numnodes = 31, engine = "Dynamic", 
#'   scheme = "MacCormack", monitor.nodes = c(1, 11, 21, 31), 
#'   monitor.times = seq(1, length(times), by = 10))
#' # mixed boundary conditions (sudden gate closure) using Lax scheme
#' lax = route_wave(0.00008, 0.013, 1, 9.81, 6.1, 1.5, 
#'   initial.condition = 126, boundary.condition = rep(5.79, 2001), 
#'   downstream.condition = rep(0, 2001), timestep = 1, spacestep = 10, 
#'   numnodes = 501, monitor.nodes = c(1, 151, 251, 301, 501), 
#'   monitor.times = c(1, 501, 1001, 1501, 2001), 
#'   engine="Dynamic", scheme="Lax", boundary.type="yQ")
#' # extract data for a monitored point
#' require(dplyr)
#' filter(lax, monitor.type == "node", node == 151)
#' }
#' @importFrom utils stack
#' @export
route_wave = function(So, n, Cm, g, B, SS, 
  initial.condition, boundary.condition, downstream.condition, 
  timestep, spacestep, numnodes, monitor.nodes, monitor.times, 
  engine=c('Dynamic', 'Kinematic'), scheme = c("MacCormack", "Lax"),
  boundary.type=c("QQ", "Qy", "yQ", "yy")){ 
  engine = match.arg(engine, c('Dynamic', 'Kinematic'))
  boundary.type = match.arg(boundary.type, c("QQ", "Qy", "yQ", "yy"))
  # ensure boundary conditions and initial conditions are included in 
  # monitoring data---boundary condition included by default
  monitor.nodes = unique(as.integer(c(1, monitor.nodes)))
  monitor.times = unique(as.integer(c(1, monitor.times)))
  # do some checks
  if(any(monitor.nodes > numnodes))
    stop("Monitor.nodes index out of bounds. Downstream boundary node index is ", numnodes, ".")
  if(any(monitor.times > length(boundary.condition)))
    stop("Monitor.times index out of bounds. Final time step is ", length(boundary.condition), ".")
  if(!missing(downstream.condition))
    if((length(boundary.condition) != length(downstream.condition)))
      stop("lengths of boundary.condition and downstream.condition do not match.")
  if(engine == 'Kinematic'){
    reslist = kinematic_wave(So, n, Cm, g, B, SS, numnodes, boundary.condition,
      initial.condition, timestep, spacestep, monitor.nodes - 1L, monitor.times - 1L)
    scheme = NULL
  } else {
    scheme = match.arg(scheme, c("MacCormack", "Lax"))
    if(scheme == "MacCormack"){
      reslist = characteristic_wave(So, n, Cm, g, B, SS, numnodes, 
        boundary.condition, downstream.condition, initial.condition, 
        timestep, spacestep, monitor.nodes - 1L, monitor.times - 1L, 
        boundary.type)
    } else{
      reslist = diffusive_wave(So, n, Cm, g, B, SS, numnodes, 
        boundary.condition, downstream.condition, initial.condition, 
        timestep, spacestep, monitor.nodes - 1L, monitor.times - 1L, 
        boundary.type)
    }
  }
  # extract monitor data
  stackmat = function(m){  
    d = cbind(rid = rep(rownames(m), ncol(m)), stack(as.data.frame(m))[2:1])
    d[[1]] = as.numeric(levels(d$rid))[d$rid]
    d[[2]] = as.numeric(levels(d$ind))[d$ind]
    d
  }
  mpmat = reslist[["monitorpoints"]]
  mtmat = reslist[["monitortimes"]]
  colnames(mpmat) = as.integer(monitor.nodes)
  rownames(mpmat) = as.integer(seq(nrow(mpmat)))
  colnames(mtmat) = as.integer(seq(ncol(mtmat)))
  rownames(mtmat) = as.integer(monitor.times)
  mpdf = stackmat(mpmat)
  mtdf = stackmat(mtmat)
  # add extra data
  mpnames = c("mpdepth", "mpvelocity", "mparea")
  mtnames= c("mtdepth", "mtvelocity", "mtarea")
  addnames = c("depth", "velocity", "area")
  for(nm in mpnames){
    thismat = reslist[[nm]]
    colnames(thismat) = as.integer(monitor.nodes)
    rownames(thismat) = as.integer(seq(nrow(thismat)))
    thisdf = stackmat(thismat)
    mpdf = cbind(mpdf, thisdf[, 3])	
  }
  for(nm in mtnames){
    thismat = reslist[[nm]]
    colnames(thismat) = as.integer(seq(ncol(thismat)))
    rownames(thismat) = as.integer(monitor.times)
    thisdf = stackmat(thismat)
    mtdf = cbind(mtdf, thisdf[, 3])	
  }
  names(mpdf) = c("step", "node", "flow", addnames)  
  names(mtdf) = c("step", "node", "flow", addnames)  
  mpdf['monitor.type'] = "node"
  mtdf['monitor.type'] = "timestep"
  allres = rbind(mpdf, mtdf)
  allres["distance"] = (allres$node - 1)*spacestep
  allres["time"] = (allres$step - 1)*timestep
  retnames = c("step", "node", "time", "distance", "flow", 
    addnames, "monitor.type")  
  retobj = allres[retnames]
  class(retobj) = c("rivr", "data.frame")  
  attr(retobj, "call") = match.call()
  attr(retobj, "simtype") = "usf"
  attr(retobj, "modspec") = list(dx = spacestep, dt = timestep, 
    nodes.monitored = monitor.nodes, 
    steps.monitored = monitor.times)
  attr(retobj, "engine") = paste(c(engine, scheme), collapse = ":")
  attr(retobj, "channel.geometry") = list(So = So, n = n, B = B, SS = SS)
  return(retobj)
} 
