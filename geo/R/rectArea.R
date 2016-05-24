#' Given rectangle code return area in square kilometers or nautical miles
#' 
#' Rectangle area is returned through a call to
#' \code{\link{rectPeri}}-functions and \code{\link{geoarea}}.
#' 
#' @name rectArea
#' @aliases rA srA mrA drA
#' @param r,sr,mr,dr rectangle code, as in \code{\link{rect2deg}} and
#' \code{\link{deg2rect}}.
#' @param scale \code{nmi, km}, default \code{nmi} returns values in square
#' nautical miles for all except \code{drA} returns area in square kilometers.
#' @param dlat,dlon Dimensions of latitude and longitude given in minutes and
#' degrees for \code{mrPeri} and \code{drPeri}, respectively.
#' @return Rectangle area in square nautical miles or kilometers.
#' @note Unit \code{nmi} is used for historical/acoustical (sA) reasons.
#' @seealso \code{\link{rectPeri}}, \code{\link{deg2rect}},
#' \code{\link{geoarea}}.
#' @keywords arith manip
#' @examples
#' 
#'   srA(7121)
#'   srA(7121, "km")
#'   srA(7121, "km")/1.852^2
#'   srA(7121, "km")
#'   rA(712)
#'   srA(7121) + srA(7122) + srA(7123) + srA(7124)
#'

#' @export rA
#' @rdname rectArea
rA <-
function (r, scale = "nmi") 
{
  if(!(scale == "nmi" | scale == "km")) 
    stop("Unit square (nautical) miles or kilometers only")
  A <- sapply(r, function(x) geoarea(rPeri(x)))
  if(scale == "nmi") {
    A/1.852^2
  }
  else {
    A
  }
}

#' @export srA
#' @rdname rectArea
srA <-
function (sr, scale = "nmi") 
{
  if(!(scale == "nmi" | scale == "km")) 
    stop("Unit square nautical miles or kilometers only")
  A <- sapply(sr, function(x) geoarea(srPeri(x)))
  if(scale == "nmi") {
    A/1.852^2
  }
  else {
    A
  }
}

#' @export mrA
#' @rdname rectArea
mrA <-
function (mr, dlat = 5, dlon = 10, scale = "nmi") 
{
  if(!(scale == "nmi" | scale == "km")) 
    stop("Unit square nautical miles or kilometers only")
  A <- sapply(mr, function(x, dlat, dlon) 
    geoarea(mrPeri(x, dlat = dlat, dlon = dlon)), 
      dlat = dlat, dlon = dlon)
  if(scale == "nmi") A/1.852^2
    else A
}

#' @export drA
#' @rdname rectArea
drA <-
function (dr, dlat = 1, dlon = 2, scale = "km") 
{
  if(!(scale == "nmi" | scale == "km")) 
    stop("Unit square nautical miles or kilometers only")
  A <- sapply(dr, function(x, dlat, dlon) 
    geoarea(drPeri(x, dlat = dlat, dlon = dlon)),
      dlat = dlat, dlon = dlon)
  if(scale == "nmi") A/1.852^2
    else A
}

