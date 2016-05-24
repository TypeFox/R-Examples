#' As ("trip", other-classes)
#'
#' Coercing \code{trip} objects to other classes.
#'
#' @name as.Other
# aliases as.SpatialLinesDataFrame.trip
# section Methods:
#
# ##\describe{
#
# ##\item{coerce}{\code{signature(from="trip", to="SpatialLinesDataFrame")}}
# ##}
setAs("trip", "SpatialLinesDataFrame", function(from) {
  split.from <- split(from, from[[getTORnames(from)[2]]])
  sdf <- suppressWarnings(summary(from))
  df <- data.frame(tripID=sdf$tripID, tripStart=sdf$tmins,
                   tripEnd=sdf$tmaxs,
                   tripDur=as.vector(sdf$tripDurationSeconds),
                   row.names=sdf$tripID)
  lns <- vector("list", nrow(df))
  for (i in 1:length(lns)) {
    lns[[i]] <- Lines(list(Line(coordinates(split.from[[i]]))),
                      ID=sdf$tripID[i])
  }
  SpatialLinesDataFrame(SpatialLines(lns,
                                     proj4string=CRS(proj4string(from))),
                        df)
})



setAs("trip", "ltraj", function(from) {
  if(!requireNamespace("adehabitatLT")) stop("adhabitatLT not available")
  tor <- getTORnames(from)
  crds <- coordinates(from)
  adehabitatLT::as.ltraj(as.data.frame(crds), date=from[[tor[1]]],
                         id=from[[tor[2]]], typeII=TRUE, slsp="remove")
})




## @importClassesFrom maptools owin ppp psp
#' @importFrom spatstat as.ppp
#' @importFrom maptools as.ppp.SpatialPointsDataFrame
#' @param X \code{trip} object.
#' @param fatal Logical value, see Details of \code{\link[spatstat]{as.ppp}}
#' @return ppp object
#' @rdname as.Other
#' @method as.ppp trip
#' @examples
#' \dontrun{
#'   ## Continuing the example from '?trip-methods:
#' utils::example("trip-methods", package="trip",
#'            ask=FALSE, echo=FALSE)
#'  as(tr, "ppp")
#' }
as.ppp.trip <- function(X, ..., fatal) {
  as.ppp.SpatialPointsDataFrame(X)
}
setAs("trip", "ppp", function(from) as.ppp.trip(from))

#' @export
#' @importFrom spatstat as.psp
#' @param x \code{trip} object
#' @param from see \code{\link[spatstat]{as.psp}} for that method.
#' @param to See \code{\link[spatstat]{as.psp}}.
#' @return psp object
#' @rdname as.Other
#' @method as.psp trip
#' @examples
#' \dontrun{
#'  ## Continuing the example from '?trip-methods:
#' utils::example("trip-methods", package="trip",
#'            ask=FALSE, echo=FALSE)
#'  as.psp.trip(tr)
#' }
as.psp.trip <- function(x, ..., from, to) {
  split.X <- split(x, x[[getTORnames(x)[2]]])
  ow <- owin(bbox(x)[1,], bbox(x)[2,])
  as.psp.trip1 <- function(this, ow=NULL) {
    if (is.null(ow)) ow <- owin(bbox(this)[1,], bbox(this)[2,])
    tor <- getTORnames(this)
    cc <- coordinates(this)
    xs <- coordinates(this)[, 1]
    ys <- coordinates(this)[, 2]
    dt <- diff(unclass(this[[tor[1]]]))
    psp(xs[-length(xs)], ys[-length(ys)],
        xs[-1], ys[-1], window=ow, marks=dt)
  }
  do.call("superimpose", lapply(split.X, as.psp.trip1, ow=ow))
}
setAs("trip", "psp", function(from) as.psp.trip(from))



#' Break a trip into its component line segments
#'
#' Function to create a SpatialLinesDataFrame from a trip object, resulting in
#' a line segment for each implicit segment along the tracks. The object stores
#' the start and end times, duration and the ID of the segment.
#'
#' @param ... reserved for future methods
#' @return SpatialLinesDataFrame
#' @examples
#' ## Continuing the example from '?trip-methods:
#' utils::example("trip-methods", package="trip",
#'            ask=FALSE, echo=FALSE)
#' spldf <- explode(tr)
#' summary(tr)
#' @return SpatialLinesDataFrame object with each individual line segment identified by start/end time and trip ID
#' @rdname as.Other
#' @export explode
explode <- function(x, ...) {
  tor <- getTORnames(x)
  id <- x[[tor[2]]]
  xs <- split(x, id)
  df <- do.call("rbind",
                lapply(xs, function(x) {
                  n <- nrow(x)
                  tms <- x[[tor[1]]]
                  data.frame(starttime = tms[-n], endtime = tms[-1], timedur = diff(unclass(tms)), id = x[[tor[2]]][-1])
                }
                )
  )
  Linelist <- vector("list", nrow(df))
  cnt <- 0
  for (i in seq_along(xs)) {
    this.x <- xs[[i]]
    this.coords <- coordinates(this.x)
    for (j in seq_len(nrow(this.x)-1)) {
      cnt <- cnt + 1
      Linelist[[cnt]] <- Lines(list(Line(this.coords[j:(j+1), ])), rownames(df)[cnt])
    }

  }
  splines <- SpatialLines(Linelist, proj4string = CRS(proj4string(x)))
  SpatialLinesDataFrame(splines, df)
}


# setMethod("lines", signature(x="trip"),
#           function(x,
#                    col=hsv(seq(0, 0.9, length = length(unique(x[[getTORnames(x)[2]]]))),
#                            0.8, 0.95),
#                    ...) {
#             x <- .explode(x)
#             col <- heat_hcl(25, h = c(0, -100), l = c(55, 40), c = c(40, 80), power = 3)
#             times <- x$time
#             val <- scl(unclass(times)) * (length(col)-1) + 1
#
#             plot(x,  col=col[val], add=TRUE, ...)
#
#           })
