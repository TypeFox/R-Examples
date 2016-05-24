#' Continuous position scales (T, L & R).
#' 
#' @rdname scale_tern_continuous
#' @inheritParams ggplot2::scale_x_continuous
#' @param ... not used
#' @export
scale_T_continuous <- function(name         = waiver(), 
                               limits       = c(0,1), 
                               breaks       = getBreaks(limits,TRUE),
                               minor_breaks = getBreaks(limits,FALSE),
                               labels       = 100*breaks,
                               expand       = waiver(),...) {
  sc <- continuous_scale(
    c("T"),
    "tern_T", identity, name = name, breaks = breaks,
    minor_breaks = minor_breaks, labels = labels, limits = limits,
    expand = waiver(), oob = censor, na.value = NA_real_, trans = "identity",
    guide = "none"
  )
  # TODO: Fix this hack. We're reassigning the parent ggproto object, but this
  # object should in the first place be created with the correct parent.
  sc$super <- ScaleContinuousPosition
  class(sc) <- class(ScaleContinuousPosition)
  sc
}


#' @name scale_tern_continuous
#' @rdname scale_tern_continuous
#' @export
scale_L_continuous <- function(name         = waiver(), 
                               limits       = c(0,1), 
                               breaks       = getBreaks(limits,TRUE),
                               minor_breaks = getBreaks(limits,FALSE),
                               labels       = 100*breaks,
                               expand       = waiver(),...) {
  sc <- continuous_scale(
    c("L"),
    "tern_L", identity, name = name, breaks = breaks,
    minor_breaks = minor_breaks, labels = labels, limits = limits,
    expand = waiver(), oob = censor, na.value = NA_real_, trans = "identity",
    guide = "none"
  )
  # TODO: Fix this hack. We're reassigning the parent ggproto object, but this
  # object should in the first place be created with the correct parent.
  sc$super <- ScaleContinuousPosition
  class(sc) <- class(ScaleContinuousPosition)
  sc
}

#' @name scale_tern_continuous
#' @rdname scale_tern_continuous
#' @export
scale_R_continuous <- function(name         = waiver(), 
                               limits       = c(0,1), 
                               breaks       = getBreaks(limits,TRUE),
                               minor_breaks = getBreaks(limits,FALSE),
                               labels       = 100*breaks,
                               expand       = waiver(),...) {
  sc <- continuous_scale(
    c("R"),
    "tern_R", identity, name = name, breaks = breaks,
    minor_breaks = minor_breaks, labels = labels, limits = limits,
    expand = waiver(), oob = censor, na.value = NA_real_, trans = "identity",
    guide = "none"
  )
  # TODO: Fix this hack. We're reassigning the parent ggproto object, but this
  # object should in the first place be created with the correct parent.
  sc$super <- ScaleContinuousPosition
  class(sc) <- class(ScaleContinuousPosition)
  sc
}
