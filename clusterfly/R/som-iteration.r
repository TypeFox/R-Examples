som_iterate <- function(df, grid, nsteps = 100, stepsize = 10, alpha = 0.05, radius = NULL) {
  if (is.null(radius)) {
    radius <- c(quantile(unit.distances(grid, FALSE), 0.67), 0)
  }

  # Alpha decreases linearly
  alpha_step <- function(i) alpha - (alpha - 0.01) * (i + c(-1, 0)) / nsteps
  # Radius decraeses linearly to 1 by 1/3 of steps
  radius_step <- function(i) {
    r <- radius[1] - (radius[1] - radius[2]) * 3.0 * (i - 1) / nsteps;
    ifelse(r < radius[2], max(0.5, radius[2]), r)
  }

  print(radius_step(1:nsteps))

  i <- 1
  fit <- kohonen::som(df, grid, rlen = stepsize, alpha = alpha_step(i), radius = radius_step(i))

  step <- function() {
    i <<- i + 1
    fit <<- kohonen::som(df, grid, rlen = stepsize, init = fit$codes,
      alpha = alpha_step(i), radius=radius_step(i), keep.data = TRUE
    )

    fit
  }

  structure(c(list(fit), replicate(nsteps - 1, step(), simplify=FALSE)), class="somiter")
}

#' @export
summary.somiter <- function(object, ...) {
  interesting <- function(fit) {
    df <- data.frame(
      alpha_start = fit$alpha[1], alpha_end = fit$alpha[2],
      radius = fit$radius[1],
      change_mean = mean(fit$changes),
      dist_mean = mean(fit$distances),
      dist_sd = sd(fit$distances)
    )
    df$codes <- list(fit$codes)
    df$map <- list(fit$grid$pts[fit$unit.classif, ])
    df
  }
  df <- do.call("rbind", lapply(object, interesting))
  df$step <- 1:nrow(df)
  class(df) <- c("somitersum", class(df))
  rownames(df) <- paste("step", 1:nrow(df), sep="")
  df
}

#' @importFrom reshape2 melt
#' @importFrom RGtk2 gSignalConnect ==.RGtkObject
#' @export
ggobi.somiter <- function(data, extra = NULL, ...) {

  g <- ggobi(data[[1]], extra=extra)
  all_fits <- summary(data)
  fits <- fits[setdiff(names(fits), c("codes", "map"))]

  jittering <- jitter(all_fits[[1, "map"]]) - all_fits[[1, "map"]]

  distances <- melt(sapply(data, function(x) x$distances))
  names(distances) <- c("oid", "step", "value")
  distances$oid <- factor(distances$oid)
  oid <- NULL # stupid hack for R CMD check
  ggobi_longitudinal(distances, step, oid, g = g)

  ggobi_longitudinal(fits, step, g = g)
  d <- display(g["fits"], vars = list(X = "step", Y = "dist_mean"))
  edges(d) <- g["fits-edges"]

  gSignalConnect(g, "identify-point", function(gg, plot, id, data) {
    if (id == -1 || !"==.RGtkObject"(data, gg$fits)) return()
    id <- id + 1

    codes <- all_fits[[id, "codes"]]
    gg$df[gg$df$net == TRUE, 1:ncol(codes)] <- codes

    map <- all_fits[[id, "map"]] + jittering
    gg$df[gg$df$net == FALSE, c("map1", "map2")] <- map

    # dist <- all_fits[[id, "dist"]]
    # gg$df[gg$df$net == FALSE, c("distance")] <- as.matrix(dist)
  })

  invisible(g)
}
