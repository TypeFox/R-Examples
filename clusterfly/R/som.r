#' Visualise Kohonen self organising maps with GGobi
#' Displays both data, and map in original high-d space.
#'
#' Map variables added as map1 and map2.  Plot these to
#' get traditional SOM plot.  Tour over all other variables to
#' see how well the map fits the original data.
#'
#' @param data SOM object
#' @param ... ignored
#' @method ggobi som
#' @keywords cluster dynamic
#' @export
#' @examples
#' \dontrun{
#' d.music <- read.csv("http://www.ggobi.org/book/data/music-all.csv")
#'
#' music <- rescaler(d.music)[complete.cases(d.music), 1:10]
#' music.som <- som::som(music[,-(1:3)], 6, 6, neigh="bubble", rlen=1000)
#' ggobi(music.som)
#' }
#' \dontrun{
#' d.music <- read.csv("http://www.ggobi.org/book/data/music-all.csv")
#'
#' music <- rescaler(d.music)[complete.cases(d.music), 1:10]
#' music.hex <- kohonen::som(music[,-(1:3)], grid = somgrid(3, 3, "hexagonal"), rlen=1000)
#' music.rect <- kohonen::som(music[,-(1:3)], grid = somgrid(6, 6, "rectangular"), rlen=1000)
#' ggobi(music.rect)
#' }
ggobi.som <- function(data, ...) {
  som <- data
  original <- data.frame(
    som$data,
    map1 = jitter(som$visual$x) + 1,
    map2 = jitter(som$visual$y) + 1,
    net = factor(FALSE)
  )

  xs <- som$xdim
  ys <- som$ydim

  net <- som$code
  colnames(net) <- colnames(som$data)
  net <- cbind(net, expand.grid(map1=1:xs, map2=1:ys), net=factor(TRUE))
  rownames(net) <- paste("net", 1:nrow(net), sep="")
  names(net) <- names(original)

  df <- rbind(original, net)

  g <- ggobi(df)
  glyph_colour(g[1]) <- c(1,3)[df$net]
  shadowed(g[1]) <- c(FALSE,TRUE)[df$net]
  d <- displays(g)[[1]]
  variables(d) <- list(X = "map1", Y = "map2")

  # Add net edges
  netlines <- make_rect_net(xs, ys)
  edges(g) <- netlines
  glyph_colour(g[2]) <- 3
  edges(d) <- g[2]

  invisible(g)
}

#' @export
#' @importFrom plyr rbind.fill
ggobi.kohonen <- function(data, extra = NULL, ...) {

  som <- data

  original <- data.frame(
    som$data,
    map1 = jitter(som$grid$pts[som$unit.classif, 1]),
    map2 = jitter(som$grid$pts[som$unit.classif, 2]),
    distance = som$distance,
    net = factor(FALSE),
    oid = factor(1:nrow(som$data))
  )
  if (!is.null(extra)) original <- cbind(original, extra)

  net <- data.frame(
    map1 = som$grid$pts[, 1],
    map2 = som$grid$pts[, 2],
    som$codes,
    net = factor(TRUE),
    oid = factor(paste("net"), 1:nrow(som$grid$pts))
  )
  rownames(net) <- paste("net", 1:nrow(net), sep="")

  df <- rbind.fill(original, net)

  g <- ggobi(df)
  glyph_colour(g[1]) <- c(1,3)[df$net]
  shadowed(g[1]) <- c(FALSE,TRUE)[df$net]
  d <- displays(g)[[1]]
  variables(d) <- list(X = "map1", Y = "map2")

  if (som$grid$topo == "rectangular") {
    netlines <- make_rect_net(som$grid$xdim, som$grid$ydim)
    edges(g) <- netlines
    glyph_colour(g[2]) <- 3
    edges(d) <- g[2]
  }

  invisible(g)
}

make_rect_net <- function(xs, ys) {
  netlines <- with(expand.grid(y=1:(xs-1), x=1:(ys)), rbind(
    cbind((x - 1) * xs + y, (x - 1)     * xs + y + 1),
    cbind((x - 1) * xs + y,   x         * xs + y)
  ))
  netlines <- rbind(netlines, cbind(1:(ys-1) * xs, 2:ys * xs))
  netlines <- apply(netlines, 2, function(x) paste("net", x, sep=""))
  netlines
}

