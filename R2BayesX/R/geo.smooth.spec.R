geo.smooth.spec <- function(object, dir, prg, data, type)
{
  if(!is.list(object$xt))
    object$xt <- list(object$xt)
  map.name <- help.map.name(deparse(substitute(object, env = .GlobalEnv), 
    backtick = TRUE, width.cutoff = 500L))
  if(!is.null(object$xt$map.name))
    map.name <- object$xt$map.name
  map.name <- rmf(map.name)
  map <- object$xt$map
  if(is.null(map)) {
    if(!is.null(object$xt$polys))
      map <- object$xt$polys
    if(!is.null(object$xt$penalty))
      map <- object$xt$penalty
  }
  if(is.null(map)) {
    if(!is.list(object$xt[[1L]]))
      map <- object$xt
    else {
      map <- NULL
      for(i in 1L:length(object$xt))
        if(inherits(object$xt[[i]], "bnd") || inherits(object$xt[[i]], "list"))
          map <- object$xt[[i]]
    }
    if(is.null(map)) {
      map <- object$xt
      if(inherits(map, "SpatialPolygons"))
        map <- sp2bnd(map)
      if(is.null(map) || (!is.list(map) && !inherits(map, "bnd")))
        stop("need to supply a bnd file object in argument xt!")
    }
  }
  if(inherits(map, "SpatialPolygons"))
    map <- sp2bnd(map)
  if(!inherits(map, "bnd"))
    class(map) <- "bnd"
  counter <- NULL
  ok <- TRUE
  if(length(map) < 2L && is.null(map[[1L]]))
    stop("map is missing!")
  if(!missing(dir)) {
    files <- list.files(dir)
    while(ok) {
      mapfile <- paste(map.name, counter, ".bnd", sep = "")
      if(any(grepl(mapfile, files))) {
        if(is.null(counter))
          counter <- 0L
        counter <- counter + 1L
      } else
        ok <- FALSE
    }
    mapfile <- file.path(dir, mapfile)
    prgfile <- file.path(dir, prg)
    cat("map", map.name, "\n", file = prgfile, append = TRUE)
    if(!any(is.na(poly.names <- as.integer(names(map))))) {
      poly.names <- sort(poly.names)
      poly.names <- as.character(poly.names)
    } else poly.names <- sort(names(map))
    map <- map[poly.names]
    class(map) <- "bnd"
    write.bnd(map = map, file = mapfile, replace = TRUE)
    cmd <- paste(map.name, ".infile using ", mapfile, "\n", sep = "")
    cat(cmd, file = prgfile, append = TRUE)
  }
  term <- if(length(object$term) > 2) paste(object$term[2L:1L], collapse = "*") else object$term
  if(length(object$p.order) == 1L) 
    m <- rep(object$p.order, 2L)
  else 
    m <- object$p.order
  m[is.na(m)] <- 2L
  object$p.order <- m
  object$p.order[1L] <- object$p.order[1L] + 1L
#  if(object$p.order[2L] > 1L) {
#    object$p.order[2L] <- 1L
#    if(type == "geospline")
#      warning("only random walks of order 1 supported for geosplines, set to default!")
#  }
  if(object$bs.dim < 0L) {
    if(type == "geokriging")
      object$bs.dim <- as.integer(length(map) / 2)
    else
      object$bs.dim <- 7L
  } else {
    if(object$bs.dim >= length(map))
      stop("basis dimension is larger than existing polygons in bnd object!")
  }
  if(type != "geokriging")
    nrknots <- object$bs.dim - object$p.order[1L] + 1L
  else
    nrknots <- object$bs.dim
  if(type == "geokriging") {
    if(!is.null(object$xt$full)) {
      term <- paste(term, "(", type, ",map=", map.name, ",full", sep = "")
      object$xt$full <- NULL
    } else term <- paste(term, "(", type, ",map=", map.name, ",nrknots=", nrknots, sep = "")
  } else {
    term <- paste(term, "(", type, ",map=", map.name,
      ",nrknots=", nrknots, ",degree=", object$p.order[1L], sep = "")
  }
  term <- paste(do.xt(term, object, c("map", "polys", "penalty", "map.name")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

