bayesx.construct.mrf.smooth.spec <- bayesx.construct.spatial.smooth.spec <- function(object, dir, prg, data)
{
  if(missing(prg))
    prg <- " "
  if(missing(dir))
    dirok <- FALSE
  else
    dirok <- TRUE
  if(is.null(object$xt))
    stop("need to supply a map object in argument xt!")  
  map.name <- help.map.name(deparse(substitute(object, env = .GlobalEnv), 
    backtick = TRUE, width.cutoff = 500L))
  if(!is.null(object$xt$map.name))
    map.name <- object$xt$map.name
  if(!is.list(object$xt))
    object$xt <- list(object$xt)
  map.name <- rmf(gsub("\\s", "", paste(map.name, sep = "", collapse = "")))

  map <- object$xt$map
  if(is.null(map)) {
    if(!is.null(object$xt$polys))
      map <- object$xt$polys
    if(!is.null(object$xt$penalty))
      map <- object$xt$penalty
  }
  if(is.null(map))
    map <- object$xt$gra
  if(is.null(map)) {
    if(!is.list(object$xt[[1L]])) {
      if(inherits(object$xt[[1L]], "gra"))
        map <- object$xt[[1L]]
      else
        map <- object$xt
    } else map <- object$xt[[1L]]
    if(is.null(map)) {
      map <- object$xt
      if(inherits(map, "SpatialPolygons"))
        map <- sp2bnd(map)
      if(is.null(map) || (!is.list(map) && !inherits(map, "bnd") || !inherits(map, "gra")))
        stop("need to supply a bnd or graph file object in argument xt!")
    }
  }
  if(is(map, "nb"))
    map <- nb2gra(map)
  if(inherits(map, "SpatialPolygons"))
    map <- sp2bnd(map)
  if(!inherits(map, "bnd") && !inherits(map, "gra")) {
    if(is.list(map))
      class(map) <- "bnd"
    else
      class(map) <- "gra"
  }
  if(dirok) {
    counter <- NULL
    ok <- TRUE
    files <- list.files(dir)
    while(ok) {
      classm <- class(map)
      if(length(classm) > 1L)
        if("list" %in% classm)
          class(map) <- classm[classm != "list"]
      mapfile <- paste(map.name, counter, ".", class(map), sep = "")[1]
      if(any(grepl(mapfile, files))) {
        if(is.null(counter))
          counter <- 0L
        counter <- counter + 1L
      } else ok <- FALSE
    }
    mapfile <- file.path(dir, mapfile)
    prgfile <- file.path(dir, prg)
    prgok <- file.exists(prgfile)
  } else prgok <- FALSE
  if(prgok)
    cat("map", map.name, "\n", file = prgfile, append = TRUE)
  if(dirok) {
    if(inherits(map, "bnd")) {
      if(!any(is.na(poly.names <- as.integer(names(map))))) {
        poly.names <- sort(poly.names)
        poly.names <- as.character(poly.names)
      } else poly.names <- sort(names(map))
      map <- map[poly.names]
      class(map) <- "bnd"
      write.bnd(map = map, file = mapfile, replace = TRUE)
      cmd <- paste(map.name, ".infile using ", mapfile, "\n", sep = "")
    } else {
      if(!is.character(map)) {
        dx <- as.character(unique(data[[object$term]]))
        cnm <- colnames(map)
        if(!all(dx %in% cnm))
          stop(paste("not all regions specified in variable", object$term, "in adjacency matrix!"))
        write.gra(map = map, file = mapfile, replace = TRUE)
        cmd <- paste(map.name, ".infile, graph using ", mapfile, "\n", sep = "")
      } else {
        stopifnot(is.character(map))
        pos <- regexpr("\\.([[:alnum:]]+)$", map)
        fext <- ifelse(pos > -1L, substring(map, pos + 1L), "")
        if(fext == "gra")
          cmd <- paste(map.name, ".infile, graph using ", path.expand(map), "\n", sep = "")
        else
          cmd <- paste(map.name, ".infile using ", path.expand(map), "\n", sep = "")
      }
    }
  }
  if(prgok)
    cat(cmd, file = prgfile, append = TRUE)
  term <- object$term
  term <- paste(term, "(spatial,map=", map.name, sep = "")
  term <- paste(do.xt(term, object, c("map", "polys", "penalty", "map.name")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

