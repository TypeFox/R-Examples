

# Hauptfunktion, die dem Anwender praesentiert wird:
# dots = arguments to be passed to plot.routes()

movie <- function(x, ...) {
    UseMethod("movie")
}

movie.routes <- function(x, show.progress = TRUE, save = TRUE,
                         width = 1024, height = 768, ...){

  from <- min(x[x$lsk == "S",]$time)
  to <- max(x[x$lsk == "L",]$time)

  stopifnot(from < to)

  dat.tm <- seq(from = from, to = to, by = 60)

  n <- length(dat.tm)
  # create progress bar
  if(show.progress)
      pb <- txtProgressBar(title = "progress bar", min = 0, max = n, width = 300)

  shapeEU <- NULL
  spCaps <- NULL
  load(system.file("maps/shapeEU.Rda", package = "MUCflights"))
  load(system.file("maps/spCaps.Rda", package = "MUCflights"))
  stopifnot(!is.null(shapeEU))
  stopifnot(!is.null(spCaps))

  dir <- NULL
  if ( save )
      dir <- tempdir()

  for ( i in 1:n ) {
    if ( save )
        jpeg(file.path(dir, sprintf("frame%s.jpg", formatC(i, width = 4, flag = "0"))),
             width = width, height = height) # FIXME: 4

    plot(x, time = dat.tm[i], borders = shapeEU, capitals = spCaps, ...)

    if ( show.progress )
        setTxtProgressBar(pb, i)

    if ( save )
        dev.off()
  }

  if ( show.progress )
      close(pb)

  dir
}

ffmpeg <- function(dir, ffmpeg = "ffmpeg") {
    wd <- setwd(dir)
    cmd <- paste(ffmpeg, 
        "-f image2 -i frame%04d.jpg -r 12 -target dvd -s xga movie.mpg")
    system(cmd)
    setwd(wd)
    file.path(dir, "movie.mpg")
}

