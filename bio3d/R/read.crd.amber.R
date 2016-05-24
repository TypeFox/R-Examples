"read.crd.amber" <- function(file, ...) {
  if (missing(file)) {
    stop("read.prmtop: please specify a crd 'file' for reading")
  }

  toread <- file.exists(file)
  if (!toread) {
    stop("No input crd file found: check filename")
  }

  trim <- function (s) {
    ## Remove leading and traling
    ## spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-NA
    s
  }

  parse.line <- function(line, fmt) {
    tmp <- seq(1, as.numeric(fmt[1])*as.numeric(fmt[2]), by=as.numeric(fmt[2]))
    substring(line, tmp, c(tmp[2:length(tmp)]-1, nchar(line)))
  }

  ## Read and parse file
  raw.lines <- readLines(file)

  name <- trim(raw.lines[1])
  ##num.atoms <- as.numeric(trim(raw.lines[2]))
  info <- unlist(strsplit(trim(raw.lines[2]), "  "))
  num.atoms <- as.numeric(info[1])

  simtime <- NULL
  if(length(info)>1)
    simtime=as.numeric(info[2])
  
  num.crdlines <- ceiling(num.atoms*3 / 6)


  if(length(raw.lines) > num.crdlines*2) {
    vel=TRUE
    boxline.ind <- (num.crdlines*2)+3
  }
  else {
    vel <- FALSE
    boxline.ind <- num.crdlines+3
  }

  ## parse coordinates
  fmt <- c(6, 12, 0, "a")
  tmplines <- raw.lines[3:(num.crdlines+2)]
  crds <- trim(unlist(lapply(tmplines, parse.line, fmt)))
  crds=as.numeric(crds)
  crds=crds[!is.na(crds)]
  
  ## parse velocities
  if(vel) {
    fmt <- c(6, 12, 0, "a")
    tmplines <- raw.lines[(num.crdlines+3):((num.crdlines*2)+2)]
    vels <- trim(unlist(lapply(tmplines, parse.line, fmt)))
    vels=as.numeric(vels)
    vels=vels[!is.na(vels)]
  }
  else {
    vels <- NULL
  }

  boxline <- raw.lines[boxline.ind]
  if(!is.na(boxline))
    box <- as.numeric(trim(parse.line(boxline, fmt)))
  else
    box <- NULL

  out <- list(xyz=as.xyz(crds), velocities=vels, time=simtime, natoms=num.atoms, box=box)
  class(out) <- c("amber", "crd")
  return(out)
}
