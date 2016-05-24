addResLab <- function(x, ...)
  UseMethod("addResLab")

addResLab.atoms <- function(x, at.centre = TRUE, col = "black", ...){
  if(at.centre){
    labels <- paste0(x$resname, x$resid)
    labels <- labels[!duplicated(labels)]
    xyz <- centres(x)
  }
  else{
    labels <- paste0(x$resname, x$resid)
    xyz <- coords(x)
  }
  txt.ids <- rgl.texts(xyz, text=labels, col = col, ...)
  invisible(txt.ids)
}

addResLab.pdb <- function(x, at.centre = TRUE, col = "black", ...)
  addResLab.atoms(x$atoms, at.centre = at.centre, col = col, ...)

addEleLab <- function(x, ...)
  UseMethod("addEleLab")

addEleLab.atoms <- function(x, eleid = FALSE, col = "black", ...){
  if(eleid){
    labels <- paste(x$elename, x$eleid, sep=":")
  }
  else{
    labels <- x$elename    
  }
  xyz <- coords(x)
  txt.ids <- rgl.texts(xyz, text=labels, col = col, ...)
  invisible(txt.ids)
}

addEleLab.pdb <- function(x, eleid = FALSE, col = "black", ...)
  addEleLab.atoms(x$atoms, eleid = eleid, col = col, ...)

info3d <- function(...)
  UseMethod("info3d")

info3d.atoms <- function(x, id = rgl.ids(), col = "black", verbose = TRUE, adj = 0, ...){
  cat("Presse ESC when you have finish your selections.\n")
  ids <- id[id$type!="text",]$id
  if(basis(x) == "abc")
    stop("'x' must contain Cartesian coordinates")
  xyz <- coords(x, basis = "xyz")
  
  sph.ids <- NULL
  txt.ids <- NULL
  clean <- function(){
    if(!is.null(sph.ids))
      rgl.pop(id = sph.ids)
    if(!is.null(txt.ids))
      rgl.pop(id = txt.ids)
  }
  on.exit(clean())
  
  dist <- 0.0015
  sel <- NULL
  while(TRUE){
    f <- rgl.select3d(button = "right", ...)
    if(is.null(f)) 
      break
    e <- environment(f)
    rgl.info <- lapply(ids,
                       function(id){
                         verts <- rgl.attrib(id,"vertices")
                         info <- verts
                         return(info)
                       })
    info <- do.call(rbind, rgl.info)
    verts <- info[,1:3]
    hits <- f(verts)
    if(!any(hits)){
      wincoords <- rgl.user2window(verts, projection = e$proj)
      hits <- (0 <= wincoords[,3]) && (wincoords[,3] <= 1)
      if(any(hits)){
        dists <- (wincoords[, 1] - e$llx)^2 +
          (wincoords[, 2] - e$lly)^2
        hits <- (dists < dist) & (dists == min(dists))
      }
    }
    if(any(hits)){
      verts <- unique(verts[hits, , drop = FALSE])
      sel <- rbind(sel, verts)
      res.lab <- do.call(paste0, x[,c("resid","resname")])
      ele.lab <- do.call(paste0, x[,c("eleid","elename")])
      labels <- paste(res.lab, ele.lab, sep=":")
      M <-     xyz$x1%in%round(verts[,1],digits=3)
      M <- M & xyz$x2%in%round(verts[,2],digits=3)
      M <- M & xyz$x3%in%round(verts[,3],digits=3)
      if(verbose){
        if(length(which(M))!=0)
          print(subset(x, M))
      }
      if(length(which(M))!=0){
        txt.ids <- c(txt.ids, 
                     rgl.texts(verts, text=labels[M], col = col,, adj = adj, ...))
      }
    }
    else{
      clean()
      sph.ids <- NULL
      txt.ids <- NULL
      sel <- NULL
    }
  }
}

info3d.pdb <- function(x, id = rgl.ids(), col = "black", verbose = TRUE, adj = 0, ...)
  info3d.atoms(x$atoms, id = id, col = col, verbose = verbose, adj = adj, ...)
