measure <- function(...)
  UseMethod("measure")

measure.default <- function(id = rgl.ids(), verbose = TRUE, ...){
  cat("Presse ESC when you have finish your selections.\n")
  ids <- id[id$type!="text",]$id

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
        radii <- rgl.attrib(id,"radii")
        if(nrow(radii)==0)
          radii <- rep(0, nrow(verts))
        radii <- radii + 0.2
        info <- cbind(verts, radii)
        return(info)
      })
    info <- do.call(rbind, rgl.info)
    verts <- info[,1:3]
    radii <- info[,4]
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
      radii <- unique(radii[hits,   drop = FALSE])
      sel <- rbind(sel, verts)
      if(nrow(sel) > 4){
        clean()
        sph.ids <- NULL
        txt.ids <- NULL
        sel <- NULL
      }
      else{
        if(verbose)
          print(verts)
        
        sph.ids <- c(sph.ids,
                     spheres3d(verts, color = "red", alpha = 0.3, radius = radii))
        if(nrow(sel)==2){
          sel.coords <- coords(sel, basis = "xyz")
          B <- bond(sel.coords, 1, 2)
          if(verbose) print(B)
          B <- round(B, digits = 4)
          txt.ids <- rgl.texts(centres(sel.coords), text = B, col = "red", cex = 1.5, adj = 0)
        }
        if(nrow(sel)==3){
          sel.coords <- coords(sel, basis = "xyz")
          A <- angle(coords(sel.coords, basis = "xyz"), 1, 2, 3)
          if(verbose) print(A)
          A <- round(A, digits = 2)
          rgl.pop(id = txt.ids)
          txt.ids <- rgl.texts(centres(sel.coords), text = A, col = "red", cex = 1.5, adj = 0)
        }
        if(nrow(sel)==4){
          sel.coords <- coords(sel, basis = "xyz")
          D <- dihedral(coords(sel.coords, basis = "xyz"), 1, 2, 3, 4)
          if(verbose) print(D)
          D <- round(D, digits = 2)
          rgl.pop(id = txt.ids)
          txt.ids <- rgl.texts(centres(sel.coords), text = D, col = "red", cex = 1.5, adj = 0)
        }
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

measure.coords <- function(x, id = rgl.ids(), verbose = TRUE, ...){
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
                         radii <- rgl.attrib(id,"radii")
                         if(nrow(radii)==0)
                           radii <- rep(0, nrow(verts))
                         radii <- radii + 0.2
                         info <- cbind(verts, radii)
                         return(info)
                       })
    info <- do.call(rbind, rgl.info)
    verts <- info[,1:3]
    radii <- info[,4]
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
      radii <- unique(radii[hits,   drop = FALSE])
      sel <- rbind(sel, verts)
      if(nrow(sel) > 4){
        clean()
        sph.ids <- NULL
        txt.ids <- NULL
        sel <- NULL
      }
      else{
        if(verbose){
          M <-     xyz$x1%in%round(verts[,1],digits=3)
          M <- M & xyz$x2%in%round(verts[,2],digits=3)
          M <- M & xyz$x3%in%round(verts[,3],digits=3)
          if(length(which(M))!=0)
            print(subset(x, M))
        }
        
        sph.ids <- c(sph.ids,
                     spheres3d(verts, color = "red", alpha = 0.3, radius = radii))
        if(nrow(sel)==2){
          sel.coords <- coords(sel, basis = "xyz")
          B <- bond(sel.coords, 1, 2)
          if(verbose) print(B)
          B <- round(B, digits = 4)
          txt.ids <- rgl.texts(centres(sel.coords), text = B, col = "red", cex = 1.5, adj = 0)
        }
        if(nrow(sel)==3){
          sel.coords <- coords(sel, basis = "xyz")
          A <- angle(coords(sel.coords, basis = "xyz"), 1, 2, 3)
          if(verbose) print(A)
          A <- round(A, digits = 2)
          rgl.pop(id = txt.ids)
          txt.ids <- rgl.texts(centres(sel.coords), text = A, col = "red", cex = 1.5, adj = 0)
        }
        if(nrow(sel)==4){
          sel.coords <- coords(sel, basis = "xyz")
          D <- dihedral(coords(sel.coords, basis = "xyz"), 1, 2, 3, 4)
          if(verbose) print(D)
          D <- round(D, digits = 2)
          rgl.pop(id = txt.ids)
          txt.ids <- rgl.texts(centres(sel.coords), text = D, col = "red", cex = 1.5, adj = 0)
        }
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

measure.pdb <- function(x, id = rgl.ids(), verbose = TRUE, ...)
  measure.coords(x$atoms, id, verbose, ...)
