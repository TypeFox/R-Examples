visualize <- function(...)
  UseMethod("visualize")

visualize.coords <- function( x, elename = NULL, cryst1 = NULL, conect = NULL, mode = NULL,
                              type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                              lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                              cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                              add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  if(!is.coords(x)) stop("'x' must be an object of class coords.")
  
  if(basis(x) == "abc") x <- abc2xyz(x, cryst1)

# Unrecognized elements are considered as dummy atoms
  if(is.null(elename)){
    warning("'elename' was not specifyed. All atom have been considered as dummy atoms.")
    elename <- rep("Xx", natom(x))
  }
  
  symb <- toSymbols(elename)
  symb[is.na(symb)] <- "Xx"
  M <- match(symb,elements[,"symb"])
  
  if(is.null(col)){
    col <- elements[M,c("red","green","blue")]
    col <- do.call(rgb, col)
  }
  if(length(col)!=natom(x)) col <- rep(col, length = natom(x))
  
  if(!add){
    open3d()
    par3d(windowRect = windowRect, userMatrix=userMatrix, FOV = FOV, ...)
    bg3d(color=bg)
  }
  ids <- rgl.ids()
  par.save <- par3d(skipRedraw=TRUE)
  
  if(is.null(xyz) & is.null(cryst1))
    xyz <- TRUE
  else
    xyz <- FALSE
  
  if(is.null(abc) & is.null(cryst1))
    abc <- FALSE
  else
    abc <- TRUE
  
  if(is.null(pbc.box) & is.null(cryst1))
    pbc.box <- FALSE
  else
    pbc.box <- TRUE
  
  if(abc & is.null(cryst1)) {
    warning("Cannot find periodical boundary conditions")
    abc <- FALSE
  }
  if(pbc.box & is.null(cryst1)) {
    warning("Cannot find periodical boundary conditions")
    pbc.box <- FALSE
  }
  
  if(xyz) ids <- rbind(ids, addXYZ(lwd = lwd.xyz, cex = cex.xyz))
  if(abc) ids <- rbind(ids, addABC(cryst1, lwd = lwd.abc, cex = cex.abc))
  if(pbc.box) ids <- rbind(ids, addPBCBox(cryst1, lwd = lwd.pbc.box))
  
  if(type == "l")
  {
    if(is.null(conect)){
      warning("Undefined connectivity: Computing connectivity from coordinates...")
      conect <- conect(x)
    }
    ind <- t(conect)
    seg.id <- segments3d(
      x$x1[ind],
      x$x2[ind],
      x$x3[ind],
      color = col[ind], lwd=lwd, ...
    )
    seg.id <- data.frame(id = seg.id, type = "atom.seg")
    ids <- rbind(ids, seg.id)
  }
  
  if(type == "p"){
    pts.id <- points3d(
      x$x1,
      x$x2,
      x$x3,
      color = col, ...)
    pts.id <- data.frame(id = pts.id, type = "atom.pts")
    ids <- rbind(ids, pts.id)
  }
  
  if(type == "s"){
    if(is.character(radii[1])){
      if(! radii[1] %in% c("rcov", "rbo", "rvdw") )
        stop("'radii' must be one of 'rcov', 'rbo', 'rvdw' or a numerical vector")
      radii <- elements[M,radii[1]]
    }
    if(all(radii==0)){
      warning("All atoms are dummy atoms. 'radii' have been set to 1")
      radii <- rep(1, natom(x))
    }
    sph.id <- spheres3d(
      x$x1,
      x$x2,
      x$x3,
      color = col, radius = radii, ...)
    sph.id <- data.frame(id = sph.id, type= "atom.sph")
    ids <- rbind(ids, sph.id)
  }

  par3d(par.save)
  if(!is.null(mode)){
    if(mode == "measure"){
      measure(x)
    }
    else if(mode == "info"){
      stop("No 'info' mode for object of class 'coords'")
    }
    else{
      stop("Unrecognized visualization mode")
    }
  }
  
  invisible(ids)
  
}

visualize.data.frame <- function( x, elename = NULL, cryst1 = NULL, conect = NULL, mode = NULL,
                                  type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                                  lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                                  cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                                 add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  if(is.null(basis(x))){
    basis(x) <- "xyz"
    warning("No basis attribute were found. Coordinates are assumed to Cartesian.")
  }
  
  visualize(coords(x), elename, cryst1, conect, mode, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg,  radii, add, windowRect, FOV, userMatrix, ...)
}

visualize.matrix <- function( x, elename = NULL, cryst1 = NULL, conect = NULL, mode = NULL,
                              type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                              lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                              cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                              add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  if(is.null(basis(x))){
    basis(x) <- "xyz"
    warning("No basis attribute were found. Coordinates are assumed to Cartesian.")
  }
  
  visualize(coords(x), elename, cryst1, conect, mode, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg,  radii, add, windowRect, FOV, userMatrix, ...)
}

visualize.atoms <- function( x, cryst1 = NULL, conect = NULL, mode = NULL,
                             type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                             lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                             cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                             add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  ids <- visualize(coords(x), x$elename, cryst1, conect, mode=NULL, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg,  radii, add, windowRect, FOV, userMatrix, ...)
  
  if(!is.null(mode)){
    if(mode == "measure"){
      measure(x)
    }
    else if(mode == "info"){
      info3d(x)
    }
    else{
      stop("Unrecognized visualization mode")
    }
  }
  
  invisible(ids)
}

visualize.pdb <- function(x, mode = NULL, type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                          lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                          cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                          add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  ids <- visualize(x$atoms, cryst1 = x$cryst1, conect = x$conect, mode=NULL, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg,  radii, add, windowRect, FOV, userMatrix, ...)
  
  if(!is.null(mode)){
    if(mode == "measure"){
      measure(x)
    }
    else if(mode == "info"){
      info3d(x)
    }
    else{
      stop("Unrecognized visualization mode")
    }
  }
  
  invisible(ids)
}

visualize.character <- function(x, mode = NULL, type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                                lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                                cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                                add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  x <- read.pdb(x)
  visualize.pdb(x, mode = mode, type = type, xyz = xyz, abc = abc, pbc.box = pbc.box, lwd = lwd,
           lwd.xyz = lwd.xyz, lwd.abc = lwd.abc, lwd.pbc.box = lwd.pbc.box,
           cex.xyz = cex.xyz, cex.abc = cex.abc, col = col, bg = bg,  radii = radii,
           add = add, windowRect = windowRect, FOV = FOV, userMatrix = userMatrix, ...)
  
}
