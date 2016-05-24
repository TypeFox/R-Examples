view.cna <- function(x, pdb, layout=layout.cna(x, pdb, k=3),
                     col.sphere=NULL, 
                     col.lines="silver",
                     weights=NULL,
                     radius=table(x$communities$membership)/5,
                     alpha=1,
                     vmdfile="network.vmd", pdbfile="network.pdb",
                     launch=FALSE) {

  ## Draw a cna network in VMD

  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }
    
  if(is.null(weights)){
    weights <- igraph::E(x$community.network)$weight
    
    if(is.null(x$call$minus.log)){
      weights <- exp(-weights)
    }
    else{
      if(x$call$minus.log){
        weights <- exp(-weights)
      }
    }
  }
  
  if(is.null(col.sphere)) {
    ## Get colors from network and convert to 0:17 VMD color index
    col.sphere <- match(igraph::V(x$community.network)$color, vmd.colors())-1
  } else {
    ## Check supplied color(s) will work in VMD
    if(!all(col.sphere %in% c(0:17))) {
      warning("Input 'col.sphere' may not work properly in VMD
               - should be 0:17 color index value")
    }
  }


  ##-- VMD draw functions for sphere, lines and cone
  .vmd.sphere <- function(cent, radius=5, col="red", resolution=25) {
    ## .vmd.sphere( matrix(c(0,0,0, 1,1,1), ncol=3,byrow=T) )
    
    if(ncol(cent) != 3)
      stop("Input 'cent' should be a 3 col xyz matrix")

    n <- nrow(cent); scr <- rep(NA, n) 
    if(length(col) != n)
      
      col <- rep(col, length=n)
    
    if(length(radius) != n)
      radius <- rep(radius, length=n)
    
    for(i in 1:n) {
      scr[i] <- paste0("draw color ", col[i],
                       "\ndraw sphere {",
                       paste(cent[i,], collapse = " "),
                       "} radius ", radius[i],
                       " resolution ",resolution, "\n")
    }
    return(scr)
  }
  
  
  .vmd.lines <- function(start, end, radius=0.2, col="silver", resolution=25) {
    
    ## .vmd.lines( start=matrix(c(0,0,0), ncol=3,byrow=T),
    ##         end=matrix(c(1,1,1), ncol=3,byrow=T) )
    
    if(ncol(start) != 3)
      stop("Input 'start' and 'end' should be 3 col xyz matrices")
    n <- nrow(start); scr <- rep(NA, n)
    
    if(length(col) != n)
      col <- rep(col, length=n)
    
    if(length(radius) != n)
      radius <- rep(radius, length=n)
    
    for(i in 1:n) {
      scr[i] <- paste0("draw color ", col[i],
                       "\ndraw cylinder {",
                       paste(start[i,], collapse = " "),
                       "} {", paste(end[i,], collapse = " "),
                       "} radius ", radius[i],
                       " resolution ",resolution, "\n")
    }
    return(scr)
  }
  
  .vmd.cone <- function(start, end, radius=5, col="silver", resolution=25) {
    warning("not here yet")
  }

  ##- Set alpha if needed
  scr <- NULL
  if(alpha != 1)
    scr <- paste("material change opacity Transparent",
                 alpha,"\ndraw material Transparent\n")
  
  ##- Lets get drawing
  ##radius = V(x$community.network)$size
  ###radius = table(x$raw.communities$membership)/5
  scr <- c(scr, .vmd.sphere( layout, radius=radius, col=col.sphere))

  ## Edges
###edge.list <- unlist2(get.adjlist(x$community.network))
###start.no <- as.numeric(names(edge.list))
###end.no <- as.numeric((edge.list))
###inds <- which(end.no > start.no)
###start <- layout[start.no[inds],]
###end <- layout[end.no[inds],]
  edge.list <- igraph::get.edges(x$community.network, 1:length(igraph::E(x$community.network)))
  start <- layout[edge.list[,1],]
  end <- layout[edge.list[,2],]
  
  ###weights=E(x$community.network)$weight ##/0.2
  scr <- c(scr, .vmd.lines( start=start, end=end,
                           radius=weights, col=col.lines))

  cat(scr, file=vmdfile, sep="")

  ## Output a PDB file with chain color
  # Use the chain field to store cluster membership data for color in VMD
  ch <- vec2resno(vec=x$communities$membership, resno=pdb$atom[,"resno"])
  write.pdb(pdb, chain=LETTERS[ch], file=pdbfile)

  ## Launch option ...
  ## vmd -pdb network.pdb -e network.vmd
  if(launch) {
    cmd <- paste("vmd", pdbfile, "-e", vmdfile)

    os1 <- .Platform$OS.type
    if (os1 == "windows") {
      shell(shQuote(cmd))
    } else{
      if(Sys.info()["sysname"]=="Darwin") {
        system(paste("/Applications/VMD\\ 1.9.*app/Contents/MacOS/startup.command",pdbfile, "-e", vmdfile))
      }
      else{
        system(cmd)
      }
    }
  }
}
