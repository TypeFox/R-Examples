rivernet.read <- function(file.reachcoord,
                          file.reachattrib = NA,
                          file.nodeattrib  = NA,
                          colnames         = c(reach = "Reach_ID",
                                               node  = "Node_ID",
                                               x     = "X",
                                               y     = "Y",
                                               z     = "Z"),
                          sep              ="\t",
                          tol              = 1,
                          analyze          = FALSE,
                          verbose          = TRUE,
                          ...)
{
  # read river reaches:
  # ===================
  
  # read structure file(s) (containing columns reach, x, y, (z): 
  # ------------------------------------------------------------
  
  coord <- read.table(file.reachcoord[1],header=TRUE,sep=sep,stringsAsFactors=FALSE,...)
  if ( length(file.reachcoord) > 1 )
  {
    for ( i in 2:length(file.reachcoord) )
    {
      coord <- rbind(coord,read.table(file.reachcoord[i],header=TRUE,sep=sep,stringsAsFactors=FALSE,...))   
    }
  }
  if ( nrow(coord) < 1 )
  {
    cat("no coordinate data found","\n")
    return(NA)
  }
  
  # check presence of column names:
  # -------------------------------
  
  colnames.missing <- character(0)
  if ( is.na(colnames["reach"]) ) colnames.missing <- c(colnames.missing,"reach")
  if ( is.na(colnames["node"]) )  colnames.missing <- c(colnames.missing,"node")
  if ( is.na(colnames["x"]) )     colnames.missing <- c(colnames.missing,"x")
  if ( is.na(colnames["y"]) )     colnames.missing <- c(colnames.missing,"y") 
  if ( is.na(colnames["z"]) )     colnames.missing <- c(colnames.missing,"z") 
  if ( length(colnames.missing) > 0 )
  {
    cat("argument \"colnames\" must contain element(s) \"",paste(colnames.missing,collapse="\", \""),"\"\n",sep="")
    return(NA)
  }
  
  # check presence of columns reach, x and y:
  # -----------------------------------------
  
  col.reach <- match(colnames["reach"],colnames(coord))
  col.x     <- match(colnames["x"],colnames(coord))
  col.y     <- match(colnames["y"],colnames(coord))
  col.z     <- match(colnames["z"],colnames(coord))
  columns.missing <- character(0)
  if ( is.na(col.reach) ) columns.missing <- c(columns.missing,colnames["reach"])
  if ( is.na(col.x) ) columns.missing <- c(columns.missing,colnames["x"])
  if ( is.na(col.y) ) columns.missing <- c(columns.missing,colnames["y"])
  if ( length(columns.missing) > 0 )
  {
    cat("file \"",file.reachcoord,"\" must contain column(s) \"",paste(columns.missing,collapse="\", \""),"\"\n",sep="")
    return(NA)
  }
  
  # ensure that reach coordinates are numeric (factors must first be converted to char)
  # and identifiers are characters and extract unique identifiers:
  # -----------------------------------------------------------------------------------
  
  if ( !is.numeric(coord[,col.x]) )
  {
    coord[,col.x] <- as.numeric(as.character(coord[,col.x]))
  }
  if ( !is.numeric(coord[,col.y]) )
  {
    coord[,col.y] <- as.numeric(as.character(coord[,col.y]))
  }
  if ( !is.na(col.z) ) 
  {
    if ( !is.numeric(coord[,col.z]) )
    {
      coord[,col.z] <- as.numeric(as.character(coord[,col.z]))
    }
  }
  coord[,col.reach] <- as.character(coord[,col.reach])
  reach.ids <- sort(unique(coord[,col.reach]))
  n.reach <- length(reach.ids)
  
  # extract reach coordinates and calculate lengths:
  # ------------------------------------------------
  
  x.min   <- NA
  x.max   <- NA
  y.min   <- NA
  y.max   <- NA
  z.min   <- NA
  z.max   <- NA
  lengths <- rep(0,n.reach)
  x.start <- rep(NA,n.reach)
  x.end   <- rep(NA,n.reach)
  y.start <- rep(NA,n.reach)
  y.end   <- rep(NA,n.reach)
  z.start <- rep(NA,n.reach)
  z.end   <- rep(NA,n.reach)
  reaches <- list()
  errors  <- character(0)
  for ( i in 1:length(reach.ids) )
  {
    # initialize empty structure:
    
    reaches[[i]] <- list()
    
    # assign coordinates and update min and max:
    
    ind.reach <- coord[,col.reach] == reach.ids[i]
    n.coord   <- sum(ind.reach)
    reaches[[i]]$n <- n.coord
    reaches[[i]]$x <- coord[ind.reach,col.x]
    reaches[[i]]$y <- coord[ind.reach,col.y]
    if ( is.na(col.z) )
    {
      reaches[[i]]$z <- rep(NA,n.coord)
    }
    else
    {
      reaches[[i]]$z <- coord[ind.reach,col.z]
    }
    if ( sum(is.na(c(reaches[[i]]$x,reaches[[i]]$y))) > 0 )
    {
      errors <- c(errors,paste("x or y coordinates of reach \"",reach.ids[i],"\" contain NAs",sep=""))
    }
    x.min <- min(x.min,reaches[[i]]$x,na.rm=T)
    x.max <- max(x.max,reaches[[i]]$x,na.rm=T)
    y.min <- min(y.min,reaches[[i]]$y,na.rm=T)
    y.max <- max(y.max,reaches[[i]]$y,na.rm=T)
    z.min <- min(z.min,reaches[[i]]$z,na.rm=T)
    z.max <- max(z.max,reaches[[i]]$z,na.rm=T)
    
    # calculate length:
    
    length <- 0
    if ( n.coord > 1 )
    {
      length <- sum(sqrt(diff(reaches[[i]]$x)^2+diff(reaches[[i]]$y)^2))
    }
    lengths[i] <- length
    reaches[[i]]$length <- length
    
    # assign start and end coordinates (according to order from file):

    x.start[i] <- reaches[[i]]$x[1]
    y.start[i] <- reaches[[i]]$y[1]
    z.start[i] <- reaches[[i]]$z[1]
    x.end[i]   <- reaches[[i]]$x[n.coord]
    y.end[i]   <- reaches[[i]]$y[n.coord]
    z.end[i]   <- reaches[[i]]$z[n.coord]
  }
  if (length(errors) > 0 )
  {
    cat("Problems reading reach coordinates:\n",paste(errors,collapse="\n"),"\n")
    return(NA)
  }
  names(reaches) <- reach.ids

  # identify nodes:
  # ---------------
  
  dist2 <- (x.end-x.start)^2 + (y.end-y.start)^2  # squared distances between start and end points
  if ( sum(dist2==0) > 0 )
  {
    cat("Reach(es) has(ve) same start and end point:",paste(reach.ids[dist2==0],collapse=","),"\n")
    return(NA)
  }
  tol2 <- min(tol^2,0.01*min(dist2))  # tolerance^2; torerance is not larger than 10 % of minimum distance start-end
  
  if ( n.reach == 1 )
  {
    node.start <- 1
    node.end   <- 2
    n.node         <- 2
  }
  else
  {
    node.start <- rep(NA,n.reach)
    node.end   <- rep(NA,n.reach)
    n.node <- 0
    for ( i in 1:n.reach )
    {
      if ( is.na(node.start[i]) )
      {
        n.node <- n.node + 1
        node.start[i] <- n.node
      }
      if ( is.na(node.end[i]) )
      {
        n.node <- n.node + 1
        node.end[i] <- n.node
      }
      if ( i < n.reach )
      {
        dist2 <- (x.start[i]-x.start[(i+1):n.reach])^2 + (y.start[i]-y.start[(i+1):n.reach])^2
        if ( sum(dist2<tol2) > 0 )
        {
          offset <- which(dist2<tol2)
          node.start[i+offset] <- node.start[i] 
        }
        dist2 <- (x.start[i]-x.end[(i+1):n.reach])^2 + (y.start[i]-y.end[(i+1):n.reach])^2
        if ( sum(dist2<tol2) > 0 )
        {
          offset <- which(dist2<tol2)
          node.end[i+offset] <- node.start[i] 
        }
        dist2 <- (x.end[i]-x.start[(i+1):n.reach])^2 + (y.end[i]-y.start[(i+1):n.reach])^2
        if ( sum(dist2<tol2) > 0 )
        {
          offset <- which(dist2<tol2)
          node.start[i+offset] <- node.end[i] 
        }
        dist2 <- (x.end[i]-x.end[(i+1):n.reach])^2 + (y.end[i]-y.end[(i+1):n.reach])^2
        if ( sum(dist2<tol2) > 0 )
        {
          offset <- which(dist2<tol2)
          node.end[i+offset] <- node.end[i] 
        }
      }
    }
  }

  # construct reach attributes:
  # ---------------------------
  
  attrib.reach <- data.frame(Reach_ID   = reach.ids,
                             Reach      = 1:n.reach,
                             x_start    = x.start,
                             y_start    = y.start,
                             z_start    = z.start,
                             x_end      = x.end,
                             y_end      = y.end,
                             z_end      = z.end,
                             node_start = node.start,
                             node_end   = node.end,
                             length     = lengths)

  # construct node attributes:
  # --------------------------
  
  node.ids <- 1:n.node
  x.node <- rep(NA,n.node)
  y.node <- rep(NA,n.node)
  for ( i in 1:n.node )
  {
    x1 <- numeric(0) 
    y1 <- numeric(0)
    ind.start <- node.start == i
    if ( sum(ind.start) > 0 ) 
    {
      x1 <- x.start[ind.start]
      y1 <- y.start[ind.start]
    }
    x2 <- numeric(0) 
    y2 <- numeric(0)
    ind.end <- node.end == i
    if ( sum(ind.end) > 0 ) 
    {
      x2 <- x.end[ind.end]
      y2 <- y.end[ind.end]
    }
    x.node[i] <- mean(c(x1,x2))
    y.node[i] <- mean(c(y1,y2))
  }
  
  attrib.node <- data.frame(Node = 1:n.node,
                            x    = x.node,
                            y    = y.node)
  
  # summarize input statistics:
  # ---------------------------

  if ( verbose )
  {
    cat("Number of reaches read:             ",n.reach,"\n")
    cat("Number of nodes identified:         ",n.node,"\n")
    cat("Reach lengths:                      ",min(lengths),"-",max(lengths),"\n")
    cat("Total network length:               ",sum(lengths),"\n")
  }

  # read reach attributes:
  # ======================
  
  if ( is.na(file.reachattrib[1]) )
  {
    if ( verbose ) cat("No reach attributes provided\n")
  }
  else
  {
    # read data:
    
    attrib <- read.table(file.reachattrib[1],header=TRUE,sep=sep,stringsAsFactors=FALSE,...)
    if ( length(file.reachattrib) > 1 )
    {
      for ( i in 2:length(file.reachattrib) )
      {
        attrib <- rbind(attrib,read.table(file.reachattrib[i],header=TRUE,sep=sep,stringsAsFactors=FALSE,...))   
      }
    }
    
    # check presence of column reach:
    
    col.reach <- match(colnames["reach"],colnames(attrib))
    if ( is.na(col.reach) )
    {
      cat("File \"",file.reachattrib,"\" must contain column \"",colnames["reach"],"\"\n",sep="")
      cat("No reach attributes read\n")
    }
    else
    {
      # ensure that reach identifiers are characters:
      
      attrib[,col.reach] <- as.character(attrib[,col.reach])
      
      # sort file according to reach structure and merge with derived attributes:
      
      n.attrib <- ncol(attrib)-1
      rows <- match(names(reaches),attrib[,col.reach])
      n.noattrib <- sum(is.na(rows))
      attrib <- rbind(attrib,rep(NA,ncol(attrib)))
      rows <- ifelse(is.na(rows),nrow(attrib),rows)
      attrib.reach <- cbind(attrib.reach,attrib[rows,-col.reach])
      if ( verbose ) 
      {
        cat("Number of reach attributes read:    ",n.attrib,"\n")
        if ( n.noattrib > 0 ) cat("No attributes found for",n.noattrib,"reach(es)\n")
      }
    }
  }
  
  # read node attributes:
  # =====================
  
  # read node attribute file(s) (containing columns for x, y and additional columns for attributes: 
  
  if ( is.na(file.nodeattrib[1]) )
  {
    if ( verbose ) cat("No node attributes provided\n")
  }
  else
  {
    # read data:
    
    attrib <- read.table(file.nodeattrib[1],header=TRUE,sep=sep,stringsAsFactors=FALSE,...)
    if ( length(file.nodeattrib) > 1 )
    {
      for ( i in 2:length(file.nodeattrib) )
      {
        attrib <- rbind(attrib.node,read.table(file.nodeattrib[i],header=TRUE,sep=sep,stringsAsFactors=FALSE,...))   
      }
    }
    
    # check presence of columns x and y:
    
    columns.missing <- character(0)
    col.x <- match(colnames["x"],colnames(attrib))
    col.y <- match(colnames["y"],colnames(attrib))
    if ( is.na(col.x) ) columns.missing <- c(columns.missing,colnames["x"])
    if ( is.na(col.y) ) columns.missing <- c(columns.missing,colnames["y"])
    if ( length(columns.missing) > 0 )
    {
      cat("file \"",file.nodeattrib,"\" must contain column(s) \"",paste(columns.missing,collapse="\", \""),"\"\n",sep="")
      cat("No node attributes read")
    }
    else
    {
      
      # ensure that node coordinates are numeric (factors must first be converted to char)
      # and identifiers are characters:
      
      if ( !is.numeric(attrib[,col.x]) )
      {
        attrib[,col.x] <- as.numeric(as.character(attrib[,col.x]))
      }
      if ( !is.numeric(attrib[,col.y]) )
      {
        attrib[,col.y] <- as.numeric(as.character(attrib[,col.y]))
      }
      col.node <- match(colnames["node"],colnames(attrib.node))
      if ( !is.na(col.node) )
      {
        attrib[,col.node] <- as.character(attrib[,col.node])
      }    
      
      # match nodes:
      
      rows <- rep(NA,n.node)
      errors <- character(0)
      for ( i in 1:n.node )
      {
        ind <- which( (x.node[i]-attrib[,col.x])^2 + (y.node[i]-attrib[,col.y])^2 < tol2 )
        if ( length(ind) == 1 )
        {
          rows[i] <- ind
        }
        else
        {
          if ( length(ind) > 1 )
          {
            errors <- c(errors,paste("Multiple node attributes found for node at (",
                                     x.node[i],",",y.node[i],")",sep=""))
            rows[i] <- ind[1]
          }
        }
      }
      if (length(errors) > 0 )
      {
        cat("Problems reading node attributes:\n",paste(errors,collapse="\n"),"\n",sep="")
      }
      
      if ( sum(!is.na(rows)) == 0 )
      {
        if ( verbose ) cat("No node attributes found\n")
      }
      else
      {    
        n.attrib <- ncol(attrib)-2
        n.noattrib <- sum(is.na(rows))
        attrib <- rbind(attrib,rep(NA,ncol(attrib)))
        rows <- ifelse(is.na(rows),nrow(attrib),rows)
        attrib.node <- cbind(attrib.node,attrib[rows,])
        if ( verbose ) 
        {
          cat("Number of node attributes read:     ",n.attrib,"\n")
          if ( n.noattrib > 0 ) cat("No attributes found for",n.noattrib,"node(s)\n")
        }
      }
    }
  }
  
  # construct rivernet object:
  # ==========================
  
  rivernet <- list()
  rivernet$reaches      <- reaches
  rivernet$xlim         <- c(x.min,x.max)
  rivernet$ylim         <- c(y.min,y.max)
  rivernet$zlim         <- c(z.min,z.max)
  rivernet$htow         <- (y.max-y.min)/(x.max-x.min)
  rivernet$total.length <- sum(lengths)
  rivernet$attrib.reach <- attrib.reach
  rivernet$attrib.node  <- attrib.node
  class(rivernet)       <- "rivernet"
  
  # analyze network:
  # ================
  
  if ( analyze ) rivernet <- analyze(rivernet,verbose=verbose)
  
  return(rivernet)
}

  
  
  
  
  
  
  
  
  
  