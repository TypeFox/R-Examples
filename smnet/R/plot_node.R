
plot_node<-function(x, coords = NULL, col = 3, cex = NULL, key = key, se = F, ...)
{
  # break up beta_hat
  getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
  cov.dims  <- lapply(x$internals$X.list, getcol)    
  inds      <- unlist(cov.dims)
  cum.inds  <- cumsum(inds)
  ord       <- order(x$internals$ord)
  n.cov     <- length(x$internals$X.list)
  ests      <- vector("list", length = n.cov)
  ests[[1]] <- x$internals$beta_hat[1]
  for(i in 1:(n.cov-1)){ests[[i+1]]<-x$internals$beta_hat[(cum.inds[i]+1):(cum.inds[i+1])]} 
  fitted_segment <- ests[[n.cov]][ord]# + ests[[1]]
  if(is.null(coords)){
    # get midpoint locations for each rid, then put in order
    data          <- getSSNdata.frame(x$ssn.object, Name = "Obs")
    get_rid_midpt <- function(L) apply(L@Lines[[1]]@coords, 2, FUN = "mean") 
    midpt_rids    <- lapply(x$ssn.object@lines, FUN = get_rid_midpt)
    add.one       <- ifelse(min(as.numeric(x$internals$adjacency$rid_bid[,1])) == 0, 1, 0)
    coords        <- Reduce("rbind", midpt_rids)[as.numeric(x$internals$adjacency$rid_bid[,1]) + add.one, ]
    coords        <- coords[x$internals$ord,]
  } 
  
  if(se){
    l.covs            <- lapply(x$internals$X.list, ncol)
    component.size    <- l.covs[[length(l.covs)]]
    spatial.n         <- nrow(x$internals$adjacency[[1]])
    empty.X.list      <- lapply(l.covs, make_sparse, nrow = spatial.n)
    empty.X.list[[length(l.covs)]] <- diag.spam(spatial.n)
    new.X             <- Reduce("cbind", empty.X.list)
    left1             <- forwardsolve.spam(x$internals$U, t(new.X), transpose = T)
    left2             <- backsolve.spam(x$internals$U, left1, transpose = T)
    vec               <- x$internals$X.spam %*% left2
    fitted_segment    <- sqrt(colSums(vec*vec)*x$internals$sigma.sq)   
  }
  
  # SET UP THE COLOUR SCALE FOR PLOTTING
  nlevels      <- 25
  zlim         <- range(fitted_segment, finite = TRUE)
  brks         <- pretty(zlim, n = 10)
  nclasses     <- length(brks) - 1 
  nnodes       <- length(fitted_segment)
  colorPalette <- rev(heat.colors(nclasses))
  mar.orig     <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w            <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  mar     <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  
  # SET UP THE PLOTTING REGION
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(brks), xaxs = "i", yaxs = "i")
  
  # PLOT THE LEGEND IF REQUESTED
  if(key){
    rect(0, brks[-length(brks)], 1, brks[-1L], col = colorPalette)
    axis(4)
    box()
  }
  
  # SET UP THE PLOTTING OF THE NETWORK NODES
  nnodes    <- length(fitted_segment)
  col.nums  <- cut(as.vector(fitted_segment), breaks = brks, labels = FALSE)
  coord     <- coords
  col       <- colorPalette[col.nums]
  if(is.null(cex)) cex <- 10/log(nnodes)

  # SET UP PLOTTING DEFAULTS AND UPDATE FROM ...
  default.parameters <- list(xlim = range(coord[,1]), ylim = range(coord[,2]))
  updated.parameters <- modifyList(default.parameters, list(...))
  
  
  # PLOT THE NETWORK NODES
  mar      <- mar.orig
  mar[4L]  <- 1
  par(mar = mar)
  plot.new()
  do.call(plot.window, updated.parameters)
  for(i in 1:ncol(x$internals$adjacency[[1]])){
    nums <- which(!x$internals$adjacency[[1]][,i] == 0)
    if(length(nums) > 0){
      for(j in 1:length(nums)){
        segments(x0 = coord[i,1], y0 = coord[i,2],
                 x1 = coord[nums[j],1], y1 = coord[nums[j],2], lwd = 1)
      }
    }
  }
  do.call(title, updated.parameters)
  if(is.null(updated.parameters$xaxt)){
    xlimits <- updated.parameters$xlim
    xticks  <- format(seq(xlimits[1], xlimits[2], length.out = 5), digits = 2, zero.print = T)
    Axis(x = xlimits, at = xticks, side = 1, labels = xticks)
  }
  if(is.null(updated.parameters$yaxt)){
    ylimits <- updated.parameters$ylim
    yticks  <- format(seq(ylimits[1], ylimits[2], length.out = 5), digits = 2, zero.print = T)
    Axis(x = ylimits, at = yticks, side = 2, labels = yticks)
  }
  points(coord, pch = 21, bg =  "black", cex = as.vector(cex)+0.1, col = 1)
  points(coord, pch = 21, bg =  as.vector(col), cex = as.vector(cex), col = 1)
}