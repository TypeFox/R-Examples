###########################################
#### AUTHOR:     Arnost Komarek        ####
####             (2004)                ####
####                                   ####
#### FILE:       bayesDensity.R        ####
####                                   ####
#### FUNCTIONS:  bayesDensity          ####
####             print.bayesDensity    ####
####             plot.bayesDensity     ####
###########################################

### ======================================
### bayesDensity
### ======================================
bayesDensity <- function(dir = getwd(),
                         stgrid,
                         centgrid,
                         grid,
                         n.grid = 100,
                         skip = 0,
                         by = 1,
                         last.iter,
                         standard = TRUE,
                         center = TRUE,
                         unstandard = TRUE)
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  
  ## Check whether needed files are available
  ## * further,  whether at least first row has correct number of elements
  ## * and determine the MC sample size
  ## ==============================================================
  filesindir <- dir(dir)    ## character vector with available files
  if (!length(filesindir)) stop("Empty directory with simulated values?")  
  if (sum(!is.na(match(filesindir, "mweight.sim")))){
    mix <- read.table(paste(dir, "/mweight.sim", sep = ""), nrows = 1)
    k.max <- length(mix)
  }
  else
    stop("File with simulated values of mixture weights not found.")

  if (sum(!is.na(match(filesindir, "mmean.sim")))){
    mix <- read.table(paste(dir, "/mmean.sim", sep = ""), nrows = 1)
    kmax2 <- length(mix)
    if (k.max != kmax2) stop("Different kmax indicated by files mweight.sim and mmean.sim.")
  }
  else
    stop("File with simulated values of mixture means not found.")
   
  if (sum(!is.na(match(filesindir, "mvariance.sim")))){
    mix <- read.table(paste(dir, "/mvariance.sim", sep = ""), nrows = 1)     
    kmax2 <- length(mix)
    if (k.max != kmax2) stop("Different kmax indicated by files mweight.sim and mvariance.sim.")
  }
  else
    stop("File with simulated values of mixture variances not found.")

  ## nsimul, skip, by  
  if (sum(!is.na(match(filesindir, "mixmoment.sim")))){
    mix <- read.table(paste(dir, "/mixmoment.sim", sep = ""), header = TRUE)
    if (missing(last.iter)) M <- dim(mix)[1]
    else{
      M <- last.iter
      if (last.iter > dim(mix)[1]) M <- dim(mix)[1]
      if (last.iter <= 0)          M <- dim(mix)[1]
    }      
  }
  else
    stop("File mixmoment.sim not found.")     
    
  if (missing(skip)) skip <- 0
  else{
    if (skip > M) stop("You ask to skip more rows from the file than available.")
    if (is.na(skip) || skip < 0) skip <- 0
  }
  if (missing(by)) by <- 1
  else{
    if (is.na(by) || by <= 0) by <- 1
  }    

  lvalue <- 1 + (M - skip - 1) %/% by
#  if (missing(nwrite)) nwrite <- lvalue
#  if (nwrite > lvalue) nwrite <- lvalue
    
  k.cond <- 1:k.max
  
  ## Try to guess the grid (from first at most 20 mixtures) if not given by the user
  if (!unstandard) grid <- 0
  if (!standard) stgrid <- 0
  if (!center) centgrid <- 0
  if (missing(grid)) miss.grid <- TRUE
  else               miss.grid <- FALSE
  if (missing(centgrid)) miss.centgrid <- TRUE
  else                   miss.centgrid <- FALSE
  if (miss.grid | miss.centgrid){
    mus <- scan(paste(dir, "/mmean.sim", sep = ""), nlines = 20, skip = 1)
    sigma2s <- scan(paste(dir, "/mvariance.sim", sep = ""), nlines = 20, skip = 1)    
    mu.min <- min(mus)
    mu.max <- max(mus)
    sd.min <- sqrt(min(sigma2s))
    sd.max <- sqrt(max(sigma2s))
    sd.mean <- sqrt(mean(sigma2s))
    if (miss.grid) grid <- seq(mu.min - 2.5*sd.mean, mu.max + 2.5*sd.mean, length = n.grid)
    if (miss.centgrid) centgrid <- seq(-2.5*sd.mean, 2.5*sd.mean, length = n.grid)
  }
  if (missing(stgrid)){
    stgrid <- seq(-2.5, 2.5, length = n.grid)
  }
  if (!unstandard) grid <- 0
  if (!center) centgrid <- 0
  ngrid <- length(grid)
  nstgrid <- length(stgrid)
  ncentgrid <- length(centgrid)
  if (!unstandard) ngrid <- 0
  if (!standard) nstgrid <- 0
  if (!center) ncentgrid <- 0

  mcmc <- .C("bayesDensity", aver = double(ngrid * (1 + k.max)),           staver = double(nstgrid * (1 + k.max)),
                             centaver = double(ncentgrid * (1 + k.max)),
                             intercept = double(lvalue),                   scale = double(lvalue),
                             M.k = integer(1 + k.max),                     as.character(dir),
                             as.double(grid),                              as.double(stgrid),
                             as.double(centgrid),
                             as.integer(k.max),                            as.integer(M),
                             as.integer(skip),                             as.integer(by),
                             as.integer(ngrid),                            as.integer(nstgrid),
                             as.integer(ncentgrid),
                             err = integer(1),
             PACKAGE = thispackage)
  
  if (mcmc$err) stop("No results produced, something is wrong.")

  if (unstandard){
    mcmc$aver <- matrix(mcmc$aver, nrow = ngrid)
    mcmc$aver <- cbind(grid, mcmc$aver)
    mcmc$aver <- as.data.frame(mcmc$aver)
    rownames(mcmc$aver) <- paste(1:ngrid)
    colnames(mcmc$aver) <- c("grid", "unconditional", paste("k = ", 1:k.max, sep = ""))    
  }
  if (standard){
    mcmc$staver <- matrix(mcmc$staver, nrow = nstgrid)
    mcmc$staver <- cbind(stgrid, mcmc$staver)
    mcmc$staver <- as.data.frame(mcmc$staver)
    rownames(mcmc$staver) <- paste(1:nstgrid)
    colnames(mcmc$staver) <- c("grid", "unconditional", paste("k = ", 1:k.max, sep = ""))        
  }
  if (center){
    mcmc$centaver <- matrix(mcmc$centaver, nrow = ncentgrid)
    mcmc$centaver <- cbind(centgrid, mcmc$centaver)
    mcmc$centaver <- as.data.frame(mcmc$centaver)
    rownames(mcmc$centaver) <- paste(1:ncentgrid)
    colnames(mcmc$centaver) <- c("grid", "unconditional", paste("k = ", 1:k.max, sep = ""))        
  }
  names(mcmc$M.k) <- c("unconditional", paste(1:k.max))
  names(mcmc$intercept) <- paste(1:lvalue)
  names(mcmc$scale) <- paste(1:lvalue)    
    
  charact <- data.frame(intercept = mcmc$intercept, scale = mcmc$scale)
  rownames(charact) <- paste(1:lvalue)          

  density <- list()
  if (standard) density$standard <- mcmc$staver
  else          density$standard <- "Not asked."

  if (unstandard) density$unstandard <- mcmc$aver
  else            density$unstandard <- "Not asked."

  if (center) density$center <- mcmc$centaver
  else        density$center <- "Not asked."  
  
  attr(density, "sample.size") <- mcmc$M.k
  attr(density, "moments") <- charact
  attr(density, "k") <- data.frame(k = mix[,1])
  
  class(density) <- "bayesDensity"
  return(density)   
}  


### ======================================
### print.bayesDensity
### ======================================
print.bayesDensity <- function(x, ...)
{
  cat("\nStandardized McMC average: \n")
  print(x$standard)  

  cat("\nUnstandardized McMC average: \n")
  print(x$unstandard)  

  cat("\nCentered McMC average: \n")
  print(x$center)  
  
  return(invisible(x))
}  


### ======================================
### plot.bayesDensity
### ======================================
plot.bayesDensity <- function(x,
                              k.cond,
                              dim.plot = TRUE,
                              over = TRUE,
                              alegend = TRUE,
                              standard = TRUE,
                              center = FALSE,
                              type = "l",
                              bty = "n",
                              xlab = expression(epsilon),
                              ylab = expression(f(epsilon)),
                              lty,
                              xlim,
                              ylim,
                              xleg,
                              yleg,
                              main,
                              ...
                              )
{
  MM.k <- attr(x, "sample.size")
  kmax <- length(MM.k) - 1
  kmax.insample <- max((0:kmax)[MM.k > 0])
  
  if (missing(k.cond)) k.cond <- 0:kmax.insample
  not.sampled <- k.cond > kmax.insample
  k.cond <- k.cond[!not.sampled]

  num.plots <- length(k.cond)
  dims <- switch(num.plots, c(1, 1),
                            c(1, 2),
                            c(2, 2), c(2, 2),
                            c(2, 3), c(2, 3),
                            c(3, 3), c(3, 3), c(3, 3),
                            c(3, 4), c(3, 4), c(3, 4),
                            c(4, 4), c(4, 4), c(4, 4), c(4, 4),
                            c(5, 4), c(5, 4), c(5, 4), c(5, 4))
  if (dim.plot)        par(mfrow = dims)
  if (dim.plot & over) par(mfrow = c(1, 1))
  if (missing(lty)) lty <- 1:num.plots
  

  empty.plot <- function(mess = "Not in the sample")
  {
    plot(0:100, 0:100, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    text(50, 90, labels = mess, adj = 1)
  }    
  

  M <- MM.k[1]
  M.k <- MM.k[-1]
  nulls <- M.k == 0
  
  if (over & k.cond[1] & !sum(M.k[k.cond])){
    empty.plot(mess = "Conditional densities you asked are not in the sample.")
    return(invisible(x))
  }    

  if (standard){ mcmc <- x$standard }
  else{
    if (center){ mcmc <- x$center }
    else       { mcmc <- x$unstandard }
  }  

  if (is.character(mcmc)) stop("You have to first compute McMC averages.")
  grid <- mcmc$grid
  
  remove <- c(TRUE, FALSE, nulls)
  mcmc0 <- mcmc[!remove]
  if (missing(xlim)) xlim <- range(grid)
  if (missing(ylim)) ylim <- c(0, max(sapply(mcmc0, max, na.rm = TRUE), na.rm = TRUE))

  legg <- character(0)
  relatM <- round((M.k/M)*100, 2)
  if (over){
    i <- 1
    if (!k.cond[1]){
      dens <- mcmc$unconditional
      legg <- c(legg, paste("Uncond.,  ", "M = ", M, sep = ""))
    }
    else{
      while (!M.k[k.cond[i]]){
        legg <- c(legg, paste("k = ", k.cond[i], "   (", relatM[k.cond[i]], " %)",  sep = ""))
        i <- i + 1 
      }  
      dens <- mcmc[[paste("k = ", k.cond[i], sep = "")]]
      legg <- c(legg, paste("k = ", k.cond[i], "   (", relatM[k.cond[i]], " %)",  sep = ""))      
    }
    plot(grid, dens, type = type, bty = bty, xlab = xlab, ylab = ylab, lty = lty[i], xlim = xlim, ylim = ylim)
    if (i < length(k.cond)){
      for (j in (i+1):length(k.cond)){
        if (!M.k[k.cond[j]]){
          legg <- c(legg, paste("k = ", k.cond[j], "   (", relatM[k.cond[j]], " %)", sep = ""))
          next
        }        
        dens <- mcmc[[paste("k = ", k.cond[j], sep = "")]]
        legg <- c(legg, paste("k = ", k.cond[j], "   (", relatM[k.cond[j]], " %)", sep = ""))      
        lines(grid, dens, lty = lty[j])  
      }
    }      
    if (missing(xleg)) xleg <- min(grid)
    if (missing(yleg)) yleg <- ylim[2] - 0.1*(ylim[2] - ylim[1]) 
    if (alegend) legend(xleg, yleg, legend = legg, lty = lty, xjust = 0, yjust = 1, bty = "n")
    if (missing(main)) title(main = "McMC averages of the density")
    else               title(main = main)
  }    
  else{    ## not over
    for (k in k.cond){
      if (k == 0){
        dens <- mcmc$unconditional
        titul <- "unconditional"
        subb <- paste("M = ", M, sep = "")
        if (missing(xlim)) xlim <- range(grid)
        if (missing(ylim)) ylim <- range(dens)
        plot(grid, dens, type = type, bty = bty, xlab = xlab, ylab = ylab, lty = lty[1], xlim = xlim, ylim = ylim)
      }      
      else{
        dens <- mcmc[[paste("k = ", k, sep = "")]]
        titul <- paste("k = ", k, sep = "")
        subb <- paste("M = ", M.k[k], "   (", relatM, " %)", sep = "")
        if (!M.k[k]){
          empty.plot()
          subb <- ""
        }        
        else{
          if (missing(xlim)) xlim <- range(grid)
          if (missing(ylim)) ylim <- range(dens)          
          plot(grid, dens, type = type, bty = bty, xlab = xlab, ylab = ylab, lty = lty[1], xlim = xlim, ylim = ylim)
        }  
      }  
      title(main = titul, sub = subb)
    }
  }    

  return(invisible(x))  
}  




