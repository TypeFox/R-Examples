#' Calculate a single EOT
#' 
#' @description
#' EotCycle() calculates a single EOT and is controlled by the main eot() function
#' 
#' @param x a ratser stack used as predictor
#' @param y a RasterStack used as response. If \code{y} is \code{NULL},
#' \code{x} is used as \code{y}
#' @param n the number of EOT modes to calculate
#' @param standardised logical. If \code{FALSE} the calculated r-squared values 
#' will be multiplied by the variance
#' @param orig.var original variance of the response domain
#' @param write.out logical. If \code{TRUE} results will be written to disk 
#' using \code{path.out}
#' @param path.out the file path for writing results if \code{write.out} is \code{TRUE}.
#' Defaults to current working directory
#' @param prefix optional prefix to be used for naming of results if 
#' \code{write.out} is \code{TRUE}
#' @param type the type of the link function. Defaults to \code{'rsq'} as in original
#' proposed method from \cite{Dool2000}. If set to \code{'ioa'} index of agreement is
#' used instead
#' @param verbose logical. If \code{TRUE} some details about the 
#' calculation process will be output to the console
#' @param ... not used at the moment
#' 
#' @export EotCycle
EotCycle <- function(x, 
                     y, 
                     n = 1,
                     standardised, 
                     orig.var,
                     write.out,
                     path.out,
                     prefix,
                     type,
                     verbose,
                     ...) {
  
    ### Identification of the most explanatory pred pixel
  
  # Extract pixel entries from RasterStack objects
  x.vals <- raster::getValues(x)
  y.vals <- raster::getValues(y)
  type <- type[1]
  
  # Calculate and summarize R-squared per pred pixel
  if (verbose) {
    cat("\nCalculating linear model ...", "\n")
  }
  
  type <- type[1]
  if (type == "rsq") {
    a <- predRsquaredSum(pred_vals = x.vals, resp_vals = y.vals, 
                         standardised = standardised)
  } else {
    a <- iodaSumC(pred_vals = x.vals, resp_vals = y.vals)
  }
  
  # Identify pred pixel with highest sum of r.squared
  if (verbose) {
    cat("Locating ", n, ". EOT ...", "\n", sep = "")
  }
  
  maxxy.all <- which(a == max(a, na.rm = TRUE))
  maxxy <- maxxy.all[1]

  if (length(maxxy.all) != 1) {
    if (verbose) {
      cat("WARNING:", "\n",
          "LOCATION OF EOT AMBIGUOUS!",  "\n",
          "MULTIPLE POSSIBLE LOCATIONS DETECTED, USING ONLY THE FIRST!\n\n")
    }
  }

  if (verbose) {
    cat("Location:", raster::xyFromCell(x, maxxy), "\n", sep = " ")
  }
  
  ### Regression of most explanatory pred pixel with resp pixels
    
  ## Fit lm
  
  # lm(y.vals[i, ] ~ x.vals[maxxy, ]) with T statistics
  y.lm.param.t <- respLmParam(x.vals, y.vals, maxxy - 1) # C++ starts at 0!
  # Calculate p value from T statistics
  y.lm.param.p <- lapply(y.lm.param.t, function(i) {
    tmp <- i
    tmp[[5]] <- 2 * pt(-abs(tmp[[5]]), df = tmp[[6]])
    
    return(tmp)
  })  
  
  
  ## Rasterize lm parameters
  
  # RasterLayer template for R-squared, slope and p value
  rst.y.template <- raster::raster(nrows = raster::nrow(y), 
                                   ncols = raster::ncol(y), 
                                   xmn = raster::xmin(y), 
                                   xmx = raster::xmax(y), 
                                   ymn = raster::ymin(y), 
                                   ymx = raster::ymax(y))
  
  rst.y.r <- rst.y.rsq <- rst.y.intercept <- 
    rst.y.slp <- rst.y.p <- rst.y.template

  # RasterBrick template for residuals
  brck.y.resids <- raster::brick(nrows = raster::nrow(y), 
                                 ncols = raster::ncol(y), 
                                 xmn = raster::xmin(y), 
                                 xmx = raster::xmax(y), 
                                 ymn = raster::ymin(y), 
                                 ymx = raster::ymax(y), 
                                 nl = raster::nlayers(y))
  
  # R
  rst.y.r[] <- sapply(y.lm.param.p, "[[", 1)
  # R-squared
  rst.y.rsq[] <- sapply(y.lm.param.p, "[[", 1) ^ 2
  # Intercept
  rst.y.intercept[] <- sapply(y.lm.param.p, "[[", 2)
  # Slope
  rst.y.slp[] <- sapply(y.lm.param.p, "[[", 3)
  # P value
  rst.y.p[] <- sapply(y.lm.param.p, "[[", 5)
  # Residuals
  brck.y.resids[] <- matrix(sapply(y.lm.param.p, "[[", 4), 
                            ncol = raster::nlayers(x), byrow = TRUE)
  # EOT over time
  eot.ts <- as.numeric(raster::extract(x, maxxy)[1, ])
  
  ### Regression of most explanatory pred pixel with pred pixels
  
  # Following code is only executed when pred and resp are not equal

    ## Fit lm
    
    # lm(x.vals[i, ] ~ x.vals[maxxy, ]) with T statistics
    x.lm.param.t <- respLmParam(x.vals, x.vals, maxxy - 1) # C++ starts at 0!
    # Calculate p value from T statistics
    x.lm.param.p <- lapply(x.lm.param.t, function(i) {
      tmp <- i
      tmp[[5]] <- 2 * pt(-abs(tmp[[5]]), df = tmp[[6]])
      
      return(tmp)
    })  
    
    
    ## Rasterize lm parameters
    
    # RasterLayer template for R-squared, slope and p value
  rst.x.template <- raster::raster(nrows = raster::nrow(x), 
                                   ncols = raster::ncol(x), 
                                   xmn = raster::xmin(x), 
                                   xmx = raster::xmax(x), 
                                   ymn = raster::ymin(x), 
                                   ymx = raster::ymax(x))
    
    rst.x.r <- rst.x.rsq <- rst.x.rsq.sums <- rst.x.intercept <- 
      rst.x.slp <- rst.x.p <- rst.x.template
    
    # RasterBrick template for residuals
  brck.x.resids <- raster::brick(nrows = raster::nrow(x), 
                                 ncols = raster::ncol(x), 
                                 xmn = raster::xmin(x), 
                                 xmx = raster::xmax(x), 
                                 ymn = raster::ymin(x), 
                                 ymx = raster::ymax(x), 
                                 nl = raster::nlayers(x))
  
    # R
    rst.x.r[] <- sapply(x.lm.param.p, "[[", 1)
    # R-squared
    rst.x.rsq[] <- sapply(x.lm.param.p, "[[", 1) ^ 2
    # R-squared sums
    rst.x.rsq.sums[] <- a
    # Intercept
    rst.x.intercept[] <- sapply(x.lm.param.p, "[[", 2)
    # Slope
    rst.x.slp[] <- sapply(x.lm.param.p, "[[", 3)
    # P value
    rst.x.p[] <- sapply(x.lm.param.p, "[[", 5)
    # Residuals
    brck.x.resids[] <- matrix(sapply(x.lm.param.p, "[[", 4), 
                              ncol = raster::nlayers(x), byrow = TRUE)
  
#     #expl.var <- x[maxxy] / orig.var
#   if (!standardised) {
#     t <- mean(apply(getValues(brck.y.resids), 1, var, na.rm = TRUE), 
#               na.rm = TRUE)
#     s <- mean(apply(getValues(brck.y.resids), 2, var, na.rm = TRUE), 
#               na.rm = TRUE)
#     resid.var <- t + s
#   } else {
#     resid.var <- var(as.vector(getValues(brck.y.resids)), na.rm = TRUE)
#   }
  
  resid.var <- calcVar(brck.y.resids, standardised = standardised)
  
  cum.expl.var <- (orig.var - resid.var) / orig.var
  
  if (verbose) {
    cat("Cum. expl. variance (%):", cum.expl.var * 100, "\n", sep = " ")
  }
  
  xy <- raster::xyFromCell(x, maxxy)
  location.df <- as.data.frame(cbind(xy, paste("mode", 
                                               sprintf("%02.f", n), 
                                               sep = "_"),
                                     cum.expl.var,
                                     if (length(maxxy.all) != 1) 
                                       "ambiguous" else "ok"),
                               stringsAsFactors = FALSE)
  names(location.df) <- c("x", "y", "mode", "cum_expl_var", "comment")
  mode(location.df$x) <- "numeric"
  mode(location.df$y) <- "numeric"
  mode(location.df$cum_expl_var) <- "numeric"
    
    ### Output
    
    # Output returned by function
  out <- new('EotMode',
             mode = n,
             name = paste("mode", sprintf("%02.f", n), sep = "_"),
             eot = eot.ts,
             coords_bp = xy,
             cell_bp = maxxy,
             cum_exp_var = cum.expl.var,
             r_predictor = rst.x.r,
             rsq_predictor = rst.x.rsq,
             rsq_sums_predictor = rst.x.rsq.sums,
             int_predictor = rst.x.intercept, 
             slp_predictor = rst.x.slp,
             p_predictor = rst.x.p,
             resid_predictor = brck.x.resids,
             r_response = rst.y.r,
             rsq_response = rst.y.rsq,
             int_response = rst.y.intercept, 
             slp_response = rst.y.slp,
             p_response = rst.y.p,
             resid_response = brck.y.resids)
  
    # Output storage (optional)
    if (write.out) {
      writeEot(out, path.out = path.out, prefix = prefix, ...)
      
      df.name <- paste(prefix, "eot_locations.csv", sep = "_")
      
      if (n == 1) {        
        write.table(location.df, col.names = TRUE, 
                    paste(path.out, df.name, sep = "/"), 
                    row.names = FALSE, append = FALSE, sep = ",")
      } else {
        write.table(location.df, col.names = FALSE, 
                    paste(path.out, df.name, sep = "/"), 
                    row.names = FALSE, append = TRUE, sep = ",")
      }
      
      rm(list = c("eot.ts",
                  "maxxy",
                  "location.df",
                  "expl.var",
                  "rst.x.r",
                  "rst.x.rsq",
                  "rst.x.rsq.sums",
                  "rst.x.intercept", 
                  "rst.x.slp",
                  "rst.x.p",
                  "brck.x.resids",
                  "rst.y.r",
                  "rst.y.rsq",
                  "rst.y.intercept", 
                  "rst.y.slp",
                  "rst.y.p",
                  "brck.y.resids"))
      gc()
      
    }

  return(out)
  
}