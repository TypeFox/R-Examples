plot.lmf <-
function(x,
                     what = "total",
                     ...)
{
  # Create the DropInf function to remove any Inf leverages.
  dropInf <- function(x, h) {
    if (any(isInf <- h >= 1)) {
      warning("Not plotting observations with leverage one:\n  ", 
              paste(which(isInf), collapse = ", "), call. = FALSE)
      x[isInf] <- NaN
    }
    x
  }
  # Set the plotting device to plot one plot at the time and press enter to see next plot
  if((prod(par("mfcol")) < 4 || what == "age-year") && dev.interactive())
    devAskNewPage(TRUE)
  on.exit(devAskNewPage(FALSE))
  # Create some axis labels used in two or three plots
  # x-axis name
  xlab13 <- "Fitted values"
  # y-axis name
  ylab234 <- "Standardized residuals"
  # Create plotting functions
  # First plot: Residuals versus Fitted values
  fip <- function(fit, res, xlab = xlab13, ...){
    # y-axis limits
    ylim1 <- range(res, na.rm = TRUE)
    ylim1 <- extendrange(r = ylim1, f = 0.08)
    # plot
    plot(fit, res, xlab = xlab, ylab = "Residuals", 
         ylim = ylim1, type = "p")
    title(main = "Residuals vs Fitted", ...)
    abline(h = 0, lty = 3, col = "gray")
    panel.smooth(fit, res)
  }
  # Second plot: Normal Q-Q plot
  sp <- function(rs, ylab = ylab234, ...){
    # y-axis limits
    ylim <- range(rs, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    #plot
    qqnorm(rs, main = "", ylab = ylab, ylim = ylim)
    title(main = "Normal Q-Q", ...)
    qqline(rs, lty = 3, col = "gray50")          
  }
  # Third plot: Scale-loacation plot
  tp <- function(fit, rs, xlab = xlab13, ylab = ylab234, ...){
    # Squared standardized residuals
    sqrtabsr <- sqrt(abs(rs))
    # y-axis limits
    ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
    # y-axis name
    ylab3 <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab))))
    # plot
    plot(fit, sqrtabsr, xlab = xlab, ylab = ylab3,
         ylim = ylim)
    title(main = "Scale-Location", ...)
    abline(h = mean(sqrtabsr, na.rm = TRUE), lty = 3, col = "gray")
    panel.smooth(fit, sqrtabsr)
  }
  # Fourth plot: Leverage plot
  fop <- function(leverage, rs, r.hat, npar, ylab = ylab234, ...){
    # cook levels
    cook.levels <- c(0.5, 1)
    # Remove leverage > 1
    xx <- leverage
    xx[xx >= 1] <- NA
    # y-axis limits
    ylim <- range(rs, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    # plot
    plot(xx, rs, xlim = c(0, max(leverage, na.rm = TRUE)), 
         ylim = ylim, xlab = "Leverage", 
         ylab = ylab)
    title(main = "Residuals vs Leverage", ...)
    abline(h = 0, v = 0, lty = 3, col = "gray")
    panel.smooth(xx, rs)
    # Insert cook levels lines and legend
    usr <- par("usr")
    hh <- seq.int(min(r.hat[1L], r.hat[2L]/100), 
                  min(usr[2L], 0.999), length.out = 101)
    for (crit in cook.levels) {
      cl.h <- sqrt(crit * npar * (1 - hh)/hh)
      lines(hh, cl.h, lty = 2, col = 2)
      lines(hh, -cl.h, lty = 2, col = 2)
    }
    legend("bottomleft", legend = "Cook's distance", 
           lty = 2, col = 2, bty = "n")
    xmax <- min(0.99, usr[2L])
    ymult <- sqrt(npar * (1 - xmax)/xmax)
    aty <- c(-sqrt(rev(cook.levels)) * ymult, sqrt(cook.levels) *  ymult)
    axis(4, at = aty, labels = paste(c(rev(cook.levels),cook.levels)),
         mgp = c(0.25, 0.25, 0), las = 2, tck = 0, cex.axis = 0.75, col.axis = 2)
  }
  # Check what to plot, then plot
  if(!any(what %in% c("total", "age-year")))
    stop("Illegal 'what' argument provided.")
  if(what == "total"){
    # Extract and combine the lists containting residuals, fitted values, leverage
    # and cook's distance
    # Residuals
    res <- unlist(x$res)
    # Fitted values
    fit <- unlist(x$fit)
    # Leverage
    leverage <- unlist(x$leverage)
    # Cook's distance
    cook <- unlist(x$cook)
    # Calculate needed variables for plotting
    # Standardized residuals
    rs <- dropInf(res/(x$sigma2.d * sqrt(1 - leverage)), leverage)
    # Leverage range
    r.hat <- range(leverage, na.rm = TRUE)
    # First plot
    fip(fit = fit, res = res)
    # Second plot
    sp(rs = rs)
    # Third plot
    tp(fit = fit, rs = rs)
    # Fourth plot
    fop(leverage = leverage, rs = rs, r.hat = r.hat, npar = x$npar)
  }
  if(what == "age-year"){
    # For-loop over all ages and years
    for(i in x$uage){
      #Numeric age
      a <- which(x$uage == i)
      for(j in x$uyear){
        #Numeric year
        y <- which(x$uyear == j)
        # Extract and combine the lists containting residuals, fitted values, leverage
        # and cook's distance for the age and year in question
        # Residuals
        res <- unlist(x$res[[a]][[y]])
        # Fitted values
        fit <- unlist(x$fit[[a]][[y]])
        # Leverage
        leverage <- unlist(x$leverage[[a]][[y]])
        # Cook's distance
        cook <- unlist(x$cook[[a]][[y]])
        # Calculate needed variables for plotting
        # Standardized residuals
        rs <- dropInf(res/(x$sigma2.d * sqrt(1 - leverage)), leverage)
        # Leverage range
        r.hat <- range(leverage, na.rm = TRUE)
        # First plot
        subt <- paste(sprintf("Age %s", i), sprintf("in year %s", j))
        fip(fit = fit, res = res, sub = subt, font.sub = 3)
        # Second plot
        sp(rs = rs, sub = subt, font.sub = 3)
        # Third plot
        tp(fit = fit, rs = rs, sub = subt, font.sub = 3)
        # Fourth plot
        fop(leverage = leverage, rs = rs, r.hat = r.hat, npar = x$npar, sub = subt, font.sub = 3)
        
      }
    }
  }
}
