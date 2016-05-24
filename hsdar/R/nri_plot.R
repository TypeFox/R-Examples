setMethod("plot", signature(x = "Nri"),
          function(x,
                   coefficient = NULL,
                   predictor   = NULL,
                   xlab        = "Wavelength band 1 (nm)",
                   ylab        = "Wavelength band 2 (nm)",
                   legend      = TRUE,
                   colspace    = "hcl",
                   col         = c(10, 90, 60, 60, 10, 80),
                   digits      = 2,
                   range       = "auto",
                   constraint  = NULL,
                   upperdiag   = FALSE,
                   ...
                  )
{
  predictordefined <- TRUE
  leg_outer <- FALSE
  if (!is.logical(legend))
  {
    leg_outer <- toupper(legend) == "OUTER"
    legend <- TRUE
  }
  
  if (length(x@fwhm)!=length(x@wavelength))
    x@fwhm <- rep.int(x@fwhm[1], length(x@wavelength))
  got_col_space <- FALSE
  if (colspace == "rgb")
  {
    colfunc <- rgb
    getcol <- colorRamp(col)
    got_col_space <- TRUE
  }
  if (colspace == "hcl")
  {    
    makecolfunc <- function(col) 
    {
      if (!any(length(col) == c(6,7)))
        stop("Length of col vector must be 6 (without alpha) or 7")
      if (!(col[1] >= 0 & col[1] <= 360))
        stop("Invalid hue value")
      if (!(col[2] >= 0 & col[2] <= 360 & col[2] >= col[1]))
        stop("Invalid hue value")
      if (!(col[5] >= 0 & col[5] <= 100))
        stop("Invalid luminance value")
      if (!(col[6] >= 0 & col[6] <= 100 & col[6] >= col[5]))
        stop("Invalid luminance value")
      if (length(col) < 7)
      {
        col <- c(col, 1)
      } else {
        if (!(col[7] >= 0 & col[7] <= 1))
          stop("Invalid alpha value")
      }
      createhcl <- function(n, maxColorValue)
      {
        minmaxcol <- get("colvec")
        hcl(h = n * minmaxcol[2] + minmaxcol[1], 
            c = n * minmaxcol[4] + minmaxcol[3], 
            l = n * minmaxcol[6] + minmaxcol[5],
            alpha = minmaxcol[7],
            fixup = TRUE
           )        
      }
      assign("colvec", col, envir = environment(fun = createhcl))
      return(createhcl)
    }
    colfunc <- makecolfunc(col)
    getcol <- function(n)
      return(n)
    got_col_space <- TRUE
  }
  if (!got_col_space)
    stop("Unknown color space")
    
  if (is.character(predictor))
    predictor <- which(dimnames(x[[1]])[[1]]==predictor)
  
  if (length(x@multivariate) == 0)
  {
    coefficient <- as.matrix(new("DistMat3D", values = apply(x@nri, MARGIN = 1, FUN = mean), ncol = length(x@wavelength), nlyr = 1))
  } else {
    if (is.null(coefficient))
      coefficient <- switch(attr(x@multivariate,"function"),
                            t.test = "p.value",
                            glm = "t.value",
                            cor.test = "p.value",
                            lm = "r.squared"
                           )
    if (is.null(predictor))
    {
      predictordefined <- FALSE
      predictor <- switch(attr(x@multivariate,"function"),
                          t.test = 1,
                          glm = 2,
                          lm = 1,
                          cor.test = 1
                         )
    }
    coefficient <- x@multivariate[[which(names(x@multivariate)==coefficient)]]
    coefficient <- as.matrix(coefficient, lyr = predictor)
  }
  if (!upperdiag)
    plot(x@wavelength,x@wavelength,type="n", xlab=xlab, ylab=ylab, ...)
  
  
  if (range[1]=="auto")
  {
    minval <- min(coefficient, na.rm = TRUE)
    maxval <- max(coefficient, na.rm = TRUE)
  } else {
    minval <- range[1]
    maxval <- range[2]
  }
  
  if (!is.null(constraint))
  {
    cons_eval <- try(eval(parse(text = constraint)), silent = TRUE)
    if (inherits(cons_eval, "try-error"))
    {
      cons_eval <- try(
        eval(parse(text = paste(strsplit(constraint, "$", fixed = TRUE)[[1]][
                                length(strsplit(constraint, "$", fixed = TRUE)[[1]])],
                                sep="")), envir = x@multivariate), 
        silent = TRUE)
      if (inherits(cons_eval, "try-error"))
        stop("Could not evaluate constraint string")
    }
    cons_eval <- distMat3D(as.numeric(cons_eval), ncol(coefficient), 
                           length(cons_eval)/sum(lower.tri(coefficient)))
    if (dim(cons_eval)[3] > 1)
    {
      if (predictordefined) 
      {
        cons_eval_mat <- try(as.matrix(cons_eval, lyr = predictor), silent = TRUE)
      } else {
        conscoefficient <- sapply(c("==", "<=", ">=", ">", "<"), function(pa, te)
        {
          pa_apply <- strsplit(te, pa)
          if (length(pa_apply[[1]]) > 1)
          {
            return(pa)
          } else {
            return(NA)
          }
        }, constraint)
        conscoefficient <- conscoefficient[!is.na(conscoefficient)]
        conscoefficient <- conscoefficient[1]
        conscoefficient <- strsplit(constraint, conscoefficient)[[1]]
        if (any(names(x@multivariate) == conscoefficient[1]))
        {
          conscoefficient <- conscoefficient[1]
        } else {
          conscoefficient <- conscoefficient[2]
        }
        conslayer <- switch(conscoefficient,
                            p.value = dim(cons_eval)[3],
                            t.value = dim(cons_eval)[3],
                            z.value = dim(cons_eval)[3], 
                            std.error = dim(cons_eval)[3],
                            estimate = 1,
                            r.squared = 1,
                            1
                           )
        cons_eval_mat <- try(as.matrix(cons_eval, lyr = conslayer), silent = TRUE)
      }
      if (inherits(cons_eval_mat, "try-error"))  
        cons_eval_mat <- as.matrix(cons_eval)
    } else {
      cons_eval_mat <- as.matrix(cons_eval)
    }
    coefficient[cons_eval_mat == 0] <- minval -1    
  }
  if (upperdiag)
  {
    for (k in 1:(dim(coefficient)[1]-1))
    {
      for (i in (k+1):dim(coefficient)[2])
      {
        if (is.finite(coefficient[i,k]))
        {
          if (coefficient[i,k]>=minval & coefficient[i,k]<=maxval)
          {
            polygon(c(x@wavelength[k]-x@fwhm[k]/2,
                      x@wavelength[k]-x@fwhm[k]/2,
                      x@wavelength[k]+x@fwhm[k]/2,
                      x@wavelength[k]+x@fwhm[k]/2),
                    c(x@wavelength[i]+x@fwhm[i]/2,
                      x@wavelength[i]-x@fwhm[i]/2,
                      x@wavelength[i]-x@fwhm[i]/2,
                      x@wavelength[i]+x@fwhm[i]/2),
                    lty="blank",
                    col=colfunc(getcol((coefficient[i,k]-minval)/(maxval-minval)),
                            maxColorValue = 255)
                    )
          }
        }
      } 
    }
  } else {
    for (k in 1:(dim(coefficient)[1]-1))
    {
      for (i in (k+1):dim(coefficient)[2])
      {
        if (is.finite(coefficient[i,k]))
        {
          if (coefficient[i,k]>=minval & coefficient[i,k]<=maxval)
          {
            polygon(c(x@wavelength[i]+x@fwhm[i]/2,
                      x@wavelength[i]-x@fwhm[i]/2,
                      x@wavelength[i]-x@fwhm[i]/2,
                      x@wavelength[i]+x@fwhm[i]/2),
                    c(x@wavelength[k]-x@fwhm[k]/2,
                      x@wavelength[k]-x@fwhm[k]/2,
                      x@wavelength[k]+x@fwhm[k]/2,
                      x@wavelength[k]+x@fwhm[k]/2),
                    lty="blank",
                    col=colfunc(getcol((coefficient[i,k]-minval)/(maxval-minval)),
                            maxColorValue = 255)
                    )
          }
        }
      } 
    }
  }
  
  if (legend)
  {
    if (!upperdiag)
    {
      if (!leg_outer)
      {
        p <- par()$usr
        xcoor <- c(p[1]+0.05*(p[2]-p[1]),p[1]+0.35*(p[2]-p[1]))
        ycoor <- c(p[4]-0.085*(p[4]-p[3]),p[4]-0.05*(p[4]-p[3]))
        for (i in 1:100)
        {
          xrect <- c(xcoor[1]+(i-1)*(xcoor[2]-xcoor[1])/100,
                     xcoor[1]+i*(xcoor[2]-xcoor[1])/100)
          rect(xrect[1],ycoor[1],xrect[2],ycoor[2],
              col = colfunc(getcol(i/100), maxColorValue = 255),
              lty = "blank")
        }
        ycoor <- c(p[4]-0.12*(p[4]-p[3]),p[4]-0.05*(p[4]-p[3]))
        text(xcoor[1],ycoor[1],round(minval,digits))
        text(xcoor[2],ycoor[1],round(maxval,digits))
      } else {
        old.par <- par(no.readonly = TRUE)
        x_leg <- c(.1,.3,0)
        p <- list(fig = par()$fig, fin = par()$fin, mai = par()$mai)
        par(fig = c(p$fig[2] - ((p$fig[2]-p$fig[1])/p$fin[1]*p$mai[4]),
                    p$fig[2],
                    p$fig[3] + ((p$fig[4]-p$fig[3])/p$fin[2]*p$mai[1]), 
                    p$fig[4] - ((p$fig[4]-p$fig[3])/p$fin[2]*p$mai[3])), 
            new = TRUE, mai = c(0,0,0,0))
        plot(c(0:1),c(0:1), type = "n", xaxt = "n", yaxt = "n", bty = "n")
        xcoor <- c(x_leg[1], x_leg[1] + x_leg[2])
        ycoor <- c(0, 1)
        for (i in 1:100)
        {
          yrect <- c(ycoor[1]+(i-1)*(ycoor[2]-ycoor[1])/100,
                     ycoor[1]+i*(ycoor[2]-ycoor[1])/100)
          rect(xcoor[1],yrect[1],xcoor[2],yrect[2],
              col = colfunc(getcol(i/100), maxColorValue = 255),
              lty = "blank")
        }
        xcoor <- sum(x_leg)
        text(xcoor[1],ycoor[1],round(minval,digits), pos = 4)
        text(xcoor[1],ycoor[2],round(maxval,digits), pos = 4)
        par(old.par)
        par(mfg=old.par$mfg)   
      }
    } else {
      old.par <- par(no.readonly = TRUE)
      y_leg <- c(.1,.3,0)
      p <- list(fig = par()$fig, fin = par()$fin, mai = par()$mai)
      par(fig = c(p$fig[1] + ((p$fig[2]-p$fig[1])/p$fin[1]*p$mai[2]),
                  p$fig[2] - ((p$fig[2]-p$fig[1])/p$fin[1]*p$mai[4]),
                  p$fig[4] - ((p$fig[4]-p$fig[3])/p$fin[2]*p$mai[3]), 
                  p$fig[4]),
          new = TRUE, mai = c(0,0,0,0))
      plot(c(0:1),c(0:1), type = "n", xaxt = "n", yaxt = "n", bty = "n")
      xcoor <- c(0, 1)
      ycoor <- c(y_leg[1], y_leg[1] + y_leg[2])
      for (i in 1:100)
      {
        xrect <- c(xcoor[1]+(i-1)*(xcoor[2]-xcoor[1])/100,
                    xcoor[1]+i*(xcoor[2]-xcoor[1])/100)
        rect(xrect[1],ycoor[1],xrect[2],ycoor[2],
            col = colfunc(getcol(i/100), maxColorValue = 255),
            lty = "blank")
      }
      ycoor <- sum(y_leg)
      text(xcoor[1],ycoor[1],round(minval,digits), pos = 3)
      text(xcoor[2],ycoor[1],round(maxval,digits), pos = 3)
      par(old.par)
      par(mfg=old.par$mfg)   
    }
  }
  invisible(c(minval,maxval))
}
)
