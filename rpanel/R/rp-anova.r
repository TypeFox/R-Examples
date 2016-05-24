#     rpanel function for simple analysis of variance

rp.anova <- function(y, x, z, model = NA, model0 = NA,
                     ylab = NA, xlab = NA, zlab = NA,
                     panel = TRUE, panel.plot = TRUE, hscale = 1.3, vscale = hscale / 1.3) {

   if (!require(lattice)) stop("the lattice package is not available.")
   if (!require(denstrip)) stop("the denstrip package is not available.")

   type <- if (missing(z)) "One-way" else "Two-way"

   if (is.na(ylab)) ylab <- deparse(substitute(y))
   if (is.na(xlab)) xlab <- deparse(substitute(x))
   zlab <- if (type == "One-way") "" else deparse(substitute(z))
   xterm <- xlab
   zterm <- zlab

   if (type == c("One-way")) z <- rep(1, length(x))
   x <- as.factor(x)
   z <- as.factor(z)

   ind <- !is.na(x) & !is.na(y) & !is.na(z) 
   x <- x[ind]
   y <- y[ind]
   z <- z[ind]
   if (length(x) < 5) stop("too few non-missing data.")
 
   graphics <- "strip plot"
   jitter.x <- jitter(as.numeric(x), factor = 0.5, amount = NULL)
   
   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale

   rp.anova.draw <- function(panel) {
   	
      if (panel$type == "Two-way") {
         panel$model  <- c(panel$model11, panel$model12, panel$model13, panel$model14)
         panel$model0 <- c(panel$model01, panel$model02, panel$model03, panel$model04)
      }
      else {
         panel$model  <- c(panel$model11, panel$model12)
         panel$model0 <- c(panel$model01, panel$model02)
      }
   	  panel$model.check  <- any(panel$model)
   	  panel$model0.check <- any(panel$model0)
      if (!panel$model[1] & any(panel$model[-1])) {
         rp.messagebox("The overall mean must be included if other terms are present.")
         panel$model.check <- FALSE
      }
      if (panel$type == "Two-way" & panel$model[4] & !all(panel$model[2:3])) {
         rp.messagebox("The main effects must be included if the interaction term is present.")
         panel$model.check <- FALSE
      }
      if (any(panel$model) & !panel$model0[1] & any(panel$model0[-1])) {
         rp.messagebox("The overall mean must be included if other terms are present",
                       "in the new model.")
         panel$model0.check <- FALSE
      }
      if (any(panel$model) & panel$type == "Two-way" & panel$model0[4] & !all(panel$model0[2:3])) {
         rp.messagebox("The main effects must be included if the interaction term",
                       "is present in the new model.")
         panel$model0.check <- FALSE
      }
      
      form <- "y ~ 1"
      trms <- panel$term.names[panel$model[-1]]
      if (length(trms) > 0)
         for (i in 1:length(trms))
            form <- paste(form, trms[i], sep = " + ")
      mdl <- lm(as.formula(form), na.action = na.exclude)
      panel$df1 <- mdl$df.residual
      panel$sigma <- summary(mdl)$sigma
            
      form <- "y ~ 1"
      trms <- panel$term.names[panel$model0[-1]]
      if (length(trms) > 0)
         for (i in 1:length(trms))
            form <- paste(form, trms[i], sep = " + ")
      mdl0 <- lm(as.formula(form), na.action = na.exclude)
      panel$df0 <- mdl0$df.residual
            
      rss1        <- sum(mdl$residuals^2)
      rss0        <- sum(mdl0$residuals^2)
      panel$fstat <- (abs(rss0 - rss1) / abs(panel$df1 - panel$df0)) / panel$sigma^2

      with(panel, {
      	
         form <- if (type == "Two-way") "y ~ x | z" else "y ~ x"
         form <- as.formula(form)
         ngps <- nrow(unique(data.frame(x, z)))
         if (graphics != "boxplot") {
            print(stripplot(form, groups = x, layout = c(length(levels(z)), 1),
               ylab = ylab, xlab = xlab,
               panel = function(x, y, subscripts, groups) {
       	          # panel.grid(-1, 0)
       	          # if (type != "one sample")
                  #    for (i in 1:length(unique(x)))
                  #       panel.rect(as.numeric(unique(x)[i]) - 0.3, min(y[x == unique(x)[i]]),
                  #            as.numeric(unique(x)[i]) + 0.3, max(y[x == unique(x)[i]]),
                  #            col = "grey90", border = "grey90")
       	          # panel.abline(v = 1:(length(unique(groups)) - 1) + 0.5, 
       	          #              col = "grey", lty = 2)
       	          zgp <- levels(groups[subscripts])
       	          if (!(model0.check & model.check) | all(model == model0)) {
       	             for (k in zgp) {
       	                ind  <- which(groups[subscripts] == k)
       	                panel.denstrip(y[ind], at = match(k, zgp), width = 0.5, colmax = "grey90",
       	                               horiz = FALSE)
       	             }
       	          }
       	          else if (model.check & !all(model == model0)) {
       	          	 for (k in zgp) {
       	          	    ind  <- which(groups[subscripts] == k)
       	          	    n.se <- length(ind)
       	          	    if (n.se > 0) {
       	          	       wt   <- seq(0, 1, length = 101)
       	          	       cl   <- matrix(col2rgb("green"),  nrow = length(wt), ncol = 3, byrow = TRUE)
       	          	       bck1 <- matrix(col2rgb("grey90"), nrow = length(wt), ncol = 3, byrow = TRUE)
       	          	       bck2 <- matrix(col2rgb("white"),  nrow = length(wt), ncol = 3, byrow = TRUE)
       	                   clr  <- wt * cl + (1 - wt) * bck1
       	                   clr2 <- wt * cl + (1 - wt) * bck2
       	                   grd  <- seq(-3, 3, length = 51)
       	                   ind1 <- cut(exp(-0.5 * grd^2), nrow(clr), labels = FALSE)
       	                   clr  <- clr[ind1, ]
       	                   clr2 <- clr2[ind1, ]
       	                   fv0  <- mean(fitted(mdl0)[subscripts][ind])
       	          	       se   <- summary(mdl)$sigma * sqrt(abs(df0 - df1)) / sqrt(n.se * ngps)
       	                   lim1 <- (min(y[ind]) - fv0) / se
       	                   lim2 <- (max(y[ind]) - fv0) / se
       	          	       indb <- (grd < min(lim1, lim2)) | (grd > max(lim1, lim2))
       	          	       clr[indb, ] <- clr2[indb, ]
       	                   clr  <- rgb(clr[ , 1], clr[ , 2], clr[ , 3], maxColorValue = 255)
       	                   del  <- diff(grd)[1] / 2
       	          	       panel.denstrip(fv0 + (grd - del) * se, dnorm(grd), at = match(k, zgp),
       	          	                      width = 0.5, colmax = "green", horiz = FALSE)
       	          	       # panel.rect(which(k == zgp) - 0.15, fv0 + (grd - del) * se,
       	          	       #            which(k == zgp) + 0.15, fv0 + (grd + del) * se,
       	          	       #            col = clr, border = clr)
       	          	    }
       	          	 }
        	      }
        	      if (model.check) {
        	         panel.xyplot(jitter.x[subscripts], y, col = "blue")
       	             fv  <- tapply(fitted(mdl)[subscripts], groups[subscripts], mean)
       	             panel.lines(1:length(zgp), fv, lwd = 2, col = "grey")
       	             panel.segments(1:length(zgp) - 0.4, fv, 1:length(zgp) + 0.4, fv,
       	                            lwd = 2, col = "red")
        	      }
        	      panel.xyplot(jitter.x[subscripts], y, col = "blue")
               }
            ))
         }
         else if (graphics == "boxplot") {
            print(bwplot(form, groups = x, layout = c(length(levels(z)), 1),
               ylab = ylab, xlab = xlab,
               panel = function(x, y, subscripts, groups) {
       	          # panel.grid(-1, 0)
                  panel.bwplot(x, y, col = groups[subscripts], horizontal = FALSE)
       	          if (any(model) & model.check) {
       	             fv <- unique(cbind(groups[subscripts], fitted(mdl)[subscripts]))
       	             ind <- apply(fv, 1, function(x) !any(is.na(x)))
       	             fv <- fv[ind, ]
       	             fv <- fv[order(fv[ , 1]), ]
       	             panel.lines(fv[ , 1], fv[ , 2], lwd = 2, col = "grey")
       	             panel.segments(fv[ , 1] - 0.4, fv[ , 2], fv[ , 1] + 0.4, fv[ , 2],
       	                            lwd = 2, col = "red")
       	          }
               }
            ))
         }
      })
      panel
      }
      
   rp.anova.fplot <- function(panel) {
      with(panel, {
         if (model.check & model0.check & any(model) & any(model0) & !all(model == model0)) {
         	clr  <- paste("grey", round(seq(100, 0, length = 101)), sep = "")
         	pct  <- qf(0.99, abs(df0 - df1), min(df0, df1))
       	    xlim <- max(pct * 1.5, fstat * 1.1)
       	    grd  <- seq(0, xlim, length = 100)
       	    del  <- diff(grd)[1] / 2
       	    grd  <- grd[-1] - del
       	    ind  <- cut(df(grd, abs(df0 - df1), min(df0, df1)), length(clr), labels = FALSE)
       	    par(mar = c(1, 1, 1, 1), oma = rep(0, 4), tcl = -0.2, xaxs = "i", mgp = c(1, 0, 0))
       	    plot(c(0, xlim), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
       	    axis(1, cex.axis = 0.7)
       	    lines(par()$usr[1:2], rep(par()$usr[3], 2))
       	    rect(grd - del, 0.05, grd + del, 1, col = clr[ind], border = clr[ind])
       	    points(fstat, 0.525, col = "red", pch = 16)
       	    title(paste("p-value:", round(1 - pf(fstat, abs(df0 - df1), min(df0, df1)), 3)),
       	        cex.main = 0.8, font.main = 1)
         }
         else {
       	    par(mar = c(0, 0, 0, 0) + 0.1, bg = bgdcol)
            plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
         }
      })
      panel
   }
      
   rp.anova.redraw <- function(panel) {
      rp.tkrreplot(panel, plot)
      rp.tkrreplot(panel, fplot)
      panel
   }

   if (panel.plot & !require(tkrplot)) {
      warning("the tkrplot package is not available so panel.plot has been set to FALSE.")
      panel.plot <- FALSE
      }

      model.options <- c("overall mean")
      if (type == "One-way") {
         model.options <- c(model.options, xterm)
         term.names    <- c("x")
      }
      if (type == "Two-way") {
         model.options <- c(model.options, zterm, xterm, paste(xterm, ":", zterm))
         term.names    <- c("z", "x", "z:x")
      }
      init.model  <- model
      init.model0 <- model0
      if (any(is.na(init.model)))  init.model  <- rep(FALSE, 4)
      if (any(is.na(init.model0))) init.model0 <- rep(FALSE, 4)
      names(init.model) <- model.options
      bgdcol <- "grey85"

   if (panel) {
   	
      panel <- rp.control(paste(type, "anova"), 
                    x = x, y = y, z = z, type = type, xlab = xlab, ylab = ylab,
                    xterm = xterm, zterm = zterm, term.names = term.names, jitter.x = jitter.x, 
                    graphics = "strip plot",
                    model11 = init.model[1],  model12 = init.model[2],  model13 = init.model[3], 
                    model14 = init.model[4],  model01 = init.model0[1], model02 = init.model0[2],
                    model03 = init.model0[3], model04 = init.model0[4],
                    model.check = TRUE, model0.check = TRUE, bgdcol = bgdcol)
      
      rp.grid(panel, "controls", row = 1, column = 0, background = bgdcol)
      rp.grid(panel, "models",   grid = "controls", row = 0, column = 0, background = bgdcol)
      rp.grid(panel, "fplot",    grid = "controls", row = 0, column = 1, background = bgdcol)

      if (panel.plot) {
         rp.grid(panel, "dataplot", row = 0, column = 0, background = "white")
         rp.tkrplot(panel, plot,  rp.anova.draw,  hscale = hscale, vscale = vscale, 
                   grid = "dataplot", row = 0, column = 0, background = "white")
      	 rp.text(panel, "", grid = "fplot", row = 0, column = 0, background = bgdcol)
      	 rp.text(panel, "", grid = "fplot", row = 1, column = 0, background = bgdcol)
         rp.tkrplot(panel, fplot, rp.anova.fplot, hscale = hscale * 0.7, vscale = vscale * 0.2, 
                   grid = "fplot", row = 2, column = 0, background = bgdcol)
      	 rp.text(panel, "", grid = "fplot", row = 3, column = 0, background = bgdcol)
         action.fn <- rp.anova.redraw
      }
      else
         action.fn <- rp.anova.draw

      rp.text(panel, "        Model", grid = "models", row = 0, column = 1, background = bgdcol)
      rp.text(panel,       "current", grid = "models", row = 1, column = 0, background = bgdcol)
      rp.text(panel,           "new", grid = "models", row = 1, column = 2, background = bgdcol)
      rp.checkbox(panel, model11, action.fn, labels = "", initval = init.model[1],
            grid = "models", row = 2, column = 0, background = bgdcol)
      rp.checkbox(panel, model12, action.fn, labels = "", initval = init.model[2],
            grid = "models", row = 3, column = 0, background = bgdcol)
      for (i in 1:length(model.options)) rp.text(panel, model.options[i], 
            grid = "models", row = i + 1, column = 1, background = bgdcol)
      rp.checkbox(panel, model01, action.fn, labels = "", initval = init.model0[1],
            grid = "models", row = 2, column = 2, background = bgdcol)
      rp.checkbox(panel, model02, action.fn, labels = "", initval = init.model0[2],
            grid = "models", row = 3, column = 2, background = bgdcol)

      if (type == "Two-way") {
         rp.checkbox(panel, model13, action.fn, labels = "", initval = init.model[3],
               grid = "models", row = 4, column = 0, background = bgdcol)
         rp.checkbox(panel, model14, action.fn, labels = "", initval = init.model[4],
               grid = "models", row = 5, column = 0, background = bgdcol)
         rp.checkbox(panel, model03, action.fn, labels = "", initval = init.model0[3],
               grid = "models", row = 4, column = 2, background = bgdcol)
         rp.checkbox(panel, model04, action.fn, labels = "", initval = init.model0[4],
               grid = "models", row = 5, column = 2, background = bgdcol)
      }
      rp.do(panel, action.fn)
   }
   else {
      panel <- list(x = x, y = y, z = z, type = type, xlab = xlab, ylab = ylab,
                    xterm = xterm, zterm = zterm, term.names = term.names, jitter.x = jitter.x, 
                    graphics = "strip plot",
                    model11 = init.model[1],  model12 = init.model[2],  model13 = init.model[3], 
                    model14 = init.model[4],  model01 = init.model0[1], model02 = init.model0[2],
                    model03 = init.model0[3], model04 = init.model0[4],
                    model.check = TRUE, model0.check = TRUE, bgdcol = bgdcol)
      rp.anova.draw(panel)
   }
      
   invisible()
   
}
