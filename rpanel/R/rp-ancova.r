#   ANCOVA

rp.ancova <- function(x, y, group, panel = TRUE, panel.plot = TRUE,
                      model = NA, model0 = NA, xlab, ylab, glab,
                      hscale = NA, vscale = hscale, style = "new") {
                    
   if(missing(x) || missing(y) || missing(group)){
      stop("rp.ancova requires x, y, and group.")
      }

   xterm <- deparse(substitute(x))
   yterm <- deparse(substitute(y))
   zterm <- deparse(substitute(group))
   if (missing(xlab)) xlab <- xterm
   if (missing(ylab)) ylab <- yterm
   if (missing(glab)) glab <- zterm
      
   group <- factor(group)
   ind   <- !is.na(x + y + as.numeric(group))
   x     <- x[ind]
   y     <- y[ind]
   group <- group[ind]
         
   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale
      
   if (style == "new") rp.ancova.new(x, y, group, panel = panel, panel.plot = panel.plot,
                          model = model, model0 = model0, xlab = xlab, ylab = ylab, glab = glab,
                          xterm = xterm, yterm = yterm, zterm = zterm,
                          hscale = hscale, vscale = vscale)
   else                rp.ancova.old(x, y, group, panel = panel, panel.plot = panel.plot,
                          model = "None", xlab = xlab, ylab = ylab,
                          hscale = hscale, vscale = vscale)
   invisible()
}
                  
rp.ancova.new <- function(x, y, group, panel = TRUE, panel.plot = TRUE,
                          model = NA, model0 = NA,
                          xlab = deparse(substitute(x)), ylab = deparse(substitute(y)),
                          glab = deparse(substitute(group)),
                          xterm = xterm, yterm = yterm, zterm = zterm,
                          hscale = NA, vscale = hscale) {

   rp.ancova.draw <- function(panel) {
   	
      panel$model  <- c(panel$model11, panel$model12, panel$model13, panel$model14)
      panel$model0 <- c(panel$model01, panel$model02, panel$model03, panel$model04)
   	  panel$model.check  <- any(panel$model)
   	  panel$model0.check <- any(panel$model0)
      if (!panel$model[1] & any(panel$model[-1])) {
         rp.messagebox("The overall mean must be included if other terms are present.")
         panel$model[1] <- TRUE
         panel$model.check <- FALSE
      }
      if (panel$model[4] & !all(panel$model[2:3])) {
         rp.messagebox("The main effects must be included if the interaction term is present.")
         panel$model.check <- FALSE
      }
      if (!panel$model0[1] & any(panel$model0[-1])) {
         rp.messagebox("The overall mean must be included if other terms are present",
                       "in the new model.")
         panel$model0.check <- FALSE
      }
      if (any(panel$model) & panel$model0[4] & !all(panel$model0[2:3])) {
         rp.messagebox("The main effects must be included if the interaction term",
                       "is present in the new model.")
         panel$model0.check <- FALSE
      }

      x <- panel$x
      y <- panel$y
      z <- panel$z
      
      form <- "y ~ 1"
      trms <- panel$term.names[panel$model[-1]]
      if (length(trms) > 0)
         for (i in 1:length(trms))
            form <- paste(form, trms[i], sep = " + ")
      mdl1        <- lm(as.formula(form), na.action = na.exclude)
      panel$mdl1  <- mdl1
      panel$df1   <- mdl1$df.residual
      panel$sigma <- summary(mdl1)$sigma
            
      form <- "y ~ 1"
      trms <- panel$term.names[panel$model0[-1]]
      if (length(trms) > 0)
         for (i in 1:length(trms))
            form <- paste(form, trms[i], sep = " + ")
      mdl0       <- lm(as.formula(form), na.action = na.exclude)
      panel$mdl0 <- mdl0
      panel$df0  <- mdl0$df.residual
      
      rss1        <- sum(mdl1$residuals^2)
      rss0        <- sum(mdl0$residuals^2)
      panel$fstat <- (abs(rss0 - rss1) / abs(panel$df1 - panel$df0)) / panel$sigma^2
                  
      with(panel, {
         n.groups <- length(levels(z))
         par(mar = c(5, 4, 2, 2) + 0.1)
         plot(x, y, type = "n", xlab = xlab, ylab = ylab)
         for (i in 1:n.groups)
            points(x[z == levels(z)[i]],
                   y[z == levels(z)[i]], col = i, pch = i)
         ind   <- (!is.na(x) & !is.na(y) & !is.na(z))
         x <- x[ind]
         y <- y[ind]
         z <- z[ind]
         if (model.check) {
            if (all(model == c(TRUE, FALSE, FALSE, FALSE)))
               abline(h = coef(mdl1))
            else if (all(model == c(TRUE, TRUE, FALSE, FALSE)))
               abline(coef(mdl1))
            else if (any(model)) {
               for (i in 1:n.groups) {
                  ind  <- (z == levels(z)[i])
                  xgp  <- x[ind]
                  fgp  <- fitted(mdl1)[ind]
                  ind1 <- order(xgp)
                  lines(xgp[ind1[range(ind1)]], fgp[ind1[range(ind1)]], col = i, lty = i, lwd = 2)
               }
            }
         }
      })
      panel
   }
      
   rp.ancova.fplot <- function(panel) {
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
            # denstrip(grd, df(grd, abs(df0 - df1), min(df0, df1)), 0.5, 0.9, colmax = "black")
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
      
   rp.ancova.redraw <- function(panel) {
      rp.tkrreplot(panel, plot)
      rp.tkrreplot(panel, fplot)
      panel
   }

   if (panel.plot & !require(tkrplot)) {
      warning("the tkrplot package is not available so panel.plot has been set to FALSE.")
      panel.plot <- FALSE
      }
   	
   model.options <- c("overall mean", xterm, zterm, paste(xterm, ":", zterm))
   term.names    <- c("x", "z", "x:z")
   init.model    <- model
   init.model0   <- model0
   if (any(is.na(init.model)))  init.model  <- rep(FALSE, 4)
   if (any(is.na(init.model0))) init.model0 <- rep(FALSE, 4)
   names(init.model) <- model.options
   bgdcol <- "grey85"

   if (panel) {
      panel <- rp.control("Analysis of covariance", 
                    x = x, y = y, z = factor(group), xlab = xlab, ylab = ylab,
                    xterm = xterm, zterm = zterm, term.names = term.names, 
                    model11 = init.model[1],  model12 = init.model[2],  model13 = init.model[3], 
                    model14 = init.model[4],  model01 = init.model0[1], model02 = init.model0[2],
                    model03 = init.model0[3], model04 = init.model0[4],
                    model.check = TRUE, model0.check = TRUE, bgdcol = bgdcol)
      
      rp.grid(panel, "controls", row = 1, column = 0, background = bgdcol)
      rp.grid(panel, "models",   grid = "controls", row = 0, column = 0, background = bgdcol)
      rp.grid(panel, "fplot",    grid = "controls", row = 0, column = 1, background = bgdcol)

      if (panel.plot) {
         rp.grid(panel, "dataplot", row = 0, column = 0, background = "white")
         rp.tkrplot(panel, plot,  rp.ancova.draw,  hscale = hscale, vscale = vscale, 
                   grid = "dataplot", row = 0, column = 0, background = "white")
      	 rp.text(panel, "", grid = "fplot", row = 0, column = 0, background = bgdcol)
      	 rp.text(panel, "", grid = "fplot", row = 1, column = 0, background = bgdcol)
         rp.tkrplot(panel, fplot, rp.ancova.fplot, hscale = hscale * 0.7, vscale = vscale * 0.2, 
                   grid = "fplot", row = 2, column = 0, background = bgdcol)
      	 rp.text(panel, "", grid = "fplot", row = 3, column = 0, background = bgdcol)
         action.fn <- rp.ancova.redraw
      }
      else
         action.fn <- rp.ancova.draw

      rp.text(panel, "        Model", grid = "models", row = 0, column = 1, background = bgdcol)
      rp.text(panel,       "current", grid = "models", row = 1, column = 0, background = bgdcol)
      rp.text(panel,           "new", grid = "models", row = 1, column = 2, background = bgdcol)
      rp.checkbox(panel, model11, action.fn, labels = "", initval = init.model[1],
            grid = "models", row = 2, column = 0, name = "model11", background = bgdcol)
      rp.checkbox(panel, model12, action.fn, labels = "", initval = init.model[2],
            grid = "models", row = 3, column = 0, name = "model12", background = bgdcol)
      for (i in 1:length(model.options)) rp.text(panel, model.options[i], 
            grid = "models", row = i + 1, column = 1, background = bgdcol)
      rp.checkbox(panel, model01, action.fn, labels = "", initval = init.model0[1],
            grid = "models", row = 2, column = 2, name = "model01", background = bgdcol)
      rp.checkbox(panel, model02, action.fn, labels = "", initval = init.model0[1],
            grid = "models", row = 3, column = 2, name = "model02", background = bgdcol)

      # rp.checkbox(panel, model, action.fn, rep("", length(model.options)), title = "current", 
      #      initval = init.model, grid = "controls", row = 1, column = 0, background = bgdcol)
      rp.checkbox(panel, model13, action.fn, labels = "", initval = init.model[3],
            grid = "models", row = 4, column = 0, name = "model13", background = bgdcol)
      rp.checkbox(panel, model14, action.fn, labels = "", initval = init.model[4],
            grid = "models", row = 5, column = 0, name = "model14", background = bgdcol)
      # rp.grid(panel, "models", row = 1, column = 1, grid = "controls")
      rp.checkbox(panel, model03, action.fn, labels = "", initval = init.model0[1],
            grid = "models", row = 4, column = 2, name = "model03", background = bgdcol)
      rp.checkbox(panel, model04, action.fn, labels = "", initval = init.model0[1],
            grid = "models", row = 5, column = 2, name = "model04", background = bgdcol)
      rp.do(panel, action.fn)
   }
   else {
      panel <- list(x = x, y = y, z = factor(group), xlab = xlab, ylab = ylab,
                    xterm = xterm, zterm = zterm, term.names = term.names, 
                    model11 = init.model[1],  model12 = init.model[2],  model13 = init.model[3], 
                    model14 = init.model[4],  model01 = init.model0[1], model02 = init.model0[2],
                    model03 = init.model0[3], model04 = init.model0[4],
                    model.check = TRUE, model0.check = TRUE, bgdcol = bgdcol)
      rp.ancova.draw(panel)
   }
      
   invisible()
   
}

#-----------------------------------------------------------------------------

rp.ancova.old <- function(x, y, group, panel = TRUE, panel.plot = TRUE, model = "None",
                          xlab = deparse(substitute(x)), ylab = deparse(substitute(y)),
                          hscale = NA, vscale = hscale) {

if (any(is.na(model))) model <- "None"

rp.ancova.draw <- function(panel) {
   with(panel, {
      group    <- factor(group)
      n.groups <- length(levels(group))
      plot(x, y, type = "n", xlab = xlab, ylab = ylab)
      for (i in 1:n.groups)
         points(x[group == levels(group)[i]],
                y[group == levels(group)[i]], col = i, pch = i)
      ind   <- (!is.na(x) & !is.na(y) & !is.na(group))
      x <- x[ind]
      y <- y[ind]
      group <- group[ind]
      if      (model == "Single mean")     lm.model <- lm(y ~ 1)
      else if (model == "Single line")     lm.model <- lm(y ~ x)
      else if (model == "Parallel lines")  lm.model <- lm(y ~ group + x)
      else if (model == "Different lines") lm.model <- lm(y ~ group * x)
      title.text <- paste("Model:", model)
      if (model == "Single mean")
         abline(h = coef(lm.model))
      else if (model == "Single line")
      abline(coef(lm.model))
      else if (!(model == "None")) {
         if (model == "Parallel lines") {
           pval <- drop1(lm.model, test = "F")[["Pr(>F)"]][2]
           pval <- round(pval, 3)
           title.text <- paste(title.text, "\n", "Test of equal groups:", pval)
         }
         if (model == "Different lines") {
           pval <- drop1(lm.model, test = "F")[["Pr(>F)"]][2]
           pval <- round(pval, 3)
           title.text <- paste(title.text, "\n", "Test of parallelism:", pval)
         }
         for (i in 1:n.groups) {
            ind  <- (group == levels(group)[i])
            xgp  <- x[ind]
            fgp  <- fitted(lm.model)[ind]
            ind1 <- order(xgp)
            lines(xgp[ind1[range(ind1)]], fgp[ind1[range(ind1)]], col = i, lty = i, lwd = 2)
         }
      }
      title(title.text, cex.main = 1)
   })
   panel
}

rp.ancova.redraw <- function(panel) {
   rp.tkrreplot(panel, plot)
   panel
   }

   if (panel.plot & !require(tkrplot)) {
      warning("the tkrplot package is not available so panel.plot has been set to FALSE.")
      panel.plot <- FALSE
      }

   if (panel) {
      panel <- rp.control("One-way ancova", y = y, x = x, group = group,
                                 xlab = xlab, ylab = ylab)
      if (panel.plot) {
         rp.tkrplot(panel, plot, rp.ancova.draw, pos = "right",
                    hscale = hscale, vscale = vscale)
         action.fn <- rp.ancova.redraw
         }
      else
         action.fn <- rp.ancova.draw
      rp.radiogroup(panel, model,
         c("None", "Single mean", "Single line", "Parallel lines", "Different lines"),
         action = action.fn)
      rp.do(panel, action.fn)
      }
   else {
      panel <- list(x = x, y = y, group = group, xlab = xlab, ylab = ylab, model = model)   
      rp.ancova.draw(panel)
      invisible()
      }
   invisible()
}
