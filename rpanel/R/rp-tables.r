#   Statistical tables

rp.tables <- function(panel = TRUE, panel.plot = TRUE, hscale = NA, vscale = hscale,
                      distribution = "normal", degf1 = 5, degf2 = 30,
                      observed.value = " ", observed.value.showing = !is.na(observed.value),
                      probability = 0.05, tail.probability, tail.direction) {
	
   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale
      
   xobs.prob <- c(observed.value, probability)
   if (missing(tail.probability) | is.na(as.numeric(observed.value))) tail.probability <- "none"
   tail.area <- tail.probability   
   if (missing(tail.direction))
   	 tail.direction <- if (distribution %in% c("normal", "t")) "two-sided" else "upper"

tables.draw <- function(tables) {
  with(tables, {
    xobs <- xobs.prob[1]
    xobs <- if (!is.na(xobs)) as.numeric(xobs) else NA
    extend <- if (is.na(xobs)) 0 else xobs * 1.1
    prob <- as.numeric(xobs.prob[2])
    ngrid <- 100
    if (distribution == "normal") {
      xrange <- c(-4, 4)
      extend <- if (is.na(xobs)) 0 else xobs * 1.1
      x      <- seq(min(xrange[1], extend), max(xrange[2], extend), length = ngrid)
      dens   <- dnorm(x)
      ylim   <- c(0, 0.4)
      pval   <- pnorm(xobs)
      pshade <- min(pval, 1 - pval)
      qts    <- qnorm(c(prob, 1 - prob, prob/2, 1 - prob/2, pshade, 1 - pshade))
    }
    if (distribution == "t") {
      xrange <- c(-4, 4)
      x      <- seq(min(xrange[1], extend), max(xrange[2], extend), length = ngrid)
      dens   <- dt(x, degf1)
      ylim   <- c(0, 0.4)
      pval   <- pt(xobs, degf1)
      pshade <- min(pval, 1 - pval)
      qts    <- qt(c(prob, 1 - prob, prob/2, 1 - prob/2, pshade, 1 - pshade), degf1)
#     xp   <- if (two.sided) qt(c(1 - prob/2, prob/2), degf)     else qt(1 - prob, degf)
#     pval <- if (two.sided) 2 * (1 - pt(abs(xobs), degf)) else 1 - pt(xobs, degf)
    }
    if (distribution == "chi-squared") {
      xr1    <- max(degf1 + 3 * sqrt(2 * degf1), qchisq(1 - prob, degf1) + 1)
      xrange <- c(0.01, xr1)
      x      <- seq(xrange[1], max(xrange[2], extend), length = ngrid)
      dens   <- dchisq(x, degf1)
      if (degf1 <= 2) ylim   <- c(0, max(dens))
        else          ylim   <- c(0, 1.1 * max(dens))
      pval   <- pchisq(xobs, degf1)
      pshade <- min(pval, 1 - pval)
      qts    <- qchisq(c(prob, 1 - prob, prob/2, 1 - prob/2, pshade, 1 - pshade), degf1)
    }
    if (distribution == "F") {
      xrange <- c(0.01, max(10, qf(1 - prob, degf1, degf2) + 1))
      x      <- seq(xrange[1], max(xrange[2], extend), length = ngrid)
      dens   <- df(x, degf1, degf2)
      if (degf1 <= 2) ylim   <- c(0, max(dens))
        else          ylim   <- c(0, 1.1 * max(dens))
      pval   <- pf(xobs, degf1, degf2)
      pshade <- min(pval, 1 - pval)
      qts    <- qf(c(prob, 1 - prob, prob/2, 1 - prob/2, pshade, 1 - pshade), degf1, degf2)
    }
    plot(x, dens, type = "l", ylim = ylim, ylab = paste(distribution, "density"))
    abline(h = 0, lty = 3)
    if (distribution == "normal")       title.text <- "Normal distribution"
    if (distribution == "t")            title.text <- paste("t(", degf1,") distribution", sep = "")
    if (distribution == "chi-squared")  title.text <- paste("Chi-squared(", degf1,
                                                            ") distribution", sep = "")
    if (distribution == "F")            title.text <- paste("F(", degf1,",", degf2,
                                                            ") distribution", sep = "")
    xobs <- as.numeric(xobs)
    if (observed.value.showing & !is.na(xobs)) {
      lines(rep(xobs, 2), c(0, ylim[2] * 0.95), lty = 2)
      text(xobs, ylim[2] * 0.97, signif(xobs, 4))
      title.text <- paste(title.text, ";  xobs =", round(xobs, 2))
    }

    if (tail.area != "none") {
      if (tail.area == "fixed probability") {
        if      (tail.direction == "lower") shade <- c(qts[1], NA)
        else if (tail.direction == "upper") shade <- c(NA, qts[2])
        else                                shade <- qts[3:4]
        title.text <- paste(title.text, ";  p =", round(prob, 3))
      }
      else {
        if      (tail.direction == "lower") shade <- c(xobs, NA)
        else if (tail.direction == "upper") {
          shade <- c(NA, xobs)
          pval  <- 1 - pval
        }
        else {
          shade <- qts[5:6]
          pval  <- 2 * min(pval, 1 - pval)
        }
        title.text <- paste(title.text, ";  pval =", round(pval, 3))
      }
      if (!is.na(shade[1])) {
        ind  <- max((1:ngrid)[x < shade[1]])
        intp <- (shade[1] - x[ind]) / (x[ind + 1] - x[ind])
        dend <- (1 - intp) * dens[ind] + intp * dens[ind + 1]
        dt   <- c(dens[x < shade[1]], dend)
        xt   <- c(x[x < shade[1]], shade[1])
        polygon(c(xt, rev(xt)), c(dt, rep(0, length(xt))),
          density = -1, col = "red", border = "red")
      }
      if (!is.na(shade[2])) {
        ind  <- min((1:ngrid)[x > shade[2]])
        intp <- (shade[2] - x[ind - 1]) / (x[ind] - x[ind - 1])
        dend <- (1 - intp) * dens[ind - 1] + intp * dens[ind]
        dt <- c(dend, dens[x > shade[2]])
        xt <- c(shade[2], x[x > shade[2]])
        polygon(c(xt, rev(xt)), c(dt, rep(0, length(xt))),
          density = -1, col = "red", border = "red")
      }
    }
    title(title.text, cex.main = 1)
  })
  tables
  }

tables.redraw <- function(object) {
  rp.tkrreplot(object, plot)
  object
  }

  if (panel) {
  tables.panel <- rp.control("Distributions",
                  xobs.prob = xobs.prob, distribution = distribution, degf1 = degf1, degf2 = degf2,
                  observed.value.showing = observed.value.showing,
                  tail.area = tail.area, tail.direction = tail.direction)
  rp.grid(tables.panel, "controls", row = 0, column = 0)
  rp.grid(tables.panel, "plot",     row = 1, column = 0, background = "white")
  if (panel.plot && require(tkrplot)) {
    rp.tkrplot(tables.panel, plot, plotfun = tables.draw,
                    hscale = hscale, vscale = vscale, row = 1, column = 0, columnspan = 4,
                    grid = "plot", sticky = "ew", background = "white")
    action.fn <- tables.redraw
    }
  else
    action.fn <- tables.draw
  rp.radiogroup(tables.panel, distribution, c("normal", "t", "chi-squared", "F"),
                  title = "Distribution", action = action.fn,
                  grid = "controls", row = 0, column = 0, sticky = "ns")
  rp.grid(tables.panel, "dfgrid", grid = "controls", row = 0, column = 1, sticky = "ns")
  rp.doublebutton(tables.panel, degf1, 1, range = c(1, NA),
                  title = "df1", action = action.fn, grid = "dfgrid", row = 0, column = 0)
  rp.doublebutton(tables.panel, degf2, 1, range = c(1, NA),
                  title = "df2", action = action.fn, grid = "dfgrid", row = 0, column = 1)
  rp.checkbox(tables.panel, observed.value.showing,
                  labels = "Show observed value", action = action.fn,
                  grid = "dfgrid", row = 1, column = 0, columnspan = 2)                    
  rp.textentry(tables.panel, xobs.prob, action.fn, 
                  c("Observed value", "Fixed probability, p"), title = "",
                  grid = "dfgrid", row = 2, column = 0, columnspan = 2, sticky = "ns")                    
  rp.radiogroup(tables.panel, tail.area,
                  c("none", "from observed value", "fixed probability"),
                  title = "Tail probability", action = action.fn,
                  grid = "controls", row = 0, column = 2, sticky = "ns")
  rp.radiogroup(tables.panel, tail.direction,
                  c("lower", "upper", "two-sided"),
                  title = "Tail direction", action = action.fn,
                  grid = "controls", row = 0, column = 3, sticky = "ns")
  rp.do(tables.panel, action.fn)
  }
  else {
    panel <- list(xobs.prob = xobs.prob, distribution = distribution,
                  degf1 = degf1, degf2 = degf2,
                  observed.value.showing = observed.value.showing,
                  tail.area = tail.area, tail.direction = tail.direction)
    tables.draw(panel)
  }
  
  invisible()
  }
