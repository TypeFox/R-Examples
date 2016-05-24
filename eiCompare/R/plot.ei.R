plot.ei <- function (x, ...) 
{

  # ei package hidden functions
  .tomog <- function (ei.object, title = "Tomography Plot with the Data", 
                      lci = T) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    bounds <- bounds1(x, t, n)
    bbounds <- cbind(bounds[, 1], bounds[, 2])
    wbounds <- cbind(bounds[, 4], bounds[, 3])
    n <- dim(bounds)[1]
    length <- NA
    for (i in 1:n) {
      length[i] <- sqrt(abs(bbounds[i, 1] - bbounds[i, 2])^2 + 
                          abs(wbounds[i, 1] - wbounds[i, 2])^2)
    }
    scale <- ((length - min(length))/(max(length) - min(length))) * 
      100
    plot(c(100, 200), xlim = c(0, 1), ylim = c(0, 1), col = "white", 
         ylab = "betaW", xlab = "betaB", xaxs = "i", yaxs = "i", 
         main = title)
    if (lci == T) {
      for (i in 1:n) {
        lines(bbounds[i, ], wbounds[i, ], col = hcl(h = 30, 
                                                    c = 100, l = scale[i], alpha = 1))
      }
    }
    else {
      for (i in 1:n) {
        lines(bbounds[i, ], wbounds[i, ], col = "yellow")
      }
    }
  }
  
  # bounds1
  bounds1 <- function (x, t, n) 
  {
    homindx <- NULL
    tx <- NULL
    tomx <- NULL
    LbetaB <- NULL
    UbetaB <- NULL
    LbetaW <- NULL
    UbetaW <- NULL
    omx = 1 - x
    Nb = x * n
    Nw = omx * n
    p = length(x)
    homoindx <- ifelse(x == 0, 1, 0)
    homoindx <- ifelse(x == 1, 2, homoindx)
    tx <- as.matrix(t/x)
    tomx = as.matrix(t/omx)
    tomxx <- as.matrix(tx - (omx/x))
    txx <- as.matrix(tomx - x/(1 - x))
    LbetaB <- apply(tomxx, 1, function(x) max(0, x))
    UbetaB <- apply(tx, 1, function(x) min(x, 1))
    LbetaW <- apply(txx, 1, function(x) max(0, x))
    UbetaW <- apply(tomx, 1, function(x) min(x, 1))
    bl <- homoindx == 2
    LbetaB[bl] = t[bl]
    UbetaB[bl] = t[bl]
    LbetaW[bl] = NA
    UbetaW[bl] = NA
    wh <- homoindx == 1
    LbetaB[wh] = NA
    UbetaB[wh] = NA
    LbetaW[wh] = t[wh]
    UbetaW[wh] = t[wh]
    return(cbind(LbetaB, UbetaB, LbetaW, UbetaW))
  }
  
  .tomogd <- function (x, t, n, title, lci = T) {
    bounds <- bounds1(x, t, n)
    bbounds <- cbind(bounds[, 1], bounds[, 2])
    wbounds <- cbind(bounds[, 4], bounds[, 3])
    n <- dim(bounds)[1]
    length <- NA
    for (i in 1:n) {
        length[i] <- sqrt(abs(bbounds[i, 1] - bbounds[i, 2])^2 + 
            abs(wbounds[i, 1] - wbounds[i, 2])^2)
    }
    scale <- ((length - min(length))/(max(length) - min(length))) * 
        100
    plot(c(100, 200), xlim = c(0, 1), ylim = c(0, 1), col = "white", 
        ylab = "betaW", xlab = "betaB", xaxs = "i", yaxs = "i", 
        main = title)
    if (lci == T) {
        for (i in 1:n) {
            lines(bbounds[i, ], wbounds[i, ], col = hcl(h = 30, 
                c = 100, l = scale[i], alpha = 1))
        }
    }
    else {
        for (i in 1:n) {
            lines(bbounds[i, ], wbounds[i, ], col = "yellow")
        }
    }
   }
  
  .repar <- function (Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw) 
  {
    sb = exp(sb0)
    sw = exp(sw0)
    bb = Bb0 * (0.25 + sb^2) + 0.5 + as.matrix(apply(Zb, 2, function(x) x - 
                                                       mean(x))) %*% as.matrix(Bb0v)
    bw = Bw0 * (0.25 + sw^2) + 0.5 + as.matrix(apply(Zw, 2, function(x) x - 
                                                       mean(x))) %*% as.matrix(Bw0v)
    rho = (exp(2 * rho0) - 1)/(exp(2 * rho0) + 1)
    return(c(t(bb), t(bw), sb, sw, rho))
  }
  
  .tomog3 <- function(bb, bw, sb, sw, rho) {
    lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2), scale = c(sb, 
                                                                 sw), centre = c(mean(bb), mean(bw)), level = 0.914), 
          col = "blue", lwd = 4)
    lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2), scale = c(sb, 
                                                                 sw), centre = c(mean(bb), mean(bw)), level = 0.35), 
          col = "red", lwd = 4)
    points(mean(bb), mean(bw), col = "pink", pch = 15)
  }
  
  .tomogl <- function (ei.object, lci = TRUE) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    Zb <- ei.object$Zb
    Zw <- ei.object$Zw
    phi <- ei.object$phi
    .tomogd(x, t, n, "Tomography Plot with ML Contours", lci = lci)
    numb <- dim(Zb)[2]
    numw <- dim(Zw)[2]
    Bb0 <- phi[1]
    Bw0 <- phi[2]
    sb0 <- phi[3]
    sw0 <- phi[4]
    rho0 <- phi[5]
    Bb0v <- phi[6:(5 + numb)]
    Bw0v <- phi[(6 + numb):length(phi)]
    vars <- .repar(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, 
                   Zw)
    bb <- vars[1:length(x)]
    bw <- vars[(length(x) + 1):(2 * length(x))]
    sb <- vars[2 * length(x) + 1]
    sw <- vars[2 * length(x) + 2]
    rho <- vars[2 * length(x) + 3]
    
    .tomog3(bb, bw, sb, sw, rho)
  }
  
  .aggs <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betab <- ei.object$betabs[ok, ]
      betaw <- ei.object$betaws[ok, ]
      omx <- 1 - x
      Nb <- n * x
      Nw <- n * omx
      Bbgg <- vector(mode = "numeric", length = dim(betab)[2])
      for (i in 1:dim(betab)[2]) {
        Bbgg[i] <- weighted.mean(betab[, i], Nb)
      }
      Bwgg <- vector(mode = "numeric", length = dim(betaw)[2])
      for (i in 1:dim(betaw)[2]) {
        Bwgg[i] <- weighted.mean(betaw[, i], Nw)
      }
      return(cbind(Bbgg, Bwgg))
    }
  }
  
  .tomog80CI <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betabs <- ei.object$betabs[ok, ]
      betaws <- ei.object$betaws[ok, ]
      .tomogd(x, t, n, "Tomography Plot with 80% CIs", lci = F)
      betabcd <- apply(betabs, 1, function(x) quantile(x, probs = c(0.1, 
                                                                    0.9)))
      betawcd <- apply(betaws, 1, function(x) quantile(x, probs = c(0.1, 
                                                                    0.9)))
      n <- dim(betabcd)[2]
      for (i in 1:n) {
        lines(betabcd[, i], sort(betawcd[, i], decreasing = T), 
              col = "red", lwd = 3)
      }
    }
  }

  .tomog95CI <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betabs <- ei.object$betabs[ok, ]
      betaws <- ei.object$betaws[ok, ]
      .tomogd(x, t, n, "Tomography Plot with 95% CIs", lci = F)
      betabcd <- apply(betabs, 1, function(x) quantile(x, probs = c(0.025, 
                                                                    0.975)))
      betawcd <- apply(betaws, 1, function(x) quantile(x, probs = c(0.025, 
                                                                    0.975)))
      n <- dim(betabcd)[2]
      for (i in 1:n) {
        lines(betabcd[, i], sort(betawcd[, i], decreasing = T), 
              col = "red", lwd = 3)
      }
    }
  }
    
  .tomogE <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betabs <- ei.object$betabs[ok, ]
      betaws <- ei.object$betaws[ok, ]
      .tomogd(x, t, n, "Tomography Plot with Mean Posterior Betabs and\nBetaws")
      betabm <- apply(betabs, 1, mean)
      betawm <- apply(betaws, 1, mean)
      points(betabm, betawm, col = "red", pch = 19)
    }
  }
    
  .tomogP2 <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      x <- ei.object$x
      t <- ei.object$t
      n <- ei.object$n
      psi <- ei.object$psi
      .tomogd(x, t, n, "Tomography Plot with Contours Based on Posterior")
      bbp <- psi[, 1:length(x)]
      bwp <- psi[, (length(x) + 1):(2 * length(x))]
      sbp <- psi[, 2 * length(x) + 1]
      swp <- psi[, 2 * length(x) + 2]
      rhop <- psi[, 2 * length(x) + 3]
      points(mean(bbp), mean(bwp), col = "red", pch = 19)
      .tomog3(bbp, bwp, mean(sbp), mean(swp), mean(rhop))
    }
  }
    
  .betabd <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab)
      betabs <- ei.object$betabs[ok, ]
      betabm <- apply(betabs, 1, mean)
      plot(density(betabm), xlim = c(0, 1), col = "green", 
           xlab = "betaB", ylab = "density across precincts, f(betaB)", 
           main = "Density of\nbetaB")
      vb <- as.vector(betabm)
      for (i in 1:length(vb)) {
        lines(c(vb[i], vb[i]), c(0, 0.25))
      }
    }
  }
    
  .betawd <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betaw)
      betaws <- ei.object$betaws[ok, ]
      betawm <- apply(betaws, 1, mean)
      plot(density(betawm), xlim = c(0, 1), col = "green", 
           xlab = "betaW", ylab = "density across precincts, f(betaW)", 
           main = "Density of\nbetaW")
      vw <- as.vector(betawm)
      for (i in 1:length(vw)) {
        lines(c(vw[i], vw[i]), c(0, 0.25))
      }
    }
  } 
  
  .xt <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    plot(x, t, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", 
         main = "X and T\nScatterplot", ylab = "T", xlab = "X", 
         pch = 20)
  }
    
  .xtc <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    circ <- 0.04
    plot(x, t, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", 
         main = "X and T\nScatterplot with Population Density Circles", 
         ylab = "T", xlab = "X", pch = 20)
    minn <- min(n)
    maxn <- max(n)
    for (i in 1:length(x)) {
      radius = (n[i] - minn + 1)/(1 + maxn - minn)
      draw.circle(x[i], t[i], radius * circ)
    }
  }
    
  .xtfit <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betabs <- ei.object$betabs[ok, ]
      betaws <- ei.object$betaws[ok, ]
      low <- 0.1
      up <- 0.9
      circ <- 0.04
      plot(x, t, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
           yaxs = "i", main = "X and T\nScatterplot with E(T|X) and 80% CIs", 
           ylab = "T", xlab = "X", pch = 20)
      minn <- min(n)
      maxn <- max(n)
      for (i in 1:length(x)) {
        radius = (n[i] - minn + 1)/(1 + maxn - minn)
        draw.circle(x[i], t[i], radius * circ)
      }
      x <- seq(0, 1, by = 0.01)
      betabs <- as.vector(betabs)
      betaws <- as.vector(betaws)
      t <- matrix(ncol = length(x), nrow = length(betabs))
      for (i in 1:length(x)) {
        t[, i] <- betabs * x[i] + betaws * (1 - x[i])
      }
      et <- apply(t, 2, mean)
      lines(x, et, col = "yellow")
      lwr <- apply(t, 2, function(x) quantile(x, probs = c(low)))
      upr <- apply(t, 2, function(x) quantile(x, probs = c(up)))
      lines(x, lwr, col = "red")
      lines(x, upr, col = "red")
    }
  }
    
  .xtfitg <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betabs <- ei.object$betabs[ok, ]
      betaws <- ei.object$betaws[ok, ]
      low <- 0.1
      up <- 0.9
      circ <- 0.04
      plot(x, t, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
           yaxs = "i", main = "X and T\nScatterplot with E(T|X), 80% CIs, and Goodman", 
           ylab = "T", xlab = "X", pch = 20)
      minn <- min(n)
      maxn <- max(n)
      for (i in 1:length(x)) {
        radius = (n[i] - minn + 1)/(1 + maxn - minn)
        draw.circle(x[i], t[i], radius * circ)
      }
      x <- seq(0, 1, by = 0.01)
      betabs <- as.vector(betabs)
      betaws <- as.vector(betaws)
      t <- matrix(ncol = length(x), nrow = length(betabs))
      for (i in 1:length(x)) {
        t[, i] <- betabs * x[i] + betaws * (1 - x[i])
      }
      et <- apply(t, 2, mean)
      lines(x, et, col = "yellow")
      lwr <- apply(t, 2, function(x) quantile(x, probs = c(low)))
      upr <- apply(t, 2, function(x) quantile(x, probs = c(up)))
      lines(x, lwr, col = "red")
      lines(x, upr, col = "red")
      t <- ei.object$t
      x <- ei.object$x
      lm.fit <- lm(t ~ x)
      abline(lm.fit, col = "green")
    }
  }
    
  .estsims <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This plot function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      betabs <- ei.object$betabs[ok, ]
      betaws <- ei.object$betaws[ok, ]
      colors = runif(length(betabs), 26, 51)
      plot(betabs, betaws, xlim = c(0, 1), ylim = c(0, 1), 
           xaxs = "i", yaxs = "i", main = "Simulations of betaW and betaB", 
           ylab = "betaW\nsimulations", xlab = "betaB simulations", 
           pch = 20, col = colors, lty = 2, cex = 0.25)
    }
  }

  .boundXb <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    truebb <- ei.object$truth[, 1]
    bounds <- bounds1(x, t, n)
    plot(x, truebb, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
         yaxs = "i", main = "Aggregation Bias for betaB", ylab = "True betab", 
         xlab = "X", pch = 20)
    for (i in 1:length(x)) {
      lines(c(x[i], x[i]), c(bounds[, 1][i], bounds[, 2][i]))
    }
    lm.xb <- lm(truebb ~ x)
    abline(lm.xb, lty = 2)
  }
    
  .boundXw <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    truebw <- ei.object$truth[, 2]
    bounds <- bounds1(x, t, n)
    plot(x, truebw, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
         yaxs = "i", main = "Aggregation Bias for betaW", ylab = "True betaw", 
         xlab = "X", pch = 20)
    for (i in 1:length(x)) {
      lines(c(x[i], x[i]), c(bounds[, 3][i], bounds[, 4][i]))
    }
    lm.xw <- lm(truebw ~ x)
    abline(lm.xw, lty = 2)
  }
  
  .CI80b <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      betab <- ei.object$betabs[ok, ]
      lwr <- vector(mode = "numeric", length = length(ei.object$x))
      upr <- vector(mode = "numeric", length = length(ei.object$x))
      lwr[ok] <- apply(betab, 1, function(x) quantile(x, probs = c(0.1)))
      lwr[!ok] <- NA
      upr[ok] <- apply(betab, 1, function(x) quantile(x, probs = c(0.9)))
      upr[!ok] <- NA
      return(cbind(lwr, upr))
    }
  }
  
  .CI80w <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      betaw <- ei.object$betaws[ok, ]
      lwr <- vector(mode = "numeric", length = length(ei.object$x))
      upr <- vector(mode = "numeric", length = length(ei.object$x))
      lwr[ok] <- apply(betaw, 1, function(x) quantile(x, probs = c(0.1)))
      lwr[!ok] <- NA
      upr[ok] <- apply(betaw, 1, function(x) quantile(x, probs = c(0.9)))
      upr[!ok] <- NA
      return(cbind(lwr, upr))
    }
  }
  
  .truthfn <- function (ei.object) 
  {
    n <- ei.object$n
    x <- ei.object$x
    omx <- 1 - x
    truebb <- ei.object$truth[, 1]
    truebw <- ei.object$truth[, 2]
    betabs <- ei.object$betabs
    betaws <- ei.object$betaws
    betab <- ei.object$betab
    betaw <- ei.object$betaw
    truthbb <- sum(truebb * n)/sum(n)
    truthbw <- sum(truebw * n)/sum(n)
    circ = 0.04
    par(mfrow = c(2, 2))
    ag <- .aggs(ei.object)
    plot(density(ag[, 1]), xlim = c(0, 1), ylim = c(0, max(density(ag[, 
                                                                      1])$y) + 1), yaxs = "i", xaxs = "i", main = "Density of Bb Posterior & Truth", 
         xlab = "Bb", ylab = "Density")
    lines(c(truthbb, truthbb), c(0, 0.25 * (max(density(ag[, 
                                                           1])$y) + 1)), lwd = 3)
    plot(density(ag[, 2]), xlim = c(0, 1), ylim = c(0, max(density(ag[, 
                                                                      2])$y) + 1), yaxs = "i", xaxs = "i", main = "Density of Bw Posterior & Truth", 
         xlab = "Bw", ylab = "Density")
    lines(c(truthbw, truthbw), c(0, 0.25 * (max(density(ag[, 
                                                           2])$y) + 1)), lwd = 3)
    plot(betab, truebb, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
         yaxs = "i", xlab = "Estimated betab", cex = 0.1, ylab = "True betab")
    minn <- min(x * n)
    maxn <- max(x * n)
    for (i in 1:length(betab)) {
      radius = (n[i] * x[i] - minn + 1)/(1 + maxn - minn)
      draw.circle(betab[i], truebb[i], radius * circ)
    }
    ci80b = .CI80b(ei.object)
    low = mean(abs(ci80b[, 1] - betab))
    high = mean(abs(ci80b[, 2] - betab))
    abline(0, 1)
    lines(c(0, 1), c(-low, 1 - low), lty = 2)
    lines(c(0, 1), c(high, 1 + high), lty = 2)
    plot(betaw, truebw, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
         yaxs = "i", xlab = "Estimated\nbetaw", ylab = "True betaw", 
         cex = 0.1)
    minn <- min(omx * n)
    maxn <- max(omx * n)
    for (i in 1:length(betaw)) {
      radius = (omx[i] * n[i] - minn + 1)/(1 + maxn - minn)
      draw.circle(betaw[i], truebw[i], radius * circ)
    }
    ci80w = .CI80w(ei.object)
    low = mean(abs(ci80w[, 1] - betaw))
    high = mean(abs(ci80w[, 2] - betaw))
    abline(0, 1)
    lines(c(0, 1), c(-low, 1 - low), lty = 2)
    lines(c(0, 1), c(high, 1 + high), lty = 2)
  }
  
  .bndplot <- function (dbuf) 
  {
    form <- dbuf$formula
    total <- dbuf$total
    data <- dbuf$data
    n <- nrow(data)
    print(n)
    covariate <- NA
    rows <- c(all.names(form)[6:(length(all.names(form)))])
    names = rows
    cols <- c(all.names(form)[3])
    if (sum(data[, rows][, 1] < 1.1) == length(data[, rows][, 
                                                            1])) {
      data <- round(data * data[, total])
    }
    bnds <- bounds(form, data = data, rows = rows, column = cols, 
                   threshold = 0)
    dv <- data[, all.names(form)[3]]
    bndsoth <- bnds$bounds[[3]]
    oth <- data[, all.names(form)[length(all.names(form))]]
    bndsx <- bnds$bounds[[1]]
    xcat <- data[, all.names(form)[6]]
    bndsy <- bnds$bounds[[2]]
    ycat <- data[, all.names(form)[7]]
    minx <- bndsx[, 1]
    miny <- bndsy[, 1]
    minoth <- bndsoth[, 1]
    maxx <- bndsx[, 2]
    maxy <- bndsy[, 2]
    maxoth <- bndsoth[, 2]
    newdv <- dv - (minx * xcat)
    newtot <- oth + ycat
    t <- newdv/newtot
    y <- ycat/newtot
    lby <- cbind(miny, (t - maxoth * oth/newtot)/(y))
    lby[, 2] <- ifelse(y == 0, 0, lby[, 2])
    lowy <- apply(lby, 1, max)
    hby <- cbind((t - minoth * oth/newtot)/y, maxy)
    highy <- apply(hby, 1, min)
    newtot <- oth + xcat
    newdv <- dv - (miny * ycat)
    x <- xcat/newtot
    t <- newdv/newtot
    lbx <- cbind(minx, (t - maxoth * oth/newtot)/x)
    lbx[, 2] <- ifelse(x == 0, 0, lbx[, 2])
    lowx <- apply(lbx, 1, max)
    hbx <- cbind((t - minoth * oth/newtot)/x, maxx)
    highx <- apply(hbx, 1, min)
    hstr <- cbind(minx, highy)
    hend <- cbind(highx, miny)
    lstr <- cbind(minx, lowy)
    lend <- cbind(lowx, miny)
    if (!is.na(covariate)) {
      redg <- data[, covariate]/(oth - data[, covariate])/max(data[, 
                                                                   covariate]/(oth - data[, covariate]))
      blug <- 1 - redg
    }
    if (is.na(covariate)) {
      redg <- rep(0.5, length(minx))
      blug <- rep(0.5, length(minx))
    }
    xl <- paste("Percent", names[1], "Vote Democrat")
    yl <- paste("Percent", names[2], "Vote Democrat")
    mn <- paste("Tomography Plot in a 2x3 Table (", names[3], 
                " Other Category)", sep = "")
    plot(c(0, 0), xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", 
         yaxs = "i", xlab = xl, ylab = yl, col = "white", main = mn)
    ok <- !is.na(hstr[, 2]) & !is.na(hend[, 1])
    exp1 <- hstr[ok, 2] >= maxy[ok]
    exp2 <- hend[ok, 1] >= maxx[ok]
    exp3 <- lstr[ok, 2] <= miny[ok]
    exp4 <- lend[ok, 1] <= minx[ok]
    hstr <- hstr[ok, ]
    hend <- hend[ok, ]
    lstr <- lstr[ok, ]
    lend <- lend[ok, ]
    dv <- dv[ok]
    ycat <- ycat[ok]
    oth <- oth[ok]
    minoth <- minoth[ok]
    xcat <- xcat[ok]
    maxy <- maxy[ok]
    maxx <- maxx[ok]
    redg <- redg[ok]
    blug <- blug[ok]
    for (i in 1:dim(hstr)[1]) {
      if ((exp1[i] + exp2[i] + exp3[i] + exp4[i]) == 0) {
        xaxs <- c(hstr[i, 1], lstr[i, 1], lend[i, 1], hend[i, 
                                                           1])
        yaxs <- c(hstr[i, 2], lstr[i, 2], lend[i, 2], hend[i, 
                                                           2])
        side1 <- c(xaxs, xaxs[1])
        side2 <- c(yaxs, yaxs[1])
        c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
        c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 
                                                        1)])
        area <- abs(c1 - c2)/2
        if (area > 0 & !is.nan(area)) {
          alpha = min(0.5/(area * (n)), 1)
        }
        if (area == 0 | is.nan(area)) {
          alpha = 0.05
        }
        border = alpha
        polygon(xaxs, yaxs, col = rgb(redg[i], 0, blug[i], 
                                      alpha = alpha), border = rgb(redg[i], 0, blug[i], 
                                                                   alpha = 1), lty = 2)
      }
      if ((exp1[i] == 1) & (exp2[i]) == 0) {
        cut <- (dv[i] - (oth[i]) * minoth[i])/xcat[i] - maxy[i] * 
          ycat[i]/xcat[i]
        kink1x <- c(cut)
        kink1y <- c(maxy[i])
      }
      if ((exp2[i] == 1) & (exp1[i]) == 0) {
        cut <- (dv[i] - (oth[i]) * minoth[i])/ycat[i] - maxx[i] * 
          xcat[i]/ycat[i]
        kink1x <- c(maxx[i])
        kink1y <- c(cut)
      }
      if ((exp2[i] == 1 & exp1[i] == 1)) {
        cut <- (dv[i] - (oth[i]) * minoth[i])/ycat[i] - maxx[i] * 
          xcat[i]/ycat[i]
        cut2 <- (dv[i] - (oth[i]) * minoth[i])/xcat[i] - 
          maxy[i] * ycat[i]/xcat[i]
        kink1x <- c(maxx[i], cut2)
        kink1y <- c(cut, maxy[i])
      }
      if ((exp3[i] == 1) & (exp4[i]) == 0) {
        cut <- (dv[i] - (oth[i]) * maxoth[i])/xcat[i] - miny[i] * 
          ycat[i]/xcat[i]
        kink2x <- c(cut)
        kink2y <- c(miny[i])
      }
      if ((exp4[i] == 1) & (exp3[i]) == 0) {
        cut <- (dv[i] - (oth[i]) * maxoth[i])/ycat[i] - minx[i] * 
          xcat[i]/ycat[i]
        kink2x <- c(minx[i])
        kink2y <- c(cut)
      }
      if ((exp3[i] == 1 & exp4[i] == 1)) {
        cut <- (dv[i] - (oth[i]) * maxoth[i])/ycat[i] - minx[i] * 
          xcat[i]/ycat[i]
        cut2 <- (dv[i] - (oth[i]) * maxoth[i])/xcat[i] - 
          miny[i] * ycat[i]/xcat[i]
        kink2x <- c(minx[i], cut2)
        kink2y <- c(cut, miny[i])
      }
      if ((exp3[i] + exp4[i]) == 0 & (exp1[i] + exp2[i] + exp3[i] + 
                                      exp4[i]) != 0) {
        xaxs <- c(hstr[i, 1], lstr[i, 1], lend[i, 1], hend[i, 
                                                           1], kink1x)
        xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
        yaxs <- c(hstr[i, 2], lstr[i, 2], lend[i, 2], hend[i, 
                                                           2], kink1y)
        yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
        side1 <- c(xaxs, xaxs[1])
        side2 <- c(yaxs, yaxs[1])
        c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
        c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 
                                                        1)])
        area <- abs(c1 - c2)/2
        if (area > 0 & !is.nan(area)) {
          alpha = min(0.5/(area * (n)), 1)
        }
        if (area == 0 | is.nan(area)) {
          alpha = 0.05
        }
        border = alpha
        polygon(xaxs, yaxs, col = rgb(redg[i], 0, blug[i], 
                                      alpha = alpha), border = rgb(redg[i], 0, blug[i], 
                                                                   alpha = 1), lty = 2)
      }
      if ((exp1[i] + exp2[i]) == 0 & (exp1[i] + exp2[i] + exp3[i] + 
                                      exp4[i]) != 0) {
        xaxs <- c(hstr[i, 1], lstr[i, 1], kink2x, lend[i, 
                                                       1], hend[i, 1])
        xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
        yaxs <- c(hstr[i, 2], lstr[i, 2], kink2y, lend[i, 
                                                       2], hend[i, 2])
        yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
        side1 <- c(xaxs, xaxs[1])
        side2 <- c(yaxs, yaxs[1])
        c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
        c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 
                                                        1)])
        area <- abs(c1 - c2)/2
        if (area > 0 & !is.nan(area)) {
          alpha = min(0.5/(area * (n)), 1)
        }
        if (area == 0 | is.nan(area)) {
          alpha = 0.05
        }
        border = alpha
        polygon(xaxs, yaxs, col = rgb(redg[i], 0, blug[i], 
                                      alpha = alpha), border = rgb(redg[i], 0, blug[i], 
                                                                   alpha = 1), lty = 2)
      }
      if ((exp1[i] + exp2[i]) != 0 & (exp3[i] + exp4[i]) != 
          0) {
        xaxs <- c(hstr[i, 1], lstr[i, 1], kink2x, lend[i, 
                                                       1], hend[i, 1], kink1x)
        xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
        yaxs <- c(hstr[i, 2], lstr[i, 2], kink2y, lend[i, 
                                                       2], hend[i, 2], kink1y)
        yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
        side1 <- c(xaxs, xaxs[1])
        side2 <- c(yaxs, yaxs[1])
        c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
        c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 
                                                        1)])
        area <- abs(c1 - c2)/2
        if (area > 0 & !is.nan(area)) {
          alpha = min(0.5/(area * (n)), 1)
        }
        if (area == 0 | is.nan(area)) {
          alpha = 0.05
        }
        border = alpha
        polygon(xaxs, yaxs, col = rgb(redg[i], 0, blug[i], 
                                      alpha = alpha), border = rgb(redg[i], 0, blug[i], 
                                                                   alpha = 1), lty = 2)
      }
    }
  }
   
  .tomogonce <- function (input, last.input, ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    bounds <- bounds1(x, t, n)
    bbounds <- cbind(bounds[, 1], bounds[, 2])
    wbounds <- cbind(bounds[, 4], bounds[, 3])
    n <- dim(bounds)[1]
    par(mfrow = c(1, 1))
    plot(c(100, 200), xlim = c(0, 1), ylim = c(0, 1), col = "white", 
         ylab = "betaW", xlab = "betaB", xaxs = "i", yaxs = "i", 
         main = "Tomography Plot")
    if (input == "") {
      last.input <- as.integer(last.input)
      input <- last.input + 1
      for (i in 1:n) {
        lines(bbounds[i, ], wbounds[i, ], col = "yellow")
      }
      lines(bbounds[input, ], wbounds[input, ], col = "black")
    }
    else {
      input <- as.integer(input)
      for (i in 1:n) {
        lines(bbounds[i, ], wbounds[i, ], col = "yellow")
      }
      lines(bbounds[input, ], wbounds[input, ], col = "black")
    }
    return(input)
  }
   
  .getinput <- function () 
  {
    readline("Hit <enter> for next observation, enter observation number, or hit <s> to stop: ")
  }
  
  .postonce <- function (input, last.input, betab, betaw, betabs, betaws) 
  {
    if (input == "") {
      last.input <- as.integer(last.input)
      input <- last.input + 1
    }
    else {
      input <- as.integer(input)
    }
    par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
    plot(density(betab[input, ]), xlim = c(0, 1), ylim = c(0, 
                                                           max(density(betab[input, ])$y) + 1), yaxs = "i", xaxs = "i", 
         main = "Posterior Distribution of betaB", xlab = "Bb", 
         ylab = "Density")
    lines(c(0, 0.25 * (max(density(betab[input, ])$y) + 1)), 
          lwd = 3)
    plot(density(betaw[input, ]), xlim = c(0, 1), ylim = c(0, 
                                                           max(density(betaw[input, ])$y) + 1), yaxs = "i", xaxs = "i", 
         main = "Posterior Distribution of betaW", xlab = "Bw", 
         ylab = "Density")
    lines(c(0, 0.25 * (max(density(betaw[input, ])$y) + 1)), 
          lwd = 3)
    colors = runif(length(betabs), 26, 51)
    plot(betabs[input, ], betaws[input, ], xlim = c(0, 1), ylim = c(0, 
                                                                    1), xaxs = "i", yaxs = "i", main = "Simulations of betaW and betaB", 
         ylab = "betaW simulations", xlab = "betaB simulations", 
         pch = 20, col = colors, lty = 2, cex = 0.25)
    mtext(sprintf("Plots for Observation %d", input), line = 0.5, 
          outer = TRUE)
    return(input)
  }
  
  .movieD <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    bounds <- bounds1(x, t, n)
    bbounds <- cbind(bounds[, 1], bounds[, 2])
    wbounds <- cbind(bounds[, 4], bounds[, 3])
    n <- dim(bounds)[1]
    plot(c(100, 200), xlim = c(0, 1), ylim = c(0, 1), col = "white", 
         ylab = "betaW", xlab = "betaB", xaxs = "i", yaxs = "i", 
         main = "Tomography Plot")
    input <- 0
    while (input != "s") {
      input <- .tomogonce(input, last.input, ei.object)
      last.input <- input
      input <- .getinput()
    }
  }
  
  .movie <- function (ei.object) 
  {
    ok <- !is.na(ei.object$betab)
    betabs <- ei.object$betabs[ok, ]
    ok <- !is.na(ei.object$betaw)
    betaws <- ei.object$betaws[ok, ]
    betab <- ei.object$betabs
    betaw <- ei.object$betaws
    input <- 1
    last.input <- 0
    while (input != "s") {
      input <- .postonce(input, last.input, betab, betaw, betabs, 
                         betaws)
      last.input <- input
      input <- .getinput()
    }
  }
    
  ei.object <- x
  lci.function.list <- list(tomogD = .tomog, tomog = .tomogl)
  function.list <- list(tomogCI = .tomog80CI, tomogCI95 = .tomog95CI, 
                        tomogE = .tomogE, tomogP = .tomogP2, betab = .betabd, 
                        betaw = .betawd, xt = .xt, xtc = .xtc, xtfit = .xtfit, 
                        xtfitg = .xtfitg, estsims = .estsims, boundXb = .boundXb, 
                        boundXw = .boundXw, truth = .truthfn, eiRxCtomog = .bndplot, 
                        movieD = .movieD, movie = .movie)
  arguments <- list(...)
  if ("lci" %in% names(arguments)) {
    lci <- arguments$lci
    arguments$lci <- NULL
  }
  else {
    lci <- TRUE
  }
  results <- list()
  if (length(arguments) != 1) {
    row = ceiling(length(arguments)/2)
    par(mfrow = c(row, 2))
  }
  for (arg in arguments) {
    if (arg %in% names(function.list)) {
      results[[arg]] <- function.list[[arg]](ei.object = ei.object)
    }
    else if (arg %in% names(lci.function.list)) {
      results[[arg]] <- lci.function.list[[arg]](ei.object = ei.object, 
                                                 lci = lci)
    }
    else results[[arg]] <- NA
  }
}
