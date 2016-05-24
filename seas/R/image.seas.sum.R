"image.seas.sum" <-
function(x, var, norm = "days", start = 1, rep = 0, zlim, alim,
         palette = colorRampPalette(c("white", "blue"))(64),
         year.filter, power, contour = TRUE, show.median, main, ...) {
  old.par <- par(no.readonly=TRUE); on.exit(par(old.par))
  orig <- as.character(substitute(x))[[1]]
  x <- seas.sum.check(x, orig, var, norm, year.filter)
  var <- x$var
  if (missing(show.median))
    show.median <- if (norm == "days" &&
                      start == 1 && rep == 0) TRUE else FALSE
  if (show.median && (norm != "days" || start != 1 || rep != 0)) {
    warning(gettextf("option %s only works if %s (all are default)",
                     sQuote("show.median = TRUE"),
                     sQuote("norm=\"days\", start=1, rep=0")))
    show.median <- FALSE
  }
  fun <- sprintf("%s/%s", var, norm)
  if (missing(power) || power == 1)
    power <- FALSE
  if (power)
    fun <- sprintf("(%s)^%s", fun, round(power, 1))
  n.years <- length(x$years)
  n.bins <- length(x$bins)
  num <- n.bins + rep
  # .b suffix is a matrix of sums in the width of the bin for each year
  var.b <- x$seas[,,var, drop=TRUE]
  norm.b <- x$norm[,,var, drop=TRUE]
  seas <- var.b / norm.b
  seas[is.finite(var.b) & norm.b==0] <- 0  # avoid div/0 -> Inf
  if (missing(zlim))
    zlim <- c(0, max(seas, na.rm=TRUE))
  else if (length(zlim) == 1)
    zlim <- c(0, zlim)
  seas[seas < zlim[1]] <- zlim[1]  # trim lower values to range
  seas[seas > zlim[2]] <- zlim[2]  # trim upper
  if (power) {
    seas.p <- seas^power
    zlim.p <- zlim^power
  } else {
    seas.p <- seas
    zlim.p <- zlim
  }
  dn <- dimnames(seas.p)
  sel <- 1:num
  if (rep != 0 || start != 1) {
    sel <- (start:(start - 1 + num) - 1) %% n.bins + 1
    dn[[2]] <- dn[[2]][sel]
    seas.p <- array(seas.p[,sel], dim=c(n.years, num), dimnames=dn)
  }
  if (missing(main))
    main <- .seastitle(id=x$id, orig=orig, name=x$name,
                       fun=fun, range=x$year.range)
  main.height <- if (is.null(main) || is.na(main)) 0.1 else 1.0
  if (show.median) {
    var.a <- x$ann[,var]
    nm <- seas.norm(x, var=var, fun="median")
    nm.quan <- nm$quantile[1, var]
    pr <- 0.8  # percent of window is the main plot on left
    nf <- layout(matrix(c(1, 1, 1, 2, 5, 4, 3, 6, 4), nrow=3, byrow=TRUE),
                 widths=c(pr, 1 - pr, lcm(1.6)), heights=lcm(main.height))
    mar2 <- c(1.0, 4.1, 2.6, 0.0)
    mar3 <- c(2.6, 4.1, 2.0, 0.0)
    mar4 <- c(1.0, 0.0, 2.6, 0.6)
    mar5 <- c(2.6, 0.0, 2.0, 0.6)
  } else {
    nf <- layout(matrix(c(1, 1, 2, 4, 3, 4), nrow=3, byrow=TRUE),
                 widths=c(1, lcm(1.6)), heights=lcm(main.height))
    mar2 <- c(1.0, 4.1, 2.6, 0.6)
    mar3 <- c(2.6, 4.1, 2.0, 0.6)
  }
  par(mar=rep(0, 4))  # title: no margins
  lwd <- par("lwd")
  frame()
  text(0.5, 0.5, main, cex=par("cex.main"), font=par("font.main"))
  par(mar=mar2)
  if (x$start.day == 1) {
    ylab2 <- gettext("year")
    year.lab <- dn[[1]]
  } else {
    ylab2 <- gettext("annum")
    year.lab <- sub("_", "\n", dn[[1]])
  }
  image(1:num, 1:length(x$years), t(seas.p), col=palette, xlab=NA, ylab=ylab2,
        xaxt="n", yaxt="n", zlim=zlim.p, bty="n")
  .seasmonthgrid(x$width, x$bin.lengths, start, rep, x$start.day)
  box()
  axis(1, 1:num, dn[[2]], lwd=lwd)  # x-axis
  axis(2, 1:n.years, year.lab, lwd=lwd)  # y-axis
  op.na <- getOption("seas.na")
  if (!is.null(op.na$pch)) {
    na.xy <- which(is.na(t(seas.p)), TRUE)
    if (length(na.xy) > 0)
      points((1:num)[na.xy[,1]], (1:n.years)[na.xy[,2]],
             pch=op.na$pch, col=op.na$col)
  }
  seas.s <- apply(seas, 2, sort, na.last=FALSE)
  if (power)
    seas.p <- seas.s^power
  else
    seas.p <- seas.s
  if (rep != 0 || start != 1) {
    seas.p <- array(seas.p[,sel], dim=c(n.years, num), dimnames=dn)
    seas.s <- array(seas.s[,sel], dim=c(n.years, num), dimnames=dn)
  }
  sam.q <- (1:n.years - 1) / (n.years - 1) * 100
  par(mar=mar3)
  xlab3 <- .seasxlab(x$width, x$start.day)
  ylab3 <- gettext("sample quantiles (%)")
  image(1:num, sam.q, t(seas.p), xaxt="n", yaxt="n", bty="n",
        zlim=zlim.p, col=palette, xlab=xlab3, ylab=ylab3)
  .seasmonthgrid(x$width, x$bin.lengths, start, rep, x$start.day, FALSE)
  box()
  if (contour) {
    contour(1:num, sam.q, t(seas.s), add=TRUE)
  }
  if (show.median) {
    op.median <- getOption("seas.median")
    op.mean <- getOption("seas.mean")
    abline(h=nm.quan * 100, col=op.median$col,
           lwd=op.median$lwd * lwd, lty=op.median$lty)
  }
  axis(2, lwd=lwd)
  axis(3, 1:num, NA, lwd=lwd)
  title(xlab=xlab3, line=0.5)
  par(mar=c(2.6, 0.0, 2.6, 3.6), xaxs="i", yaxs="i", bty="n")
  frame()
  plot.window(c(0, 1), zlim)
  lev <- seq(zlim.p[1], zlim.p[2], length.out=length(palette) + 1)
  if (power)
    lev <- lev^(1/power)
  rect(0, lev[-length(lev)], 1, lev[-1], col=palette, border=NA)
  box(bty="o", yaxt="s")
  par(yaxt="s")
  #if (power) {
  #  l <- pretty(zlim, 10)
  #  ll <- l^power
  #  axis(4, ll * zlim.p[2] / max(ll[-length(ll)]), l)
  #} else axis(4)
  axis(4, lwd=lwd)
  seas.units <- if (is.null(x$units[[var]])) NULL else gettextf("%s/day", x$units[[var]])
  mtext(.seasylab(x$var, x$long.name[[var]], seas.units),
        4, 2.1, cex=par("cex.axis") * 0.8)
  if (show.median) {
    nm.mean <- mean(var.a, na.rm=TRUE)
    nm.median <- nm$ann[1, var]
    ann.a <- apply(seas.s, 2, quantile, probs=sam.q/100, na.rm=TRUE)
    days.a <- x$days
    ann.at <- rowSums(ann.a * days.a, na.rm=TRUE)
    sam.q <- c(sam.q, sam.q[length(sam.q)])
    ann.at <- c(ann.at, ann.at[length(ann.at)])
    if (missing(alim))  # weight'd mean between theoretical max and real max
      alim <- c(0, max(var.a, na.rm=TRUE) * 0.9 + max(ann.at, na.rm=TRUE) * 0.1)
    else if (length(alim) == 1)
      alim <- c(0, alim)
    var.a <- c(var.a, var.a[length(var.a)])
    if (is.na(var.a[length(var.a)]) && length(var.a) > 2)
      var.a[length(var.a) - 1] <- var.a[length(var.a) - 2]
    par(mar=mar4, bg="white", bty="o")
    plot(var.a, 1:(n.years + 1), type="S", xaxt="n", yaxt="n",
         xaxs="i", yaxs="i", xlim=alim, ylim=c(1, n.years + 1))
    abline(v=c(nm.mean, nm.median), col=c(op.mean$col, op.median$col),
           lwd=c(op.mean$lwd, op.median$lwd) * lwd,
           lty=c(op.mean$lty, op.median$lty))
    axis(1, labels=FALSE, lwd=lwd)
    mtext(.seasylab(orig, getOption("seas.label")$ann, x$units[[var]]),
          1, line=par("mgp")[2], cex=par("cex.axis") * 0.8)
    par(mar=mar5, bg="white")
    sam.qp <- seq(0, 100, length.out=length(ann.at))
    plot(ann.at, sam.qp, type="S", yaxt="n", xaxt="n", xaxs="i",
         yaxs="i", xlim=alim, ylim=c(0, 100))
    abline(v=c(nm.mean, nm.median), col=c(op.mean$col, op.median$col),
           lwd=c(op.mean$lwd, op.median$lwd) * lwd,
           lty=c(op.mean$lty, op.median$lty))
    abline(h=(nm.quan + (0.5 - nm.quan) / n.years) * 100, col=op.median$col,
           lwd=op.median$lwd * lwd, lty=op.median$lty)
    axis(1, lwd=lwd)
    axis(3, labels=FALSE, lwd=lwd)
  }
}
