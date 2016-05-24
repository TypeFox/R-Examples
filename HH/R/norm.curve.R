"norm.setup" <-
  function(xlim=c(-2.5,2.5),
           ylim=c(0,.4)/se,
           mean=0,
           main=main.calc,
           se=sd/sqrt(n), sd=1, n=1,
           df.t=NULL,
           Use.obs.mean=TRUE,
           ...) {
    if (se <= 0) stop("Standard deviation must be positive")
    main.calc <-
      if (is.null(df.t) || df.t==Inf)  ## normal
        ifelse(!(((missing(se) && missing(sd)) || se==1) &&
                 missing(n) && mean==0),
               paste("normal density:  se =", round(se,3)),
               "Standard Normal Density N(0,1)")
      else { ## t distribution
        if (length(df.t) != 1) stop("df.t must have length 1")
        ifelse(!(missing(se) && missing(sd) && missing(n) && mean==0),
               paste("t density:  se = ", round(se,3), ", df = ", df.t, sep=""),
               paste("standard t density, df =", df.t))
      }
    plot(xlim=xlim, ylim=ylim,
         x=xlim, y=ylim,
         yaxt="n", type="n",
         las=1,
         xlab="",
         ylab=ifelse(is.null(df.t) || df.t==Inf, "f(z)", "f(t)"),
         main=main)
    axis(4, las=1)
}

"norm.curve" <-
function(mean=0, se=sd/sqrt(n),
         critical.values=mean + se*c(-1, 1)*z.975,
         z=if(se==0) 0 else do.call("seq", as.list(c((par()$usr[1:2]-mean)/se, length=109))),
         shade, col="blue",
         axis.name=ifelse(is.null(df.t) || df.t==Inf, "z", "t"),
         second.axis.label.line=3,
         sd=1, n=1,
         df.t=NULL,
         axis.name.expr=axis.name,
         Use.obs.mean=TRUE,
         col.label=col,
         hypoth.or.conf="Hypoth",
         col.conf.arrow=par("col"),
         col.conf.label=par("col"),
         col.crit=ifelse(hypoth.or.conf=="Hypoth", 'blue', col.conf.arrow),
         cex.crit=1.2,
         polygon.density=-1,
         polygon.lwd=4,
         col.border=ifelse(is.na(polygon.density), FALSE, col),
         ...) {
    old.err <- par(err=-1) ## prevent "out of bounds" errors from displaying
    on.exit(par(old.err))
  ## Valid values for shade are "right", "left", "inside", "outside",
  ## "none".  Default is "right" for one-sided critical.values and
  ## "outside" for two-sided critical values.  "none" is used to
  ## redraw an outline of the curve that would otherwise be obscured
  ## by a solid color from the shaded area of another curve.

  cex.center <- ifelse(hypoth.or.conf=="Hypoth", 1, cex.crit)
  z.975 <- if (is.null(df.t) || df.t==Inf) qnorm(.975) else qt(.975, df.t)

  if (missing(shade))
    shade <- switch(length(critical.values)+1,
                    "none",
                    "right",
                    "outside",
                    stop("Specify no more than 2 critical values."))

  cex.small <- par()$cex*.7

  z.critical.values <- if(se==0) c(0,0) else (critical.values-mean)/se

  dfunction <-  function(z, df.t=NULL)
    if (is.null(df.t) || df.t==Inf) dnorm(z) else dt(z, df.t)
  pfunction <-  function(z, df.t=NULL)
    if (is.null(df.t) || df.t==Inf) pnorm(z) else pt(z, df.t)

  x.z <- mean + z*se
  lines(y=dfunction(z, df.t)/se, x=x.z)
  zvals <- trunc(range(z))
  zvals <- seq(zvals[1], zvals[2], 1)
  if (axis.name=="z" || axis.name=="t") {
    if (!(missing(se) && missing(sd) && missing(n) && mean==0)) {
      axis.list <- list(1, at=mean+se*zvals, labels=zvals, tick=FALSE,
                        line=2)
      if.R(r=axis.list$cex.axis <- cex.small,
           s=axis.list$cex      <- cex.small)
      do.call("axis", axis.list)
    }
  }
  else {
    axis.list <- list(1, at=mean+se*zvals, labels=zvals, tick=FALSE,
                      line=5)
      if.R(r=axis.list$cex.axis <- cex.small,
           s=axis.list$cex      <- cex.small)
      do.call("axis", axis.list)
  }
  y.ticks <- pretty(par()$usr[3:4]*se)
  if (se != 0) axis(2, at=y.ticks/se, labels=y.ticks, las=1)
  if (!(missing(se) && missing(sd) && missing(n) && mean==0)) {
    mtext(xpd=NA,side=4,
          text=ifelse(is.null(df.t) || df.t==Inf,
            if.R(r=expression(g( ~ bar(x) ~ "") == f(( ~ bar(x)-mu[~i])/sigma[~bar(x)])/sigma[~bar(x)]),
                 s="g(xbar) = f((xbar-mu_i)/se) / se"),
            if.R(r=expression(g( ~ bar(x) ~ "" ) == f(( ~ bar(x)-mu[~i])/s[~bar(x)])/s[~bar(x)]),
                 s="g(xbar) = f((xbar-mu_i)/se) / se")),
          line=second.axis.label.line, cex=par()$cex,
          col=ifelse(second.axis.label.line==3, 1, col))
    mtext(xpd=NA,side=2,
          text=ifelse(is.null(df.t) || df.t==Inf,
            "f(z)",
            "f(t)"),
          line=second.axis.label.line, cex=par()$cex,
          col=ifelse(second.axis.label.line==3, 1, col))
  }
  critical.one <- TRUE
  if (length(critical.values)==1) {
    if (shade=="right") {
      x <- seq(z.critical.values, max(z), length=51)*se + mean
      shaded.area <- 1-pfunction(z.critical.values, df.t)
    }
    else {
      x <- seq(min(z), z.critical.values, length=51)*se + mean
      shaded.area <- pfunction(z.critical.values, df.t)
    }
  }
  if (length(critical.values)==2) {
    if (shade=="outside") {
      critical.one <- FALSE
      x1 <- seq(min(z), z.critical.values[1], length=51)*se + mean
      x2 <- seq(z.critical.values[2], max(z), length=51)*se + mean
      shaded.area <- 1-diff(pfunction(z.critical.values, df.t))
    }
    else { ## shade == "inside"
      x <- seq(z.critical.values[1], z.critical.values[2], length=51)*se + mean
      shaded.area <- diff(pfunction(z.critical.values, df.t))
    }
  }
  if (shade != "none") {
    if (critical.one) {
      polygon(x=c(x[1], x, x[length(x)]),
              y=c(0, dfunction((x-mean)/se, df.t)/se, 0),
              col=col, density=polygon.density, lwd=polygon.lwd,
              border=col.border,
              angle=if (shade=="left") 45 else -45)
      if (hypoth.or.conf=="Conf") {
        y.arrow <- par()$usr[3] - diff(par()$usr[3:4])*.06
        arrows(x[1], y.arrow, x[length(x)], y.arrow,
               length=.05, code=3, col=col.conf.arrow, lwd=2, xpd=TRUE)
      }
    }
    else {
      polygon(x=c(x1[1], x1, x1[length(x1)]),
              y=c(0, dfunction((x1-mean)/se, df.t)/se, 0),
              col=col, density=polygon.density, lwd=polygon.lwd,
              border=col.border,
              angle=45)
      polygon(x=c(x2[1], x2, x2[length(x2)]),
              y=c(0, dfunction((x2-mean)/se, df.t)/se, 0),
              col=col, density=polygon.density, lwd=polygon.lwd,
              border=col.border,
              angle=-45)
    }
  }

  axis(1, at=critical.values, tck=-.09, labels=FALSE)
  left.margin <- .15*diff(par()$usr[1:2])
  if (axis.name=="z" || axis.name=="t") {
    if (length(critical.values)>0)
    axis(1, at=critical.values, tick=FALSE, line=3,
         labels=round(z.critical.values, 3))
    if (!(missing(se) && missing(sd) && missing(n) && mean==0))
      mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=3, text=axis.name, cex=cex.small)
    else
      mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=1, text=axis.name, cex=cex.small)
    mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=4, text=axis.name, cex=cex.small)
    if (shade != "none") {
      mtext(xpd=NA,side=1, at=par()$usr[2]+.9*left.margin, line=1,
            text="shaded area", cex=par()$cex)
      pval <- format(round(shaded.area,4), digits=4, nsmall=4, scientific=FALSE)
      if (hypoth.or.conf=="Hypoth")
        mtext(xpd=NA,side=1, at=par()$usr[2]+1.45*left.margin, line=4, adj=1,
              text=if.R(r=substitute(list(alpha * " = " * group("",list(p),"")), list(p=pval)),
                s=paste("alpha =", pval)),
              cex=par()$cex, col=col.label)
      else {
        mtext(xpd=NA,side=1, at=par()$usr[2]+left.margin*(.08/.15), line=4,
              text=if.R(r=substitute(list("Conf Level= " * group("",list(p),"")), list(p=pval)),
                s=paste("Conf Level =", pval)),
              cex=par()$cex, col=col.conf.label)
      }
    }
    if (!(missing(se) && missing(sd) && missing(n) && mean==0)) {
      axis(1, at=critical.values, tick=FALSE, line=1,
         labels=round(critical.values, 3), col.axis=col.crit, cex.axis=cex.crit)
      mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=1,
            text=if.R(r=expression(bar(x)), s="xbar"), cex=cex.small)
      mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=2,
            text=if.R(r=expression(bar(x)), s="xbar"), cex=cex.small)
      mtext(xpd=NA,side=3, at=par()$usr[1]-left.margin, line=.5,
            text=if (Use.obs.mean)
            if.R(r=expression(mu ~~ bar(x)), s="mu   mean")
            else if.R(r=expression(mu), s="mu"),
            cex=par()$cex)
    }
    else
      mtext(xpd=NA,side=3, at=par()$usr[1]-left.margin, line=.5,
            text=if (Use.obs.mean)
            if.R(r=expression(mu ~~ bar(x)), s="mu   mean")
            else if.R(r=expression(mu), s="mu"),
            cex=par()$cex)

  }
  else { ## (axis.name=="z1" || axis.name=="t1")
    axis(1, at=critical.values, tick=FALSE, line=6,
         labels=round(z.critical.values, 3))
    mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=6,
          text=axis.name.expr, cex=cex.small)
    mtext(xpd=NA,side=1, at=par()$usr[1]-left.margin, line=7,
          text=axis.name.expr, cex=cex.small)
    if (shade != "none") {
      pval <- format(round(shaded.area,4), digits=4, nsmall=4, scientific=FALSE)
      mtext(xpd=NA,side=1, at=par()$usr[2]+1.45*left.margin, line=7, adj=1,
            text=if.R(
              r=substitute(list(beta * " = " * group("",list(p),"")), list(p=pval)),
              s=paste("beta =", pval)),
            cex=par()$cex, col=col.label)
    }
  }
  axis(1, at=mean, tck=-.04, labels="")
  axis(3, at=mean, tck=-.02, labels="", xpd=TRUE)
  axis(3, at=mean, tick=FALSE, xpd=TRUE, line=-.5, labels=round(mean, 3),
                    cex.axis=cex.center)
  abline(h=0, v=mean)
  invisible(NULL)
}

norm.observed <-
  function(xbar, t.xbar, t.xbar.H1=NULL,
           col="green",
           p.val=NULL, p.val.x=par()$usr[2] + left.margin,
           t.or.z=ifelse(is.null(deg.free) || deg.free==Inf, "z", "t"),
           t.or.z.position=par()$usr[1]-left.margin,
           cex.small=par()$cex*.7, col.label=col,
           xbar.negt=NULL, cex.large=par()$cex,
           left.margin=.15*diff(par()$usr[1:2]),
           sided="", deg.free=NULL) {
  abline(v=xbar, col=col, lty=5, lwd=3)
  xbar2 <- c(xbar, xbar)
  arrows(xbar2, par()$usr[3:4]+c(-.01,.01), xbar2, par()$usr[3:4],
         xpd=TRUE, col=col.label, length=.1)
  axis(side=1, at=xbar, labels=FALSE, col=col)
  axis(side=3, at=xbar, labels=FALSE, col=col)
  mtext(xpd=NA,side=3, text=round(xbar,3), at=xbar, line=.5,
        cex=cex.large, col=col.label)
  mtext(xpd=NA,side=1, text=round(t.xbar,3), at=xbar,
        line=5, cex=par()$cex, col=col.label)
  mtext(xpd=NA,side=1, text=t.or.z, at=t.or.z.position,
        line=5, cex=cex.small)
  if (!is.null(xbar.negt))
    mtext(xpd=NA,side=1, text=round(-t.xbar,3), at=xbar.negt,
          line=5, cex=par()$cex, col=col.label)
  if (!is.null(p.val)) {
    pval <- format(round(p.val,4), digits=4, nsmall=4, scientific=FALSE)
    mtext(xpd=NA,side=1, at=par()$usr[2]+1.45*left.margin, line=5, adj=1,
          text=if.R(r=substitute(list(p * " = " * group("",list(pv),"")), list(pv=pval)),
            s=paste("p =", pval)),
          cex=par()$cex, col=col.label)
  }
  if (!is.null(t.xbar.H1)) {
    mtext(xpd=NA,side=1, text=round(t.xbar.H1,3), at=xbar,
          line=8, cex=par()$cex, col=col.label)
    t.or.z.expr <-
      if (t.or.z=="z") if.R(r=expression(z[1]), s="z1")
      else if.R(r=expression(t[1]), s="t1")
    mtext(xpd=NA,side=1, text=t.or.z.expr, at=t.or.z.position,
          line=8, cex=cex.small)
  }

}

norm.outline <- function(dfunction, left, right, mu.H0, se, deg.free=NULL,
                         col.mean="green") {
        x.p <- seq(left, right, length=109)
        z.p.H0 <- (x.p-mu.H0)/se
        y.p <-
          if (dfunction=="dnorm")
            dnorm(z.p.H0)
          else
            dt(z.p.H0, df=deg.free)
        lines(y=c(0, y.p/se, 0, 0),
              x=c(left, x.p, right, left),
              col=col.mean, lwd=3)
      }



normal.and.t.dist<-
  function(
           mu.H0          = 0,
           mu.H1          = NA,
           obs.mean       = 0,
           std.dev        = 1,
           n              = NA,
           deg.freedom    = NA,
           alpha.left     = alpha.right,
           alpha.right    = .05,
           Use.mu.H1      = FALSE,
           Use.obs.mean   = FALSE,
           Use.alpha.left = FALSE,
           Use.alpha.right= TRUE,
           hypoth.or.conf = 'Hypoth',
           xmin           = NA,
           xmax           = NA,
           gxbar.min      = NA,
           gxbar.max      = NA,
           cex.crit       = 1.2,
           polygon.density= -1,
           polygon.lwd    = 4,
           col.mean       = 'limegreen',
           col.mean.label = 'limegreen',
           col.alpha      = 'blue',
           col.alpha.label= 'blue',
           col.beta       = 'red',
           col.beta.label = 'red',
           col.conf       = 'palegreen',
           col.conf.arrow = 'darkgreen',
           col.conf.label = 'darkgreen'
           )
  {
    pfunction <-  function(z, df.t=NULL)
      if (is.null(df.t) || df.t==Inf) pnorm(z) else pt(z, df.t)

    qfunction <-  function(p, df.t=NULL)
      if (is.null(df.t) || df.t==Inf) qnorm(p) else qt(p, df.t)

    dfunction <-  function(z, df.t=NULL)
      if (is.null(df.t) || df.t==Inf) dnorm(z) else dt(z, df.t)

    is.na.or.blank <- function(x) is.na(x) || x==""

    old.par <- par(oma=c(4,0,2,5), mar=c(7,7,4,2)+.1)
    deg.free <- if(is.na.or.blank(deg.freedom))  NULL else deg.freedom
    dfunction.name <- if (is.null(deg.free) || deg.free==Inf) "dnorm" else "dt"
    normal <- is.na.or.blank(deg.freedom)
    standard <- is.na.or.blank(n) && (is.na.or.blank(std.dev) || std.dev==1) && mu.H0==0
    if (!(is.na.or.blank(std.dev)) && std.dev <= 0) stop("Standard deviation must be positive")
    standard.normal <- standard && normal

    n.conf <- if (is.na.or.blank(n))       1 else n
    sd     <- if (is.na.or.blank(std.dev)) 1 else std.dev
    se     <- sd/sqrt(n.conf)

    center <- if (hypoth.or.conf=="Hypoth") mu.H0 else obs.mean
    crit.val.z <- qfunction(1-alpha.right, deg.free)
    crit.val <- center + crit.val.z * se
    crit.val.left.z <- qfunction(1-alpha.left, deg.free)
    crit.val.left <- center - crit.val.left.z * se

    cv <- c(crit.val.left[Use.alpha.left], crit.val[Use.alpha.right])

    if (is.na.or.blank(xmin)) xmin <-
      if (hypoth.or.conf=='Hypoth')
        min(mu.H0-3*se, mu.H1-2.5*se, obs.mean-.5*se, na.rm = TRUE)
      else ## 'Conf'
        min(crit.val.left-.5*se, obs.mean-3*se, na.rm = TRUE)

    if (is.na.or.blank(xmax)) xmax <-
      if (hypoth.or.conf=='Hypoth')
        max(mu.H0+3*se, mu.H1+2.5*se, obs.mean+.5*se, na.rm = TRUE)
      else ## 'Conf'
        max(crit.val+.5*se, obs.mean+3*se, na.rm = TRUE)

    if (is.na.or.blank(gxbar.min)) gxbar.min <- 0
    if (is.na.or.blank(gxbar.max)) gxbar.max <- dfunction(0, df.t=deg.free) / se

    conf.level.fract <- 1
    if (Use.alpha.left) conf.level.fract <- conf.level.fract - alpha.left
    if (Use.alpha.right) conf.level.fract <- conf.level.fract - alpha.right

    ## main titles for graph
    if.R(r={
      standard.normal.main <- "Standard Normal Density N(0,1)"

      normal.main <-
        substitute(list("normal density:  " *
                        sigma[bar(bold(x))] * " = " * group("",list(se),"")) * ", " *
                        " n = " * group("",n.conf,"") ,
                   list(se=round(se,3), n.conf=n.conf))

      normal.conf.main <-
        substitute(list(conf.level * "%  Normal Confidence Limits:  " *
                        sigma[bar(bold(x))] * " = " * group("",list(se),"") * ", " *
                        "n = " * group("",n.conf,"")),
                   list(se=round(se,3), n.conf=n.conf,
                        conf.level=100*conf.level.fract))

      standard.t.main <-
        substitute(list("Standard t Density, " *  nu * " = " * group("",list(df),"")),
                   list(df=deg.free))

      t.dist.main <-
        substitute(list("t density:  " *
                        "s"[bar(bold(x))] * " = " * group("",list(se),"") * ", " *
                        "n =" * group("",n.conf,"") * ", " *
                        nu * " = " * group("",list(df),"")),
                   list(se=round(se,3), n.conf=n.conf, df=deg.free))

      t.conf.main <-
        substitute(list(conf.level * "%  t Confidence Limits:  " *
                        s[bar(bold(x))] * " = " * group("",list(se),"") * ", " *
                        "n = " * group("",n.conf,"") * ", " *
                        nu * " = " * group("",list(df),"")),
                   list(se=round(se,3), n.conf=n.conf, df=deg.free,
                        conf.level=100*conf.level.fract))

    },s={
      standard.normal.main <- "Standard Normal Density N(0,1)"

      normal.main <- paste("normal density:  se =", round(se,3),
                           ",  n = ", n.conf, sep="")

      normal.conf.main <- paste(round(100*conf.level.fract),
                                "%  Normal Confidence Limits:  se = ", round(se,3),
                                ",  n = ", n.conf, sep="")

      standard.t.main <- paste("Standard t Density, df =", deg.free)

      t.dist.main <- paste("t density:  se =", round(se,3),
                                ",  n = ", n.conf, ",  df = ", deg.free, sep="")

      t.conf.main <- paste(round(100*conf.level.fract),
                                "%  t Confidence Limits:  se = ", round(se,3),
                                ",  n = ", n.conf, ",  df = ", deg.free, sep="")
    })

    if (normal) {
      if (hypoth.or.conf=="Hypoth") {
        if (standard.normal)
          norm.setup(mean=mu.H0, xlim=c(xmin,xmax),
                     ylim=c(gxbar.min, gxbar.max), main=standard.normal.main)
        else
          norm.setup(mean=mu.H0, xlim=c(xmin,xmax),
                     ylim=c(gxbar.min, gxbar.max), main=normal.main,
                     se=se, n=n.conf)
      } else { ## Conf
        norm.setup(mean=obs.mean, xlim=c(xmin,xmax),
                   ylim=c(gxbar.min, gxbar.max),
                   sd=sd,
                   se=se, n=n.conf,
                   main=normal.conf.main)
                   ##paste(
                   ##  "Normal Confidence Limits:  se =", round(se,3),
                   ##  "n =", n.conf))
      }
    } else { ## t
      if (hypoth.or.conf=="Hypoth") {
        if (standard)
          norm.setup(mean=mu.H0, df.t=deg.free,
                     xlim=c(xmin,xmax), ylim=c(gxbar.min, gxbar.max),
                     main=standard.t.main)
        else
          norm.setup(mean=mu.H0, n=n.conf, se=se, df.t=deg.free,
                     xlim=c(xmin,xmax), ylim=c(gxbar.min, gxbar.max),
                     main=t.dist.main)
      } else { ## Conf
        norm.setup(mean=obs.mean, n=n.conf, se=se, df.t=deg.free,
                   xlim=c(xmin,xmax), ylim=c(gxbar.min, gxbar.max),
                   main=t.conf.main)
                   ## paste(
                   ##  "t Confidence Limits:  se =", round(se,3),
                   ##  "n =", n.conf, "df =", deg.free))
      }
    }

    if ( Use.alpha.left &&  Use.alpha.right) {cv.shade <- "outside"; cv.altshade <- "inside"}
    if (!Use.alpha.left &&  Use.alpha.right) {cv.shade <- "right"  ; cv.altshade <- "left"  }
    if ( Use.alpha.left && !Use.alpha.right) {cv.shade <- "left"   ; cv.altshade <- "right" }
    if (!Use.alpha.left && !Use.alpha.right) {cv.shade <- "none"   ; cv.altshade <- "none"  }

    if (hypoth.or.conf=="Hypoth" && Use.mu.H1) {
        t.or.z.expr <-
            if (is.null(deg.free)) if.R(r=expression(z[1]), s="z1")
            else if.R(r=expression(t[1]), s="t1")
        if (standard.normal) {
            norm.curve(critical.values=cv, mean=mu.H1,
                       col=col.beta, shade=cv.altshade,
                       axis.name=if(is.null(deg.free)) 'z1' else 't1',
                       axis.name.expr=t.or.z.expr,
                       Use.obs.mean=Use.obs.mean, col.label=col.beta.label,
                       polygon.density=polygon.density, polygon.lwd=polygon.lwd)
        } else {
            norm.curve(critical.values=cv, se=se, df.t=deg.free, mean=mu.H1,
                       col=col.beta, shade=cv.altshade,
                       axis.name=if(is.null(deg.free)) 'z1' else 't1',
                       axis.name.expr=t.or.z.expr,
                       Use.obs.mean=Use.obs.mean, col.label=col.beta.label,
                       polygon.density=polygon.density, polygon.lwd=polygon.lwd)
        }
    }

    if (hypoth.or.conf=="Hypoth") {
      if (standard) {
      if (standard.normal) {
        norm.curve(critical.values=cv, mean=mu.H0,
                   col=col.alpha, shade=cv.shade,
                   Use.obs.mean=Use.obs.mean, col.label=col.alpha.label,
                   polygon.density=polygon.density, polygon.lwd=polygon.lwd)
      }
      else
        norm.curve(critical.values=cv, df.t=deg.free, mean=mu.H0,
                   col=col.alpha, shade=cv.shade,
                   Use.obs.mean=Use.obs.mean, col.label=col.alpha.label,
                   polygon.density=polygon.density, polygon.lwd=polygon.lwd)
      }
      else {
        norm.curve(critical.values=cv, se=se, df.t=deg.free, mean=mu.H0,
                   col=col.alpha, shade=cv.shade,
                   Use.obs.mean=Use.obs.mean, col.label=col.alpha.label,
                   polygon.density=polygon.density, polygon.lwd=polygon.lwd)
      }
    } else { ## "Conf"
        ## if (standard.normal) {
        ##   norm.curve(critical.values=cv, mean=obs.mean,
        ##              col=col.alpha, shade=cv.shade,
        ##              Use.obs.mean=Use.obs.mean, col.label=col.alpha.label,
        ##              polygon.density=polygon.density, polygon.lwd=polygon.lwd)
        ## } else {
          norm.curve(critical.values=cv,
                     se=se, df.t=deg.free, mean=obs.mean,
                     col=col.conf, shade=cv.altshade,
                     Use.obs.mean=Use.obs.mean, col.label=col.alpha.label,
                     hypoth.or.conf="Conf",
                     col.conf.arrow=col.conf.arrow,
                     col.conf.label=col.conf.label,
                     cex.crit=cex.crit,
                     polygon.density=polygon.density, polygon.lwd=polygon.lwd)
        ##}
    }

    if (!Use.obs.mean) obs.mean <- center
    obs.mean.H0.z <- if((obs.mean-center)==0) 0 else (obs.mean-center)/se
                     ## yes, exact 0 stays 0 in any scale.
    obs.mean.H0.p.val <- pfunction(obs.mean.H0.z, df.t=deg.free)

    if ( Use.alpha.left &&  Use.alpha.right) {
      obs.mean.H0.p.val <- ifelse(obs.mean.H0.z >= 0,
                                  2*(1-obs.mean.H0.p.val),
                                  2*obs.mean.H0.p.val)
      side <- "two-sided"
      if (hypoth.or.conf=="Hypoth" && Use.obs.mean) {
        obs.mean.z.pos <- abs(obs.mean.H0.z)
        obs.mean.x.pos <- mu.H0 + obs.mean.z.pos*se
        obs.mean.x.neg <- mu.H0 - obs.mean.z.pos*se
        ## right side
        norm.outline(dfunction.name, obs.mean.x.pos, par()$usr[2], mu.H0, se, deg.free, col.mean)
        ## left side
        norm.outline(dfunction.name, par()$usr[1], obs.mean.x.neg, mu.H0, se, deg.free, col.mean)
      }
      else {
        obs.mean.z.pos <- 0
        obs.mean.x.pos <- 0
        obs.mean.x.neg <- 0
      }
    }
    if (!Use.alpha.left &&  Use.alpha.right) {
      obs.mean.H0.p.val <- 1-obs.mean.H0.p.val
    obs.mean.z.pos <- obs.mean.H0.z
    obs.mean.x.pos <- obs.mean
    obs.mean.x.neg <- 0 ## place holder
      side <- "right-sided"
      if (hypoth.or.conf=="Hypoth" && Use.obs.mean) {
        norm.outline(dfunction.name, obs.mean, par()$usr[2], mu.H0, se, deg.free, col.mean)
      }
    }
    if ( Use.alpha.left && !Use.alpha.right) {
      obs.mean.H0.p.val <- obs.mean.H0.p.val ## redundant but legible
      obs.mean.z.pos <- abs(obs.mean.H0.z)
      obs.mean.x.pos <- 0 ## place holder
      obs.mean.x.neg <- obs.mean
      side <- "left-sided"
      if (hypoth.or.conf=="Hypoth" && Use.obs.mean) {
        norm.outline(dfunction.name, par()$usr[1], obs.mean, mu.H0, se, deg.free, col.mean)
      }
    }
    if (!Use.alpha.left && !Use.alpha.right) {
      obs.mean.H0.p.val <- 0
      obs.mean.z.pos <- 0
      obs.mean.x.pos <- 0
      obs.mean.x.neg <- 0
      side <- ""
    }

    obs.mean.H1.z <- (obs.mean-mu.H1)/se
    if (hypoth.or.conf=="Hypoth" && Use.obs.mean) {
      left.margin <- .15*diff(par()$usr[1:2])
      norm.observed(obs.mean, obs.mean.H0.z,
                    if (Use.mu.H1) obs.mean.H1.z else NULL,
                    col=col.mean,
                    p.val=obs.mean.H0.p.val, p.val.x=par()$usr[2]+ left.margin,
                    t.or.z=ifelse(is.null(deg.free) || deg.free==Inf, "z", "t"),
                    t.or.z.position=par()$usr[1]-left.margin,
                    cex.small=par()$cex*.7, col.label=col.mean.label,
                    xbar.negt=if (side=="two-sided") {
                      if (obs.mean==obs.mean.x.pos) obs.mean.x.neg else obs.mean.x.pos
                    }
                    else NULL,
                    cex.large=cex.crit)
    }
    crit.val.H1 <- if(Use.mu.H1) (crit.val-mu.H1)/se else ""
    crit.val.H1.left <- if(Use.mu.H1) (crit.val.left-mu.H1)/se else ""
    beta.left <- if (Use.mu.H1) pfunction(crit.val.H1, deg.free) else ""
    beta.right <- if (Use.mu.H1) 1-pfunction(crit.val.H1.left, deg.free) else ""
    beta.middle <- if (Use.mu.H1) diff(pfunction(c(crit.val.H1.left, crit.val.H1), deg.free)) else ""

    ## par(old.par) ## oma and mar are potentially useful
    invisible(list(beta.left=beta.left,
                   beta.middle=beta.middle,
                   beta.right=beta.right,
                   crit.val=crit.val,
                   crit.val.H1=crit.val.H1,
                   crit.val.H1.left=crit.val.H1.left,
                   crit.val.left=crit.val.left,
                   crit.val.left.z=-crit.val.left.z,
                   crit.val.z=crit.val.z,
                   obs.mean.H0.p.val=obs.mean.H0.p.val,
                   obs.mean.H0.side=side,
                   obs.mean.H0.z=obs.mean.H0.z,
                   obs.mean.H1.z=obs.mean.H1.z,
                   obs.mean.x.neg=obs.mean.x.neg,
                   obs.mean.x.pos=obs.mean.x.pos,
                   obs.mean.z.pos= obs.mean.z.pos,
                   standard=standard,
                   standard.error=se,
                   standard.normal=standard.normal))
  }

## source("~/HH-R.package/HH/R/norm.curve.R")
