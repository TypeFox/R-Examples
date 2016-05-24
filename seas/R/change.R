"change" <-
function(x1, x2, var1, var2 = var1, width = "mon",
         cent = "mean", sprd = "sd", disc = FALSE, inter = FALSE,
         p.cut = 0.3, start.day = 1, calendar) {
  orig1 <- as.character(substitute(x1))[[1]]
  sc1 <- seas.df.check(x1, orig1, var1)
  orig2 <- as.character(substitute(x2))[[1]]
  sc2 <- seas.df.check(x2, orig2, var2)
  if (is.function(cent))
    cent <- as.character(substitute(cent))[[1]]
  if (is.function(sprd))
    sprd <- as.character(substitute(sprd))[[1]]
  if (missing(calendar))
    calendar <- c(sc1$calendar, sc2$calendar)
  else if (length(calendar) == 1)
    calendar[2] <- calendar[1]
  calendar <- c(calendar[1], calendar[2])
  l <- list()
  l$orig <- list(orig1, orig2)
  l$var <- list(var1, var2)
  l$id <- list(sc1$id, sc2$id)
  l$name <- list(sc1$name, sc2$name)
  l$year.range <- list(sc1$year.range, sc2$year.range)
  l$long.name <- list(sc1$long.name, sc2$long.name)
  l$units <- list(sc1$units, sc2$units)
  l$width <- width
  l$cent.fun <- cent
  l$sprd.fun <- sprd
  l$disc <- disc
  if (!disc) {  # temperature, etc. ... generally normally distributed vars
                # s1, s2: seasonal groups
    x1$fact <- mkseas(x1, width, start.day, calendar[1])
    x2$fact <- mkseas(x2, width, start.day, calendar[2])
    s1 <- split(x1[,var1], x1$fact)
    if (length(var1) > 1)
      s1 <- lapply(s1, unlist, recursive=FALSE, use.names=FALSE)
    s2 <- split(x2[,var2], x2$fact)
    if (length(var2) > 1)
      s2 <- lapply(s2, unlist, recursive=FALSE, use.names=FALSE)
    cent1 <- eval(call("sapply", s1, cent, na.rm=TRUE))
    cent2 <- eval(call("sapply", s2, cent, na.rm=TRUE))
    l$abs <- cent2 - cent1
    if (min(x1[,var1], na.rm=TRUE) >= 0) {  # only for absolute values
      abs.vals <- TRUE  # e.g, precip, etc.
      l$rel <- cent2 / cent1
    } else abs.vals <- FALSE  # doesn't make sense for negative values
    sprd1 <- eval(call("sapply", s1, sprd, na.rm=TRUE))
    sprd2 <- eval(call("sapply", s2, sprd, na.rm=TRUE))
    l$sprd.rel <- sprd2 / sprd1
    if (length(var1) > 1)
      a1 <- unlist(x1[,var1], recursive=FALSE, use.names=FALSE)
    else
      a1 <- x1[,var1]
    if (length(var2) > 1)
      a2 <- unlist(x2[,var2], recursive=FALSE, use.names=FALSE)
    else
      a2 <- x2[,var2]
    ann.cent1 <- eval(call(cent, a1, na.rm=TRUE))
    ann.cent2 <- eval(call(cent, a2, na.rm=TRUE))
    l$ann.abs <- ann.cent2 - ann.cent1
    if (abs.vals)
      l$ann.rel <- ann.cent2 / ann.cent1
    ann.sprd1 <- eval(call(sprd, a1, na.rm=TRUE))
    ann.sprd2 <- eval(call(sprd, a2, na.rm=TRUE))
    l$ann.sprd.rel <- ann.sprd2 / ann.sprd1
  } else {  # discontinuous variable, such as precipitation
    if (length(p.cut) == 1)
      p.cut[2] <- p.cut[1]
    l$p.cut <- p.cut
    ss1 <- seas.sum(x1, var1, width, start.day, a.cut=p.cut[1])
    snm1 <- anm1 <- vnm1 <- matrix(nrow=length(var1), ncol=length(ss1$bins),
                                   dimnames=list(var1, ss1$bins))
    for (v in var1) {
      nm1 <- seas.norm(ss1, var=v, fun=cent)$seas
      snm1[v,] <- nm1[,v]
      anm1[v,] <- nm1[,"active"]
      vnm1[v,] <- seas.norm(ss1, var=v, fun=sprd)$seas[,v]
    }
    dayspy1 <- eval(call(cent, ss1$ann$days, na.rm=TRUE))
    ss2 <- seas.sum(x2, var2, width, start.day, a.cut=p.cut[2])
    snm2 <- anm2 <- vnm2 <- matrix(nrow=length(var2), ncol=length(ss2$bins),
                                   dimnames=list(var2, ss2$bins))
    for (v in var2) {
      nm2 <- seas.norm(ss2, var=v, fun=cent)$seas
      snm2[v,] <- nm2[,v]
      anm2[v,] <- nm2[,"active"]
      vnm2[v,] <- seas.norm(ss2, var=v, fun=sprd)$seas[,v]
    }
    dayspy2 <- eval(call(cent, ss2$ann$days, na.rm=TRUE))
    cent1 <- apply(snm1, 2, cent)
    cent2 <- apply(snm2, 2, cent)
    l$abs <- cent2 - cent1
    l$rel <- cent2 / cent1
    anm1 <- apply(anm1, 2, cent)
    anm2 <- apply(anm2, 2, cent)
    l$active.abs <- anm2 - anm1
    l$active.rel <- anm2 / anm1
    sprd1 <- apply(vnm1, 2, cent)  # average of the variations
    sprd2 <- apply(vnm2, 2, cent)
    l$sprd.rel <- sprd2 / sprd1
    if (length(var1) > 1)
      a1 <- unlist(ss1$ann[,var1], recursive=FALSE, use.names=FALSE)
    else
      a1 <- ss1$ann[,var1]
    if (length(var2) > 1)
      a2 <- unlist(ss2$ann[,var2], recursive=FALSE, use.names=FALSE)
    else
      a2 <- ss2$ann[,var2]
    ann.cent1 <- eval(call(cent, a1, na.rm=TRUE)) / dayspy1
    ann.cent2 <- eval(call(cent, a2, na.rm=TRUE)) / dayspy2
    l$ann.abs <- ann.cent2 - ann.cent1
    l$ann.rel <- ann.cent2 / ann.cent1
    ann.active1 <- eval(call(cent, ss1$ann[,"active"], na.rm=TRUE)) / dayspy1
    ann.active2 <- eval(call(cent, ss2$ann[,"active"], na.rm=TRUE)) / dayspy2
    l$ann.active.abs <- ann.active2 - ann.active1
    l$ann.active.rel <- ann.active2 / ann.active1
    ann.sprd1 <- eval(call(sprd, a1, na.rm=TRUE)) / dayspy1
    ann.sprd2 <- eval(call(sprd, a2, na.rm=TRUE)) / dayspy2
    l$ann.sprd.rel <- ann.sprd2 / ann.sprd1
    # evaluate interarrivals, or wet- and dry-spells
    if (inter) {
      for (v in var1) {
        ir <- interarrival(x1, var=v, p.cut=p.cut[1])
        if (v == var1[1])
          ir1 <- ir
        else
          ir1 <- rbind(ir1, ir)
      }
      ir1$fact <- mkseas(ir1, width, start.day)
      wet1 <- eval(call("tapply", ir1$wet, ir1$fact, cent, na.rm=TRUE))
      dry1 <- eval(call("tapply", ir1$dry, ir1$fact, cent, na.rm=TRUE))
      for (v in var2) {
        ir <- interarrival(x2, var=v, p.cut=p.cut[2])
        if (v == var2[1])
          ir2 <- ir
        else
          ir2 <- rbind(ir2, ir)
      }
      ir2$fact <- mkseas(ir2, width, start.day)
      wet2 <- eval(call("tapply", ir2$wet, ir2$fact, cent, na.rm=TRUE))
      dry2 <- eval(call("tapply", ir2$dry, ir2$fact, cent, na.rm=TRUE))
      l$wet.abs <- wet2 - wet1
      l$wet.rel <- wet2 / wet1
      l$dry.abs <- dry2 - dry1
      l$dry.rel <- dry2 / dry1
      ann.wet1 <- eval(call(cent, ir1$wet, na.rm=TRUE)) / dayspy1
      ann.wet2 <- eval(call(cent, ir2$wet, na.rm=TRUE)) / dayspy2
      ann.dry1 <- eval(call(cent, ir1$dry, na.rm=TRUE)) / dayspy1
      ann.dry2 <- eval(call(cent, ir2$dry, na.rm=TRUE)) / dayspy2
      l$ann.wet.abs <- ann.wet2 - ann.wet1
      l$ann.wet.rel <- ann.wet2 / ann.wet1
      l$ann.dry.abs <- ann.dry2 - ann.dry1
      l$ann.dry.rel <- ann.dry2 / ann.dry1
    }
  }
  class(l) <- "seas.change"
  l
}
