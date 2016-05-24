#a) plot 3-d varios
#b) plot upscaled point variograms

errorBar <- function(x, y, upper, lower=upper, length=0.1,...){
 if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
 stop("vectors must be same length")
 arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

checkVario.rtop = function(object,  acor = 1, log = "xy", cloud = FALSE, gDist = TRUE, params = list(), ...) {
  params = getRtopParams(object$params, newPar = params, ...)
  dots = list(...)
  askpar = par("ask")
  if (dev.interactive()) par("ask" = TRUE) else par("ask" = FALSE)
  variogramModel = object$variogramModel
# variogram preferred as sampleVariogram, if not variogramCloud (if exisiting, NULL otherwise)  
  sampleVariogram = object$variogram
  if (is.null(sampleVariogram)) sampleVariogram = object$variogramCloud
  observations = object$observations
  formulaString = object$formulaString
  amul = object$params$amul
  varFit = object$varFit
  abins = adfunc(NULL, observations, amul)
  observations$acl = findInterval(observations$area, abins)
  observations$n = 1
  obsvar = aggregate(observations@data[,as.character(formulaString[[2]])], 
           by = list(acl = observations$acl), FUN = var)
  obsvar$area = aggregate(observations$area, by = list(acl = observations$acl), FUN = mean)[,2]*acor
  obsvar$n = aggregate(observations$n, by = list(acl = observations$acl), FUN = sum)[,2]
  obsvar$n = obsvar$n/max(obsvar$n)*20
  obsvar = obsvar[!is.na(obsvar$x),]
  plot(obsvar$area, obsvar$x, xlab = "area", ylab = "variance", 
      cex = sqrt(obsvar$n), pch = 16, log = log, ylim = c(min(obsvar$x)/5, max(obsvar$x)*5))

  print(paste("cloud is ",cloud))
  if (cloud | is(sampleVariogram, "rtopVariogramCloud")) {
    print("Creating cloud variogram; this might take some time")
    if (!is(sampleVariogram, "rtopVariogramCloud")) {
      if (!("variogramCloud" %in% names(object))) 
        object$variogramCloud = rtopVariogram(observations, formulaString, params, cloud = TRUE)
      clvar = object$variogramCloud
    } else clvar = sampleVariogram
    if (gDist) {
      if (!("gdistObs" %in% names(object))) {
        if (!("dObs" %in% names(object))) 
          object$dObs = rtopDisc(observations, params = params)
        dObs = object$dObs
        object$gDistObs = gDist(dObs, dObs, params = params)
      } 
      gdists = object$gDistObs
      gDiag = diag(gdists)
      clvar$gdist = 0
      for (ic in 1:dim(clvar)[1]) {
        ia = clvar$acl1[ic]
        ib = clvar$acl2[ic]
        clvar$dist[ic] = gdists[ia,ib] - 0.5*(gDiag[ia] + gDiag[ib])
      }
    } 
    clvar$np = clvar$ord  
    if (!"identify" %in% names(dots) | !dev.interactive()) dots$identify = FALSE
    print(plot(clvar, xlab = "distance", unlist(dots)))      
  }
  if (!is.null(varFit) & is(sampleVariogram, "rtopVariogram")) {
    gammar = varFit[,c("np", "gamma", "gammar")]
    gammar$nnp = sqrt(gammar$np)/max(sqrt(gammar$np))*20
    gmax = max(gammar[,c("gamma", "gammar")])
    gmin = quantile(c(gammar$gammar, gammar$gamma), 0.05)
    nnp = 0 # dummy variable to avoid check warning
    plot(gammar ~ gamma, gammar, xlim = c(ifelse(length(grep("x",log)) >0, gmin,0), gmax), 
            ylim = c(ifelse(length(grep("x",log)) >0, gmin,0),gmax), cex = sqrt(nnp),
      xlab = "gamma", ylab = "gamma regularized", log = log)
    abline(0,1)  
  } else if (!is.null(varFit) & is(sampleVariogram, "rtopVariogramCloud")) {
    gammar = varFit[,c("np", "gamma", "gammar")]
    gammar = gammar[order(gammar$gamma),]
    ng = dim(gammar)[1]  
    groups = ifelse(ng > 200, 20, ng/10)
    npg = ng/groups
    gammar$group = c(1:ng) %/% npg
    ngammar = aggregate(list(gamma = gammar$gamma, gammar = gammar$gammar), 
          by = list(gammar$group), FUN = mean)
    ngammar = cbind(ngammar,
       aggregate(list(gammav = gammar$gamma, gammarv = gammar$gammar), 
       by = list(gammar$group), FUN = var))
    xmax = max(c(ngammar$gamma, ngammar$gammar))
    plot(gammar ~ gamma, ngammar, 
      xlab = "regularized gamma", ylab = "gamma", xlim = c(0,xmax), ylim = c(0,xmax), pch = 16)
    errorBar(ngammar$gammar, ngammar$gamma, upper = sqrt(ngammar$gammav))
    abline(0,1)
  }
  if (is.null(variogramModel)) {
    if (is.null(sampleVariogram)) sampleVariogram = rtopVariogram(observations)
    checkVario(sampleVariogram, observations, params = object$params, log = log, ...)    
  } else {
    if (is.null(sampleVariogram)) {
      object$checkVario = checkVario(object$variogramModel, observations = object$observations, 
           params = object$params, acor = acor, log = log, ...)
    } else {
      object$checkVario = checkVario(object$variogramModel, sampleVariogram = sampleVariogram, 
          observations = object$observations, params = object$params, acor = acor, 
          log = log, ...) 
    }
  }
  par(ask = askpar)
  invisible(object)
}






checkVario.rtopVariogramModel = function(object, 
           sampleVariogram = NULL, observations = NULL, areas = NULL, dists = NULL, acomp = NULL, 
           params = list(), compVars = list(), acor = 1, log = "xy", legx = NULL, legy = NULL, 
           plotNugg = TRUE, ...) {
variogramModel = object
params = getRtopParams(params, ...)
askpar = par("ask")
if (dev.interactive()) par("ask" = TRUE) else par("ask" = FALSE)

if (is.null(areas)) areas = params$amul
if (is.null(dists)) dists = params$dmul

if (length(areas) == 1) areas = adfunc(sampleVariogram, observations, areas)
if (length(dists) == 1) dists = dfunc(sampleVariogram, observations, dists)
 
Srl = list()
icomb = 0
polylist = list()
aavg = areas[1:(length(areas)-1)]
dists = c(0,dists)
adists = dists[1:length(dists)]
for (iarea in 1:(length(areas)-1)) {
  area =  mean(c(areas[iarea], areas[iarea + 1]))
  aavg[iarea] = area
  for (idist in 1:(length(dists))) {
    icomb = icomb + 1
    ddist = ifelse(idist == 1, 0, mean(c(dists[idist],dists[idist-1])))
    if (iarea == 1) adists[idist] = ddist
    cs = sqrt(area)/2
    x1 = ddist -cs
    x2 = ddist + cs
    y1 = -cs
    y2 = cs                                           
    boun = cbind(x = c(x1,x2,x2,x1,x1),y = c(y1,y1,y2,y2,y1))
    polyBoun = Polygon(boun)
    Srl[[icomb]] = Polygons(list(polyBoun),ID = as.character(icomb))
  }
}

polys = SpatialPolygons(Srl)
vmats = list()
iplot = 0
na = length(areas)
if (is.null(acomp) | length(acomp) == 1) {
  if (is.null(acomp)) acomp = 5
  if (!is.null(sampleVariogram) & is(sampleVariogram,"rtopVariogram")) {
    samp = aggregate(sampleVariogram$np, 
        by = list(acl1 = sampleVariogram$acl1, acl2 = sampleVariogram$acl2), sum)
    if (acomp > dim(samp)[1]) acomp = dim(samp)[1]
    acomp = samp[order(samp$x, decreasing = TRUE)[1:acomp],1:2]
  } else {
    aacomp = expand.grid(a1 = c(1:(na-1)),a2 = c(1:(na-1)))
    aacomp = aacomp[aacomp$a1 >= aacomp$a2,]
    if (acomp > dim(aacomp)[1]) acomp = dim(aacomp)[1]
    acomp = aacomp[sample(dim(aacomp)[1],acomp),]
  }
} else {
  acomp = acomp[acomp$acl1 < length(areas) & acomp$acl2 < length(areas),]
}
vmats = matrix(0, ncol = length(dists), nrow = dim(acomp)[1])
for (iplot in 1:dim(acomp)[1]) {
  i1 = acomp[iplot,2]
  i2 = acomp[iplot,1]
  ld = length(adists)
  poly1 = polys[unique(c((i1-1)*ld+1,((i2-1)*ld+1):(i2*ld)))]
  lobject = createRtopObject(SpatialPolygonsDataFrame(poly1, 
       data = data.frame(obs = c(1:length(poly1))), match.ID = FALSE), params = params,
       formulaString = obs~1)
  lobject$variogramModel = variogramModel
  nadists = adists
  if (i1 != i2) nadists = c(0, nadists)
  overlapObs = findVarioOverlap(data.frame(a1 = poly1[1]@polygons[[1]]@area,
       a2 = poly1[2]@polygons[[1]]@area, dist = nadists)) 
  lobject$overlapObs =  t(matrix(rep(overlapObs, ld+(i1 != i2)),ncol = ld+(i1 != i2)))
  vmat = varMat(lobject, cv = TRUE)$varMatObs
#  vmat = varMat(poly1,variogramModel = variogramModel, params = params)
#  vmat = vmat-diag(vmat)
#    vmats[iplot,] = vmat[1,2:ld]
#
  if (i1 == i2) {
    vmats[iplot,2:ld] = vmat[1,2:ld]
  } else {
    vmats[iplot,] = vmat[1,2:(ld+1)]
  }   
}


pvar = apply(as.matrix(adists),1, varioEx, variogramModel = variogramModel)+ 
   ifelse(plotNugg, nuggEx(1,variogramModel)*acor, 0)
ymin = max(min(vmats[vmats > 0]),min(sampleVariogram$gamma))
ymax = max(pvar)
if (inherits(sampleVariogram, "rtopVariogramCloud")) {
  xmin = min(sampleVariogram$dist)/1.3
} else {
  xmin = min(sampleVariogram$dist[sampleVariogram$np > 2]/1.3)
}
xmax = max(adists)
if (acor != 1) {
  Rver = R.Version()
  if (as.numeric(Rver$major)*100 + as.numeric(Rver$minor) >= 214) {
    xTicks = axTicks(1,c(xmin,xmax,3),usr = c(log10(xmin),log10(xmax)), 
             log=TRUE, nintLog = Inf)
  } else xTicks = axTicks(1,c(xmin,xmax,3),usr = c(log10(xmin),log10(xmax)), log=TRUE)
  xlabs = xTicks*sqrt(acor)
} else {
  xTicks = NULL
  xlabs = TRUE
}
plot(adists,pvar,ylim = c(ymin, ymax), xlim = c(xmin, xmax), log = log, 
    type="l", col = "black", lwd = 2, ylab = "gamma", 
    xlab = "distance", xaxt = "n")
axis(1,at = xTicks, labels = xlabs)

legende = list(text = "point", col = c("black"), lty = c(1), pch = 16)

bcols = c("red", "blue", "green", "orange", "brown", "violet", "yellow", "salmon")
cols1 = bcols[1:length(areas)]
cols2 = bcols[1:dim(acomp)[1]]
for (iplot in 1:dim(acomp)[1]) {
  i1 = acomp[iplot,2]
  i2 = acomp[iplot,1]
  ld = length(adists)
  if (i1 == i2) {
    lt = 1
    lcol = cols1[i1]
  } else {
    lt = 2
    lcol = cols2[iplot]
  }
  lines(adists,vmats[iplot,1:ld],lty = lt, lwd = 2, col = lcol)
  legende$text = c(legende$text, paste(aavg[i1]*acor, "vs", aavg[i2]*acor))
  legende$col = c(legende$col, lcol)
  legende$lty = c(legende$lty, lt)
  if (!is.null(sampleVariogram) & is(sampleVariogram, "rtopVariogram")) {
    ppts = sampleVariogram[sampleVariogram$acl2 == i1 & sampleVariogram$acl1 == i2,]
    lpch = 16+lt
    np = 0 # Dummy variable to avoid check warning
    points(gamma ~ dist,ppts,col = lcol, pch = lpch, cex = sqrt(sqrt(np/max(sampleVariogram$np)*10)) )
    legende$pch = c(legende$pch, lpch)
  }
}
if (length(compVars) > 0) {
  for (ic in 1: length(compVars)) {
    cvar = compVars[ic]
    clines = variogramLine(cvar[[1]], dist_vector = adists)
    lines(clines, lty = 3, lwd = 2, col = cols2[ic])
    legende$text = c(legende$text,names(cvar))
    legende$col = c(legende$col, cols2[ic])
    legende$lty = c(legende$lty, 3)
    legende$pch = c(legende$pch, 26)
  }
}
if (is.null(legx)) legx = ifelse(length(grep("x", log)) > 0, max(adists)/log(xmax/xmin,5), max(adists)*0.7)
if (is.null(legy)) legy = ifelse(length(grep("y", log)) > 0, sqrt(ymin*ymax/1.5), ymax*0.35)
warn = options("warn")
options(warn = -1)
legend(legx, legy, legende$text, col = legende$col, 
    lty = legende$lty, lwd = rep(2, length(legende$pch)), pch = legende$pch, merge = TRUE)

checkVarioRes = list(vmats = rbind(vmats, pvar), acomp = acomp)
par(ask = askpar)
options(warn = warn$warn)
invisible(checkVarioRes)
}


checkVario.rtopVariogram = function(object, ...) {

}

checkVario.rtopVariogramCloud = function(object, ...) {


}

