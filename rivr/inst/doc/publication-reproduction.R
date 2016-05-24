## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(echo = TRUE, eval = FALSE)

## ----the-code------------------------------------------------------------
#  ## ---- load-libs -------------------------------------------------------
#  library(knitr)
#  opts_chunk$set(fig.width = 6, fig.height = 4, fig.path = "fig/",
#    warning = FALSE, message = FALSE, error = FALSE, results = "asis",
#    out.width = '\\textwidth', cache = TRUE, cache.path = "cache/")
#  library(rivr)
#  library(ggplot2)
#  library(scales)
#  library(RColorBrewer)
#  library(dplyr)
#  library(xtable)
#  
#  ## ---- startparams -----------------------------------------------------
#  plotopts = list(theme_bw(), xlab(expression(Distance~from~control~section~~(ft))))
#  g = 32.2
#  Cm = 1.486
#  slope = 0.001
#  mannings = 0.045
#  flow = 250
#  width = 100
#  sideslope = 0
#  # calculate control depth as 1ft above the normal depth
#  depth.m1 = round(1 + normal_depth(slope, mannings, flow, 2, Cm,
#    width, sideslope), 2)
#  depth.m2 = round(1.1*critical_depth(flow, 1, g, width, sideslope), 2)
#  rivdist = 3000
#  
#  ## ---- restest ---------------------------------------------------------
#  # Test sensitivity of step size
#  stepsizes = c(500, 100, 50, 10, 1)
#  resolution.test.m1 = list()
#  resolution.test.m2 = list()
#  for(r in stepsizes){
#    resolution.test.m1[[paste0('dx=', r)]] = compute_profile(slope, mannings,
#      flow, depth.m1, Cm, g, width, sideslope, stepdist = r,
#      totaldist = rivdist)
#    resolution.test.m2[[paste0('dx=', r)]] = compute_profile(slope, mannings,
#      flow, depth.m2, Cm, g, width, sideslope, stepdist = r,
#      totaldist = rivdist)
#  }
#  resolution.plot = NULL
#  for(lbl in names(resolution.test.m1)){
#    f1 = resolution.test.m1[[lbl]]
#    f2 = resolution.test.m2[[lbl]]
#    f1['type'] = 'M1'
#    f2['type'] = 'M2'
#    f = rbind(f1, f2)
#    f['run'] = lbl
#    f['res'] = as.numeric(substr(lbl, 4, 10L))
#    resolution.plot = rbind(resolution.plot, f)
#  }
#  resolution.plot['res'] = factor(resolution.plot$res, levels = rev(sort(stepsizes)))
#  resolution.plot['run'] = factor(resolution.plot$run,
#    levels = names(resolution.test.m1))
#  resolution.plot['type'] = factor(resolution.plot$type, levels = c('M1', 'M2'))
#  stepsize = 10
#  
#  ## ----resfig ----------------------------------------------------------
#  ggplot(resolution.plot, aes(x = x, y = y + z)) + plotopts +
#  geom_line(aes(color = res, linetype = res), size = 1) +
#  facet_wrap(~type) + ylab(expression(River~stage~~(ft))) +
#  scale_linetype_manual("resolution (ft)", labels = levels(resolution.plot$res),
#    values = c("solid", "longdash", "dashed", "dotdash", "dotted")) +
#  scale_color_manual("resolution (ft)", values = brewer.pal(7, "YlGnBu")[3:7],
#    labels = levels(resolution.plot$res))
#  
#  ## ----roughness -------------------------------------------------------
#  roughness.test.m1 = list()
#  roughness.test.m2 = list()
#  roughnesses = 0.045*seq(0.5, 1.5, length = 5)
#  for(n in roughnesses){
#    roughness.test.m1[[paste0('n=',n)]] = compute_profile(slope, n, flow,
#      depth.m1, Cm, g, width, sideslope, stepdist=stepsize, totaldist=rivdist)
#    roughness.test.m2[[paste0('n=',n)]] = compute_profile(slope, n, flow,
#      depth.m2, Cm, g, width, sideslope, stepdist=stepsize, totaldist=rivdist)
#  }
#  roughness.plot = NULL
#  for(lbl in names(roughness.test.m1)){
#    f1 = roughness.test.m1[[lbl]]
#    f2 = roughness.test.m2[[lbl]]
#    f1['type'] = 'M1'
#    f2['type'] = 'M2'
#    f = rbind(f1, f2)
#    f['run'] = lbl
#    f['mannings'] = as.numeric(substr(lbl, 4, 10L))
#    roughness.plot = rbind(roughness.plot, f)
#  }
#  roughness.plot['run'] = factor(roughness.plot$run,
#    levels = names(roughness.test.m1))
#  roughness.plot['mannings'] = factor(roughness.plot$mannings,
#    levels = sort(roughnesses))
#  roughness.plot['type'] = factor(roughness.plot$type, levels=c('M1', 'M2'))
#  
#  ## ---- roughplot -------------------------------------------------------
#  ggplot(roughness.plot, aes(x = x, y = y + z)) +
#  geom_line(aes(linetype = mannings, color = mannings), size = 1) +  plotopts +
#  facet_wrap(~ type) + ylab(expression(River~stage~~(ft))) +
#  scale_linetype_manual("Bed roughness", labels = levels(roughness.plot$mannings),
#    values = c("solid", "longdash", "dashed", "dotdash", "dotted")) +
#  scale_color_manual("Bed roughness", labels = levels(roughness.plot$mannings),
#    values = brewer.pal(6, "Oranges")[2:6])
#  
#  ## ---- relroughplot ----------------------------------------------------
#  roughness.rel.m1 = list()
#  roughness.rel.m2 = list()
#  for(n in roughnesses){
#    thisyn = normal_depth(slope, n, flow, 2, Cm, width, sideslope)
#    roughness.rel.m1[[paste0('n=',n)]] = compute_profile(slope, n, flow,
#      1.25*thisyn, Cm, g, width, sideslope, stepdist=stepsize, totaldist=rivdist)
#    roughness.rel.m2[[paste0('n=',n)]] = compute_profile(slope, n, flow,
#      0.75*thisyn, Cm, g, width, sideslope, stepdist=stepsize, totaldist=rivdist)
#    roughness.rel.m1[[paste0('n=',n)]]['pd.yn'] =
#      (roughness.rel.m1[[paste0('n=',n)]]$y - thisyn)/thisyn
#    roughness.rel.m2[[paste0('n=',n)]]['pd.yn'] =
#      (roughness.rel.m2[[paste0('n=',n)]]$y - thisyn)/thisyn
#  }
#  roughness.rel.plot = NULL
#  for(lbl in names(roughness.rel.m1)){
#    f1 = roughness.rel.m1[[lbl]]
#    f2 = roughness.rel.m2[[lbl]]
#    f1['type'] = 'M1'
#    f2['type'] = 'M2'
#    f = rbind(f1, f2)
#    f['run'] = lbl
#    f['mannings'] = as.numeric(substr(lbl, 4, 10L))
#    roughness.rel.plot = rbind(roughness.rel.plot, f)
#  }
#  roughness.rel.plot['mannings'] = factor(roughness.plot$mannings,
#    levels = sort(roughnesses))
#  roughness.rel.plot['run'] = factor(roughness.rel.plot$run,
#    levels = names(roughness.rel.m1))
#  roughness.rel.plot['type'] = factor(roughness.rel.plot$type,
#    levels = c('M1', 'M2'))
#  ggplot(roughness.rel.plot, aes(x = x, y = pd.yn)) +
#  geom_line(aes(linetype = mannings, color = mannings), size = 1) +
#  plotopts + facet_wrap(~type, scales = 'free_y') +
#  scale_y_continuous(expression(Percent~difference~from~normal~depth), labels = percent) +
#  scale_linetype_manual("Bed roughness", labels = levels(roughness.rel.plot$mannings),
#    values = c("solid", "longdash", "dashed", "dotdash", "dotted")) +
#  scale_color_manual("Bed roughness", labels = levels(roughness.rel.plot$mannings),
#    values = brewer.pal(6, "Oranges")[2:6])
#  
#  
#  ## ---- loadup-kwm ------------------------------------------------------
#  oldscipen = options('scipen')
#  options(scipen = 1000)
#  
#  slope = 0.001
#  extent = 150000
#  mannings = 0.045
#  B = 100
#  SS = 0
#  Cm = 1.486
#  g = 32.2
#  iflow = 250
#  
#  # keep Courant number at 0.7 to balance temporal and spatial resolution
#  idepth = normal_depth(slope, mannings, 250, 10, 1.485, 100, SS)
#  iarea = channel_geom(idepth, B, SS)["A"]
#  cn = 0.7*iarea/iflow
#  # define upstream boundary condition assuming a timestep in seconds
#  bcfunc = function(x)
#    ifelse(x < 9000, 250 + (750/pi)*(1 - cos(pi*x/(60*75))), 250)
#  bctime = 76000
#  xnodes = extent/c(25, 50, 125, 250, 500, 1000, 2000, 5000) + 1
#  myxres = extent/(xnodes - 1)
#  mytres = cn*myxres
#  plotopts = list(theme_bw())
#  
#  # run the model
#  modtime = list()
#  modelresults = list()
#  for(i in seq(length(xnodes))){
#    numnodes = xnodes[i]
#    xstep = myxres[i]
#    tstep = mytres[i]
#    bc = bcfunc(seq(0, bctime, by=tstep))
#    mp = c(1, as.integer(50000/xstep + 1), as.integer(100000/xstep + 1), numnodes)
#    mt = as.integer(round(seq(1, length(bc), length.out=10)))[seq(7)]
#    if(xstep == min(myxres)){
#      mt = as.integer(round(seq(1, length(bc), length.out=125)))
#      mtslice = mt[round(seq(1, length(mt), length.out=7))]
#    }
#    lbl = paste0('dx=',xstep)
#    modtime[[lbl]] <- system.time(modelresults[[lbl]] <-
#      route_wave(slope, mannings, Cm, g, B, SS, iflow, bc,
#  	timestep=tstep, spacestep=xstep, numnodes=numnodes,
#  	  monitor.nodes=mp, monitor.times=mt, engine="Kinematic"))
#    modelresults[[lbl]]['deltax'] = xstep
#    modelresults[[lbl]]['deltat'] = tstep
#    modelresults[[lbl]]['computationtime'] = modtime[[lbl]][[3]]
#    modelresults[[lbl]]['label'] = lbl
#  }
#  allresults = do.call(rbind.data.frame, modelresults)
#  row.names(allresults) = NULL
#  names(allresults)[which(names(allresults)=='flow')] = "Q"
#  names(allresults)[which(names(allresults)=='distance')] = "x"
#  names(allresults)[which(names(allresults)=='time')] = "t"
#  
#  ## ---- floodhydrograph -------------------------------------------------
#  floodhydro = data.frame(time = seq(15000), flow=bcfunc(seq(15000)))
#  ggplot(floodhydro, aes(x=time/60, y=flow)) + geom_line(color="darkblue",
#    size = 1) + plotopts +
#  xlab(expression(time~~(minutes))) + ylab(expression(flow~~(ft^3~~s^-1)))
#  
#  ## ---- routeplot-kinematic ---------------------------------------------
#  # still plot
#  surfaces = group_by(filter(allresults, monitor.type=='timestep',
#    deltax==min(deltax)))
#  ggplot(filter(surfaces, step %in% mtslice[2:5]),
#    aes(x=x, y=Q, linetype=factor(round(t/60)), color=factor(round(t/60)))) +
#  geom_line(size = 1) + plotopts +
#  scale_color_manual("time\n(minutes)", values = brewer.pal(7, "YlGnBu")[3:6]) +
#  scale_linetype_manual("time\n(minutes)", values = c("solid",
#    "dashed", "dotdash", "dotted")) +
#  xlab(expression(distance~~downstream~~(ft))) +
#  ylab(expression(flow~~(ft^3~~s^-1)))
#  # animation
#  routeylims = c(min(surfaces$Q), max(surfaces$Q))
#  routexlims = c(min(surfaces$x), max(surfaces$x))
#  thisdt = unique(surfaces$deltat)
#  for(n in unique(surfaces$step))
#    print(
#      ggplot(filter(surfaces, step==n), aes(x=x, y=Q)) +
#      geom_line(color = "darkblue", size = 1) + plotopts +
#      scale_linetype_discrete("time\n(minutes)") +
#      scale_x_continuous(expression(distance~~downstream~~(ft)),
#        limits=routexlims) +
#  	scale_y_continuous(expression(flow~~(ft^3~~s^-1)), limits=routeylims) +
#  	ggtitle(paste(format(round((n - 1)*thisdt/60), width = 4), "minutes"))
#    )
#  
#  ## ---- xsectionplot-kinematic ------------------------------------------
#  xsections = filter(allresults, monitor.type == "node",
#    x %in% c(0,50000))[c("x","t","Q", "deltax")]
#  xsections['deltax'] = factor(paste(xsections$deltax, "feet"),
#    levels=paste(sort(unique(xsections$deltax)), "feet"))
#  ggplot(xsections, aes(x=t/60, y=Q, linetype=factor(format(x, big.mark=",",
#    trim=TRUE), levels = format(sort(unique(x)), big.mark=",", trim=TRUE)),
#    color = deltax)) +
#  geom_path(size = 1) + facet_grid(deltax~.) + plotopts +
#  scale_y_continuous(expression(flow~~(ft^3~~s^-1)), breaks=c(300, 500,700)) +
#  scale_linetype_manual("downstream\ndistance (ft)", values = c("solid",
#    "dashed", "dotted")) +
#  scale_x_continuous(expression(time~~(min)), limits=c(0,600)) +
#  scale_color_manual(values = c("#7570B3", "#8DA0CB", "#66C2A5", "#A6D854",
#    "#E5C494", "#E6AB02", "#FC8D62", "#E78AC3"), guide = FALSE)
#  
#  ## ---- peakstable-kinematic --------------------------------------------
#  peaks = as.data.frame(summarize(group_by(filter(allresults, x==50000,
#    monitor.type=='node'), computationtime, deltax, deltat, monitor.type),
#    peak.flow=max(Q), time.to.peak=t[which(Q==peak.flow)]))
#  real.peak = max(floodhydro$flow)
#  peaks['flow.percent.error'] = (peaks$peak.flow - real.peak)/real.peak
#  tbld = as.data.frame(peaks[c("deltax", "deltat", "flow.percent.error",
#    "computationtime")])
#  tbld["deltax"] = round(tbld$deltax)
#  tbld["deltat"] = round(tbld$deltat, 2)
#  tbld["flow.percent.error"] = round(100*tbld$flow.percent.error, 2)
#  names(tbld) = c("$\\Delta x$ (ft)", "$\\Delta t$ (s)", "\\% error (peak flow)",
#    "cost (s)")
#  ptable = xtable(tbld[order(tbld[,1]),], digits = c(0,0,2,2,2))
#  print(ptable, include.rownames = FALSE, sanitize.colnames.function = identity,
#    sanitize.text.function = function(x){x}, floating=FALSE, hline.after=NULL,
#    add.to.row=list(pos=list(-1,0, nrow(ptable)),
#    command=c('\\toprule ', '\\midrule ', '\\bottomrule ')),
#    format.args = list(big.mark = ","))
#  
#  ## ---- loadup-dwm ------------------------------------------------------
#  oldscipen = options('scipen')
#  options(scipen = 1000)
#  
#  slope = 0.001
#  extent = 150000
#  mannings = 0.045
#  B = 100
#  SS = 0
#  g = 32.2
#  Cm = 1.486
#  iflow = 250
#  
#  # define upstream boundary condition assuming a timestep in seconds
#  idepth = normal_depth(slope, mannings, 250, 10, Cm, 100, SS)
#  iarea = channel_geom(idepth, B, SS)[["A"]]
#  # keep Courant number at 0.06 to balance temporal and spatial resolution
#  cn = 0.06*iarea/250
#  bcfunc = function(x)
#    ifelse(x < 9000, 250 + (750/pi)*(1 - cos(pi*x/(60*75))), 250)
#  bctime = 76000
#  xnodes = extent/c(50, 125, 250, 500, 1000) + 1
#  myxres = extent/(xnodes - 1)
#  mytres = cn*myxres
#  plotopts = list(theme_bw())
#  
#  # run the model
#  modtime = list()
#  modelresults = list()
#  for(i in seq(length(xnodes))){
#    numnodes = xnodes[i]
#    xstep = myxres[i]
#    tstep = mytres[i]
#    bc = bcfunc(seq(0, bctime, by=tstep))
#    dc = rep(-1, length(bc))
#    mp = c(1, as.integer(c(1000, 50000, 100000, 149000)/xstep + 1), numnodes)
#    mt = as.integer(round(seq(1, length(bc), length.out=10)))[seq(7)]
#    if(xstep == min(myxres)){
#  	mt = as.integer(round(seq(1, length(bc), length.out=125)))	
#      mtslice = mt[round(seq(1, length(mt), length.out=7))]
#    }
#    lbl = paste0('dx=',xstep)
#    modtime[[lbl]] <- system.time(modelresults[[lbl]] <-
#      route_wave(slope, mannings, Cm, g, B, SS, iflow, bc, dc, timestep=tstep,
#        spacestep=xstep, numnodes=numnodes, monitor.nodes=mp, monitor.times=mt,
#  	  engine="Dynamic", boundary.type="QQ"))
#    modelresults[[lbl]]['deltax'] = xstep
#    modelresults[[lbl]]['deltat'] = tstep
#    modelresults[[lbl]]['computationtime'] = modtime[[lbl]][[3]]
#    modelresults[[lbl]]['label'] = lbl
#  }
#  # combine results
#  allresults = do.call(rbind.data.frame, modelresults)
#  row.names(allresults) = NULL
#  names(allresults)[which(names(allresults)=='flow')] = "Q"
#  names(allresults)[which(names(allresults)=='distance')] = "x"
#  names(allresults)[which(names(allresults)=='time')] = "t"
#  allresults['t'] = (allresults$step - 1)*allresults$deltat
#  allresults['x'] = (allresults$node - 1)*allresults$deltax
#  
#  ## ---- routeplot-characteristic ----------------------------------------
#  # still plot
#  surfaces = filter(allresults, monitor.type=='timestep',
#    deltax==min(deltax))
#  ggplot(filter(surfaces, step %in% mtslice[2:5]),
#    aes(x=x, y=Q, linetype=factor(round(t/60)), color=factor(round(t/60)))) +
#  geom_line(size = 1) + plotopts +
#  scale_color_manual("time\n(minutes)", values = brewer.pal(7, "YlGnBu")[3:6]) +
#  scale_linetype_manual("time\n(minutes)", values = c("solid",
#    "dashed", "dotdash", "dotted")) +
#    xlab(expression(distance~~downstream~~(ft))) +
#    ylab(expression(flow~~(ft^3~~s^-1)))
#  # animation
#  routeylims = c(min(surfaces$Q), max(surfaces$Q))
#  routexlims = c(min(surfaces$x), max(surfaces$x))
#  thisdt = unique(surfaces$deltat)
#  for(n in unique(surfaces$step))
#    print(
#      ggplot(filter(surfaces, step==n), aes(x=x, y=Q)) +
#        geom_line(color = "darkblue", size = 1) + plotopts +
#        scale_linetype_discrete("time\n(minutes)") +
#        scale_x_continuous(expression(distance~~downstream~~(ft)),
#          limits=routexlims) +
#  	  scale_y_continuous(expression(flow~~(ft^3~~s^-1)), limits = routeylims) +
#  	  ggtitle(paste(format(round((n - 1)*thisdt/60), width = 4), "minutes"))
#    )
#  
#  ## ---- xsectionplot-characteristic -------------------------------------
#  xsections = filter(allresults, monitor.type == "node",
#    x %in% c(0, 50000, 150000))[c("x","t","Q", "deltax")]
#  xsections['deltax'] = factor(paste(xsections$deltax, "feet"),
#    levels=paste(sort(unique(xsections$deltax)), "feet"))
#  data(waterolympics)
#  realdat = NULL
#  for(x in unique(xsections$deltax)){
#    thisdat = waterolympics[waterolympics$t > 200*60,]
#    thisdat['deltax'] = x
#    realdat = rbind(realdat, thisdat)
#  }
#  realdat['deltax'] = factor(realdat$deltax)
#  ggplot(xsections, aes(x=t/60, y=Q, linetype=factor(format(x, big.mark=",",
#    trim=TRUE), levels = format(sort(unique(x)), big.mark=",", trim=TRUE)),
#    color = deltax)) +
#  geom_path(size = 1) + geom_point(data = realdat[realdat$t > 200*60,], color = "black") +
#  facet_grid(deltax~.) +scale_y_continuous(expression(flow~~(ft^3~~s^-1)),
#    breaks=c(300, 500,700)) +
#  scale_linetype_manual("downstream\ndistance (ft)", values = c("solid",
#    "dashed", "dotted")) + plotopts +
#  scale_x_continuous(expression(time~~(min)), limits=c(0,1250)) +
#  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#A6D854", "#E5C494",
#    "#E6AB02"), guide = FALSE)
#  
#  ## ---- peakstable-characteristic ---------------------------------------
#  peaks = as.data.frame(summarize(group_by(filter(allresults, x==50000,
#    monitor.type=='node'), computationtime, deltax, deltat, monitor.type),
#    peak.flow=max(Q), time.to.peak=t[which(Q==peak.flow)]))
#  real.peak.flow = max(waterolympics$Q)
#  real.time.to.peak = waterolympics$t[which(waterolympics$Q==real.peak.flow)][1]
#  peaks['flow.percent.error'] =
#    (peaks$peak.flow - real.peak.flow)/real.peak.flow
#  peaks['time.percent.error'] =
#    (peaks$time.to.peak - real.time.to.peak)/real.time.to.peak
#  tbld = as.data.frame(peaks[c("deltax", "deltat", "flow.percent.error",
#    "time.percent.error", "computationtime")])
#  tbld["flow.percent.error"] = round(100*tbld$flow.percent.error, 2)
#  tbld["time.percent.error"] = round(100*tbld$time.percent.error, 2)
#  names(tbld) = c("$\\Delta x$ (ft)", "$\\Delta t$ (s)", "\\% error (peak flow)",
#    "\\% error (time to peak)", "cost (s)")
#  ptable = xtable(tbld[order(tbld[,1]),], digits = c(0,0,2,2,2,2))
#  print(ptable, include.rownames=FALSE, sanitize.colnames.function=identity,
#    sanitize.text.function = function(x){x}, floating=FALSE, hline.after=NULL,
#    add.to.row=list(pos=list(-1,0, nrow(ptable)),
#    command=c('\\toprule ', '\\midrule ', '\\bottomrule ')),
#    format.args = list(big.mark = ","))
#  
#  ## ---- loadup-boundaries -----------------------------------------------
#  oldscipen = options('scipen')
#  options(scipen = 1000)
#  
#  plotopts = list(theme_bw())
#  slope = 0.00008
#  extent = 5000
#  mannings = 0.013
#  B = 6.1
#  SS = 1.5
#  g = 9.81
#  Cm = 1
#  
#  ic = 126
#  id = 5.79
#  ia = channel_geom(id, B, SS)[["A"]]
#  CN = 0.9
#  
#  dx = 10
#  dt = round(dx*CN/(ic/ia + sqrt(id*g)), 2)
#  numnodes = extent/dx + 1
#  
#  bctime = 2000
#  bc = rep(id, round(bctime/dt) + 1)
#  dc = rep(0, length(bc))
#  dt = round(bctime/(length(bc) - 1), 2)
#  
#  CN = dt*(ic/ia + sqrt(id*g))/dx
#  
#  mp = c(1, as.integer(c(1500, 2500, 3000, 5000)/dx + 1))
#  mt = as.integer(round(seq(0, length(bc)-1, by=25))) + 1L
#  mtslice = c(1, as.integer(c(500, 1000, 1500, 2000)/dt + 1))
#  
#  dclose.lax = route_wave(slope, mannings, Cm, g, B, SS, ic, bc, dc,
#    timestep=dt, spacestep=dx, numnodes=numnodes, monitor.nodes=mp,
#    monitor.times=mt, engine="Dynamic", scheme="Lax", boundary.type="yQ")
#  dclose.mac = route_wave(slope, mannings, Cm, g, B, SS, ic, bc, dc,
#    timestep=dt, spacestep=dx, numnodes=numnodes, monitor.nodes=mp,
#    monitor.times=mt, engine="Dynamic", scheme="MacCormack", boundary.type="yQ")
#  dclose.lax['scheme'] = "Lax diffusive"
#  dclose.mac['scheme'] = "MacCormack predictor-corrector"
#  dclose = rbind(dclose.lax, dclose.mac)
#  row.names(dclose) = NULL
#  dclose["CN"] = dt*(dclose$flow/dclose$area + sqrt(dclose$depth*g))/dx
#  
#  ## ---- through-time ----------------------------------------------------
#  dclose.times = filter(dclose, monitor.type=="timestep")
#  ggplot(filter(dclose.times, step %in% mtslice),
#    aes(x=distance, y=depth, linetype = factor(round(time)),
#      color = factor(round(time)))) + geom_path(size = 1) +
#  scale_y_continuous(expression(depth~~(ft))) + facet_wrap(~scheme) +
#  scale_x_continuous(expression(distance~~downstream~~(ft))) + plotopts +
#  scale_color_manual("time (s)", values = brewer.pal(7, "YlGnBu")[3:7]) +
#  scale_linetype_manual("time (s)", values = c("solid", "longdash", "dashed",
#    "dotdash", "dotted"))
#  # animation
#  routeylims = c(min(dclose.times$depth), max(dclose.times$depth))
#  routexlims = c(min(dclose.times$distance), max(dclose.times$distance))
#  for(n in unique(dclose.times$step))
#    print(
#      ggplot(filter(dclose.times, step==n), aes(x=distance, y=depth)) +
#      geom_line(color = "darkblue", size = 1) + plotopts +
#      scale_linetype_discrete("time\n(minutes)") +
#      scale_x_continuous(expression(distance~~downstream~~(ft)),
#        limits=routexlims) +
#  	scale_y_continuous(expression(depth~~(ft)), limits=routeylims) +
#  	ggtitle(paste(format(round((n - 1)*dt), width = 4), "seconds")) +
#      facet_wrap(~scheme)
#    )

