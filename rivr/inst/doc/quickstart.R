## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, error=FALSE, dev='svg', fig.show='hold')

## ----ynyc----------------------------------------------------------------
require(rivr)
flow = 250; mannings = 0.045 ; Cm = 1.486; gravity = 32.2
width = 100; slope = 0.001; sideslope = 0

yn = normal_depth(slope, mannings, flow, yopt = 2, Cm, width, sideslope)
yc = critical_depth(flow, yopt = 2, gravity, width, sideslope)

print(c(normal.depth = yn, critical.depth = yc))

## ----gvf-----------------------------------------------------------------
flow = 250; mannings = 0.045 ; Cm = 1.486; gravity = 32.2
width = 100; slope = 0.001; sideslope = 0

gvf = compute_profile(slope, mannings, flow, y0 = 2.7, Cm, gravity, width, 
  sideslope, stepdist=50, totaldist=3000)

## ----plot-gvf------------------------------------------------------------
require(ggplot2)
ggplot(gvf, aes(x = x, y = y + z)) + geom_line(color='blue')
# or try the default plot method
# plot(gvf)

## ----uf------------------------------------------------------------------
baseflow = 250; mannings = 0.045 ; Cm = 1.486; gravity = 32.2
width = 100; slope = 0.001; sideslope = 0

numnodes = 301; xresolution = 250; tresolution = 10; 
times = seq(0, 30000, by = tresolution)
wave = ifelse(times >= 9000, baseflow,
  baseflow + (750/pi)*(1 - cos(pi*times/(60*75))))
downstream = rep(-1, length(wave))

mn = c(1, 101, 201)
mt = c(501, 1501, 3001)
uf = route_wave(slope, mannings, Cm, gravity, width, sideslope, 
  baseflow, wave, downstream, tresolution, xresolution, numnodes, 
  mn, mt, "Dynamic", "MacCormack", "QQ") 

## ----uf-access-----------------------------------------------------------
require(dplyr)
uf.nodes = filter(uf, monitor.type == "node")
ggplot(uf.nodes, aes(x=time, y=flow, color=factor(distance))) + geom_line()
uf.times = filter(uf, monitor.type == "timestep")
ggplot(uf.times, aes(x=distance, y=flow, color=factor(time))) + geom_line()
# or try the default plot method
# plot(uf)

