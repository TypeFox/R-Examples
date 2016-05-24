### R code from vignette source 'landscape_simulation.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: landscape_simulation.Rnw:23-28
###################################################
library(MetaLandSim)

rl <- rland.graph(mapsize = 1000, dist_m = 60, 
                  areaM = 0.5, areaSD = 0.2, Npatch =70, 
                  disp = 100, plotG = TRUE)


###################################################
### code chunk number 2: landscape_simulation.Rnw:33-40
###################################################
library(MetaLandSim)

#The occupation of a landscape is simulated by:
sp_t0 <- species.graph(rl=rl, method="percentage", parm=50, 
                       nsew="none", plotG=TRUE)

names(sp_t0)


###################################################
### code chunk number 3: landscape_simulation.Rnw:46-65
###################################################
data(param1)

sp_t1 <- spom(
sp_t0,
kern="op1",
conn="op1",
colnz="op1",
ext="op1",
param_df=param1,
beta1=NULL,
b=1,
c1=NULL,
c2=NULL,
z=NULL,
R=NULL
)

#Which has the following elements:
names(sp_t1)


###################################################
### code chunk number 4: landscape_simulation.Rnw:74-109
###################################################

#Loading species parameters

data(param1)

#Simulating occupation in dynamic landscape

it1 <- iterate.graph(
iter = 2, 
mapsize = 1000, 
dist_m = 30, 
areaM = 0.5, 
areaSD= 0.1,
Npatch = 200, 
disp = 800, 
span = 100, 
par1 = "stoc", 
par2 = 2,
par3 = 2, 
method = "percentage",
parm = 50, 
nsew = "none", 
succ = "none", 
param_df = param1, 
kern = "op1", 
conn = "op1", 
colnz = "op1",
ext = "op1", 
b = 1, 
graph  = FALSE
)

#This file is composed by the following elements: 

names(it1)


