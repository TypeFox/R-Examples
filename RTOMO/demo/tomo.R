

#####  load up the Mount St. Helens tomography
data(HELMOD)

###  load up the map information (GEOmap structure)

data(HELMAP)

#####  show whole model layers 1-15

SHOWTOMO(HELMOD, MAP=HELMAP,  bkgr="beige", I=1, J=15)

readline("To Continue Hit Enter Key\n")


################   show 4 layers
par(mfrow=c(2,2))
for(i in 6:9)
{
##  i = 4
FANCY.TOMO(HELMOD, i, MAP=HELMAP, bkgr="beige")
}

readline("To Continue Hit Enter Key\n")



par(mfrow=c(1,1))


FANCY.TOMO(HELMOD, 7, MAP=HELMAP, bkgr="beige")
readline("To Continue Hit Enter Key\n")


L=list()
L$x=c(3.4332423635852,24.5907738623029)
L$y=c(11.1625031193420,13.2100710944510)

points(L)
lines(L)

readline("To Continue Hit Enter Key\n")

dev.new()

### get(getOption("device"))()

TOMOXSEC(HELMOD, L$x[1],  L$y[1],  L$x[2],  L$y[2], zmax=20, depth = c(-25, 0),
    COL = tomo.colors(100), LIM = NULL, STA = NULL, PLOT = TRUE)



###################
