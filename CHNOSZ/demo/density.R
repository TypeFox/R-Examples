# make T-P diagram for H2O, colored according to density

# IAPWS95 or SUPCRT92
method <- "IAPWS95"
# low or high T,P range
TPrange <- "low"

blue <- "blue"
if(TPrange=="low") {
  T <- seq(300, 800, 10)
  #P <- seq(1, 600, 11.98)
  P <- seq(0, 600, 12)
  bias <- 1.68
} else {
  # upper T,P limit for SUPCRT92: 2250 degC, 30000 bar
  T <- seq(273, 2523, 50)
  P <- seq(1, 30000, length.out=50)
  # to attempt to match the colors using the different methods
  # (ranges are different because IAPWS95 reports higher density in
  #  the high-P, low-T region, where SUPCRT92 doesn't give output)
  if(method=="IAPWS95") bias <- 2.2
  else if(method=="SUPCRT92") {
    bias <- 2.1
    blue <- "#0d0dff"
  }
}
TP <- expand.grid(T=T, P=P)
if(method=="IAPWS95") {
  # the following should trigger parallel calculations
  # if nrow(TP) (5751 for TPrange="low") is >= thermo$opt$paramin (default 1000)
  rho <- palply("TP", 1:nrow(TP), function(i){CHNOSZ::rho.IAPWS95(TP$T[i], TP$P[i])})
} else if(method=="SUPCRT92") {
  rho <- water.SUPCRT92("rho", TP$T, TP$P)
  # water.SUPCRT92 returns 0 when the density can't be calculated
  rho[rho==0] <- NA
}
rho.num <- unlist(rho)
rho.mat <- matrix(rho.num, nrow=length(T), ncol=length(P))
# blueest for most dense, reddest for least dense
# bias is adjusted to white for the critical density
ncol <- 500
col <- colorRampPalette(c("red", "white", blue), bias=bias)(ncol)
# first make a background image (for debugging - 
# will be visible only if some density calculations fail)
fill.mat <- matrix(0, nrow=length(T), ncol=length(P))
image(T, P, fill.mat, col="black", xlab=axis.label("T", "K"), ylab=axis.label("P"), useRaster=TRUE, yaxt="n")
axis(2, at=c(1, seq(100, 600, 100)))
# now plot densities
image(T, P, rho.mat, col=col, add=TRUE, useRaster=TRUE)
# add a title and calculate saturation line
if(method=="IAPWS95") {
  title(main=expression("Density of"~H[2]*O~"inverted from IAPWS-95 equations"))
##  title(main=expression("Line calculated using auxiliary equations for saturation"), line=0.8)
  Psat <- convert(WP02.auxiliary("P.sigma", T), "bar")
} else if(method=="SUPCRT92") {
  title(main=expression("Density of"~H[2]*O~"calculated using SUPCRT92"))
  Psat <- water.SUPCRT92("Psat", T, "Psat")[,1]
}
### plot saturation line
##lines(T, Psat, lwd=6)
##lines(T, Psat, lwd=3, col="gold")
# add a color key
if(TPrange=="low") {
  x <- c(355, 395, 402)
  y <- c(170, 520)
} else if(TPrange=="high") {
  x <- c(600, 780, 800)
  y <- c(10000, 25000)
}
ykey <- seq(y[1], y[2], length.out=ncol+1)
for(i in 1:ncol) rect(x[1], ykey[i], x[2], ykey[i+1], col=col[i], border=NA)
rect(x[1], ykey[1], x[2], rev(ykey)[1])
# label the extrema
rrange <- range(rho.num, na.rm=TRUE)
text(x[3], ykey[1], as.expression(substitute(x~kg/m^3, list(x=round(rrange[1], 4)))), adj=0, col="white")
text(x[3], rev(ykey)[1], as.expression(substitute(x~kg/m^3, list(x=round(rrange[2], 4)))), adj=0, col="white")
# label the critical density
rlevels <- seq(rrange[1], rrange[2], length.out=ncol+1)
rho.critical <- 322
icrit <- which.min(abs(rlevels-rho.critical))
text(x[3], ykey[icrit], as.expression(substitute(x~kg/m^3~group("(", rho[c], ")"), list(x=rho.critical))), adj=0, col="white")

#if(method=="IAPWS95") {
#  # the saturation line is very accurate but not quite perfect;
#  # we can show whether it is on the liquid or vapor side
#  ina <- is.na(P.sigma)
#  rho.sigma <- rho.IAPWS95(T[!ina], P.sigma[!ina])
#  col <- rep("blue", length(rho.sigma))
#  col[rho.sigma < 322] <- "red"
#  points(T[!ina], P.sigma[!ina], col=col, pch=20)
#}

