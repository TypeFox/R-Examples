## Eh-pH diagrams for copper-water-glycine
## After Fig. 2 of Aksu and Doyle, 2001
## (Aksu, S. and Doyle, F. M., 2001. Electrochemistry of copper in aqueous glycine 
## solutions. J. Electrochem. Soc., 148, B51-B57. doi:10.1149/1.1344532)

# add some new species to thermo$obigt
m1 <- makeup(info(c("Cu+", "glycinate", "glycinate")), sum=TRUE)
mod.obigt(name="Cu(Gly)2-", formula=as.chemical.formula(m1))
m2 <- makeup(info(c("Cu+2", "glycinate", "H+")), sum=TRUE)
mod.obigt(name="HCu(Gly)+2", formula=as.chemical.formula(m2))
# Gibbs energies from A&D's Table 1 and Table II
Cu_s <- c("copper", "cuprite", "tenorite")
Gly <- c("glycinium", "glycine", "glycinate")
Cu_aq <- c("Cu+", "Cu+2", "CuO2-2", "HCuO2-")
CuGly <- c("Cu(Gly)+", "Cu(Gly)2", "Cu(Gly)2-", "HCu(Gly)+2")
names <- c(Cu_s, Gly, Cu_aq, CuGly)
G <- c(
  convert(c(0, -146, -129.7,
  -384.061, -370.647, -314.833,
  49.98, 65.49, -183.6, -258.5, -298.2)*1000, "cal"),
  convert(c(15.64, 10.1, 2.92), "G"))
# run updates in order so later species take account of prev. species' values
getG <- function(x) info(info(x))$G
for(i in 1:length(G)) {
  myG <- G[i]
  if(i==12) myG <- myG + getG("Cu+2") + 2*getG("glycinate")
  if(i==13) myG <- myG + getG("Cu+") + 2*getG("glycinate")
  if(i==14) myG <- myG + getG("Cu(Gly)+")
  mod.obigt(names[i], G=myG)
}  

# in Fig. 2b, total log activities of Cu (Cu_T) and glycine (L_T) are -4 and -1
basis(c("Cu+2", "H2O", "H+", "e-", "glycinium", "CO2"), c(999, 0, 999, 999, -1, 999))
# add solids and aqueous species
species(Cu_s)
species(c(Cu_aq, CuGly), -4)
names <- c(Cu_s, Cu_aq, CuGly)
# mosaic diagram with to speciate glycine as a function of pH
m <- mosaic(bases=Gly, pH=c(0, 16, 300), Eh=c(-0.6, 1.0, 300))
fill <- c(rep("lightgrey", 3), rep("white", 4), rep("lightblue", 4))
d <- diagram(m$A.species, fill=fill, names=NULL, tplot=FALSE, xaxs="i", yaxs="i")
# to make the labels look nicer
names <- names[sort(unique(as.numeric(d$predominant)))]
for(i in 1:length(names)) {
  if(i %in% 1:3) lab <- names[i] else lab <- expr.species(names[i])
  # some manual adjustment so labels don't collide
  srt <- dy <- dx <- 0
  if(names[i]=="tenorite") dy <- -0.1
  if(names[i]=="CuO2-2") dy <- -0.1
  if(names[i]=="HCu(Gly)+2") srt <- 90
  if(names[i]=="HCu(Gly)+2") dx <- -0.2
  if(names[i]=="Cu(Gly)+") srt <- 90
  text(d$lx[i]+dx, d$ly[i]+dy, lab, srt=srt)
}

# add glycine ionization lines
d <- diagram(m$A.bases, add=TRUE, col="darkblue", dotted=c(2, 3), names=NULL)
text(d$lx, -0.5, Gly, col="darkblue")

# add water lines and title and re-draw a box around the plot
# because the filling of fields masks it
water.lines()
box()
mtitle(expression("Copper-water-glycine at 25"~degree*"C and 1 bar",
  "After Aksu and Doyle, 2001 (Fig. 2b)"), line=0.5)

# done!
data(thermo)
