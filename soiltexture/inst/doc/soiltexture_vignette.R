### R code from vignette source 'soiltexture_vignette.Rnw'

###################################################
### code chunk number 1: soiltexture_vignette.Rnw:99-111
###################################################
# Set a few Sweave options:
options( 
    width       = 65,  # width of R output
    prompt      = " ", # Sign preceding R input in R-GUI
    continue    = " "  # same, but after 2nd line
)   # 

# The working directory:
# setwd("C:/_RTOOLS/SWEAVE_WORK/SOIL_TEXTURES/rforge/pkg/soiltexture/inst/doc/INOUT") 

# And load the xtable package:
library( "xtable" )


###################################################
### code chunk number 2: soiltexture_vignette.Rnw:130-153
###################################################
old.wd  <- getwd() 

# setwd("C:/_RTOOLS/SWEAVE_WORK/SOIL_TEXTURES/rforge/pkg/soiltexture/inst/doc/INOUT") 

# if( !("soiltexture" %in%  as.character( installed.packages()[,1] )) ) 
# {   #
#     suppressMessages( 
#         install.packages( 
#             pkgs  = "soiltexture"  
#             # repos = "http://R-Forge.R-project.org" 
#         )   #
#     )   #
# }   #

suppressPackageStartupMessages( library( "soiltexture" ) )

# library( 
#     package        = "soiltexture", 
#     character.only = TRUE, 
#     quietly        = TRUE 
# )   

# setwd(old.wd) 


###################################################
### code chunk number 3: COVERFIG
###################################################
TT.plot(class.p.bg.col=T,class.sys="USDA.TT",main=NA)


###################################################
### code chunk number 4: soiltexture_vignette.Rnw:367-407
###################################################
bornes <- c(0,2,20,50,200,2e3,20e3)
noms   <- c("Cl","FiSi","CoSi","FiSa","CoSa","Gr","St")
txt.b  <- expression( 0*mu*m, 2*mu*m, 20*mu*m, 50*mu*m, 200*mu*m, 2*'mm', 2*'cm')

tmp <- data.frame(bornes,noms) # ,txt.b
#tmp$"txt.b" <- as.character(tmp$"txt.b")

par(  "mar"=c(4,1,1,1)+0.1  )  #  c(bottom, left, top, right)

plot( 
	x		= tmp$"bornes"[-1],  
	y		= rep(1,dim(tmp[-1,])[1]),  
	type	= "n",  
	main	= "",  
	xlab	= "Soil particule sizes",  
	ylab	= "",  
	yaxt	= "n",  xaxt = "n",  
	log		= "x",  
	xlim	= c(0.2,75e3), 
	bty		= "n", 
	cex.lab = 2  
)	#

abline(v=tmp$"bornes",lty=3,lwd=c(2,4,2,4,2,4,2))
abline(h=par("usr")[3:4],lty=1,lwd=4)

mtext( 
    text    = txt.b[-1], 
    side    = 1, 
    line    = rep( 
        c(0.5,1.25), 
        (dim(tmp)[1]-1)/2
    ),  #
    at  = tmp$"bornes"[-1], 
    cex = 2  
)   #

xtxt <- (tmp$"bornes"[1:(length(tmp$"bornes"))]+c(tmp$"bornes"[2:length(tmp$"bornes")],75e3))/2

text(x=xtxt,y=rep(1,length(xtxt)),labels=tmp$"noms",cex=2) 


###################################################
### code chunk number 5: soiltexture_vignette.Rnw:496-505
###################################################
TT.plot( 
    class.sys   = "none", 
    tri.data    = data.frame( 
        "CLAY"  = 45, 
        "SILT"  = 38, 
        "SAND"  = 17 
    ),  #
    main        = NA  
)   #


###################################################
### code chunk number 6: soiltexture_vignette.Rnw:543-552
###################################################
TT.plot( 
    class.sys   = "HYPRES.TT", 
    tri.data    = data.frame( 
        "CLAY"  = 45, 
        "SILT"  = 38, 
        "SAND"  = 17 
    ),  
    main        = NA  
)   


###################################################
### code chunk number 7: soiltexture_vignette.Rnw:557-558
###################################################
library( "xtable" ) 


###################################################
### code chunk number 8: soiltexture_vignette.Rnw:561-567
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "HYPRES.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the European system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 9: soiltexture_vignette.Rnw:640-641 (eval = FALSE)
###################################################
## install.packages( pkgs = "soiltexture" ) 


###################################################
### code chunk number 10: soiltexture_vignette.Rnw:650-654 (eval = FALSE)
###################################################
## install.packages( 
##     pkgs  = "soiltexture", 
##     repos = "http://R-Forge.R-project.org" 
## )   


###################################################
### code chunk number 11: soiltexture_vignette.Rnw:661-662
###################################################
library( soiltexture ) 


###################################################
### code chunk number 12: soiltexture_vignette.Rnw:670-672 (eval = FALSE)
###################################################
## detach( "package:soiltexture" ) 
## remove.packages( "soiltexture" ) 


###################################################
### code chunk number 13: soiltexture_vignette.Rnw:771-772
###################################################
TT.plot( class.sys = "none" ) 


###################################################
### code chunk number 14: soiltexture_vignette.Rnw:797-798
###################################################
TT.plot( class.sys = "USDA.TT" ) 


###################################################
### code chunk number 15: soiltexture_vignette.Rnw:812-818
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "USDA.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the USDA system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 16: soiltexture_vignette.Rnw:836-837
###################################################
TT.plot( class.sys = "USDA1911" ) 


###################################################
### code chunk number 17: soiltexture_vignette.Rnw:851-857
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "USDA1911" ) 
xtable( 
    x       = tex.tbl[, -3 ],  #
    caption = "Texture classes of the Whitney 1911 system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 18: soiltexture_vignette.Rnw:886-887
###################################################
TT.plot( class.sys = "HYPRES.TT" ) 


###################################################
### code chunk number 19: soiltexture_vignette.Rnw:907-913
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "HYPRES.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the European system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 20: soiltexture_vignette.Rnw:933-934
###################################################
TT.plot( class.sys = "FR.AISNE.TT" ) 


###################################################
### code chunk number 21: soiltexture_vignette.Rnw:948-954
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "FR.AISNE.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the French 'Aisne' system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 22: soiltexture_vignette.Rnw:973-974
###################################################
TT.plot( class.sys = "FR.GEPPA.TT" ) 


###################################################
### code chunk number 23: soiltexture_vignette.Rnw:985-991
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "FR.GEPPA.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the French 'GEPPA' system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 24: soiltexture_vignette.Rnw:1014-1015
###################################################
TT.plot( class.sys = "DE.BK94.TT" ) 


###################################################
### code chunk number 25: soiltexture_vignette.Rnw:1026-1032
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "DE.BK94.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the German system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 26: soiltexture_vignette.Rnw:1057-1058
###################################################
TT.plot( class.sys = "DE.SEA74.TT" ) 


###################################################
### code chunk number 27: soiltexture_vignette.Rnw:1063-1064
###################################################
plLim <- TT.get("DE.SEA74.TT")[["base.css.ps.lim"]][3]


###################################################
### code chunk number 28: soiltexture_vignette.Rnw:1075-1081
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "DE.SEA74.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the German SEA 1974 system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 29: soiltexture_vignette.Rnw:1096-1102
###################################################
TT.plot( 
    class.sys = "DE.SEA74.TT", 
    blr.clock   = rep(T,3), 
    tlr.an      = rep(60,3), 
    blr.tx      = c("SAND","CLAY","SILT"), 
)   #


###################################################
### code chunk number 30: soiltexture_vignette.Rnw:1118-1119
###################################################
TT.plot( class.sys = "DE.TGL85.TT" ) 


###################################################
### code chunk number 31: soiltexture_vignette.Rnw:1124-1125
###################################################
plLim <- TT.get("DE.TGL85.TT")[["base.css.ps.lim"]][3]


###################################################
### code chunk number 32: soiltexture_vignette.Rnw:1136-1142
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "DE.TGL85.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the German TGL 1985 system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 33: soiltexture_vignette.Rnw:1153-1159
###################################################
TT.plot( 
    class.sys = "DE.TGL85.TT", 
    blr.clock   = rep(T,3), 
    tlr.an      = rep(60,3), 
    blr.tx      = c("SAND","CLAY","SILT"), 
)   #


###################################################
### code chunk number 34: soiltexture_vignette.Rnw:1174-1175
###################################################
TT.plot( class.sys = "UK.SSEW.TT" ) 


###################################################
### code chunk number 35: soiltexture_vignette.Rnw:1185-1191
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "UK.SSEW.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the UK system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 36: soiltexture_vignette.Rnw:1208-1209
###################################################
TT.plot( class.sys = "AU2.TT" ) 


###################################################
### code chunk number 37: soiltexture_vignette.Rnw:1220-1226
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "AU2.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Australian system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 38: soiltexture_vignette.Rnw:1248-1249
###################################################
TT.plot( class.sys = "BE.TT" ) 


###################################################
### code chunk number 39: soiltexture_vignette.Rnw:1264-1270
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "BE.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Belgian system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 40: soiltexture_vignette.Rnw:1288-1289
###################################################
TT.plot( class.sys = "CA.EN.TT" ) 


###################################################
### code chunk number 41: soiltexture_vignette.Rnw:1296-1297
###################################################
TT.plot( class.sys = "CA.FR.TT" ) 


###################################################
### code chunk number 42: soiltexture_vignette.Rnw:1310-1316
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "CA.EN.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Canadian (en) system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 43: soiltexture_vignette.Rnw:1323-1329
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "CA.FR.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Canadian (fr) system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 44: soiltexture_vignette.Rnw:1355-1356
###################################################
TT.plot( class.sys = "ISSS.TT" ) 


###################################################
### code chunk number 45: soiltexture_vignette.Rnw:1369-1375
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "ISSS.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the ISSS system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 46: soiltexture_vignette.Rnw:1393-1394
###################################################
TT.plot( class.sys = "ROM.TT" ) 


###################################################
### code chunk number 47: soiltexture_vignette.Rnw:1407-1413
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "ROM.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Romanian system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 48: soiltexture_vignette.Rnw:1424-1430
###################################################
TT.plot( 
    class.sys = "ROM.TT", 
    blr.clock   = c(F,T,NA), 
    tlr.an      = c(45,90,45), 
    blr.tx      = c("SILT","CLAY","SAND"), 
)   #


###################################################
### code chunk number 49: soiltexture_vignette.Rnw:1451-1458
###################################################
test <- try( TT.plot( class.sys = "PL.TT" ) ) 

#   In case the posih triangle was not loaded at startup
if( "try-error" %in% class(test) ){ 
    plot(1,1,type="n",) 
    text(1,1,label="Plotting failed. Polish triangle not loaded")
}   


###################################################
### code chunk number 50: soiltexture_vignette.Rnw:1482-1484
###################################################
test <- try( plLim <- TT.get("PL.TT")[["base.css.ps.lim"]][3] ) 
if( "try-error" %in% class(test) ){ plLim <- NA_real_ }


###################################################
### code chunk number 51: soiltexture_vignette.Rnw:1507-1513
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "PL.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Polish system / triangle", 
    label   = NULL  
)   #


###################################################
### code chunk number 52: soiltexture_vignette.Rnw:1539-1540
###################################################
TT.plot( class.sys = "BRASIL.TT" ) 


###################################################
### code chunk number 53: soiltexture_vignette.Rnw:1553-1559
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "BRASIL.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Brazilian system (1996)", 
    label   = NULL  
)   #


###################################################
### code chunk number 54: soiltexture_vignette.Rnw:1577-1578
###################################################
TT.plot( class.sys = "SiBCS13.TT" ) 


###################################################
### code chunk number 55: soiltexture_vignette.Rnw:1591-1597
###################################################
tex.tbl <- TT.classes.tbl( class.sys = "SiBCS13.TT" ) 
xtable( 
    x       = tex.tbl[,-3],  #
    caption = "Texture classes of the Brazilian system (2013)", 
    label   = NULL  
)   #


###################################################
### code chunk number 56: soiltexture_vignette.Rnw:1620-1637
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles
TT.plot( 
    class.sys       = "USDA.TT", 
    class.p.bg.col  = TRUE
)   #

TT.plot( 
    class.sys       = "HYPRES.TT", 
    class.p.bg.col  = TRUE
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 57: soiltexture_vignette.Rnw:1645-1662
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles
TT.plot( 
    class.sys       = "FR.AISNE.TT", 
    class.p.bg.col  = TRUE
)   #

TT.plot( 
    class.sys       = "FR.GEPPA.TT", 
    class.p.bg.col  = TRUE
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 58: soiltexture_vignette.Rnw:1670-1687
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles
TT.plot( 
    class.sys       = "UK.SSEW.TT", 
    class.p.bg.col  = TRUE
)   #

TT.plot( 
    class.sys       = "DE.BK94.TT", 
    class.p.bg.col  = TRUE
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 59: soiltexture_vignette.Rnw:1694-1711
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles
TT.plot( 
    class.sys       = "AU2.TT", 
    class.p.bg.col  = TRUE
)   #

TT.plot( 
    class.sys       = "BE.TT", 
    class.p.bg.col  = TRUE
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 60: soiltexture_vignette.Rnw:1719-1736
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles
TT.plot( 
    class.sys       = "CA.EN.TT", 
    class.p.bg.col  = TRUE
)   #

TT.plot( 
    class.sys       = "CA.FR.TT", 
    class.p.bg.col  = TRUE
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 61: soiltexture_vignette.Rnw:1750-1754
###################################################
TT.plot( 
    class.sys       = "HYPRES.TT", 
    class.p.bg.col  = c("red","green","blue","pink","purple") 
)   #


###################################################
### code chunk number 62: soiltexture_vignette.Rnw:1775-1792
###################################################
# First plot the USDA texture triangle, and retrieve its 
#   geometrical features, silently outputted by TT.plot 
geo <- TT.plot( 
    class.sys   = "USDA.TT", 
    main        = "USDA and French Aisne triangles, overplotted"  
)   # 

# Then overplot the French Aisne texture triangle, 
#   and customise the colors so triangles are well distinct.
TT.classes(
    geo             = geo, 
    class.sys       = "FR.AISNE.TT", 
    # Additional "graphical" options
    class.line.col  = "red", 
    class.lab.col   = "red", 
    lwd.axis        = 2  
)   #


###################################################
### code chunk number 63: soiltexture_vignette.Rnw:1810-1827
###################################################
# First plot the USDA texture triangle, and retrieve its 
#   geometrical features, silently outputted by TT.plot 
geo <- TT.plot( 
    class.sys   = "FR.AISNE.TT", 
    main        = "French Aisne and GEPPA triangles, overplotted"  
)   # 

# Then overplot the French Aisne texture triangle, 
#   and customise the colors so triangles are well distinct.
TT.classes(
    geo             = geo, 
    class.sys       = "FR.GEPPA.TT", 
    # Additional "graphical" options
    class.line.col  = "red", 
    class.lab.col   = "red", 
    lwd.axis        = 2  
)   #


###################################################
### code chunk number 64: soiltexture_vignette.Rnw:1844-1854
###################################################
# Create a dummy data frame of soil textures:
my.text <- data.frame( 
    "CLAY"  = c(05,60,15,05,25,05,25,45,65,75,13,47), 
    "SILT"  = c(05,08,15,25,55,85,65,45,15,15,17,43), 
    "SAND"  = c(90,32,70,70,20,10,10,10,20,10,70,10), 
    "OC"    = c(20,14,15,05,12,15,07,21,25,30,05,28)  
)   #

# Display the table:
my.text


###################################################
### code chunk number 65: soiltexture_vignette.Rnw:1861-1866
###################################################
TT.plot( 
    class.sys   = "HYPRES.TT", 
    tri.data    = my.text, 
    main        = "Soil texture data" 
)   #


###################################################
### code chunk number 66: soiltexture_vignette.Rnw:1884-1890
###################################################
TT.plot( 
    class.sys   = "none", 
    tri.data    = my.text, 
    z.name      = "OC", 
    main        = "Soil texture triangle and OC bubble plot" 
)   #


###################################################
### code chunk number 67: soiltexture_vignette.Rnw:1912-1913
###################################################
rand.text	<- TT.dataset(n=100,seed.val=1980042401)


###################################################
### code chunk number 68: soiltexture_vignette.Rnw:1917-1923
###################################################
TT.plot( 
    class.sys   = "none", 
    tri.data    = rand.text, 
    z.name      = "Z", 
    main        = "Soil texture triangle and Z bubble plot" 
)   #


###################################################
### code chunk number 69: soiltexture_vignette.Rnw:1934-1983
###################################################
TT.plot( 
    class.sys   = "none", 
    tri.data    = my.text, 
    z.name      = "OC", 
    main        = "Soil texture triangle and OC bubble plot" 
)   #

# Recompute some internal values:
z.cex.range <- TT.get("z.cex.range") 
def.pch     <- par("pch") 
def.col     <- par("col")
def.cex     <- TT.get("cex") 
oc.str      <- TT.str( 
    my.text[,"OC"], 
    z.cex.range[1], 
    z.cex.range[2]
)   #

# The legend:
legend( 
    x           = 80, 
    y           = 90, 
    title       = 
        expression( bold('OC [g.kg'^-1 ~ ']') ), 
    legend      = formatC( 
        c( 
            min( my.text[,"OC"] ), 
            quantile(my.text[,"OC"] ,probs=c(25,50,75)/100), 
            max( my.text[,"OC"] ) 
        ), 
        format  = "f", 
        digits  = 1, 
        width   = 4, 
        flag    = "0" 
    ),  #
    pt.lwd      = 4, 
    col         = def.col, 
    pt.cex      = c( 
            min( oc.str ), 
            quantile(oc.str ,probs=c(25,50,75)/100), 
            max( oc.str ) 
    ),  #, 
    pch         = def.pch, 
    bty         = "o", 
    bg          = NA, 
    #box.col    = NA, # Uncomment this to remove the legend box
    text.col    = "black", 
    cex         = def.cex  
)   #


###################################################
### code chunk number 70: soiltexture_vignette.Rnw:2030-2049
###################################################
geo <- TT.geo.get() 
#
iwd.res <- TT.iwd( 
    geo         = geo, 
    tri.data    = rand.text, 
    z.name      = "Z", 
)   #
#
TT.image( 
    x       = iwd.res, 
    geo     = geo, 
    main    = "Soil texture triangle and Z heatmap" 
)   # 
#
TT.plot( 
    geo         = geo, 
    grid.show   = FALSE, 
    add         = TRUE  #  <<-- important 
)   #


###################################################
### code chunk number 71: soiltexture_vignette.Rnw:2079-2097
###################################################
TT.image( 
    x       = iwd.res, 
    geo     = geo, 
    main    = "Soil texture triangle and Z heatmap" 
)   # 
#
TT.contour( 
    x       = iwd.res, 
    geo     = geo, 
    add     = TRUE, #  <<-- important
    lwd     = 2  
)   # 
#
TT.plot( 
    geo         = geo, 
    grid.show   = FALSE, 
    add         = TRUE  #  <<-- important
)   #


###################################################
### code chunk number 72: soiltexture_vignette.Rnw:2132-2154
###################################################
geo <- TT.geo.get()  
#
kde.res <- TT.kde2d( 
    geo         = geo, 
    tri.data    = rand.text  
)   #
#
TT.contour( 
    x       = kde.res, 
    geo     = geo, 
    main    = "Probability density estimate of the texture data", 
    lwd     = 2, 
    col     = "red"  
)   # 
#
TT.plot( 
    tri.data    = rand.text, 
    geo         = geo, 
    grid.show   = FALSE, 
    add         = TRUE, #  <<-- important 
    col         = "gray"
)   #


###################################################
### code chunk number 73: soiltexture_vignette.Rnw:2200-2222
###################################################
geo <- TT.geo.get() 
#
maha <- TT.mahalanobis( 
    geo         = geo, 
    tri.data    = rand.text  
)   #
#
TT.contour( 
    x       = maha, 
    geo     = geo, 
    main    = "Texture data Mahalanobis distance", 
    lwd     = 2, 
    col     = "blue"  
)   # 
#
TT.plot( 
    tri.data    = rand.text, 
    geo         = geo, 
    grid.show   = FALSE, 
    add         = TRUE, #  <<-- important 
    col         = "gray"
)   #


###################################################
### code chunk number 74: soiltexture_vignette.Rnw:2250-2274
###################################################
geo <- TT.geo.get() 
#
maha <- TT.mahalanobis( 
    geo         = geo, 
    tri.data    = rand.text, 
    alr         = TRUE  #  <<-- important 
)   #
#
TT.contour( 
    x       = maha, 
    geo     = geo, 
    main    = "Texture data Mahalanobis distance", 
    lwd     = 2, 
    col     = "blue", 
    levels  = c(0.5,1,2,4,8)  #  <<-- manually set. Otherwise 
)   #                                 ugly plot
#
TT.plot( 
    tri.data    = rand.text, 
    geo         = geo, 
    grid.show   = FALSE, 
    add         = TRUE,  #  <<-- important 
    col         = "gray"
)   #


###################################################
### code chunk number 75: soiltexture_vignette.Rnw:2303-2318
###################################################
# Display the USDA texture triangle:
geo     <- TT.plot(class.sys="USDA.TT") 

# Create some custom labels:
labelz  <- letters[1:dim(my.text)[1]] 
labelz 

# Display the text
TT.text( 
    tri.data    = my.text, 
    geo         = geo, 
    labels      = labelz, 
    font        = 2, 
    col         = "blue"  
)   #


###################################################
### code chunk number 76: soiltexture_vignette.Rnw:2356-2357
###################################################
TT.data.test( tri.data = rand.text ) 


###################################################
### code chunk number 77: soiltexture_vignette.Rnw:2380-2391
###################################################
res <- TT.normalise.sum( tri.data = rand.text ) 
#
# With output of the residuals:
res <- TT.normalise.sum( 
    tri.data    = rand.text, 
    residuals   = TRUE  #  <<-- default = FALSE 
)   #
#
colnames( rand.text )
colnames( res )  #  "Z" has been dropped
max( res[ , "residuals" ] ) 


###################################################
### code chunk number 78: soiltexture_vignette.Rnw:2419-2423
###################################################
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT"  
)   #


###################################################
### code chunk number 79: soiltexture_vignette.Rnw:2432-2436
###################################################
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "USDA.TT"  
)   #


###################################################
### code chunk number 80: soiltexture_vignette.Rnw:2448-2453
###################################################
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT", 
    PiC.type    = "l" 
)   #


###################################################
### code chunk number 81: soiltexture_vignette.Rnw:2464-2469
###################################################
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT", 
    PiC.type    = "t" 
)   #


###################################################
### code chunk number 82: soiltexture_vignette.Rnw:2480-2486
###################################################
TT.points.in.classes( 
    tri.data    = my.text[1:5,], 
    class.sys   = "HYPRES.TT", 
    PiC.type    = "t", 
    collapse    = ";"
)   #


###################################################
### code chunk number 83: soiltexture_vignette.Rnw:2563-2654
###################################################
tmp.cex <- 1.5
old.par <- par(no.readonly = TRUE)
par(cex=tmp.cex,cex.axis=tmp.cex,cex.lab=tmp.cex,cex.main=tmp.cex)

tmp.text <- data.frame( "CLAY" = 20, "SILT" = 15, "SAND" = 65 ) 

plot( 
    x       = TT.dia2phi( c(2,20,2000) ), 
    y       = cumsum( unlist(tmp.text[1,]) ), 
    ylim    = c(0,100), 
    xlim    = TT.dia2phi( c(1,2000) ), 
    xaxt    = "n", 
    xlab    =  
expression( 'Particle size['~ mu * 'm] (log'[2] * 'scale)' ), 
    ylab    = "Cumulated particle size distribution [%]", 
    bty     = "n", 
    type    = "b", 
    main    = 
"Principle of particle size log-linear transformation", 
    cex     = tmp.cex  
)   #

lines( 
    spline( 
        y   = rev(cumsum( unlist(tmp.text[1,]) )), 
        x   = TT.dia2phi( c(2000,20,2))
    ),  # 
    col = "green"  
)   #

segments( 
    x0  = TT.dia2phi( c(2,20,2000) ), 
    x1  = TT.dia2phi( c(2,20,2000) ), 
    y0  = rep(0,3), 
    y1  = cumsum( unlist(tmp.text[1,]) ), 
    col = "red"  
)   #

new.tmp.text <- TT.text.transf( 
    tri.data        = tmp.text,  
    base.css.ps.lim = c(0,2,50,2000),  
    dat.css.ps.lim  = c(0,2,20,2000)   
)   #

new.silt.c <- cumsum( unlist(new.tmp.text[1,]) )[2]

arrows( 
    x0  = TT.dia2phi( c(50,50) ), 
    x1  = TT.dia2phi( c(50,1) ), 
    y0  = c(0,new.silt.c), 
    y1  = c(new.silt.c,new.silt.c), 
    col = "blue"  
)   #

text( 
    x       = TT.dia2phi( c(2,20,2000) ), 
    y       = cumsum( unlist(tmp.text[1,]) ), 
    pos     = 2, 
    offset  = 1, 
    labels  = c("Clay","Silt","Sand"), 
    col     = "red", 
    cex     = tmp.cex  
)   #

text( 
    x       = TT.dia2phi( c(50) ), 
    y       = new.silt.c, 
    pos     = 4, 
    offset  = 1, 
    labels  = "new Silt", 
    col     = "blue", 
    cex     = tmp.cex  
)   #

axis( 
    side    = 1, 
    at      = TT.dia2phi( c(2,20,50,2000) ), 
    labels  = c(2,20,50,2000) 
)   #

text( 
    x       = TT.dia2phi( 500 ), 
    y       = 65, 
    #pos    = 4, 
    #offset = 1, 
    labels  = "real distribution?", 
    col     = "green", 
    cex     = tmp.cex  
)   #

par(old.par) 


###################################################
### code chunk number 84: soiltexture_vignette.Rnw:2666-2667
###################################################
my.text[1:5,]   


###################################################
### code chunk number 85: soiltexture_vignette.Rnw:2678-2683
###################################################
TT.text.transf( 
	tri.data        = my.text[1:5,],  
	base.css.ps.lim = c(0,2,50,2000),  
	dat.css.ps.lim  = c(0,2,63,2000)   
)   #


###################################################
### code chunk number 86: soiltexture_vignette.Rnw:2691-2697
###################################################
# Copy the data.frame
my.text.fr  <- my.text 
# Curent columns names:
colnames(my.text.fr) 
# New columns names: 
colnames(my.text.fr) <- c("ARGILE","LIMON","SABLE","CO") 


###################################################
### code chunk number 87: soiltexture_vignette.Rnw:2704-2710
###################################################
TT.text.transf( 
    tri.data        = my.text.fr[1:5,],  
    base.css.ps.lim = c(0,2,50,2000),  
    dat.css.ps.lim  = c(0,2,63,2000),  
    css.names       = c("ARGILE","LIMON","SABLE")   
)   #


###################################################
### code chunk number 88: soiltexture_vignette.Rnw:2749-2760
###################################################
# Create a random fraction between 0 and 1
r.frac <- runif(n=dim(my.text)[1]) 
#
my.text4 <- cbind( 
    "CLAY"          = my.text[,"CLAY"], 
    "FINE_SILT"     = my.text[,"SILT"] * r.frac, 
    "COARSE_SILT"   = my.text[,"SILT"] * (1-r.frac), 
    "SAND"          = my.text[,"SAND"]  
)   #
#
my.text4[1:5,] 


###################################################
### code chunk number 89: soiltexture_vignette.Rnw:2769-2774
###################################################
TT.text.transf.X( 
    tri.data        = my.text4[1:5,], 
    base.ps.lim = c(0,2,20,50,2000),  
    dat.ps.lim  = c(0,2,20,63,2000)   
)   #


###################################################
### code chunk number 90: soiltexture_vignette.Rnw:2788-2793
###################################################
TT.text.transf.X( 
    tri.data        = my.text4[1:5,], 
    base.ps.lim = c(0,2,50,2000),  
    dat.ps.lim  = c(0,2,20,63,2000)   
)   #


###################################################
### code chunk number 91: soiltexture_vignette.Rnw:2807-2824
###################################################
# First, plot the data without transformation:
geo <- TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    tri.data    = my.text, 
    col         = "red", 
    main        = "Transformed and untransformed data"
)   #

# Then, re-plot them with transformation:
TT.points( 
    tri.data        = my.text, 
    geo             = geo, 
    dat.css.ps.lim  = c(0,2,63,2000),  
    css.transf      = TRUE, 
    col             = "blue", 
    pch             = 3  
)   #


###################################################
### code chunk number 92: soiltexture_vignette.Rnw:2847-2867
###################################################
# Not transformed
geo <- TT.plot( 
    class.sys   = "UK.SSEW.TT", 
    base.css.ps.lim = c(0,2,50,2000), 
    main        = 
        "Dummy transformation of the UK texture triangle"  
)   # 

# Transformed
TT.classes(
    geo             = geo, 
    class.sys       = "UK.SSEW.TT", 
    css.transf      = TRUE, 
    # Additional "graphical" options
    class.line.col  = "red", 
    class.lab.col   = "red", 
    lwd.axis        = 2, 
    class.lab.show  = "none", 
    class.lty       = 2 
)   #


###################################################
### code chunk number 93: soiltexture_vignette.Rnw:2883-2901
###################################################
# No transformation needed or stated
geo <- TT.plot( 
    class.sys   = "USDA.TT", 
    main        = 
        "USDA and transformed UK triangle, overplotted"  
)   # 

# Transformed
TT.classes(
    geo             = geo, 
    class.sys       = "UK.SSEW.TT", 
    css.transf      = TRUE,  #  <<-- important
    # Additional "graphical" options
    class.line.col  = "blue", 
    class.lab.col   = "blue", 
    lwd.axis        = 2, 
    class.lty       = 2 
)   #


###################################################
### code chunk number 94: soiltexture_vignette.Rnw:2912-2931
###################################################
# Untransformed
geo <- TT.plot( 
    class.sys   = "USDA.TT", 
    main        = 
        "(Dummy) transformation of the USDA texture triangle"  
)   # 

# Transformed
TT.classes(
    geo             = geo, 
    class.sys       = "USDA.TT", 
    tri.css.ps.lim  = c(0,2,20,2000), 
    css.transf      = TRUE,  #  <<-- important
    # Additional "graphical" options
    class.line.col  = "blue", 
    class.lab.col   = "blue", 
    lwd.axis        = 2, 
    class.lty       = 2 
)   #


###################################################
### code chunk number 95: soiltexture_vignette.Rnw:2942-2960
###################################################
geo <- TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    blr.tx      = c("SAND","CLAY","SILT"), 
    main        = 
        "(Dummy) transformation of the GEPPA texture triangle"  
)   # 

TT.classes(
    geo             = geo, 
    class.sys       = "FR.GEPPA.TT", 
    tri.css.ps.lim  = c(0,2,20,2000), 
    css.transf      = TRUE,  #  <<-- important
    # Additional "graphical" options
    class.line.col  = "blue", 
    class.lab.col   = "blue", 
    lwd.axis        = 2, 
    class.lty       = 2 
)   #


###################################################
### code chunk number 96: soiltexture_vignette.Rnw:2970-2990
###################################################
# Not transformed
geo <- TT.plot( 
    class.sys       = "FR.GEPPA.TT", 
    blr.tx          = c("SAND","CLAY","SILT"), 
    base.css.ps.lim  = c(0,2,20,2000), 
    main        = 
        "(Dummy) transformation of the GEPPA texture triangle"  
)   # 

# Transformed
TT.classes(
    geo             = geo, 
    class.sys       = "FR.GEPPA.TT", 
    css.transf      = TRUE,  #  <<-- important
    # Additional "graphical" options
    class.line.col  = "blue", 
    class.lab.col   = "blue", 
    lwd.axis        = 2, 
    class.lty       = 2 
)   #


###################################################
### code chunk number 97: soiltexture_vignette.Rnw:3006-3012
###################################################
TT.points.in.classes( 
    tri.data        = my.text[1:5,], 
    class.sys       = "USDA.TT", 
    dat.css.ps.lim  = c(0,2,20,2000), 
    css.transf      = TRUE   #  <<-- important
)   #


###################################################
### code chunk number 98: soiltexture_vignette.Rnw:3017-3024
###################################################
TT.plot( 
    class.sys       = "USDA.TT", 
    tri.data        = my.text, 
    dat.css.ps.lim  = c(0,2,20,2000), 
    css.transf      = TRUE,  #  <<-- important
    col             = "red"  
)   # 


###################################################
### code chunk number 99: soiltexture_vignette.Rnw:3032-3039
###################################################
TT.points.in.classes( 
    tri.data        = my.text[1:5,], 
    class.sys       = "USDA.TT", 
    dat.css.ps.lim  = c(0,2,20,2000), 
    base.css.ps.lim = c(0,2,20,2000), 
    css.transf      = TRUE  
)   #


###################################################
### code chunk number 100: soiltexture_vignette.Rnw:3045-3053
###################################################
TT.plot( 
    class.sys       = "USDA.TT", 
    tri.data        = my.text, 
    dat.css.ps.lim  = c(0,2,20,2000), 
    base.css.ps.lim = c(0,2,20,2000), 
    css.transf      = TRUE, 
    col             = "red"  
)   # 


###################################################
### code chunk number 101: soiltexture_vignette.Rnw:3097-3111
###################################################
# Create a new function, in fact the copy of TT.text.transf()
TT.text.transf2 <- TT.text.transf

# Imagine some changes in TT.text.transf2...

# Use your new function (will give identical results)
TT.points.in.classes( 
    tri.data        = my.text[1:5,], 
    class.sys       = "USDA.TT", 
    dat.css.ps.lim  = c(0,2,20,2000), 
    base.css.ps.lim = c(0,2,20,2000), 
    css.transf      = TRUE, 
    text.transf.fun = "TT.text.transf2"  #  <<-- important
)   #


###################################################
### code chunk number 102: soiltexture_vignette.Rnw:3118-3129
###################################################
TT.plot( 
    class.sys       = "USDA.TT", 
    tri.data        = my.text, 
    dat.css.ps.lim  = c(0,2,20,2000), 
    base.css.ps.lim = c(0,2,20,2000), 
    css.transf      = TRUE, 
    col             = "red", 
    text.transf.fun = "TT.text.transf2", #  <<-- important
    main            = 
        "Test of a (dummy) new transformation function"
)   # 


###################################################
### code chunk number 103: soiltexture_vignette.Rnw:3184-3189
###################################################
TT.plot( 
    class.sys   = "USDA.TT", 
    tlr.an      = c(45,90,45), 
    main        = "Re-projected USDA triangle (angles)"  
)   # 


###################################################
### code chunk number 104: soiltexture_vignette.Rnw:3214-3219
###################################################
TT.plot( 
    class.sys   = "FR.AISNE.TT", 
    blr.tx      = c("CLAY","SILT","SAND"), 
    main        = "Re-projected French Aisne triangle (axis)"  
)   # 


###################################################
### code chunk number 105: soiltexture_vignette.Rnw:3248-3253
###################################################
TT.plot( 
    class.sys   = "HYPRES.TT", 
    blr.clock   = c(FALSE,TRUE,NA), 
    main        = "Re-projected European triangle (axis directions)"  
)   # 


###################################################
### code chunk number 106: soiltexture_vignette.Rnw:3274-3281
###################################################
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    tlr.an      = c(60,60,60), 
    blr.tx      = c("SAND","CLAY","SILT"), 
    blr.clock   = c(TRUE,TRUE,TRUE), 
    main        = "Fully re-projected GEPPA triangle"  
)   # 


###################################################
### code chunk number 107: soiltexture_vignette.Rnw:3299-3315
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles with different geometries:
TT.plot( class.sys = "USDA.TT" ) 

TT.plot( 
    class.sys   = "USDA.TT", 
    blr.tx      = c("SILT","SAND","CLAY"), 
    blr.clock   = c(FALSE,FALSE,FALSE), 
    main        = "USDA triangle with a different geometry"  
)   # 

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 108: soiltexture_vignette.Rnw:3339-3356
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles with different languages:
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "fr" 
)   #

TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "de" 
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 109: soiltexture_vignette.Rnw:3364-3381
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles with different languages:
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "es" 
)   #

TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "it" 
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 110: soiltexture_vignette.Rnw:3391-3408
###################################################
# Set a 2 by 2 plot matrix:
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles with different languages:
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "nl" 
)   #

TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "fl" 
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 111: soiltexture_vignette.Rnw:3416-3434
###################################################
# Set a 2 by 2 plot matrix (for size):
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles with different languages:
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "se" 
)   #

# Plot the triangles with different languages:
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "ro" 
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 112: soiltexture_vignette.Rnw:3448-3460
###################################################
# Set a 2 by 2 plot matrix (for size):
old.par <- par(no.readonly=T)
par("mfcol" = c(1,2),"mfrow"=c(1,2)) 

# Plot the triangles with different languages:
TT.plot( 
    class.sys   = "FR.GEPPA.TT", 
    lang        = "en" 
)   #

# Back to old parameters:
par(old.par)


###################################################
### code chunk number 113: soiltexture_vignette.Rnw:3485-3490
###################################################
TT.plot( 
    tri.data    = my.text.fr, 
    class.sys   = "HYPRES.TT", 
    css.names   = c("ARGILE","LIMON","SABLE") 
)   #


###################################################
### code chunk number 114: soiltexture_vignette.Rnw:3514-3522
###################################################
TT.plot( 
    tri.data    = my.text.fr, 
    class.sys   = "HYPRES.TT", 
    css.names   = c("ARGILE","LIMON","SABLE"), 
    css.lab     = c("l'argile [%]","Le limon [%]","Le sable [%]"), 
    main        = 
        "A texture triangle with (dummy) custom axis names"  
)   #


###################################################
### code chunk number 115: soiltexture_vignette.Rnw:3530-3542
###################################################
TT.plot( 
    tri.data    = my.text.fr, 
    class.sys   = "HYPRES.TT", 
    css.names   = c("ARGILE","LIMON","SABLE"), 
    css.lab     = expression( 
        bold(sqrt('Argile'^2)~'[%]'), 
        bold(sqrt('Limon'^2)~'[%]'), 
        bold(sqrt('Sable'^2)~'[%]')
    ),  #
    main        = 
        "A texture triangle with (dummy) custom axis names"  
)   #


###################################################
### code chunk number 116: soiltexture_vignette.Rnw:3566-3577
###################################################
# Fisrt, retrieve all the data about 
#   the USDA texture triangle
tmp <- TT.get("USDA.TT") 

# It is not displayed here because it is to big
#   The list names are:
names(tmp) 

# If we drop "tt.points" and "tt.polygons", that will be 
#   presented later, the list size is more reasonable
tmp[ !names(tmp) %in% c("tt.points","tt.polygons") ]


###################################################
### code chunk number 117: soiltexture_vignette.Rnw:3595-3603
###################################################
# Retrieve and save the table:
tmp2 <- TT.classes.tbl( class.sys = "HYPRES.TT" ) 

# Display the first part:
tmp2[,1:2] 

# Then display the last column (and the 1st again):
tmp2[,c(1,3)] 


###################################################
### code chunk number 118: soiltexture_vignette.Rnw:3618-3619
###################################################
TT.vertices.tbl( class.sys = "HYPRES.TT" ) 


###################################################
### code chunk number 119: soiltexture_vignette.Rnw:3635-3647
###################################################
geo <- TT.plot( 
    class.sys   = "HYPRES.TT", 
    main        = "Vertices numbers. USDA texture triangle"
)   # 

TT.vertices.plot( 
    geo         = geo, 
    class.sys   = "HYPRES.TT", 
    col         = "red", 
    cex         = 2, 
    font        = 2  
)   #


###################################################
### code chunk number 120: soiltexture_vignette.Rnw:3676-3689
###################################################
# Step 1 
HYPRES63 <- TT.get("HYPRES.TT") 
#
# Visualize the data that will be modified
HYPRES63[[ "base.css.ps.lim" ]] 
HYPRES63[[ "tri.css.ps.lim" ]] 
#
# Step 2 
HYPRES63[[ "base.css.ps.lim" ]][3] <- 63 
HYPRES63[[ "tri.css.ps.lim" ]][3]  <- 63 
#
# Step 3: Load the new texture triangle
TT.add( "HYPRES63.TT" = HYPRES63 ) 


###################################################
### code chunk number 121: soiltexture_vignette.Rnw:3697-3701
###################################################
TT.plot( 
    class.sys   = "HYPRES63.TT", 
    main        = "Modified European soil texture triangle"
)   # 


###################################################
### code chunk number 122: soiltexture_vignette.Rnw:3712-3732
###################################################
# Get the definition of the HYPRES texture triangle
HYPRES <- TT.get( "HYPRES.TT" ) 
#
# Check its class (list) 
class( HYPRES ) 
#
# Check its parameters names 
names( HYPRES ) 
#
# Check its parameters class 
for( i in 1:length(HYPRES) )
{   
    print( 
        paste( 
            names( HYPRES )[i], 
            class( HYPRES[[i]] ), 
            sep = ": "
        )   
    )    
}   


