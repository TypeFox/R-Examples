

#######  load in the data:

data(cosogeol)
data(cosomap)
     data(faults)
     data(hiways)
     data(owens)

###  set the projection
proj = cosomap$PROJ
     plotGEOmapXY(cosomap, PROJ=proj,  add=FALSE, ann=FALSE, axes=FALSE)


##########  add the geology as it is stored

     plotGEOmapXY(cosogeol, PROJ=proj,  add=TRUE, ann=FALSE, axes=FALSE)


#########  this does not look too good.
########   the color palette is not very attractive or organized.
####  change the colors

####  set up colors from geotouch
XMCOL=setXMCOL()

####  assign colors
newcol = XMCOL[cosogeol$STROKES$col+1]
cosocolnums = cosogeol$STROKES$col
cosogeol$STROKES$col = newcol

#####   get unique names of units
ss = strsplit(cosogeol$STROKES$nam, split="_")     

geo = unlist(lapply(ss  , FUN="getmem", mem=1))
UGEO = unique(geo)


mgeo = match( geo, UGEO )

gcol = paste(sep=".", geo, cosogeol$STROKES$col)


ucol = unique(gcol)



ucol = unique(gcol)

N = length(ucol)


spucol = strsplit(ucol,split="\\.")     

       
names = unlist(lapply(spucol  , FUN="getmem", mem=1))

shades = unlist(lapply(spucol  , FUN="getmem", mem=2))

ORDN = order(names)

geoLEGEND(names[ORDN], shades[ORDN], .28, .14, 16, 6)


