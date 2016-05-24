opar <- par(ask = dev.interactive(orNone = TRUE))

##################################################
##
##
## Introduction


## loads the data
data(puechabon)
locs <- puechabon$locs
head(locs)
kasc <- puechabon$kasc
kasc
## kasc is a multilayer map of class "kasc"

image(kasc)

## Single maps of class asc can be extracted from this map:
(elev <- getkasc(kasc,"Elevation")) ## elev is numeric
(aspe <- getkasc(kasc,"Aspect"))    ## aspe is a factor

## Note that a kasc object is created using maps of class 'asc'
ka2 <- as.kasc(list(tutu=elev, toto=aspe))
image(ka2)
## maps of class asc can be imported from a GIS using import.asc.

## Now see the map of elevation
image(elev)
plot(elev)
opar<-par(mar = c(0,0,3,0), bg = "slategray")
persp(elev, scale = FALSE,
      box = FALSE, border = NA, shade = 0.75,
      col = "darkolivegreen3", expand = 2, theta = -60, phi = 30,
      main = "The topography of Puechabon")
par(opar)

## Come back to the map "kasc"
if (require(tkrplot)) {
    explore.kasc(kasc)
}

##################################################
##
##
## A multivariate analysis of the multi-layer map


## First transform the kasc in data frame
## (Mapped pixel in rows X mapped variables in columns)
litab <- kasc2df(kasc)
names(litab)

## a Hill-Smith Analysis of the table
## (a kind of mix between a principal component analysis and
## multiple correspondence analysis)
dud <- dudi.hillsmith(litab$tab, scan=F)


## map the scores of the principal components
kadud <- df2kasc(dud$li,litab$index, kasc)
image(kadud)

## Now keep the areas with elevation lower than 200 meters
mask <- elev
mask[mask>200] <- NA
glou <- setmask(kasc,mask)
image(glou)


##################################################
##
##
## Now, add the relocations:

image(elev)
points(locs[,c("X","Y")], pch=16, col="red")


## what is the value of the environmental
## variables for each relocation?
jo <- join.kasc(locs[,c("X","Y")], kasc)
head(jo)

## How many relocations in each pixel?
cp <- count.points(locs[,c("X","Y")], kasc)
image(cp)

## How many relocations in each pixel, for each animal?
cp2 <- count.points.id(locs[,c("X","Y")], locs$Name, kasc)
image(cp2)

## Define a buffer of 300 meters around the relocations
bu <- buffer(locs[,c("X","Y")], kasc, 300)
image(bu)

## Use this buffer as a mask for the multilayer maps
bumap <- setmask(kasc, bu)
image(bumap)

## Now, work with bu:
image(bu)

## How to count automatically the connex
## features on this map?
labu <- labcon(bu)
image(labu, main="One color per connex feature")
table(labu)

## Three components
## OK, now keep only the largest
labu[labu!=2] <- NA
labu[labu==2] <- 1

## The geographical attributes of the map have been lost:
class(labu)

## Copy them from another map, e.g. elev
labu <- getascattr(elev, labu)
class(labu)
image(labu)


## OK now, note that the animals are located in the
## south of the area:
image(elev)
points(locs[,c("X","Y")])

## We may also focus on this part:
xl <- c(697034, 702935)
yl <- c(3156388, 3161974)
su <- subsetmap(kasc, xlim=xl, ylim=yl)
image(su)


##################################################
##
##
## Now, work with a factor map

aspe

## aspe is a factor. See the levels
levels(aspe)

## Ok, plot it
image(aspe)

## Try someting different, reflecting
## the sunshine on the area:
image(aspe, clfac=c("purple1","green1","green4", "purple4"))
legend(696586,3166042, levels(aspe),
       c("purple1","green1","green4", "purple4"))

## Mmmmmmh... well... OK
## Now let us recode the map with dummy variables
co <- convnum(as.kasc(list(aspect=aspe)))
image(co$kasc)

## What is the distance from each category?
di <- distfacmap(aspe)
image(di)


cat("******************************************************\n",
    "The deeply commented source for this demo can be found in the file:\n",
    file.path(system.file(package = "adehabitat"), "demo", "rastermaps.R"),
    "\n******************************************************\n")
