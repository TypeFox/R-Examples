opar <- par(ask = dev.interactive(orNone = TRUE))

##############################################
##
## Design I study: several habitat variables
##


########################
##
## Data management

data(bauges)
kasc <- bauges$kasc
loc <- bauges$locs


## The data are:
image(kasc)
image(getkasc(kasc,1))
points(loc)


########################
##
## Analysis of the environment


## Now, compute a data frame containing the value of the
## habitat variables in each mapped pixel of the map
litab <- kasc2df(kasc)

## Performs a PCA:
pc <- dudi.pca(litab$tab, scannf=FALSE, nf = 2)
opar <- par(mfrow=c(2,2))
barplot(pc$eig)                         # the eigenvalue diagram
s.corcircle(pc$co)                      # Correlations between environment
                                        # variable and principal components
## maps the scores of the pixels in the geographical space
opar <- par(mar=c(0,0,0,0))
image(getkasc(df2kasc(pc$li, litab$index, kasc), 1), axes=FALSE)
box()
par(opar)

## The first axis is an axis measuring the elevation



########################
##
## Analysis of the niche of the species


## counts the points in the pixels of the map
cp <- count.points(loc, kasc)
image(cp)

## keeps only the mapped pixels
zz <- cp[litab$index]

## plots the niche of the species on these variables
histniche(litab$tab, zz)


## Compute a habitat suitability map using the Mahalanobis distances
hsm <- mahasuhab(kasc, loc)
## Note that mahasuhab returns the *squared* Mahalanobis distance.
## plot the distance:
plot(sqrt(hsm), plot.axes=points(loc, pch=16, cex=0.7))

## Another way to clarify the map: use the rank of the pixels
plot(rank(hsm, na.last="keep"), plot.axes=points(loc, pch=16, cex=0.7))


## We may use a MADIFA to see the patterns
mad <- madifa(pc, zz, scan=FALSE)
opar <- par(mfrow=c(2,2))
barplot(mad$eig)                      # the eigenvalue diagram
scatterniche(mad$li,zz)                # plots the niche
s.arrow(cor(pc$tab,mad$li))            # meaning of the axis
par(opar)
## The first axis of the MADIFA opposes the areas at high elevation,
## close to grass and far from deciduous forests


## We may use an ENFA to see the patterns
en <- enfa(pc, zz, scan=FALSE)
opar <- par(mfrow=c(2,2))
barplot(en$s)                         # the specialisation diagram
scatterniche(en$li,zz)                # plots the niche
s.arrow(cor(pc$tab,en$li))            # meaning of the axis
par(opar)

## No axis of specialisation
## Actually, all the structures of the niche-environment
## system are caused by the marginality

niche.test(kasc,loc)

## The marginality is very strong, whereas the tolerance
## (1/specialisation) is only slightly significant
## Consequently, the first axis of the MADIFA is strongly
## correlated with the marginality axis of the ENFA
plot(mad$li[,1], en$li[,1])
cor(mad$li[,1], en$li[,1])

## We focus on the first axis of the MADIFA
## We compute a reduced rank habitat suitability map
## with the first axis:
ca <- df2kasc(mad$l1,litab$index,kasc)
hsmred <- getkasc(ca,1)^2
plot(sqrt(hsmred), plot.axes=points(loc, pch=16, cex=0.7))


## Compare the reduced rank HSM with the full rank hsm
opar <- par(mfrow=c(1,2))
image(sqrt(hsmred), main="Reduced")
image(sqrt(hsm), main="Full")
par(opar)
## Reduced-rank map is clearer



## A Resource selection function (Manly et al. 2002) based on
## these results: we build a model with the variables defining
## the first axis of the MADIFA:
mod <- glm(as.numeric(zz>0)~litab$tab$Grass+
           litab$tab$Elevation+litab$tab$Deciduous, family="binomial")
anova(mod, test="Chisq")

## The Deciduous variable is highly correlated with the Grass variable.
## These two variables are redundant. We keep only "Grass" and "Elevation"
mod2 <- glm(as.numeric(zz>0)~litab$tab$Grass+
            litab$tab$Elevation, family="binomial")
anova(mod2, test="Chisq")

## As recommended by Manly et al. 2002, we take the exponential
## of the link:
pr <- exp(predict(mod2)) ## relative probability

## Build the map
hsm3 <- getkasc(df2kasc(data.frame(pr,pr), litab$index, kasc),1)
image(hsm3)
points(loc)

## There is a close relationship, though not linear, with the
## reduced rank habitat suitability map build with the MADIFA:
plot(c(hsm3), c(hsmred))




cat("************************************************************\n",
    "The deeply commented source for this demo can be found in the file:\n",
    file.path(system.file(package = "adehabitat"), "demo", "nichehs.r"),
    "\n******************************************************\n")

