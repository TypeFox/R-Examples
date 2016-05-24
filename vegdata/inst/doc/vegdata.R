## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance = TRUE,
results = "tex")
options(results = 'tex') # for pqR

## ----prep, echo=FALSE, results='hide'-----------------------------------------------------------------------
library(knitr)
options(width=110, digits=2)
opts_chunk$set(comment = "", warning = FALSE, message = TRUE, echo = TRUE, size="footnotesize")
# read_chunk("some/script/I/want/to/load.R")
tmp <- tempdir()
suppressPackageStartupMessages(library(vegdata))
options(tv_home = tmp)
dir.create(file.path(tmp, 'Species'))
dir.create(file.path(tmp, 'Popup'))
dir.create(file.path(tmp, 'Data'))
file.copy(from = file.path(path.package("vegdata"), 'tvdata', 'Popup'), to = tmp, recursive = TRUE)
file.copy(from = file.path(path.package("vegdata"), 'tvdata', 'Species'), to = tmp, recursive = TRUE)
file.copy(from = file.path(path.package("vegdata"), 'tvdata', 'Data'), to = tmp, recursive = TRUE)

## ----load, results='hide'-----------------------------------------------------------------------------------
library(vegdata)

## ----eval=TRUE----------------------------------------------------------------------------------------------
h <- tv.home()

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  options(tv_home="path_to_your_Turboveg_root_directory")

## ----dblisting, eval=FALSE----------------------------------------------------------------------------------
#  tv.db()

## -----------------------------------------------------------------------------------------------------------
tv.refl()

## ----tax----------------------------------------------------------------------------------------------------
tax('Brachythecium rutabulum')

## ----syn----------------------------------------------------------------------------------------------------
tax('Elytrigia repens')$TaxonName
syn('Elytrigia repens')

## ----childs-------------------------------------------------------------------------------------------------
childs(27, quiet=TRUE)$TaxonName
parents('ACHIMIL')

## ----db-----------------------------------------------------------------------------------------------------
db <- 'taxatest'

## ----path, results='hide', echo=FALSE-----------------------------------------------------------------------

## ----meta, eval=FALSE---------------------------------------------------------------------------------------
#  tv.metadata(db)

## ----obs----------------------------------------------------------------------------------------------------
getOption('tv_home')
obs.tax <- tv.obs(db)
# Adding species names
species <- tax('all')
obs.tax$TaxonName <-  species$TaxonName[match(obs.tax$TaxonUsageID, species$TaxonUsageID)]
head(obs.tax[,c('RELEVE_NR','TaxonUsageID','COVER_CODE','LAYER','TaxonName')])

## ----taxval, eval=TRUE--------------------------------------------------------------------------------------
obs.tax$OriginalName <- obs.tax$TaxonName
obs.taxval <- taxval(obs.tax, db=db, mono='lower', maxtaxlevel='AGG', interactive=FALSE)

## ----Taxon--------------------------------------------------------------------------------------------------
obs.taxval$OriginalName <- obs.taxval$TaxonName
obs.taxval$TaxonName <-  species$TaxonName[match(obs.taxval$TaxonUsageID, species$TaxonUsageID)]
obs.taxval[!duplicated(obs.taxval$OriginalName),c('RELEVE_NR', 'COVER_CODE', 'TaxonName', 'OriginalName')]

## ----coarsen, eval=TRUE, results='hide'---------------------------------------------------------------------
tmp <- taxval(obs.tax, refl='GermanSL 1.3', ag='adapt', rank='FAM')
tmp$newTaxon <- tax(tmp$TaxonUsageID, refl='GermanSL 1.3')$TaxonName

## ----print.coarsen------------------------------------------------------------------------------------------
head(tmp[,c('OriginalName','newTaxon')], 10)

## ----coverperc, echo=2:4, eval=TRUE-------------------------------------------------------------------------
options(width=120)
obs <- tv.obs(db)
# obs <- tv.coverperc(db, obs)
tail(obs)
options(width=110)

## ----pseudo1, eval=FALSE------------------------------------------------------------------------------------
#  data(lc.0)
#  obs <- tv.obs(db)
#  tv.veg(db, pseudo = list(lc.0, c("LAYER")), lc = "layer")

## ----lc0, warning=FALSE-------------------------------------------------------------------------------------
tmp <- tv.veg(db, tax=FALSE, pseudo = list(lc.0, "LAYER"), lc = "layer", quiet=TRUE)
names(tmp)

## ----Season-------------------------------------------------------------------------------------------------
comb <- list(data.frame(SEASON=0:4, COMB=c(0,'Spring','Summer','Autumn','Winter')),'SEASON')
names(tv.veg(db, tax=FALSE, pseudo=comb, quiet=TRUE))

## ----layer, results='hide', warning=FALSE-------------------------------------------------------------------
data(lc.1)
veg <- tv.veg(db, lc = "sum", pseudo = list(lc.1, 'LAYER'), dec = 1, check.critical = FALSE)

## ----layerdiff----------------------------------------------------------------------------------------------
veg[,1:10]

## -----------------------------------------------------------------------------------------------------------
obs.tax$TaxonUsageID[obs.tax$TaxonUsageID == 27] <- 31

## ----replace------------------------------------------------------------------------------------------------
taxon.repl <- data.frame(old=c(27), new=c(31))
obs.tax$TaxonUsageID <- replace(obs.tax$TaxonUsageID, 
    match(taxon.repl$old, obs.tax$TaxonUsageID), taxon.repl$new)

## ----comb.spec, eval=TRUE-----------------------------------------------------------------------------------
comb.species(veg, sel=c('QUERR-R','QUERR-R.Tree'))

## ----site.echo, eval=TRUE-----------------------------------------------------------------------------------
site <- tv.site('taxatest')

## ----eval=FALSE---------------------------------------------------------------------------------------------
#  tv.compRefl('taxref1', 'taxref2')

## ----elbaue, results='hide'---------------------------------------------------------------------------------
elbaue <- tv.veg('elbaue')
elbaue.env <- tv.site('elbaue')

## ----cluster------------------------------------------------------------------------------------------------
clust <- vector('integer', nrow(elbaue.env))
clust[elbaue.env$MGL < -50 & elbaue.env$SDGL < 50] <- 1		# dry sites, low deviation
clust[elbaue.env$MGL < -50 & elbaue.env$SDGL >= 50] <- 2	# dry sites, high deviation
clust[elbaue.env$MGL >= -50 & elbaue.env$SDGL >= 50] <- 3	# wet sites, high deviation
clust[elbaue.env$MGL >= -50 & elbaue.env$SDGL < 50] <- 4	# wet sites, low deviation
levels(clust) <- c('dry.ld','dry.hd', 'wet.hd','wet.ld')

## ----syntab.mupa--------------------------------------------------------------------------------------------
require(indicspecies)
st <- syntab(elbaue, clust, mupa=TRUE)
# Print Ellenberg indicator values for soil moisture and nutrient demand
traits <- tv.traits()
trait <- traits[traits$LETTERCODE %in% names(elbaue), ]
rownames(trait) <- trait$LETTERCODE
trait <- trait[,c('OEK_F', 'OEK_N')]
print(st, limit=30, trait=trait)

## ----nmds, quiet=TRUE, results='hide', eval=TRUE------------------------------------------------------------
## Data analyses
library(vegan)
veg.nmds <- metaMDS(elbaue, distance = "bray", trymax = 5, autotransform =FALSE, 
     noshare = 1, expand = TRUE, trace = 2)
eco <- tv.traits()
eco$OEK_F <- as.numeric(eco$OEK_F)
F <- isc(elbaue, eco, 'OEK_F', method = 'mean')
N <- isc(elbaue, eco, 'OEK_N', method = 'mean')
env <- envfit(veg.nmds, data.frame(F, N))

## ----nmdsplotfun, quiet=TRUE, results='hide'----------------------------------------------------------------
library(labdsv)
library(akima)
color = function(x)rev(topo.colors(x))
nmds.plot <- function(ordi, site, var1, var2, disp, plottitle =  'NMDS', env = NULL, ...) {
 lplot <- nrow(ordi$points);  lspc <- nrow(ordi$species)
 filled.contour(interp(ordi$points[, 1], ordi$points[, 2], site[, var1]), 
                ylim = c(-1, 1.1), xlim = c(-1.4, 1.4),
   color.palette = color, xlab = var1, ylab = var2, main = plottitle,
    key.title = title(main = var1, cex.main = 0.8, line = 1, xpd = NA),
    plot.axes = { axis(1);  axis(2)
      points(ordi$points[, 1], ordi$points[, 2], xlab = "", ylab = "", cex= .5, col = 2, pch = '+')
      points(ordi$species[, 1], ordi$species[, 2], xlab = "", ylab = "", cex=.2, pch = 19)
      ordisurf(ordi, site[, var2], col = 'black', choices = c(1, 2), add = TRUE)
      orditorp(ordi, display = disp, pch = " ")
      legend("topright", paste("GAM of ", var2), col = 'black', lty = 1)
      if(!is.null(env)) plot(env, col='red')
   }
  ,...)
}

## ----nmdsplot, quiet=TRUE, results='hide', eval=TRUE, warning=FALSE-----------------------------------------
nmds.plot(veg.nmds, elbaue.env, disp='species', var1="MGL", var2="SDGL", env=env, 
        plottitle = 'Elbaue floodplain dataset')

