### R code from vignette source 'ImportExport.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ImportExport.Rnw:32-45
###################################################
options(warn=2)

library(cheddar)

# Makes copy-paste much less painful
options(continue=' ')
options(width=90)
options(prompt='> ')

options(SweaveHooks = list(fig=function() par(mgp=c(2.5,1,0), 
                                              mar=c(4,4,2,1),
                                              oma=c(0,0,1,0),
                                              cex.main=0.8)))


###################################################
### code chunk number 2: ImportExport.Rnw:68-69
###################################################
stream12 <- LoadCommunity('Stream12')


###################################################
### code chunk number 3: ImportExport.Rnw:123-124
###################################################
cat(paste(paste(colnames(NPS(stream12)), collapse=' & '), '\\\\'))


###################################################
### code chunk number 4: ImportExport.Rnw:127-129
###################################################
junk <- apply(NPS(stream12), 1, 
              function(row) cat(paste(paste(replace(row, which(is.na(row)), ''), collapse=' & '), ' \\\\ \n')))


###################################################
### code chunk number 5: ImportExport.Rnw:165-166
###################################################
cat(paste(paste(colnames(TLPS(stream12)), collapse=' & '), '\\\\'))


###################################################
### code chunk number 6: ImportExport.Rnw:169-171
###################################################
junk <- apply(TLPS(stream12), 1, 
              function(row) cat(paste(paste(row, collapse=' & '), ' \\\\ \n')))


###################################################
### code chunk number 7: ImportExport.Rnw:189-190 (eval = FALSE)
###################################################
## stream12 <- LoadCommunity('c:/Stream12')


###################################################
### code chunk number 8: ImportExport.Rnw:193-196
###################################################
stream12
NumberOfNodes(stream12)
NumberOfTrophicLinks(stream12)


###################################################
### code chunk number 9: ImportExport.Rnw:199-207
###################################################
# Community properties
CPS(stream12)

# Node properties
NPS(stream12)

# Trophic links
TLPS(stream12)


###################################################
### code chunk number 10: ImportExport.Rnw:211-212
###################################################
SumBiomassByClass(stream12)


###################################################
### code chunk number 11: ImportExport.Rnw:215-216
###################################################
SumBiomassByClass(stream12, class='functional.group')


###################################################
### code chunk number 12: ImportExport.Rnw:221-222
###################################################
grassland <- LoadCollection('Grassland1994')


###################################################
### code chunk number 13: ImportExport.Rnw:253-254 (eval = FALSE)
###################################################
## grassland <- LoadCommunity('c:/Grassland1994')


###################################################
### code chunk number 14: ImportExport.Rnw:256-258
###################################################
grassland
length(grassland)


###################################################
### code chunk number 15: ImportExport.Rnw:274-302 (eval = FALSE)
###################################################
## install.packages('igraph')
## library(igraph)
## ToIgraph <- function(community, weight=NULL)
## {
##     if(is.null(TLPS(community)))
##     {
##         stop('The community has no trophic links')
##     }
##     else
##     {
##         tlps <- TLPS(community, link.properties=weight)
##         if(!is.null(weight))
##         {
##             tlps$weight <- tlps[,weight]
##         }
##         return (graph.data.frame(tlps, 
##                                  vertices=NPS(community), 
##                                  directed=TRUE))
##     }
## }
## 
## data(TL84)
## # Unweighted network
## TL84.ig <- ToIgraph(TL84)
## 
## data(Benguela)
## # Use diet fraction to weight network
## Benguela.ig <- ToIgraph(Benguela, weight='diet.fraction')


###################################################
### code chunk number 16: ImportExport.Rnw:310-329 (eval = FALSE)
###################################################
## install.packages("NetIndices")
## install.packages("foodweb", repos="http://R-Forge.R-project.org")
## 
## library(foodweb)    # Loads the dependency NetIndices
## 
## data(TL84)
## TL84.ni <- PredationMatrix(TL84, weight='diet.fraction')
## 
## # Plot the predation matrix
## foodweb::imageweb(TL84.ni)
## 
## data(Benguela)
## # Use diet fraction to weight network
## Benguela.ni <- PredationMatrix(Benguela, weight='diet.fraction')
## 
## # Compute flow-based trophic level using both Cheddar and NetIndices. Both 
## # packages should give the same results
## cbind(ni.tl=NetIndices::TrophInd(Benguela.ni)[,'TL'], 
##       cheddar.tl=FlowBasedTrophicLevel(Benguela, weight.by='diet.fraction'))


###################################################
### code chunk number 17: ImportExport.Rnw:337-355 (eval = FALSE)
###################################################
## library(cheddar)
## 
## ExportToNetwork3D <- function(community, dir, fname_root=CP(community, 'title'))
## {
##     links <- cbind(PredatorID=NodeNameIndices(community, TLPS(community)[,'consumer']),
##                    PreyID=NodeNameIndices(community, TLPS(community)[,'resource']))
##     write.table(links, file.path(dir, paste(fname_root, '_links.web', sep='')), 
##                 row.names=FALSE, sep=' ')
## 
##     species <- cbind(ID=1:NumberOfNodes(community), 
##                      CommonName=NP(community, 'node'))
##                      
##     write.table(links, file.path(dir, paste(fname_root, '_species.txt', sep='')), 
##                 row.names=FALSE, sep=' ')
## }
## 
## data(TL84)
## ExportToNetwork3D(TL84, 'c:')


###################################################
### code chunk number 18: ImportExport.Rnw:361-376 (eval = FALSE)
###################################################
## install.packages("foodweb")
## 
## library(foodweb)
## 
## data(TL84)
## 
## # Write the predation matrix to a csv file in the format required by foodweb.
## write.table(PredationMatrix(TL84), 'TL84.foodweb.csv', row.names=FALSE, 
##             col.names=FALSE, sep=',')
## 
## # foodweb::analyse.single() creates a file called 'Results-TL84.foodweb.csv' 
## # to the current directory. This file will contain some basic food-web 
## # statistics and is required by the foodweb::plotweb() function.
## foodweb::analyse.single('TL84.foodweb.csv')
## foodweb::plotweb(cols=1:7, radii=7:1)           # A 3D plot


