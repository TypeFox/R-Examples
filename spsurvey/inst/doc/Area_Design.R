### R code from vignette source 'Area_Design.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: createshape
###################################################
# Load the sp object in the data directory
data(UT_ecoregions)

# Create a shapefile
sp2shape(sp.obj=UT_ecoregions, shpfilename="UT_ecoregions")



###################################################
### code chunk number 3: att
###################################################
# Read the attribute table from the shapefile
att <- read.dbf("UT_ecoregions")



###################################################
### code chunk number 4: att
###################################################
# Display the attribute data frame
att



###################################################
### code chunk number 5: att
###################################################
# Summarize frame area by ecoregion
temp <- tapply(att$Area_ha, att$Level3_Nam, sum)
temp <- round(addmargins(temp), 0)
temp



###################################################
### code chunk number 6: Equalsites
###################################################
# Call the set.seed function so that the survey designs can be replicate
set.seed(4447864)



###################################################
### code chunk number 7: Equalsites
###################################################
# Create the design list
Equaldsgn <- list(None=list(panel=c(PanelOne=115), seltype="Equal"))



###################################################
### code chunk number 8: Equalsites
###################################################
# Select the sample
Equalsites <- grts(design=Equaldsgn,
                   DesignID="EQUAL",
                   type.frame="area",
                   src.frame="shapefile",
                   in.shape="UT_ecoregions", 
                   att.frame=att,
                   shapefile=FALSE)



###################################################
### code chunk number 9: Equalsites
###################################################
# Print the initial six lines of the survey design
head(Equalsites@data)



###################################################
### code chunk number 10: Equalsites
###################################################
# Print the survey design summary
summary(Equalsites)



###################################################
### code chunk number 11: Unequalsites
###################################################
# Create the design list
Unequaldsgn <- list(None=list(panel=c(PanelOne=115),
                              seltype="Unequal",
                              caty.n=c("Central Basin and Range"=25,
                                       "Colorado Plateaus"=25,
                                       "Mojave Basin and Range"=10,
                                       "Northern Basin and Range"=10,
                                       "Southern Rockies"=10,
                                       "Wasatch and Uinta Mountains"=25,
                                       "Wyoming Basin"=10)))



###################################################
### code chunk number 12: Unequalsites
###################################################
# Select the sample
Unequalsites <- grts(design=Unequaldsgn,
                     DesignID="UNEQUAL",
                     type.frame="area",
                     src.frame="shapefile",
                     in.shape="UT_ecoregions", 
                     att.frame=att,
                     mdcaty="Level3_Nam",									
                     shapefile=FALSE)



###################################################
### code chunk number 13: Unequalsites
###################################################
# Print the initial six lines of the survey design
head(Unequalsites@data)



###################################################
### code chunk number 14: Unequalsites
###################################################
# Print the survey design summary
summary(Unequalsites)



###################################################
### code chunk number 15: Stratsites
###################################################
# Read the shapefile
shp <- read.shape("UT_ecoregions")



###################################################
### code chunk number 16: Stratsites
###################################################
# Create the design list
Stratdsgn <- list("Central Basin and Range"=list(panel=c(PanelOne=25),
                                                 seltype="Equal"),
                  "Colorado Plateaus"=list(panel=c(PanelOne=25),
                                           seltype="Equal"),
                  "Mojave Basin and Range"=list(panel=c(PanelOne=10),
                                                seltype="Equal"),
                  "Northern Basin and Range"=list(panel=c(PanelOne=10),
                                                  seltype="Equal"),
                  "Southern Rockies"=list(panel=c(PanelOne=10),
                                          seltype="Equal"),
                  "Wasatch and Uinta Mountains"=list(panel=c(PanelOne=25),
                                                     seltype="Equal"),
                  "Wyoming Basin"=list(panel=c(PanelOne=10),
                                       seltype="Equal"))



###################################################
### code chunk number 17: Stratsites
###################################################
# Select the sample
Stratsites <- grts(design=Stratdsgn,
                   DesignID="STRATIFIED",
                   type.frame="area",
                   src.frame="sp.object",
                   sp.object=shp,
                   att.frame=att,
                   stratum="Level3_Nam",									
                   shapefile=FALSE)



###################################################
### code chunk number 18: Stratsites
###################################################
# Print the initial six lines of the survey design
head(Stratsites@data)



###################################################
### code chunk number 19: Stratsites
###################################################
# Print the survey design summary
summary(Stratsites)



###################################################
### code chunk number 20: Panelsites
###################################################
# Create the design list
Paneldsgn <- list(None=list(panel=c(Year1=50, Year2=50, Year3=50,
                                    Year4=50, Year5=50),
                            seltype="Unequal",
                            caty.n=c("Central Basin and Range"=64,
                                     "Colorado Plateaus"=63,
                                     "Mojave Basin and Range"=15,
                                     "Northern Basin and Range"=15,
                                     "Southern Rockies"=15,
                                     "Wasatch and Uinta Mountains"=63,
                                     "Wyoming Basin"=15),
                            over=100))



###################################################
### code chunk number 21: Panelsites
###################################################
# Select the sample
Panelsites <- grts(design=Paneldsgn,
                   DesignID="UNEQUAL",
                   type.frame="area",
                   src.frame="shapefile",
                   in.shape="UT_ecoregions", 
                   att.frame=att,
                   mdcaty="Level3_Nam",									
                   shapefile=FALSE)



###################################################
### code chunk number 22: Panelsites (eval = FALSE)
###################################################
## # Print the warning message
## warnings()
## 


###################################################
### code chunk number 23: Panelsites
###################################################
# Print the initial six lines of the survey design
head(Panelsites@data)



###################################################
### code chunk number 24: Panelsites
###################################################
# Print the survey design summary
summary(Panelsites)



