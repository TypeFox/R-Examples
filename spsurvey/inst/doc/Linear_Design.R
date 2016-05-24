### R code from vignette source 'Linear_Design.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: createshape
###################################################
# Load the sp object in the data directory
data(Luck_Ash_streams)

# Create a shapefile
sp2shape(sp.obj=Luck_Ash_streams, shpfilename="Luck_Ash_streams")



###################################################
### code chunk number 3: att
###################################################
# Read the attribute table from the shapefile
att <- read.dbf("Luck_Ash_streams")



###################################################
### code chunk number 4: att
###################################################
# Display the initial six lines in the attribute data frame
head(att)



###################################################
### code chunk number 5: att
###################################################
# Display number of stream segments cross-classified by the  strata and
# multidensity category variables

addmargins(table("Stream Type"=att$Per_Int, "Strahler Order"=att$Strah_Cat))



###################################################
### code chunk number 6: att
###################################################
# Summarize frame stream length by stratum and multidensity category
temp <- tapply(att$Length_km, list(att$Per_Int, att$Strah_Cat), sum)
temp <- round(addmargins(temp), 2)
names(dimnames(temp)) <- list("Stream Type", "Strahler Order")
temp



###################################################
### code chunk number 7: Equalsites
###################################################
# Call the set.seed function so that the survey designs can be replicate
set.seed(19742003)



###################################################
### code chunk number 8: Equalsites
###################################################
# Create the design list
Equaldsgn <- list(None=list(panel=c(PanelOne=50), seltype="Equal"))



###################################################
### code chunk number 9: Equalsites
###################################################
# Select the sample
Equalsites <- grts(design=Equaldsgn,
                   DesignID="EQUAL",
                   type.frame="linear",
                   src.frame="shapefile",
                   in.shape="Luck_Ash_streams",
                   att.frame=att,
                   shapefile=FALSE)



###################################################
### code chunk number 10: Equalsites
###################################################
# Print the initial six lines of the survey design
head(Equalsites@data)



###################################################
### code chunk number 11: Equalsites
###################################################
# Print the survey design summary
summary(Equalsites)



###################################################
### code chunk number 12: Stratsites
###################################################
# Create the design list
Stratdsgn <- list(Perennial=list(panel=c(PanelOne=50),
                                 seltype="Equal",
                                 over=50),
                  Intermittent=list(panel=c(PanelOne=50),
                                    seltype="Equal",
                                    over=50))



###################################################
### code chunk number 13: Stratsites
###################################################
# Select the sample
Stratsites <- grts(design=Stratdsgn,
                   DesignID="STRATIFIED",
                   type.frame="linear",
                   src.frame="shapefile",
                   in.shape="Luck_Ash_streams",
                   att.frame=att,
                   stratum="Per_Int",
                   shapefile=FALSE)



###################################################
### code chunk number 14: Stratsites
###################################################
# Print the initial six lines of the survey design
head(Stratsites@data)



###################################################
### code chunk number 15: Stratsites
###################################################
# Print the survey design summary
summary(Stratsites)



###################################################
### code chunk number 16: Unequalsites
###################################################
# Read the shapefile
shp <- read.shape("Luck_Ash_streams")



###################################################
### code chunk number 17: Unequalsites
###################################################
# Create the design list
Unequaldsgn <- list(Perennial=list(panel=c(PanelOne=75),
                                   seltype="Unequal",
                                   caty.n=c("1st"=25, "2nd"=25, "3rd+"=25),
                                   over=36),
                    Intermittent=list(panel=c(PanelOne=32),
                                      seltype="Unequal",
                                      caty.n=c("1st"=25, "2nd"=5, "3rd+"=2),
                                      over=0))



###################################################
### code chunk number 18: Unequalsites
###################################################
# Select the sample
Unequalsites <- grts(design=Unequaldsgn,
                     DesignID="UNEQUAL",
                     type.frame="linear",
                     src.frame="sp.object",
                     sp.object=shp,
                     att.frame=att,
                     stratum="Per_Int",
                     mdcaty="Strah_Cat",
                     shapefile=FALSE)



###################################################
### code chunk number 19: Unequalsites
###################################################
# Print the initial six lines of the survey design
head(Unequalsites@data)



###################################################
### code chunk number 20: Unequalsites
###################################################
# Print the survey design summary
summary(Unequalsites)



###################################################
### code chunk number 21: Panelsites
###################################################
# Create the design list
Paneldsgn <- list(Perennial=list(panel=c(Annual=16, Year1=17, Year2=17),
                                 seltype="Unequal",
                                 caty.n=c("1st"=15, "2nd"=15, "3rd+"=20),
                                 over=50),
                  Intermittent=list(panel=c(Annual=27),
                                    seltype="Unequal",
                                    caty.n=c("1st"=20, "2nd"=5, "3rd+"=2)))



###################################################
### code chunk number 22: Panelsites
###################################################
# Select the sample
Panelsites <- grts(design=Paneldsgn,
                   DesignID="UNEQUAL",
                   type.frame="linear",
                   src.frame="shapefile",
                   in.shape="Luck_Ash_streams",
                   att.frame=att,
                   stratum="Per_Int",
                   mdcaty="Strah_Cat",
                   shapefile=FALSE)



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



