### R code from vignette source 'Finite_Design.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: createshape
###################################################
# Load the sp object in the data directory
data(NE_lakes)

# Create a shapefile
sp2shape(sp.obj=NE_lakes, shpfilename="NE_lakes")



###################################################
### code chunk number 3: att
###################################################
# Read the attribute table from the shapefile
att <- read.dbf("NE_lakes")



###################################################
### code chunk number 4: att
###################################################
# Display the initial six lines in the attribute data frame
head(att)



###################################################
### code chunk number 5: att
###################################################
# Display number of lakes cross-classified by strata and multidensity
# category
addmargins(table("State"=att$State, "Lake Area Category"=att$Area_Cat))



###################################################
### code chunk number 6: Equalsites
###################################################
# Call the set.seed function so that the survey designs can be replicate
set.seed(4447864)



###################################################
### code chunk number 7: Equalsites
###################################################
# Create the design list
Equaldsgn <- list(None=list(panel=c(PanelOne=300), seltype="Equal"))



###################################################
### code chunk number 8: Equalsites
###################################################
# Select the sample
Equalsites <- grts(design=Equaldsgn,
                   DesignID="EQUAL",
                   type.frame="finite",
                   src.frame="shapefile",
                   in.shape="NE_lakes", 
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
### code chunk number 11: Stratsites
###################################################
# Create the design list
Stratdsgn <- list(CT=list(panel=c(PanelOne=125), seltype="Equal"),
                  MA=list(panel=c(PanelOne=125), seltype="Equal"),
                  RI=list(panel=c(PanelOne=50), seltype="Equal"))



###################################################
### code chunk number 12: Stratsites
###################################################
# Select the sample
Stratsites <- grts(design=Stratdsgn,
                   DesignID="STRATIFIED",
                   type.frame="finite",
                   src.frame="att.frame",
                   att.frame=att,
                   xcoord="xcoord",
                   ycoord="ycoord",
                   stratum="State",
                   shapefile=FALSE)



###################################################
### code chunk number 13: Stratsites
###################################################
# Print the initial six lines of the survey design
head(Stratsites@data)



###################################################
### code chunk number 14: Stratsites
###################################################
# Print the survey design summary
summary(Stratsites)



###################################################
### code chunk number 15: Unequalsites
###################################################
# Read the shapefile
shp <- read.shape("NE_lakes")



###################################################
### code chunk number 16: Unequalsites
###################################################
# Create the design list
Unequaldsgn <- list(None=list(panel=c(PanelOne=300),
                              seltype="Unequal",
                              caty.n=c("(0,1]"=50, "(1,5]"=120, "(5,10]"=50,
                                       "(10,50]"=50, "(50,500]"=25,
                                       "(500,1e+04]"=5),
                              over=120))



###################################################
### code chunk number 17: Unequalsites
###################################################
# Select the sample
Unequalsites <- grts(design=Unequaldsgn,
                     DesignID="UNEQUAL",
                     type.frame="finite",
                     src.frame="sp.object",
                     sp.object=shp,
                     att.frame=att,
                     mdcaty="Area_Cat",
                     shapefile=FALSE)



###################################################
### code chunk number 18: Unequalsites
###################################################
# Print the initial six lines of the survey design
head(Unequalsites@data)



###################################################
### code chunk number 19: Unequalsites
###################################################
# Print the survey design summary
summary(Unequalsites)



###################################################
### code chunk number 20: Panelsites
###################################################
# Create the design list
Paneldsgn <- list(None=list(panel=c(Annual=50, Year1=50, Year2=50, Year3=50,
                                    Year4=50, Year5=50),
                            seltype="Unequal",
                            caty.n=c("(0,1]"=50, "(1,5]"=120, "(5,10]"=50,
                                     "(10,50]"=50, "(50,500]"=25,
                                     "(500,1e+04]"=5),
                            over=100))



###################################################
### code chunk number 21: Panelsites
###################################################
# Select the sample
Panelsites <- grts(design=Paneldsgn,
                   DesignID="UNEQUAL",
                   type.frame="finite",
                   src.frame="shapefile",
                   in.shape="NE_lakes",
                   att.frame=att,
                   mdcaty="Area_Cat",
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



