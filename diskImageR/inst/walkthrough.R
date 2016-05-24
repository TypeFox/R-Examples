##############################################################
#Run this function the first time you want to run diskImageR
#(or any other package through github)
##############################################################
#install the devtools package
install.packages("devtools")

##############################################################
#Run the following functions everytime to use up diskImageR
#For all functions type ?functionName to bring up a help file
##############################################################
#load devtools
library(devtools)

#install the diskImageR package
install_github("acgerstein/diskImageR", build_vignettes = FALSE) 

#load diskImageR
library(diskImageR)

#Run the ImageJ analysis component, save the output. "newProject" shhould be changed to something of your choice (and then the same name used throughout); note that the quotation marks are required.
#To use a pop-up box interface:
IJMacro("newProject")
#To specify the appropriate directories without a popup:
IJMacro("newProject",  "/path/to/projectDir", "/path/to/projectDir/photographs/")

#Plot the result of ImageJ analysis (averaged among 72 lines draft outward from the center of the diffusion disk). Type ?plotRaw for additional parameter options.
plotRaw("newProject")

#Use maximum likelihood to fit a bilogistic and single logistic model to the data from each photograph. "clearHalo" is used to specify a picture that has a clear halo; this is used to standardize all photographs and will be most effective when photographs are taken with equal lighting without shadows. 
maxLik("newProject", clearHalo=1, RAD="all")

#Use the models to calculate resistance (20%, 50% and 80% reduction in growth = RAD20, RAD50, RAD80), tolerance (actual fraction of growth achvied above RAD relative to potential growth = FoG20, FoG50, FoG80), and sensitivity (slope at RAD50), which are saved in a .csv file. 
createDataframe("newProject", clearHalo = 1)

#[OPTIONAL] Calculate the mean and error for parameter estimates across replicate pictures
aggregateData("newProject")

#[OPTIONAL] Calculate corresponding MIC values from RAD values using either built-in data for specific species/drug combinations or supply your own RAD/MIC data to use for a standard curve
calcMIC("newProject")


