## ----setup, echo=FALSE, warning=FALSE------------------------------------
library(knitr)
library(EEM)
opts_chunk$set(fig.width=6.5, fig.height=4)

## ----readEEM, eval=FALSE-------------------------------------------------
#  library(EEM) # load library
#  
#  # read raw data files from a file
#  data <- readEEM(file)
#  data <- readEEM("sample1.txt") # read in a file
#  data <-readEEM(c("sample1.txt", "sample2.txt")) # read in two files
#  
#  # read raw data files from a folder
#  data <- readEEM(folder)
#  data <-readEEM("C:\\data") # full path. Note that the slash is doubled.
#  data <- readEEM("C:/data") # read in all files in data folder. Aside from double slashes,
#                             # a reverted slash can also be used.
#  
#  # read raw data files from the current working folder
#  setwd(choose.dir()) # set working folder interactively (only work in windows)
#  data <- readEEM(getwd()) # read raw data files in current working folder

## ----loadData------------------------------------------------------------
# load dataset
data(applejuice) 
class(applejuice) # EEM class

# get sample names
names(applejuice)

# use summary to see information about the dataset.
summary(applejuice)

## ----drawEEM-------------------------------------------------------------
# draw EEM of sample no.1
drawEEM(applejuice, n = 1) 

# draw EEM of sample no.1 with different color
drawEEM(applejuice, n = 1, color.palette = cm.colors) 

## ----drawEEM2------------------------------------------------------------
# flip the axis
drawEEM(applejuice, n = 1, flipaxis = TRUE) 

## ----drawEEMgg-----------------------------------------------------------
# draw EEM of sample no.1
drawEEMgg(applejuice, n = 1) 
 
# all functionalities in ggplot can be applied directly 
library(ggplot2)
# add grid line to the plot
drawEEMgg(applejuice, n = 1) + theme(panel.grid = element_line(color = "grey"), 
                                     panel.grid.major = element_line(colour = "grey"))

## ----delScattering1_1----------------------------------------------------
# delete scattering regions and assign them as NA
applejuice_delS <- delScattering(applejuice, rep = NA) 
drawEEM(applejuice_delS, 1)

## ----delScattering1_2----------------------------------------------------
applejuice_delS <- delScattering(applejuice, rep = NA, first = 30, second = 0, third = 0, forth = 0) 
drawEEM(applejuice_delS, 1)

## ----delScattering1_3----------------------------------------------------
applejuice_delS <- delScattering(applejuice, rep = 0, 
                                 first = 30, second = 0, third = 0, forth = 0) 

## ----delScattering2------------------------------------------------------
drawEEM(delScattering2(applejuice, NA), 1)

## ----cutEEM--------------------------------------------------------------
applejuice_delS_cut <- cutEEM(applejuice_delS, cutEX = 350:500, cutEM = 500:700)
drawEEM(applejuice_delS_cut, 1)

## ----unfold--------------------------------------------------------------
## unfold EEM into EEM_uf (matrix form with samples x variables dimension)
applejuice_delS_uf <- unfold(applejuice_delS) 

# dimension of unfolded data
dim(applejuice_delS_uf)

# take a look at unfolded data
applejuice_delS_uf[1:5 ,1:5]

## ----normalize-----------------------------------------------------------
# normalize data
applejuice_delS_uf_norm <- normalize(applejuice_delS_uf) 

# the absolute sum of each row should equal to 1
rowSums(abs(applejuice_delS_uf_norm)) 

## ----export, eval=FALSE--------------------------------------------------
#  # export as csv file
#  write.csv(applejuice_delS_uf, "applejuice.csv")

## ----pca-----------------------------------------------------------------
# perform PCA
result <- prcomp(applejuice_delS_uf_norm) # mean-centering is enabled by default

# plot scree plot
screeplot(result, npcs = 10, type = "lines", main = "Screeplot")

## ----scoreloading--------------------------------------------------------
# plot score plot 
plotScore(result, xPC = 1, yPC = 2) # pc 1 vs pc 2

# plot loading plot
plotLoading(result, ncomp = 1) # loading 1

## ----extractName---------------------------------------------------------
# extract sample name
sName <- names(applejuice) 

# country of apple production
country <- sapply(strsplit(sName, split = "-"), "[", 1) 
table(country) # counts of samples grouped by country

# cultivar of apples
cultivar <- sapply(strsplit(sName, split = "-"), "[", 2) 
table(cultivar) # counts of samples grouped by cultivar

## ----scoreg--------------------------------------------------------------
# plot score plot with grouping
plotScore(result, xPC = 1, yPC = 2,country, legendlocation = "topright")

# plot score using scatterplot matrix with grouping
plotScorem(result, ncomp = 5, country)
plotScorem(result, ncomp = 5, cultivar, cex = 1)

## ----pls, message=FALSE, warning=FALSE-----------------------------------
# load gluten data
data(gluten)
gluten_uf <- unfold(gluten) # unfold list into matrix

# delete columns with NA values
index <- colSums(is.na(gluten_uf)) == 0
gluten_uf <- gluten_uf[, index]
gluten_ratio <- as.numeric(names(gluten))

require(pls)
model <- plsr(gluten_ratio ~ gluten_uf, ncomp = 3)
plotLoading(model, ncomp = 3)
plotReg(model) 

