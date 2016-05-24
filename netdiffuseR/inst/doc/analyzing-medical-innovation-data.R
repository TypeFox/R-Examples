## ----Setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
library(knitr)
opts_chunk$set(out.width = 600, fig.align = "center", fig.width = 7, fig.height = 7) 

## ----Loading pkgs and reading the data-----------------------------------
library(netdiffuseR)
data(medInnovations)

## ----Preparing data for netdiffuseR--------------------------------------
# Creating unique ids (including for the network data)
othervars <- c("id", "toa", "city")
netvars <- names(medInnovations)[grepl("^net", names(medInnovations))]
for (i in c("id", netvars))
  medInnovations[[i]] <- medInnovations[[i]] + medInnovations$city*1000

# Leaving unsurveyed individuals with NA
surveyed <- medInnovations$id
for (i in netvars)
  medInnovations[[i]][which(!(medInnovations[[i]] %in% surveyed))] <- NA

# Reshaping data (so we have an edgelist)
medInnovations.long <- reshape(
  medInnovations[,c(othervars, netvars)], v.names= "net",
  varying = netvars,
  timevar = "level", idvar="id", direction="long")

## ----Importing data to netdiffuseR---------------------------------------
# Coersing the edgelist to an adjacency matrix. Here we are assuming that the
# network is constant through time.
graph <- with(
  medInnovations.long,
  edgelist_to_adjmat(cbind(id, net), t=18,undirected=FALSE, keep.isolates = TRUE)
)

## ------------------------------------------------------------------------
# Just to be sure. Sorting the data!
orddata <- as.numeric(as.factor(rownames(graph[[1]])))
medInnovations <- medInnovations[orddata,]

# Creating a diffnet object
diffnet <- as_diffnet(graph, medInnovations$toa, 
                      vertex.static.attrs = subset(medInnovations, select=c(-id, -toa)))

## ----Checking-the-methods------------------------------------------------
plot(diffnet, t=diffnet$meta$nper)
diffnet
summary(diffnet)

## ----graphs--------------------------------------------------------------
plot_diffnet(diffnet, slices=c(1,4,8,12,16,18))
plot_threshold(diffnet, undirected = FALSE, vertex.cex = 1/5)
plot_adopters(diffnet)
plot_hazard(diffnet)

## ----Subsetting----------------------------------------------------------
# Get cities ids so we can subset the vertices and run the test by city.
city <- diffnet$vertex.static.attrs[,"city"]

# Subsetting diffnet, notice that we can use either indices or ids to create a
# "subdiffnet". In this case we are using indices.
diffnet_city1 <- diffnet[which(city==1),]
diffnet_city2 <- diffnet[which(city==2),]
diffnet_city3 <- diffnet[which(city==3),]
diffnet_city4 <- diffnet[which(city==4),]

## ----Multithreshold------------------------------------------------------
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
plot_threshold(diffnet_city1, vertex.label = "", main="Threshold and ToA\nin City 1")
plot_threshold(diffnet_city2, vertex.label = "", main="Threshold and ToA\nin City 2")
plot_threshold(diffnet_city3, vertex.label = "", main="Threshold and ToA\nin City 3")
plot_threshold(diffnet_city4, vertex.label = "", main="Threshold and ToA\nin City 4")
par(oldpar)

## ------------------------------------------------------------------------
plot_infectsuscep(diffnet_city1, K=5, logscale = TRUE, bins=20)

## ----Stat-test-----------------------------------------------------------
# Defining the statistic
avgthr <- function(x) mean(threshold(x), na.rm = TRUE)

# Running the test by city
test1   <- struct_test(diffnet_city1, avgthr, 2000, ncpus=2, parallel="multicore")
test2   <- struct_test(diffnet_city2, avgthr, 2000, ncpus=2, parallel="multicore")
test3   <- struct_test(diffnet_city3, avgthr, 2000)
test4   <- struct_test(diffnet_city4, avgthr, 2000)

# Running the test aggregated
testall <- struct_test(diffnet, avgthr, 2000, ncpus=2, parallel="multicore")

# Printing the outcomes
test1
test2
test3
test4
testall

## ----Histograms----------------------------------------------------------
# To make it nicer, we change the parameters in using par
# (see ?par)
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,2))

# Now we use the hist method for the -diffnet_boot- class
hist(test1, main="Distribution of Statistic on rewired\nnetwork (City 1)")
hist(test2, main="Distribution of Statistic on rewired\nnetwork (City 2)")
hist(test3, main="Distribution of Statistic on rewired\nnetwork (City 3)")
hist(test4, main="Distribution of Statistic on rewired\nnetwork (City 4)")

# Returning to the previous graphical parameters
par(oldpar)

## ----Computing-exposure--------------------------------------------------
# Calculating exposure
expo <- exposure(diffnet)
head(expo)

# Netdiffuser automatically identifies whether the input is dynamic or not.
diffnet[["netexp"]] <- expo

## ------------------------------------------------------------------------
mydata <- diffnet.attrs(diffnet, as.df = TRUE)

