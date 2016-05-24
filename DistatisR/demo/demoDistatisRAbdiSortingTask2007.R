#
# How To Analyze the Results of a Sorting Task
#     in R with the package DistatisR
#     herve@utdallas.edu, www.utdallas.edu/~herve
#     
#     Library needed
#     can be loaded with the commands: 
#     install.packages('DistatisR')
#     followed by 
#     library('DistatisR)
#     Herve Abdi. November 7, 2012
#      New Version of DISTATIS
#
# Thanks to Derek Beaton (for help with prettyGraphs)
#       and Cherise Chin Fatt (for help with distatis)
# This demo shows how to analyze the data from Table 1 of
# Abdi, H., Valentin, D., Chollet, S., & Chrea, C. (2007).
#  Analyzing assessors and products in sorting tasks: DISTATIS, theory and applications. 
#  Food Quality and Preference, 18, 627-640.
# Available from www.utdallas.edu/~herve
# as paper # A59
#

#-----------------------------------------------------------------------------
# THis command is assume to have been issued prior to running the demo
# library(DistatisR)
#-----------------------------------------------------------------------------
#  1. Get the data from the 2007 sorting example
#      this is the eay they look from Table 1 of 
#      Abdi et al. (2007).
#                       Assessors
#                  1 2 3 4 5 6 7 8 9 10
# Beer        Sex  f m f f m m m m f m
#            -----------------------------                         
#Affligen          1 4 3 4 1 1 2 2 1 3
#Budweiser         4 5 2 5 2 3 1 1 4 3
#Buckler_Blonde    3 1 2 3 2 4 3 1 1 2
#Killian           4 2 3 3 1 1 1 2 1 4
#St. Landelin      1 5 3 5 2 1 1 2 1 3
#Buckler_Highland  2 3 1 1 3 5 4 4 3 1
#Fruit Defendu     1 4 3 4 1 1 2 2 2 4
#EKU28             5 2 4 2 4 2 5 3 4 5

# 1.1. Create the
#     Name of the Beers
BeerName <- c('Affligen', 'Budweiser','Buckler Blonde',
              'Killian','St. Landelin','Buckler Highland',
              'Fruit Defendu','EKU28')
# 1.2. Create the name of the Assessors 
#      (F are females, M are males)
Juges <- c('F1','M2', 'F3', 'F4', 'M5', 'M6', 'M7', 'M8', 'F9', 'M10')

# 1.3. Get the sorting data
SortData <- c(1, 4, 3, 4, 1, 1, 2, 2, 1, 3,
              4, 5, 2, 5, 2, 3, 1, 1, 4, 3,
              3, 1, 2, 3, 2, 4, 3, 1, 1, 2,
              4, 2, 3, 3, 1, 1, 1, 2, 1, 4,
              1, 5, 3, 5, 2, 1, 1, 2, 1, 3,
              2, 3, 1, 1, 3, 5, 4, 4, 3, 1,
              1, 4, 3, 4, 1, 1, 2, 2, 2, 4,
              5, 2, 4, 2, 4, 2, 5, 3, 4, 5)
# 1.4 Create a data frame
#     (alternatively we could have read a csv file)            
Sort <- matrix(SortData,ncol = 10, byrow= TRUE, dimnames = list(BeerName, Juges))

# 1.5 Example of how to read a csv filw
#Sort <- read.table("JeuxEpicesFrance.csv", header=TRUE, 
#   sep=",", na.strings="NA", dec=".", row.names=1, strip.white=TRUE)

#-----------------------------------------------------------------------------
# 2. Create the set of distance matrices (one distance matrix per assessor)
#    (ues the function DistanceFromSort)
DistanceCube <- DistanceFromSort(Sort)

#-----------------------------------------------------------------------------
# 3. Call the DISTATIS routine with the cube of distance as parameter
testDistatis <- distatis(DistanceCube)
# The factor scores for the beers are in
# testDistatis$res4Splus$F
# the factor scores for the assessors are in (RV matrice)
#  testDistatis$res4Cmat$G

#-----------------------------------------------------------------------------
# 4. Inferences on the beers obtained via bootstrap
#    here we use two different bootstraps:
#    1. Bootstrap on factors (very fast but could be too liberal 
#         when the number of assessors is very large)
#    2. Complete bootstrap obtained by computing sets of compromises
#       and projecting them (could be significantly longer because a lot
#       of computations is required)
# 
# 4.1 Get the bootstrap factor scores (with default 1000 iterations)
BootF <- BootFactorScores(testDistatis$res4Splus$PartialF)
#
# 4.2 Get the boostrap from full bootstrap (default niter = 1000)
 F_fullBoot <- BootFromCompromise(DistanceCube,niter=1000) 
 

#-----------------------------------------------------------------------------
# 5. Create the Graphics
# Test the call to the graphroutines
# Get the Factor Scores and Partial Factor Scores for the plot Routine
 LeF       <- testDistatis$res4Splus$F
 PartialFS <- testDistatis$res4Splus$PartialF
# 5.1. plot the Observations with the Bootstrapped CI from the factor scores
  PlotOfSplus <- GraphBootDistatisCpt(LeF,BootF,PartialFS,ZeTitle='Bootstrap on Factors')
# 5.2. Plot the Observations with the Bootstrapped CI from the Compromises
  PlotOfSplusCpt <- GraphBootDistatisCpt(LeF, F_fullBoot,PartialFS,ZeTitle='Full Bootstrap')
# 5.3 Plot the RV map
PlotOfRvMat <- GraphDistatisRv(testDistatis$res4Cmat$G,ZeTitle='Rv Mat')
# Et voila!


