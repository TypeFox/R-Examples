# Distatis Example. 4 Algorithms evaluate the similarity of 6 faces
# Data from Abdi et al. (2005)
#-----------------------------------------------------------------------------
# 1. Load the DistAlgo data set (available from the DistatisR package)
data(DistAlgo)
# DistAlgo is a 6*6*4 Array (face*face*Algorithm)
#-----------------------------------------------------------------------------
# 2. Call the DISTATIS routine with the array of distance (DistAlgo) as parameter
DistatisAlgo <- distatis(DistAlgo)
# The factor scores for faces are in
# DistatisAlgo$res4Splus$F
# the factor scores for the algorithms are in (RV matrice)
# DistatisAlgo$res4Cmat$G
# We do not perform Bootstrap here because N = 4 is too small 
#-----------------------------------------------------------------------------
# We are simply plotting the RV coefficient map
# and the Compromise plus partial factor scores
#
#-----------------------------------------------------------------------------
# 3. Create the Graphics
# Test the call to the graphroutines
# Get the Factor Scores and Partial Factor Scores for the plot Routine
 LeF       <- DistatisAlgo$res4Splus$F  # LeF:  Factor scores for the faces
 PartialFS <- DistatisAlgo$res4Splus$PartialF 
     # PartialFS: Partial factor scores (one per algorithm)
 LeG  <- DistatisAlgo$res4Cmat$G 
     # LeG the factor scores (derived from the Rv mat) for the algorithms
 # 3.1. Plot the Rv Map
 PlotOfRvMat <- GraphDistatisRv(LeG,ZeTitle='Algorithms: Rv Mat')
 # 3.1. PLot the compromise map
 GraphDistatisCpt(LeF,PartialFS,axis1=1,axis2=2,ZeTitle= 'Distatis-Compromise',nude=FALSE)
 # 3.2. Have the graphs without legend for further editing
  GraphDistatisCpt(LeF,PartialFS,axis1=1,axis2=2,ZeTitle= 'Distatis-Compromise',nude=TRUE)