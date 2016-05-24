# test by William Revelle
# from Jennrich, Psychometrika, 2002, solution for the Thurstone 20 box problem. # Specifying 27 elements to be 0 as discussed in that article (Table 1 at 
# page 12) and using vgQ.target as revised or vgQ.pst with a W matrix 
# and Target as specified does not yield the reported solution.  
# The solution is almost identical for the high loadings but differs slightly
# for the small loadings.  The two models have a factor congruence of .99 for
# all three factors, but do not agree completely.

# Jennrich (2002) apparently was using the oblique rotation option.  
# When running TargetQ the results are fine, or when running 
#   the vgQ.pst function with  GPFoblq.

# This a good test case for both TargetQ  
# (It could also be adapted for pst but there is already a test for it.)

require("GPArotation")
data(Thurstone)   #the 20 box problem

#solution reported in Jennrich 2002 

browne <- t(matrix(c( 
   0.013,  0.994,  0.007,
   0.991,  0.012,  0.001,
   0.018,  0.003,  0.986,
   0.772,  0.477,  0.002,
   0.003,  0.393,  0.874,
   0.409,  0.003,  0.816,
   0.548,  0.730, -0.020,
   0.023,  0.870,  0.405,
   0.799, -0.024,  0.453,
   0.664,  0.621, -0.005,
  -0.058,  0.915,  0.512,
   0.639, -0.018,  0.644,
   0.046,  0.980, -0.003,
   0.971, -0.038,  0.060,
  -0.026,  0.025,  0.965,
   0.380,  0.281,  0.726,
   0.490,  0.652,  0.286,
  -0.025,  0.971,  0.019,
   0.957,  0.061, -0.045,
   0.028,  0.000,  0.976),
 3,20,dimnames = list(c("B1", "B2", "B3"), NULL)))

#a simplified target matrix, with NAs for ? and 0 for 0s.
#   (compare to pst appproach)
Target <- t(matrix(c(
   0, NA,  0,
  NA,  0,  0,
   0,  0, NA,
  NA, NA,  0,
   0, NA, NA,
  NA,  0, NA,
  NA, NA,  0,
   0, NA, NA,
  NA,  0, NA,
  NA, NA,  0,
   0, NA, NA,
  NA,  0, NA,
   0, NA,  0,
  NA,  0,  0,
   0,  0, NA,
  NA, NA, NA,
  NA, NA, NA,
   0, NA,  0,
  NA,  0,  0,
   0,  0, NA),
 3, 20, dimnames = list(c("T1", "T2", "T3"), NULL)))

v <- targetQ(box20,Target=Target)$loadings

all.ok <- TRUE  

#slightly larger fuzz for comparison with published value.
# note max(abs(v) - abs(browne))rather than max(abs(v - browne))
# as sign change is possible
 if( 10e-4 < max(abs(v) - abs(browne))) {
    cat("Calculated value is not the same as test value in Jennrich2002. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - browne, digits=18)
    all.ok <- FALSE  
    } 

good <- t(matrix(c( 
   0.01324194563970146343, -0.99360765277094842407,  0.007265459960371034587,
   0.99121314541487770544, -0.01178320700232154961,  0.000654586020267855506,
   0.01798447315534307256, -0.00266076852016330911,  0.985581004768931734361,
   0.77198435084052174915, -0.47723548341238952730,  0.001547735983967568618,
   0.00334198654247502835, -0.39290416948063611180,  0.874043793719835537814,
   0.40934347835281348349, -0.00274610551094590233,  0.815649888720176186041,
   0.54757055519984310088, -0.72951044925148011977, -0.020211353947714422175,
   0.02292379053779741716, -0.87011712730189194609,  0.404542252780873523577,
   0.79911058029224457666,  0.02416810475294199623,  0.452727043944761764482,
   0.66393502364020362538, -0.62149665012300570055, -0.005186928343372421146,
  -0.05839790682548451350, -0.91517931889838155524,  0.511521949806932663130,
   0.63924406199386740735,  0.01841750353525576159,  0.643544196342115570886,
   0.04597086497418309547, -0.97980801598321454193, -0.002918643110053173451,
   0.97103389549392915558,  0.03847065084578840666,  0.060066450372699808913,
  -0.02622776344285615568, -0.02482060086975104718,  0.965272709232911085842,
   0.37998105522582992233, -0.28073835673932595602,  0.726047993725112084107,
   0.48985182554738604388, -0.65226812910595410866,  0.285738966726349907788,
  -0.02451057644240206557, -0.97122042802717223342,  0.019132901654980147277,
   0.95708220223038309449, -0.06086293722346142188, -0.045050942196376064786,
   0.02797903728304645954,  0.00036458752733534161,  0.976083771686937051726),
 3,20,dimnames = list(c("B1", "B2", "B3"), NULL)))

#tighter fuzz for numerical comparison with previous test value
 if( 10e-12 < max(abs(v - good))) {
    cat("Calculated value is not the same as previous test value. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - good, digits=18)
    all.ok <- FALSE  
    } 

 
cat("tests completed.\n")


if (! all.ok) stop("some tests FAILED")
