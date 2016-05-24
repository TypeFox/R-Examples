#Example from: Gradient Projection Algorithms and Software for 
#  Arbitrary Rotation Criteria in Factor Analysis.
#  by Coen A. Bernaards and Robert I. Jennrich
#  Website: http://www.stat.ucla.edu/research

 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()


fuzz <- 1e-5 # using  eps=1e-5 these tests do not do better than this
all.ok <- TRUE  

#  quartimax (orthogonal) rotation of Harman's 8 physical variables.

data("Harman", package="GPArotation")

qHarman  <- GPForth(Harman8, Tmat=diag(2), method="quartimax")
qHarman2 <- quartimax(Harman8) 

 if( fuzz < max(abs(qHarman$loadings - qHarman2$loadings))) {
    cat("Calculated value is not the same as test value in test Harman 1. Value:\n")
    print(qHarman$loadings, digits=18)
    cat("difference:\n")
    print(qHarman$loadings - qHarman2$loadings, digits=18)
    all.ok <- FALSE  
    } 

#qHarman$Th - qHarman2$Th

# with eps=1e-8
# tst <- t(matrix(c(
#  0.898754567491920398, 0.194823580226859222,
#  0.933943406208487592, 0.129748657024604030,
#  0.902131483644799892, 0.103864268239045668,
#  0.876508251941102934, 0.171284220753554678,
#  0.315572019798302239, 0.876476069451083251,
#  0.251123191235179066, 0.773488941629975613,
#  0.198007116064591759, 0.714678376605717203,
#  0.307857241091366252, 0.659334451631046314
#  ), 2, 8))

# with eps=1e-5
 tst <- t(matrix(c(
  0.898755404698461491, 0.194819718009510034,
  0.933943963768413821, 0.129744643590955028,
  0.902131929972106672, 0.103860391510923730,
  0.876508987992224209, 0.171280454135453869,
  0.315575786273609882, 0.876474713336210853,
  0.251126515144778573, 0.773487862471829213,
  0.198010187248201075, 0.714677525703678707,
  0.307860074444663512, 0.659333128670876345
   ), 2, 8))

 if( fuzz < max(abs(qHarman$loadings - tst ))) {
    cat("Calculated value is not the same as test value in test Harman 2. Value:\n")
    print(qHarman$loadings, digits=18)
    cat("difference:\n")
    print(qHarman$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 

# with eps=1e-8
# tst <- t(matrix(c(
#   0.790828307905322436, 0.612038060430562525,
#  -0.612038060430562525, 0.790828307905322214
#  ), 2, 2))

# with eps=1e-5
 tst <- t(matrix(c(
   0.790830938007507367, 0.612034662000581764,
  -0.612034662000581764, 0.790830938007507145
  ), 2, 2))
 
 if( fuzz < max(abs(qHarman$Th - tst ))) {
    cat("Calculated value is not the same as test value in test Harman 3. Value:\n")
    print(qHarman$Th, digits=18)
    cat("difference:\n")
    print(qHarman$Th - tst, digits=18)
    all.ok <- FALSE  
    } 





#  quartimin (oblique) rotation of Harman's 8 physical variables.

qminHarman  <- GPFoblq(Harman8, Tmat=diag(2), method="quartimin")
qminHarman2 <- quartimin(Harman8) 

 if( fuzz < max(abs(qminHarman$loadings - qminHarman2$loadings))) {
    cat("Calculated value is not the same as test value in test Harman 4. Value:\n")
    print(qminHarman$loadings, digits=18)
    cat("difference:\n")
    print(qminHarman$loadings - qminHarman2$loadings, digits=18)
    all.ok <- FALSE  
    } 


# with eps=1e-8
# tst <- t(matrix(c(
#   0.8918217697289939627,  0.0560146456758183961,
#   0.9536799985772628219, -0.0232460005406671701,
#   0.9291498623396581280, -0.0465027396531852502,
#   0.8766828510822184395,  0.0336582451338717017,
#   0.0136988312985193428,  0.9250013826349388069,
#  -0.0172668087945964319,  0.8212535444941218010,
#  -0.0524468998178311899,  0.7649536381341245361,
#   0.0858880630098148856,  0.6831160953442911854
#   ),2, 8))					

# with eps=1e-5
 tst <- t(matrix(c(
   0.8918219293548808047,  0.0560145122875230911,
   0.9536799846795966928, -0.0232460559140742311,
   0.9291497958388006406, -0.0465027685653178480,
   0.8766829604751505967,  0.0336581364763500201,
   0.0137008854716444972,  0.9250004106413580729,
  -0.0172649861805529957,  0.8212526839806429946,
  -0.0524452035885302342,  0.7649528396536503516,
   0.0858895830186393733,  0.6831153711863455769
   ),2, 8))				       
					       
 if( fuzz < max(abs(qminHarman$loadings - tst ))) {
    cat("Calculated value is not the same as test value in test Harman 5. Value:\n")
    print(qminHarman$loadings, digits=18)
    cat("difference:\n")
    print(qminHarman$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 


# with eps=1e-8
# tst <- t(matrix(c(
#  1.000000000000000000, 0.472747617396915065,
#  0.472747617396915065, 1.000000000000000000
#   ),2, 2))				       

# with eps=1e-5
 tst <- t(matrix(c(
  1.000000000000000222, 0.472745958387102538,
  0.472745958387102538, 1.000000000000000000
   ),2, 2))				       
					       
 if( fuzz < max(abs(qminHarman$Phi - tst ))) {
    cat("Calculated value is not the same as test value in test Harman 6. Value:\n")
    print(qminHarman$Phi, digits=18)
    cat("difference:\n")
    print(qminHarman$Phi - tst, digits=18)
    all.ok <- FALSE  
    } 


# with eps=1e-8
# tst <- t(matrix(c(
#   0.878125245495924522, 0.836723841642554422,
#  -0.478430823863515542, 0.547625065922776710
#   ),2, 2))				       

# with eps=1e-5
 tst <- t(matrix(c(
   0.878125280760480686, 0.836722770276292271,
  -0.478430759137962514, 0.547626702874473570
   ),2, 2))				       
					       
 if( fuzz < max(abs(qminHarman$Th - tst ))) {
    cat("Calculated value is not the same as test value in test Harman 7. Value:\n")
    print(qminHarman$Th, digits=18)
    cat("difference:\n")
    print(qminHarman$Th - tst, digits=18)
    all.ok <- FALSE  
    } 


cat("tests completed.\n")

if (! all.ok) stop("some tests FAILED")

