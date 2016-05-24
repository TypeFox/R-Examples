
 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()

#require("stats")  

fuzz <- 1e-6 
all.ok <- TRUE  


# test MASS 4th ed. p 322-324

  data(ability.cov)
  ability.cov
  ability.FA <- factanal(factors = 1, covmat=ability.cov)

  (ability.FA <- update(ability.FA, factors = 2))

#  ability.FA2 <- factanal(factors = 2, covmat = ability.cov)
# max(abs(ability.FA2$loadings - ability.FA$loadings))

#  summary(ability.FA) MASS ed.4 p 323 seems to be print not summary in R 2.0.1
  ability.FA

# this is default varimax rotation. There are 3rd+ digit differences with MASS
  tst <- t(matrix(c(
      0.499437829039896530, 0.54344904693111962,
      0.156070079431279873, 0.62153798991197484,
      0.205786989958578748, 0.85992588538426895,
      0.108530754440558652, 0.46776101732283504,
      0.956242470279811574, 0.18209631992182243,
      0.784768183877880943, 0.22482213687364205
      ), 2, 6))
 

 if( fuzz < max(abs(loadings(ability.FA) - tst))) {
    cat("Calculated value is not the same as test value in test 1. Value:\n")
    #print(loadings(ability.FA), digits=18) this truncates
    print(unclass(ability.FA$loadings), digits=18)
    cat("difference:\n")
    print(unclass(ability.FA$loadings) - tst, digits=18)
    all.ok <- FALSE  
    } 


  # differences with MASS here are a bit more than might be expected,
  # but there is already a difference before rotation.
  (oblirot <- oblimin(loadings(ability.FA)))

  obli2 <- factanal(factors = 2, covmat = ability.cov, rotation="oblimin")

  max(abs(loadings(oblirot) - loadings(obli2)))


# factanal(factors = 2, covmat = ability.cov, scores = Bartlett, rotation="oblimin")


 tst <- t(matrix(c(
       0.3863637969729337152,  0.4745113977203344047,
      -0.0110032278171669998,  0.6458708261423832253,
      -0.0262888675561207576,  0.8961123879025085781,
      -0.0180180060207963122,  0.4882918937716873575,
       0.9900948712271664398, -0.0370729040114848238,
       0.7905663749272058283,  0.0526099352008769991
    ), 2, 6))
 
 if( fuzz < max(abs(loadings(oblirot) - tst ))) {
    cat("Calculated value is not the same as test value in test 2. Value:\n")
    print(loadings(oblirot), digits=18)
    cat("difference:\n")
    print(loadings(oblirot) - tst, digits=18)
    all.ok <- FALSE  
    } 

cat("tests completed.\n")

if (! all.ok) stop("some tests FAILED")
