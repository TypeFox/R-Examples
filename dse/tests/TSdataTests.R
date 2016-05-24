###############################################################

####   Test TSdata methods

###############################################################

 Sys.getenv("R_LIBS")
 library()
 
 require("dse")
 search()

 Sys.info()
 DSEversion() 

###############################################################
 all.ok <- TRUE  

 cat("TSdataTests test 1 ...\n")
   z <- TSdata(output=matrix(rnorm(300), 100,3))
   nseriesOutput(z);	   ok <- nseriesOutput(z) == 3
   nseriesInput(z);	   ok <- ok & nseriesInput(z) == 0
   Tobs(outputData(z)); ok <- ok & Tobs(outputData(z)) == 100
   start(outputData(z));   ok <- ok & all(start(outputData(z)) == c(1,1))
   end(outputData(z));     ok <- ok & all(end(outputData(z)) == c(100,1))
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test 2 ...\n")
   z <- TSdata(output=matrix(rnorm(300), 100,3), input=matrix(rnorm(200), 100,2))
   nseriesInput(z);	   ok <- ok & nseriesInput(z) == 2
   nseriesOutput(z);	   ok <- nseriesOutput(z) == 3
   Tobs(inputData(z));  ok <- ok & Tobs(inputData(z)) == 100
   Tobs(outputData(z)); ok <- ok & Tobs(outputData(z)) == 100
   start(inputData(z));    ok <- ok & all(start(inputData(z)) == c(1,1))
   start(outputData(z));   ok <- ok & all(start(outputData(z)) == c(1,1))
   end(inputData(z));      ok <- ok & all(end(inputData(z)) == c(100,1))
   end(outputData(z));     ok <- ok & all(end(outputData(z)) == c(100,1))
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test 3 ...\n")
   z <- TSdata(output=matrix(rnorm(300), 100,3), input=rnorm(200))
   nseriesInput(z);	   ok <- ok & nseriesInput(z) == 1
   nseriesOutput(z);	   ok <- nseriesOutput(z) == 3
   Tobs(inputData(z));  ok <- ok & Tobs(inputData(z)) == 200
   Tobs(outputData(z)); ok <- ok & Tobs(outputData(z)) == 100
   start(inputData(z));    ok <- ok & all(start(inputData(z)) == c(1,1))
   start(outputData(z));   ok <- ok & all(start(outputData(z)) == c(1,1))
   end(inputData(z));      ok <- ok & all(end(inputData(z)) == c(200,1))
   end(outputData(z));     ok <- ok & all(end(outputData(z)) == c(100,1))
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test 4 ...\n")
   z <- TSdata(output=ts(matrix(rnorm(300), 100,3), start=c(1976,1), frequency=12),
               input=ts(rnorm(200), start=c(1976,1), frequency=12))
   nseries(inputData(z));    ok <- nseries(inputData(z)) == 1
   nseries(outputData(z));   ok <- ok & nseries(outputData(z)) == 3
   nseriesInput(z);	     ok <- ok & nseriesInput(z) == nseries(inputData(z))
   nseriesOutput(z);	     ok <- ok & nseriesOutput(z) == nseries(outputData(z))
   Tobs(inputData(z));    ok <- ok & Tobs(inputData(z)) == 200
   Tobs(outputData(z));   ok <- ok & Tobs(outputData(z)) == 100
   TobsInput(z);          ok <- ok & Tobs(inputData(z)) == TobsInput(z)
   TobsOutput(z);         ok <- ok & Tobs(outputData(z)) == TobsOutput(z)
   frequency(inputData(z));  ok <- ok & frequency(inputData(z)) == 12
   frequency(outputData(z)); ok <- ok & frequency(outputData(z)) == 12
   frequencyInput(z);        ok <- ok & frequency(inputData(z)) == frequencyInput(z)
   frequencyOutput(z);       ok <- ok & frequency(outputData(z)) == frequencyOutput(z)
   start(inputData(z));      ok <- ok & all(start(inputData(z)) == c(1976,1))
   start(outputData(z));     ok <- ok & all(start(outputData(z)) == c(1976,1))
   startInput(z);            ok <- ok & all(start(inputData(z)) == startInput(z))
   startOutput(z);           ok <- ok & all(start(outputData(z)) == startOutput(z))
   end(inputData(z));        ok <- ok & all(end(inputData(z)) == c(1992,8))
   end(outputData(z));       ok <- ok & all(end(outputData(z)) == c(1984,4))
   endInput(z);              ok <- ok & all(end(inputData(z)) == endInput(z))
   endOutput(z);             ok <- ok & all(end(outputData(z)) == endOutput(z))
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test 5 ...\n")
   z <- TSdata(output=ts(matrix(rnorm(300), 100,3), start=c(1976,1), frequency=4),
               input=ts(rnorm(200), start=c(1976,1), frequency=4))
   nseries(inputData(z));    ok <- nseries(inputData(z)) == 1
   nseries(outputData(z));   ok <- ok & nseries(outputData(z)) == 3
   nseriesInput(z);	     ok <- ok & nseriesInput(z) == nseries(inputData(z))
   nseriesOutput(z);	     ok <- ok & nseriesOutput(z) == nseries(outputData(z))
   Tobs(inputData(z));    ok <- ok & Tobs(inputData(z)) == 200
   Tobs(outputData(z));   ok <- ok & Tobs(outputData(z)) == 100
   TobsInput(z);          ok <- ok & Tobs(inputData(z)) == TobsInput(z)
   TobsOutput(z);         ok <- ok & Tobs(outputData(z)) == TobsOutput(z)
   frequency(inputData(z));  ok <- ok & frequency(inputData(z)) == 4
   frequency(outputData(z)); ok <- ok & frequency(outputData(z)) == 4
   frequencyInput(z);        ok <- ok & frequency(inputData(z)) == frequencyInput(z)
   frequencyOutput(z);       ok <- ok & frequency(outputData(z)) == frequencyOutput(z)
   start(inputData(z));      ok <- ok & all(start(inputData(z)) == c(1976,1))
   start(outputData(z));     ok <- ok & all(start(outputData(z)) == c(1976,1))
   startInput(z);            ok <- ok & all(start(inputData(z)) == startInput(z))
   startOutput(z);           ok <- ok & all(start(outputData(z)) == startOutput(z))
   end(inputData(z));        ok <- ok & all(end(inputData(z)) == c(2025,4))
   end(outputData(z));       ok <- ok & all(end(outputData(z)) == c(2000,4))
   endInput(z);              ok <- ok & all(end(inputData(z)) == endInput(z))
   endOutput(z);             ok <- ok & all(end(outputData(z)) == endOutput(z))
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test 6 ...\n")
   zz <- TSdata(output=ts(matrix(rnorm(600), 200,3), start=c(1976,1), frequency=4),
               input=ts(rnorm(100), start=c(1976,1), frequency=4))
   inputData(zz)  <- inputData(z)
   outputData(zz) <- outputData(z)
   ok <- testEqual(z,zz)
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test eg 1 ...\n")
   data("eg1.DSE.data", package="dse") 
   if (!is.TSdata(eg1.DSE.data)) {
     cat("example data not found.")
     ok <- FALSE
     }
 if (!ok) {cat("ERROR in this test.") ; all.ok <- FALSE  }

 cat("TSdataTests test eg 2 ...\n")
   # there will only be an error indicated if these fail to work
   summary(eg1.DSE.data) 
   print(eg1.DSE.data) 
   


 if (! all.ok) stop("some tests FAILED")

