if(!require("stats"))   stop("this test requires stats.")
if(!require("dse"))  stop("this test requires EvalEst.")
if(!require("setRNG"))stop("this test requires setRNG.")
 #x11()
  dir <- tempdir()
  postscript(file=paste(dir,"/lite.out.ps", sep=""),  paper="letter",
      horizontal=FALSE, onefile=TRUE)
             # width=6, height=8, pointsize=10,
   Sys.info()
   DSEversion()
   

#######################################################################

#    test functions for examples in the Brief User's Guide  Part 1

#######################################################################


guide.example.tests.part1 <- function( verbose=TRUE, synopsis=TRUE, fuzz.small=1e-14,
    fuzz.large=1e-8, graphics=TRUE)
{# test examples in Brief User's guide
 # NOTE: it was necessary to reduce fuzz from 1e-14 because of differences
 # in the results between Splus 3.2 and Splus 3.3 (C libraries were changed).
 # Differences affected lsfit (used in estVARXls) among other things.


  # If no device is active then write to postscript file 
  if (graphics)
   {if ( dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=FALSE, print.it=FALSE, append=FALSE)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
    else
      {old.par <- par()
       #on.exit(par(old.par))
      }
   }


  max.error <- NA
  if (synopsis & !verbose) cat("All Brief User Guide example part 1 tests ...")

  if (verbose) cat("Guide part 1 test 0 ... ")
  data("eg1.DSE.data.diff", package="dse")
  data("eg1.DSE.data", package="dse")
  data("egJofF.1dec93.data", package="dse")   
  if (verbose) { cat("ok\n") }
     
  if (verbose) cat("Guide part 1 test 1a ... ")

  seriesNamesInput(eg1.DSE.data)   <-  "u1"
  seriesNamesOutput(eg1.DSE.data) <-  c("y1","y2","y3")
  if (! all(c(seriesNamesInput(eg1.DSE.data) == "u1" ,
              seriesNamesOutput(eg1.DSE.data) == c("y1","y2","y3") )))
     stop("seriesNamesInput, seriesNamesOutput specification failed.")

  if (verbose) cat("Guide part 1 test 1b ... ")

  error <- abs(126943980.50000011921 - sum(outputData(eg1.DSE.data)))
  ok <- 100*fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("Guide part 1 test 2 ... ")
  model1 <- estVARXls(eg1.DSE.data, warn=FALSE)
  model2 <- estSSMittnik(eg1.DSE.data, n=14)
#  summary(model1)
#  summary(model2)
#  print(model1)
#  print(model2)
#  stability(model1)
#  stability(model2)
   if (graphics) tfplot(model1)

# In the next test
# with svd  (in Splus and previously in R)
#             sum(TSmodel(model2)$F)= -1.1078692933906153506 and 
# with La.svd sum(TSmodel(model2)$F)=  3.9469252417636165
# ACML BLAS on an amd athlon64 it is   4.277158111324035
# These differences seems fairly large, but the matrix is 14x14 and 
# the roots are almost identical

  good <-     c(15.430979953081722655, 3.9469252417636165, 2.4561249653768193468) 
 #Splus value c(15.430979953081722655, -1.1078692933906153506, 2.4561249653768193468)

  test.value <- c(sum(TSmodel(model1)$A), sum(TSmodel(model2)$F),
                         sum(roots(model2)) )
  error <- Mod(good - test.value)
  ok <- any(c(fuzz.large, 1.0,fuzz.large ) > error)
  if (!ok) {print(test.value, digits=16)
            if (is.na(max.error)) max.error <- max(error)
            else max.error <- max(error, max.error)}
  ok <-  ok & is.TSestModel(model1) & is.TSestModel(model2)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("Guide part 1 test 3 ... ")
  ar<-array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2))
  ma<-array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2))
  arma<-ARMA(A=ar,B=ma,C=NULL)
#  print(arma)
  ok <- is.TSmodel(arma) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 4 ... ")
  data.arma.sim<-simulate(arma)
  arma<-l(arma,data.arma.sim)
#  summary(arma)
  if (graphics) 
     {tfplot(data.arma.sim)
      tfplot(arma)
     }
  ok <- is.TSdata(eg1.DSE.data) & is.TSestModel(arma) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 5 ... ")
  f<-array(c(.5,.3,.2,.4),c(2,2))
  h<-array(c(1,0,0,1),c(2,2))
  k<-array(c(.5,.3,.2,.4),c(2,2))
  ss<-SS(F=f,G=NULL,H=h,K=k)
#  ss
  ok <- is.SS(ss)
  all.ok <- all.ok & ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 6 ... ")
  data.ss.sim<-simulate(ss)
  ss<-l(ss,data.ss.sim)
#  summary(ss)
  if (graphics) tfplot(ss)
  ok <- is.TSestModel(ss)
  all.ok <- all.ok & ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 7 ... ")
  ss.from.arma <- l(toSS(arma), data.arma.sim)
  arma.from.ss <- l(toARMA(ss), data.ss.sim)
#  summary(ss.from.arma)
#  summary(arma)
#  summary(arma.from.ss)
#  summary(ss)
#  stability(arma)
#  stability(ss.from.arma)
#  caution: tests on $estimates will depend on seed when data is generated.
  error <- max(Mod(c(-0.15000000000000018874 - sum(TSmodel(ss.from.arma)$F),
                         0.47999999999999998224  - sum(TSmodel(arma.from.ss)$A),
                         -1                      - sum(roots(ss.from.arma)) )))
#                        -1+0i                   - sum(roots(ss.from.arma)) )))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (synopsis) 
    {if (verbose) cat("All Brief User Guide example part 1 tests completed")
     if (all.ok) cat(" OK\n")
     else    cat(", some FAILED! max.error = ", max.error,"\n")
    }
  if (all.ok) invisible(TRUE)  else stop("FAILED")
}




   random.number.test() 
   guide.example.tests.part1(verbose=TRUE, graphics=TRUE)
#         gives  Warning cov. matrix is singular. Working on subspace

unlink(dir, recursive=TRUE)
