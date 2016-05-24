library(gdata)

options(humanReadable=FALSE)

set.seed(123456)

baseSI <- 10
powerSI <- seq(from=0, to=27, by=3)
SI0 <- (baseSI)^powerSI
k <- length(SI0) - 1
SI1 <- SI0 - SI0 / c(2, runif(n=k, min=1.01, max=5.99))
SI2 <- SI0 + SI0 / c(2, runif(n=k, min=1.01, max=5.99))

baseIEC <- 2
powerIEC <- seq(from=0, to=90, by=10)
IEC0 <- (baseIEC)^powerIEC
IEC1 <- IEC0 - IEC0 / c(2, runif(n=k, min=1.01, max=5.99))
IEC2 <- IEC0 + IEC0 / c(2, runif(n=k, min=1.01, max=5.99))

# Auto units, specify width
cbind(humanReadable(x=SI2,  standard="SI",   width=7),
      humanReadable(x=SI2,  standard="SI",   width=5),
      humanReadable(x=SI2,  standard="SI",   width=3),
      humanReadable(x=IEC2, standard="IEC",  width=7),
      humanReadable(x=IEC2, standard="IEC",  width=5),
      humanReadable(x=IEC2, standard="IEC",  width=3),
      humanReadable(x=IEC2, standard="Unix", width=7),
      humanReadable(x=IEC2, standard="Unix", width=5),
      humanReadable(x=IEC2, standard="Unix", width=3))

# Auto units, specify digits
cbind(humanReadable(x=SI2,  standard="SI",   width=NULL, digits=7),
      humanReadable(x=SI2,  standard="SI",   width=NULL, digits=3),
      humanReadable(x=SI2,  standard="SI",   width=NULL, digits=2),
      humanReadable(x=SI2,  standard="SI",   width=NULL, digits=1),
      humanReadable(x=IEC2, standard="IEC",  width=NULL, digits=7),
      humanReadable(x=IEC2, standard="IEC",  width=NULL, digits=3),
      humanReadable(x=IEC2, standard="IEC",  width=NULL, digits=2),
      humanReadable(x=IEC2, standard="IEC",  width=NULL, digits=1),
      humanReadable(x=IEC2, standard="Unix", width=NULL, digits=7),
      humanReadable(x=IEC2, standard="Unix", width=NULL, digits=3),
      humanReadable(x=IEC2, standard="Unix", width=NULL, digits=2),
      humanReadable(x=IEC2, standard="Unix", width=NULL, digits=1))

# Single unit, specify width
cbind(humanReadable(x=SI1,  units="GB",  standard="SI",   width=7),
      humanReadable(x=SI1,  units="GB",  standard="SI",   width=5),
      humanReadable(x=SI1,  units="GB",  standard="SI",   width=3),
      humanReadable(x=IEC1, units="GiB", standard="IEC",  width=7),
      humanReadable(x=IEC1, units="GiB", standard="IEC",  width=5),
      humanReadable(x=IEC1, units="GiB", standard="IEC",  width=3),
      humanReadable(x=IEC1, units="G",   standard="Unix", width=7),
      humanReadable(x=IEC1, units="G",   standard="Unix", width=5),
      humanReadable(x=IEC1, units="G",   standard="Unix", width=3)
      )

# Single unit, specify digits
cbind(humanReadable(x=SI1, units="GB", standard="SI", width=NULL, digits=7),
      humanReadable(x=SI1, units="GB", standard="SI", width=NULL, digits=3),
      humanReadable(x=SI1, units="GB", standard="SI", width=NULL, digits=2),
      humanReadable(x=SI1, units="GB", standard="SI", width=NULL, digits=1),
      humanReadable(x=IEC1, units="GiB", standard="IEC", width=NULL, digits=7),
      humanReadable(x=IEC1, units="GiB", standard="IEC", width=NULL, digits=3),
      humanReadable(x=IEC1, units="GiB", standard="IEC", width=NULL, digits=2),
      humanReadable(x=IEC1, units="GiB", standard="IEC", width=NULL, digits=1),
      humanReadable(x=IEC1, units="G", standard="Unix", width=NULL, digits=7),
      humanReadable(x=IEC1, units="G", standard="Unix", width=NULL, digits=3),
      humanReadable(x=IEC1, units="G", standard="Unix", width=NULL, digits=2),
      humanReadable(x=IEC1, units="G", standard="Unix", width=NULL, digits=1)
      )


stopifnot( is.object_sizes(as.object_sizes( 2^(1:30) ) ) )
stopifnot( format(as.object_sizes(124)) == "124 bytes")
stopifnot( format(as.object_sizes(124e8), units="auto") == "11.5 GiB")
stopifnot( format(as.object_sizes(124e8), humanReadable=TRUE) == "11.5 GiB")
stopifnot( format(as.object_sizes(124e8), units="bytes") == "1.24e+10 bytes")

tools::assertError( as.object_sizes(-1) )
tools::assertError( as.object_sizes('a') )
tools::assertError( as.object_sizes(list()) )
tools::assertError( as.object_sizes(NULL) )
tools::assertError( as.object_sizes(0+1i) )

stopifnot( format(as.object_sizes(1e40)               ) == "1e+40 bytes"     )
stopifnot( format(as.object_sizes(1e40), units="auto" ) == "8.271806e+15 YiB")
stopifnot( format(as.object_sizes(1e40), units="bytes") == "1e+40 bytes"     )
stopifnot( format(as.object_sizes(1e40), humanReadable=TRUE) == "8.271806e+15 YiB")
stopifnot( format(as.object_sizes(1e40), humanReadable=FALSE) ==  "1e+40 bytes")

options(humanReadable=TRUE)
stopifnot( format(as.object_sizes(1e40) ) == "8.271806e+15 YiB")
options(humanReadable=FALSE)
