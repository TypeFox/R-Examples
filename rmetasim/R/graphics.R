#Allan Strand 9/29/01
#
#Plot various characteristics of landscape
#both pretty lame
#
popsizedist.plot.landscape <- function(Rland)
  {
    if (is.landscape(Rland))
      {
        plot(table(landscape.populations(Rland)),
                   main="Frequency distribution of population sizes",
                   xlab=c("Population"),
                   ylab=c("Number of individuals"),
                   )
      }
  }


stgsizedist.plot.landscape <- function(Rland)
  {
    if (is.landscape(Rland))
      {
        plot(table(Rland$individuals[,1]+1),
                   main=c("Frequency distribution of demographic stage sizes"),
                   xlab=c("Demographic stage"),
                   ylab=c("Number of individuals"),
             )
      }
  }


