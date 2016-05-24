##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################


# = =============================================================== =
# =  Massive unit test to check all possible combinations of input  =
# =  parameters and make sure that the output matches the output    =
# =  from the serial version.                                       =
# = =============================================================== =

ratio <- function(d, w) sum(d$x * w)/sum(d$u * w)
trimmedmean <- function(x, d, trim=0) { return(mean(x[d], trim/length(x))) }

test.bootExample1 <- function() {

  # example from http://www.ats.ucla.edu/stat/r/library/bootstrap.htm 
  #defining the ratio function
  #using the boot function
  set.seed(1337)
  a = boot(city, ratio, R=999, stype="w")
  set.seed(1337)
  b = pboot(city, ratio, R=999, stype="w")
  # Ignore the calls having different names when testing equality.
  b$call <- a$call 
  checkEquals(a,b,"Bootstrap examples 1, weight based stype = w")

}

test.bootTrim <- function() {
 # http://www.mayin.org/ajayshah/KB/R/documents/boot.html
   set.seed(1337)
   a = boot(discoveries, trimmedmean, R=1000, trim=5)
   set.seed(1337)
	b = pboot(discoveries, trimmedmean, R=1000, trim=5)
	# Ignore the calls having different names when testing equality.
	b$call <- a$call 
   checkEquals(a,b,"Bootstrap Trim example")
}



# Stratified resampling for the difference of means.  In this
# example we will look at the difference of means between the final
# two series in the gravity data.
diff.means <- function(d, f)
{    n <- nrow(d)
     gp1 <- 1:table(as.numeric(d$series))[1]
     m1 <- sum(d[gp1,1] * f[gp1])/sum(f[gp1])
     m2 <- sum(d[-gp1,1] * f[-gp1])/sum(f[-gp1])
     ss1 <- sum(d[gp1,1]^2 * f[gp1]) - 
            (m1 *  m1 * sum(f[gp1]))
     ss2 <- sum(d[-gp1,1]^2 * f[-gp1]) - 
            (m2 *  m2 * sum(f[-gp1]))
     c(m1-m2, (ss1+ss2)/(sum(f)-2))
}
grav1 <- gravity[as.numeric(gravity[,2])>=7,]

test.grav1 <- function(){
  set.seed(37)
  a = boot(grav1, diff.means, R=999, stype="f", strata=grav1[,2])
  set.seed(37)
	b = pboot(grav1, diff.means, R=999, stype="f", strata=grav1[,2])
	# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"Bootstrap grav1 example")
}  



#  Finally a parametric bootstrap.  For this example we shall look 
#  at the air-conditioning data.  In this example our aim is to test 
#  the hypothesis that the true value of the index is 1 (i.e. that 
#  the data come from an exponential distribution) against the 
#  alternative that the data come from a gamma distribution with
#  index not equal to 1.
air.fun <- function(data)
{    ybar <- mean(data$hours)
     para <- c(log(ybar),mean(log(data$hours)))
     ll <- function(k) {
          if (k <= 0) out <- 1e200 # not NA
          else out <- lgamma(k)-k*(log(k)-1-para[1]+para[2])
         out
     }
     khat <- nlm(ll,ybar^2/var(data$hours))$estimate
     c(ybar, khat)
}


air.rg <- function(data, mle)
#  Function to generate random exponential variates.  mle will contain 
#  the mean of the original data
{    out <- data
     out$hours <- rexp(nrow(out), 1/mle)
     out
}

	# TODO. Don't expect the results from each run to be exactly equal because the random data
	# sampling is done in parallel in this case. Need to re-write this test along the lines of
	# runnit_simple.R "Test simple equals true".
#test.airparam <- function(){
#  set.seed(7)
#  a = boot(aircondit, air.fun, R=999, sim="parametric", ran.gen=air.rg, mle=mean(aircondit$hours))
#  set.seed(7)
#  b = pboot(aircondit, air.fun, R=999, sim="parametric", ran.gen=air.rg, mle=mean(aircondit$hours))
#	# Ignore the calls having different names when testing equality.
#	b$call <- a$call 
#	checkEquals(a,b,"Bootstrap parametric example")
#}


