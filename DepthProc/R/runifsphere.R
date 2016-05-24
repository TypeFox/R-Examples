#'@title Random number generation from unit sphere.
#'@export
#'  @description This function generates random numbers from p-dimensional unit sphere.
#'
#'  @param n number of random samples.
#'  @param p dimension of the unit sphere.
#'
#'  @author Daniel Kosiorowski, Mateusz Bocian, Anna Wegrzynkiewicz and Zygmunt Zawadzki from Cracow University of Economics.
#'
#'  @examples
#'  
#'  x = runifsphere(n=100)
#'  plot(x)
#'  


runifsphere = function(n, p = 2)
{
  runifsphereCPP(n,p)
}
