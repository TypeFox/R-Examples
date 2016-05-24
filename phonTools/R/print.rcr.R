# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


print.rcr <-
function (x, ...){
  cat ("\nCall:\n")
  print (x$call)
  
  cat ("\nCoefficient Means:\n")
  print (x$coefficient.means)
}
