library(datamart)

test_strhead <- function() {
  print(strhead(""))
  print(strhead("", -10))
  print(strhead("abc"))
  print(strhead("abc",-1))
  print(strhead("abc",-2))
  print(strhead("abc",-3))
  print(strhead("abc",-4))
  print(strhead("abc",1))
  print(strhead("abc",2))
  print(strhead("abc",3))
  print(strhead("abc",4))
}

test_strhead()

