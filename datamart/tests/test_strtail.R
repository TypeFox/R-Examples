library(datamart)

test_strtail <- function() {
  print(strtail(""))
  print(strtail("", -10))
  print(strtail("abc"))
  print(strtail("abc",-1))
  print(strtail("abc",-2))
  print(strtail("abc",-3))
  print(strtail("abc",-4))
  print(strtail("abc",1))
  print(strtail("abc",2))
  print(strtail("abc",3))
  print(strtail("abc",4))
}

test_strtail()

