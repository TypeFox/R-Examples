library(datamart)

test_strxxcrypt <- function() {
  print(strencrypt(""))
  print(strdecrypt(""))
  
  print(strdecrypt(strencrypt("abc")))
  print(strencrypt("abc"))
}

test_strxxcrypt()

