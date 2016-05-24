"atom2xyz" <-
function(num) {
  num3 <- num*3
  c(t(matrix(c(((num3) - 2),
               ((num3) - 1),
               (num3)), ncol=3)))
}

