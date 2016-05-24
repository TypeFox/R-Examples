x <- 1:12
matrix(x,nr=2)                     # 2 rows, entries columnwise
matrix(x,nr=3,byrow=TRUE)          # 3 rows, entries rowwise
matrix(x,nc=3,byrow=TRUE)          # 3 columns, entries rowwise
x                                  # x is unchanged
