# use as
# R -d valgrind --vanilla < valgrind.R &> res_err.txt
library(rmumps)
library(Matrix)
n=3
a=as(diag(1:n), "dgTMatrix")
am=Rmumps$new(a)
am$inv()
ai=Rmumps$new(a@i, a@j, a@x, nrow(a))
ai$inv()
