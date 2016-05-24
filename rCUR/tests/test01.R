library(rCUR)
runif(1)  #this is needed to make .Random.seed visible
rnd_save=0

load("test01.rda")
cur_sv=svd(m)

#Test for division by zero bug
#long test
#  for(cur_r in 2:400) for(cur_k in c(217,222,227,232,237)) {
#short test
for(cur_r in c(300,427)) for(cur_k in c(217,436)) {
  cur_res <- CUR(m , r=cur_r, k=cur_k, sv=cur_sv, method="ortho.top.scores", alpha=1)
}

#Test for Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
#by default it is not run
#tryCatch(
#{
#  for(cur_r in 2:400) for(cur_k in c(217,222,227,232,237)) {
#    rnd_save = .Random.seed;
#    cur_res <- CUR(m , r=cur_r, k=cur_k, sv=cur_sv, method="exact.num.random")
#  }
#}, error = function(e) e)

#somwhere in the middle it will stop with
#Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'

#to extract the offending matrix, we run CUR again with different parameters:
#.Random.seed = rnd_save
#cur_links <- CUR(m , r=cur_r, k=cur_k, sv=cur_sv, method="exact.num.random", matrix.return=F)

#the matrix is extracted:
#R=m[cur_links@R.index,]

#this sould produce the error deterministically
#try( svd(R), silent = TRUE)
#geterrmessage()
