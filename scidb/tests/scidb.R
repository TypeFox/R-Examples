# scidb array object tests
# Set the SCIDB_TEST_HOST system environment variable to the hostname or IP
# address of SciDB to run these tests. If SCIDB_TEST_HOST is not set the tests
# are not run.

check = function(a,b)
{
  print(match.call())
  stopifnot(all.equal(a,b,check.attributes=FALSE))
}

library("scidb")
host = Sys.getenv("SCIDB_TEST_HOST")
if(nchar(host)>0)
{
  scidbconnect(host)
  options(scidb.debug=TRUE)
  oplist = scidbls(type="operators")
  got = function(op)
  {
    length(grep(op,oplist))>0
  }

cat("# Upload dense matrix to SciDB\n")
  set.seed(1)
  A = matrix(rnorm(50*40),50)
  B = matrix(rnorm(40*40),40)
  X = as.scidb(A, chunkSize=c(30,40), start=c(1,0))
  Y = as.scidb(B, chunkSize=c(17,20), star=c(2,3))
  check(X[],A)
  check(Y[],B)

cat("# n-d array\n")
  x = build("k",dim=c(3,3,3),names=c("x","i","j","k"), type="double")
  check(x[][,,2] , matrix(1,3,3))

cat("# Dense matrix multiplication\n")
  check(A%*%B, (X%*%Y)[])
cat("# Transpose\n")
  check(crossprod(A), (t(X)%*%X)[])
cat("# Crossprod, tcrossprod\n")
  check(crossprod(A), crossprod(X)[])
  check(tcrossprod(B), tcrossprod(B)[])
cat("# Arithmetic on mixed R/SciDB objects\n")
  x = rnorm(40);
  check((X%*%x)[,drop=FALSE], A%*%x)
cat("# Scalar multiplication\n")
  check(2*A, (2*X)[])
cat("# Elementwise addition\n")
  check(A+A, (X+X)[])
cat("# SVD\n")
  check(svd(A)$d, as.vector(svd(X)$d[]))

cat("# Numeric array subsetting\n")
  check((X %*% X[1,,drop=TRUE])[,drop=FALSE], A %*% A[1,])
  check(X[c(2,6,16),c(11,12,25)][], A[c(2,6,16),c(12,13,26)])
  check(as.vector(diag(Y)[]), diag(B))

cat("# Filtering\n")
  W = subset(X,"val>0")
cat("# Sparse elementwise addition and scalar multiplication\n")
  D = (W + 2*W)[]
  w = W[]
  w = w + 2*w
  check(sum(D-w), 0)

cat("# Indexing by other SciDB arrays\n")
  x = build("i",c(5,3),type="double")
  a = as.scidb(c(1,5,1,5,1))
  check(x[a %>% 2, ][], matrix(c(1,1,1,3,3,3),nrow=2,byrow=TRUE))

cat("# some binary operations\n")
# XXX ADD **lots** more tests here
  a = as.scidb(Matrix::sparseMatrix(
               sample(10,100,replace=TRUE),sample(10,100,replace=TRUE),x=runif(100)))
  b = build(pi,c(10,10),names=c("a","b","c"),chunksize=c(3,2))
  check(sum((a*b)[] - a[]*b[]), 0)
  check(sum((a+b)[] - (a[]+b[])), 0)
  check(sum((a+2)[] - (a[]+2)), 0)
  a*apply(a,2,mean)
  apply(a,2,mean)*a

cat("# Make sure binary operations handle mismatched coordinates\n")
  x[1:2, 1] - x[2:3,1]

cat("# Aggregation\n")
  check( sweep(B,MARGIN=2,apply(B,2,mean)),
         sweep(Y,MARGIN=2,apply(Y,2,mean))[])


cat("# Join\n")
# We need 'subarray' here to reset the origin of Y to zero to match
# diag(Y)--diag always returns a zero indexed vector.
  check(project(bind(merge(subarray(Y),diag(Y),by.x="i",by.y="i_1"),"v","val*val_1"),"v")[], diag(B)*B)

cat("# On different dimensions\n")
  x = as.scidb(rnorm(5))
  a = as.scidb(data.frame(p=1:5),start=0)
  merge(x,a,by.x=dimensions(x),by.y=dimensions(a))

# Add many more join/merge checks here... XXX

cat("# Sparse upload, count\n")
  S = Matrix::sparseMatrix(sample(100,200,replace=TRUE),sample(100,200,replace=TRUE),x=runif(200))
  Z = as.scidb(S,start=c(1,5))
  check(count(Z), Matrix::nnzero(S))

cat("# Check that image does not fail\n")
  image(Z,plot=FALSE)

cat("# spgemm\n")
  if(got("spgemm"))
  {
    x = crossprod(Z)
  }

cat("# Misc\n")
  z = atan(tan(abs(acos(cos(asin(sin(Z)))))))
  s = atan(tan(abs(acos(cos(asin(sin(S)))))))
  check(sum(s-z[]), 0)

cat("# Labeled indices\n")
  L = c(letters,LETTERS)
  i = as.scidb(data.frame(L[1:nrow(X)]), start=1)
  j = as.scidb(data.frame(L[1:ncol(X)]), start=0)
  rownames(X) = i
  colnames(X) = j
  rownames(A) = L[1:nrow(A)]
  colnames(A) = L[1:ncol(A)]
  check(X[c("f","v","F"),c("a","A","N")][], A[c("f","v","F"),c("a","A","N")])

cat("# 4d labels and auto promotion of labels to SciDB arrays\n")
  X = build(0,dim=c(3,4,5,6))
  rownames(X) = letters[1:3]
  dimnames(X)[[2]] = letters[1:4]
  dimnames(X)[[3]] = letters[1:5]
  dimnames(X)[[4]] = letters[1:6]
  i = count(X[c("a","b"),"a",c("d","e"),"a"])
  check(i,4)

cat("# Indices not at the origin, and general ranged index check\n")
  a = scidbeval(build("random()",c(5,5),start=c(-3,-2)))
  check(count(a[-3:-2,0:2]),6)

cat("# Pseudo-uint64 support! And also simplified aggregation function syntax\n")
# for apply.
  check(sum(apply(a,2,count)),25)

cat("# Aggregation, just trying to catch errors, github issue #57\n")
  A = scidbeval(build("random()%10",c(100,100)))
  p = scidbeval(build("random()%2",100))
  aggregate(A,by=p,mean) # Aggregate by another array
  aggregate(A,by=2,mean) # Positional dimension index
  aggregate(A,FUN=mean)  # Grand aggregate
  x = build(5.5,10)
  aggregate(x, FUN=sum, window=c(0,1e6)) # GitHub issue #57
  aggregate(x,FUN=sum, by=1, variable_window=c(0,1e6))

cat("# More special index tests, sort, apply, densify giant sparse sort result.\n")
  a = build("i+j",c(5,5),type="double")
  s = sort(a)
  check(as.vector(s[s][]), sort(a[]))

cat("# GLM (requires SciDB p4 plugin)\n")
  if(got("glm"))
  {
    x = as.scidb(matrix(rnorm(5000*20),nrow=5000))
    y = as.scidb(rnorm(5000))
    M = glm.fit(x, y)
  }

cat("# slice\n")
  x = build("i+j+k",c(3,4,5),type="double")
  check(slice(x,c("i","j"),c(2,3))[0][], 5)

# na.locf
# Write me!

# hist
# Write me!

cat("# order\n")
  x = scidbeval(build("random()%100", 100, type="double", start=1))
  o1 = order(x)[]
  o2 = order(x[])
  xx = x[]
  check(xx[o1], xx[o2])

cat("# rank\n")
  check(rank(x)[,2][], rank(x[])) 


cat("# Another quantile failure case reported by Alex (bug #47)\n")
x = read.csv(file=textConnection('"tumor_type_id","sample_id","agilentg4502a_07_1_probe_id","value"
7,2742,13317,1.1024545
7,2963,8060,0.9827
7,2609,13709,-0.18572727
7,2643,13772,0.56753844
7,2629,2126,3.6668334
7,2643,10996,0.35366666
7,2594,10300,-0.534
7,2680,4252,-0.842
7,2678,17062,-1.1738
7,2765,13244,-2.1102'))
a = redimension(as.scidb(x,types=c("int64","int64","int64","double")), dim=names(x)[1:3])
check(quantile(a)[][,2], quantile(x$value))

cat("# Check that replaceNA works\n")
 replaceNA(scidb("apply(build(<x:double null>[i=1:1,1,0],null), s, string(null))"))
 check(count(replaceNA(scidb("apply(build(<x:double null>[i=1:1,1,0],null), s, string(null))"))),1)

cat("# Another merge test courtesy Alex Polyiakov\n")
  x = build(1,c(2,2))
  a = build(2,5,names=c("a","j"))
  z = merge(x,a,by="j")
  check(count(z),4)

cat("# Complicated cross_join filtering\n")
  set.seed(1)
  X = as.scidb( matrix(rnorm(20),nrow=5))
  rownames(X) = as.scidb( data.frame(letters[1:5]))
  X[c("b","a","d"), ]
gc()
  idx = rownames(X) > "b"
  check(nrow(X[idx, ]),3)

cat("# Github issue #45\n")
  x = scidbeval(build(5.5,5))
  a = scidbeval(x, eval=FALSE)
  rm(a)
  gc()
  x[]

cat("# Github issue #54\n")
  scidbremove("oh_no",error=invisible,force=TRUE)
  iquery("create array oh_no <mask:bool>[variant_id=0:1109999,10000,0,pos=0:12000000,10000000,0,chrom_id=0:0,1,0]")
  x = scidb("oh_no")
  x [, 11000000:12000000, ]
  scidbremove("oh_no",error=invisible,force=TRUE)

cat("# Github issue #52\n")
  x = scidb("build(<s:string>[i=1:3,1,0,j=1:1,1,0],'{1,1}[[(law)],[(bryan)],[(homer)]]',true)")
  a = as.scidb(x[])

cat("# Github issue #61\n")
  v = scidb("build(<val:double>[i=1:10000,1000,0], random()%10)")
  v = scidbeval(v, temp=TRUE)
  aggregate(v, FUN=count, by="val")[]

cat("# Github issue #62\n")
  a = as.scidb(matrix(rnorm(25),5))
  i = as.scidb(1:5)
  a[ i %>% 2,  between(1,3) ]

}
rm(list=ls())
gc()
