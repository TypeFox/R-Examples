#### Testing the  Exponential Condition Number computations

library(expm)

mSource <- function(file, ...) source(system.file(file, ...,
						  package = "expm", mustWork=TRUE),
				      keep.source = FALSE)
mSource("test-tools.R")## -> assertError(), rMat()
## getting examples where we know  expm(.) "exactly":
mSource("demo", "exact-fn.R")

M <- xct10$m
eC <- list(expmCondF = 566.582631819923,
           expmCond1 = 137.455837652872)
C1 <- expmCond(M, "exact")
(C2 <- expmCond(M, "1.est", expm=FALSE))
(C3.  <- expmCond(M, "F.est", abstol = 0.1)[[1]])
(C3.1 <- expmCond(M, "F.est", abstol = 0.01, reltol = 1e-12)[[1]])

stopifnot(all.equal(C1[1:2], eC, tolerance = 1e-14),
	  all.equal(C2	, eC$expmCond1, tolerance = 1e-14),
	  all.equal(C3. , eC$expmCondF, tolerance = 1e-14, check.attributes = FALSE),
	  all.equal(C3.1, eC$expmCondF, tolerance = 1e-14, check.attributes = FALSE))

cat('Time elapsed: ', (p1 <- proc.time()),'\n') # for ``statistical reasons''

## cat('Time elapsed: ',(p2 <- proc.time())-p1,'\n') # for ``statistical reasons''

## cat('Time elapsed: ',(p3 <- proc.time())-p2,'\n') # for ``statistical reasons''
