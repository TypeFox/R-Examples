notests <- FALSE
if(notests) q(save = "no")
stopifnot(require(FAiR))

## FAiR_PACE_by_RGENOUD
Sigma <- diag(5)
Sigma[1,] <- Sigma[,1] <- c(1,.49, .14, .48, .22)
Sigma[2,] <- Sigma[,2] <- c(.49, 1, .07, .42, .14)
Sigma[3,] <- Sigma[,3] <- c(.14, .07, 1, .48, .58)
Sigma[4,] <- Sigma[,4] <- c(.48, .42, .48, 1, .60)
Sigma[5,] <- Sigma[,5] <- c(.22, .14, .58, .60, 1)

rownames(Sigma) <- colnames(Sigma) <- paste("Y", 1:5, sep = ".")
test_PASS <- FAiR:::FAiR_PACE_by_RGENOUD(Sigma, 2)
stopifnot(isTRUE(all.equal(c(0.5, 0.49, 0.50, 0.72, 0.68), test_PASS, 
		check.attributes = FALSE)))

## FAiR_Landahl
test_PASS <- crossprod(FAiR:::FAiR_Landahl(3))
stopifnot(isTRUE(all.equal(diag(3), test_PASS, check.attributes = FALSE)))

## FAiR_make_Tmat
I <- diag(3)
test_PASS <- crossprod(FAiR:::FAiR_make_Tmat(c(I[-3,], I[3,])))
stopifnot(isTRUE(all.equal(diag(3), test_PASS, check.attributes = FALSE)))

## FAiR_Tmat2par
I <- diag(3)
test_PASS <- FAiR:::FAiR_make_Tmat(FAiR:::FAiR_Tmat2par(I))
stopifnot(isTRUE(all.equal(I, test_PASS, check.attributes = FALSE)))

## FAiR_check_Reiersol
beta <- kronecker(c(1,1,1), diag(3))
test_PASS <- FAiR:::FAiR_check_Reiersol(beta)
stopifnot(isTRUE(all.equal(-1, test_PASS, check.attributes = FALSE)))

test_FAIL <- FAiR:::FAiR_check_Reiersol(beta[-c(1,4),])
stopifnot(!isTRUE(all.equal(-1, test_FAIL, check.attributes = FALSE)))

test_FAIL <- FAiR:::FAiR_check_Reiersol(beta[,-1])
stopifnot(!isTRUE(all.equal(-1, test_FAIL, check.attributes = FALSE)))

test_FAIL <- FAiR:::FAiR_check_Reiersol(beta[-c(1,2,7,8),])
stopifnot(!isTRUE(all.equal(-1, test_FAIL, check.attributes = FALSE)))

## FAiR_check_Howe
beta <- kronecker(c(1,1,1), diag(3))
test_PASS <- FAiR:::FAiR_check_Howe(beta)
stopifnot(isTRUE(all.equal(-1, test_PASS, check.attributes = FALSE)))

test_FAIL <- FAiR:::FAiR_check_Howe(beta[-c(1,2,7,8),])
stopifnot(!isTRUE(all.equal(-1, test_FAIL, check.attributes = FALSE)))

