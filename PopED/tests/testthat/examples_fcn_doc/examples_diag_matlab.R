
diag_matlab(3)
diag_matlab(c(1,2,3))
diag_matlab(cbind(1,2,3))
diag_matlab(rbind(1,2,3))

diag_matlab(matrix(c(1, 2, 3),6,6))

# here is where the R default does something different
diag(cbind(1,2,3))
diag(rbind(1,2,3))

