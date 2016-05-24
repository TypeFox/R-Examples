# y-hat values:
H %*% y
#
# beta-hat values:
#
solve(t(X) %*% X) %*% t(X) %*% y
