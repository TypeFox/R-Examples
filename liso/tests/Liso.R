library(liso)
set.seed(121)


X = matrix(rnorm(50*40),40)
Y = X[,1]+ X[,2]+X[,3] + 0.5*rnorm(40)

cvobj = cv.liso(X,Y)
fitobj = liso.backfit(X,Y, lambda = cvobj$optimlam)
plot(Y, X * fitobj)

