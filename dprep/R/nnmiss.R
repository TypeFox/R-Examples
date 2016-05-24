nnmiss <-
function(x, xmiss, ismiss,xnom, K=1) 
{
#x:submatrix of complete rows from original matrix
#xmiss: a row with a missing value
#ismiss: vector that indicates whether a value in xmiss is missing or not
#xnom: vector with indexes of nominal variables
#require(StatMatch)
#Find distance between xmiss (not NA) and each row of x
#print(xmiss)
#print(x)
dd=StatMatch::gower.dist(xmiss,x)
#print(dd)
#order of the rows of x according to their closeness to xmiss
  od <- order(dd)[seq(K)]
#print(od)
#if column of ismiss is nominal, find mode otherwise find mean of KNN
  
  ismiss.nom=ismiss[]&xnom[]
  ismiss.con=ismiss[]&!xnom[]
#print(xnom)
#print(ismiss)
#print(ismiss.nom)
colnom=seq(1:dim(x)[2])[ismiss.nom]
colcon=seq(1:dim(x)[2])[ismiss.con]
#print(seq(1:dim(x)[2])[ismiss.con])
xmiss[ismiss.nom] <- apply(as.data.frame(x[od,colnom]),2, function(x)moda(x)[1] )
  xmiss[ismiss.con] <- apply(as.matrix(x[od, colcon]),2,mean)
  xmiss
}
