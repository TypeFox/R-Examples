cover<-function(x,sobj){				

# Count the number of distinct events

k<-sum(sobj$n.event>0)

# Initialize the matrix of a's & b's

abmat<-matrix(0,k+1,2)

# Take x as the current critical value

cv<-x

abmat<-exact(sobj,cv)

cp<-noe(k,abmat[k:1,1],abmat[(k+1):2,2])

cp}
