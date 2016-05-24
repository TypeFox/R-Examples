library( fields)
options(echo=FALSE)
test.for.zero.flag<- 1

n<- 50
x<- matrix( runif(n*2), n, 2)
A<- Exp.cov( x,x,theta=.2)
B<- Exp.cov( x,x, theta=.5)

fields.diagonalize(A,B)-> look
fields.diagonalize2(A,B, verbose=FALSE)-> look2

test.for.zero( look$D, look2$D, tol=1E-8,tag="eigenvalues of both versions")

G1<- look$G
G2<-look2$G
a1<- sign( G1[1,])
a2<- sign(G2[1,])
a<- a1*a2


lambda<- .8
test.for.zero( solve( A + lambda* B), G2%*%diag( 1/(1+ lambda*look2$D))%*%t(G2), tag="inverse A+lambda*B", tol=1e-8 )
test.for.zero( solve( A + lambda* B), G1%*%diag( 1/(1+ lambda*look$D))%*%t(G1), tag="inverse A+lambda*B", tol=1e-8 )
test.for.zero( G2%*%diag( 1/(1+ lambda*look2$D))%*%t(G2) ,
                    G1%*%diag( 1/(1+ lambda*look$D))%*%t(G1), tag="inverse A+lambda*B" , tol=1e-8)

options( echo=TRUE)
cat("all done testing both versions of simultaneous diagonalization ", fill=TRUE)

