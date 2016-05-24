library(kappalab)

## some alternatives
a <- c(18,11,18,11,11)
b <- c(18,18,11,11,11)
c <- c(11,11,18,18,11)
d <- c(18,11,11,11,18)
e <- c(11,11,18,11,18)
    
## preference threshold relative
## to the preorder of the alternatives
delta.C <- 1

## corresponing Choquet preorder constraint matrix 
Acp <- rbind(c(d,a,delta.C),
             c(a,e,delta.C),
             c(e,b,delta.C),
             c(b,c,delta.C)
            )

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,-delta.I),
             c(3,4,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
delta.I <- 0.01
Aii <- rbind(c(1,2,-0.2,-0.15))


## a minimum variance 2-additive solution
min.var <- mini.var.capa.ident(5,2,A.Choquet.preorder = Acp)              
m <- min.var$solution
m

## the resulting global evaluations
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))

## the Shapley value
Shapley.value(m)

## a minimum variance 3-additive more constrainted solution
min.var2 <- mini.var.capa.ident(5,3,
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp)
m <- min.var2$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
Shapley.value(m)

## a minimum variance 5-additive more constrainted solution
min.var3 <- mini.var.capa.ident(5,5,
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp,
                                   A.Shapley.interval = Asi,
                                   A.interaction.preorder = Aip,
                                   A.interaction.interval = Aii)

m <- min.var3$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
summary(m)



