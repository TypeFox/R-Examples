library("planor")
 haies.fnames <- c("A1","A2","B1","B2","C1","C2","P1","P2","Q1","Q2")
 haies.key1 <- planor.designkey(factors=haies.fnames,
 nlevels=rep(2,10), nunits=2^5, model=~A1*A2+B1*B2+C1*C2+P1*P2+Q1*Q2,
 base=~A1+B1+C1+P1+Q1,max.sol=1) 
planor.designkey(factors=haies.fnames, nlevels=rep(2,10), nunits=2^5, model=~A1*A2+B1*B2+C1*C2+P1*P2+Q1*Q2+A1*B1*C1*P1*Q1, base=~A1+B1+C1+P1+Q1,max.sol=1) 
