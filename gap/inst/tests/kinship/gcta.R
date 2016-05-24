# 12-12-2013 MRC-Epid JHZ

p <- matrix(0,N,4)
for(i in 1:N) p[i,] <- with(ped51[i,],c(i,i,qt,bt))
write(t(p),file="51.txt",4,sep="\t")
NN <- rep(51, N * (N + 1)/2)
WriteGRM(51,p[,1:2],NN,k2)
one <- ReadGRM(51)
grm <- one$grm
WriteGRMBin(51,grm,NN,p[,1:2])
two <- ReadGRMBin(51,TRUE)
sum(one$GRM-two$GRM)
