block_solve <-
function(Xscale,Yscale,samplesize,ii,sortcol,sortcol2)
{
n <- samplesize
TXT <- NULL
XYmain <- matrix(0,length(ii),1)
TXTmain <- matrix(0,length(ii),length(ii))
for (i in 1:max(sortcol))
{
if(sum(sortcol[ii] == i) > 0) 
{
XY <- NULL
XT <- NULL
XT <- Xscale[,sortcol2[ii][which(sortcol[ii]==i)]]
XY <- t(XT)%*%Yscale[,i]
ind <- c(1:length(ii))[which(sortcol[ii] == i)]
XYmain[ind,] <- XY
if(qr(t(XT)%*%XT)$rank == ncol(t(XT)%*%XT))
 {
  TXTmain[ind,ind] <- solve(t(XT)%*%XT)
 }
else
 {
  TXTmain[ind,ind] <- MPinverse(t(XT)%*%XT)
 }
}
}
result <- TXTmain%*%XYmain
return(result)
}
