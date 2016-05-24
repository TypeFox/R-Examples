knngow <-
function (train,test,k) 
{#require(StatMatch)
 p=dim(train)[2]
 ntest=dim(test)[1]
 ntrain=dim(train)[1]
 classes=rep(0,ntest)
 if(ntest==ntrain)
{     for(i in 1:ntest)
     { tempo=order(gower.dist(test[i,-p],train[-i,-p]))[1:k]
       classes[i]=moda(train[tempo,p])[1]}
 }
     else{
  for(i in 1:ntest)
{ tempo=order(StatMatch::gower.dist(test[i,-p],train[,-p]))[1:k] 
classes[i]=moda(train[tempo,p])[1]}
}
classes
}
