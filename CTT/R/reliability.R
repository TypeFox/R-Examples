`reliability` <-
function(items, itemal=TRUE, NA.Delete=TRUE){

if(!all(apply(items,c(1,2),is.numeric))) { items <- apply(items,c(1,2),as.numeric)
          warning("Data is not numeric. Data has been coerced to be numeric.")}

if(NA.Delete==FALSE) { items[is.na(items)] <- 0
                       warning("Missing values or NA values are converted to zeros.")} 
     

items <- na.omit(items)
s <- apply(items,2,var)
N <- ncol(items)

X <- rowSums(items)
alpha <- (N/(N-1))*(1 - sum(s)/var(X))

if(itemal){
  alphad <- array(dim=N)
  pbis <- array(dim=N)
  bis <- array(dim=N)
  for(i in 1:N){
    Xd <- rowSums(items[,-i])
    pvalu <- colMeans(items)
    alphad[i] <- ((N-1)/(N-2))*(1 - sum(s[-i])/var(Xd))
    pbis[i] <- cor(items[,i],Xd)
	bis[i] <- polyserial(Xd, items[,i])
    out <- list(nItem=N,nPerson=nrow(items),alpha=alpha, scaleMean=mean(X), scaleSD=sd(X),
                alphaIfDeleted=alphad, pBis=pbis, bis=bis, itemMean=pvalu)
  }
} 
else out <- list(nItem=N,nPerson=nrow(items),alpha=alpha, scaleMean=mean(X), scaleSD=sd(X))
class(out) <- "reliability"
out
}
