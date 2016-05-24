kNN_MLE <-
function(X, Y, kmax=ceiling(length(Y)*0.5), plot=TRUE){
d<-numeric(kmax)
Q <- length(unique(Y))
if (Q<2)
    stop("error - fewer than 2 levels in Y")
y <- factor(Y)
for (k in 1:kmax){
    z <- nnc(X = X, Y = y, k=k)
    d[k] <- NA
    if (Q == 2)
        d[k]<- deviance(stats::glm.fit(x = z, y = y, family = binomial(link = "logit"))) else
        	capture.output(d[k] <- nnet::multinom(y~z)$deviance)
}
kHat <- which.min(d)
etaOpt <-misclassificationrate(y, class::knn.cv(X, y, k = kHat))
etaOpt <- paste0(round(100*etaOpt, 2),"%")
#plot and/or return
if (plot) {
    plot(d, xlab="k", ylab="deviance")
    points(kHat, d[kHat], col="red", pch=16, cex=0.5)
    title(sub=bquote(hat(k) == .(kHat)))
    title(main=paste("MLE misclassification rate =", etaOpt))
    }
 kHat
}

nnc <-
function (X, Y, k) {
    n <- length(Y)
    if (nrow(X) != n)
        stop("nrow(X) != length(Y)")
    if (k <= 0) stop("error: k <= 0")
    if (k >= n) stop("error: k >= n")
    classes <- unique(Y)
    Q <- length(classes)
    if (Q < 2)
        stop("invalid, need Q >= 2, Q = number of classes")
    y <- numeric(n)
    if (Q==2){
        if (!all(y %in% c(-1, 1))){
            ind <- Y==classes[1]
            y[ind]<- -1
            y[!ind]<- 1
            }
        ans <- class::knn.cv(train = X, cl = as.factor(y), k = k, prob = TRUE)
        pr<- attr(ans, "prob") #proportion of votes for winning class
        yfit<-as.numeric(as.character(ans))
        numVoteWinner <- pr * k
        numVoteLoser <- k - numVoteWinner
        ConsensusProportion <- (numVoteWinner - numVoteLoser)/k
        z <- yfit * ConsensusProportion
        }
    else {
        ind1 <- Y==classes[1]
        y[ind1] <- -1
        z <- matrix(numeric(n*(Q-1)), nrow=n)
        for (j in 2:Q){
            indk <- Y==classes[j]
            indOther <- !(ind1|indk)
            y[indk] <- 1
            y[indOther] <- -1
            zA <- nnc(X, y, k)
            y[indOther] <- 1
            zB <- nnc(X, y, k)
            z[,j-1] <- (zA+zB)/2
            }
    }
    z
}

`kNN_LOOCV` <-
  function(X, y, kmax=ceiling(length(y)*0.5), plot=FALSE){
    stopifnot (nrow(X)==length(y), is.factor(y))
    ks<-seq(1,kmax,2)
    nks<-length(ks)
    m<-numeric(nks)
    for (i in 1:nks)
      m[i]<-misclassificationrate(y, class::knn.cv(X, y, k = ks[i]))
    names(m)<-ks
    indOpt <- which.min(m)
    etaOpt <- paste0(round(100*m[indOpt],2),"%")
    kopt<-as.numeric(names(sort(m)[1]))
    if (plot) {
      plot(ks,m, xlab="k", ylab="error rate")
      points(kopt, m[indOpt], col="red", pch=16, cex=0.5)
      title(main=paste("LOOCV misclassification rate =", etaOpt))
      title(sub=bquote(hat(k) == .(kopt)))
    }
    ks[which.min(m)]
   }








