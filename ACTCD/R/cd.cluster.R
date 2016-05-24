cd.cluster <-
function(Y, Q, method = c("HACA","Kmeans"), Kmeans.centers = NULL, 
                       Kmeans.itermax = 10,Kmeans.nstart = 1,
                       HACA.link = c("complete", "ward", "single","average", 
                      "mcquitty", "median", "centroid"),HACA.cut = NULL) {

  
  cluster.method <- match.arg(method)
  HACA.link <- match.arg(HACA.link)
  
#------------------------Input check-------------------------#
input.check(Y,Q,cluster.method=cluster.method)

#-----------------------Basic variables----------------------#
#N:number of examinees
#J:number of items
#K:number of attributes
#M:number of ideal attribute patterns, which is equal to 2^K
N <- dim(Y)[1]
J <- dim(Y)[2]
K <- dim(Q)[2]
M <- 2^K
#------------------------------------------------------------#

# a vector of summed scores for items that measure each of the K attributes
w <- as.matrix(Y)%*%as.matrix(Q)

#========================KMEANS==============================#
if (cluster.method=="Kmeans") {
  if (is.null(Kmeans.centers)){
    Kmeans.centers <- M
  }
kmns <- kmeans(w, Kmeans.centers, Kmeans.itermax, Kmeans.nstart)
Kmeans.class <- kmns$cluster
Kmeans.size <- kmns$size
Kmeans.mean.w <- kmns$centers
Kmeans.wss.w <- kmns$withinss
Kmeans.class <- kmns$cluster
Kmeans.sqmwss.w <- NULL
Kmeans.mean.y <- NULL

for (i in 1:Kmeans.centers)
{
        c <- which(Kmeans.class==i)
        if (length(c)>1)
        {
                #the mean of total score for each cluster
                mean.y <- mean(apply(Y[c,],1,sum))
                #the average within-cluster sum of squares
                sqmwss.w <- sqrt(mean(apply((w[c,]-matrix(rep(Kmeans.mean.w[i,], length(c)), length(c),K,byrow=TRUE))^2,1,sum)))
        }
        if (length(c)==1)
        {
                mean.y <- sum(Y[c,])
                sqmwss.w <- 0
                wss.w <- 0
        }
        if (length(c)==0)
        {
                mean.y <- NA
                sqmwss.w <- NA
                wss.w <- NA
        }
        Kmeans.mean.y <- c(Kmeans.mean.y, mean.y)
        Kmeans.sqmwss.w <- c(Kmeans.sqmwss.w, sqmwss.w)
}

output <- list(W=w,size=Kmeans.size, mean.w=Kmeans.mean.w,
wss.w=Kmeans.wss.w, sqmwss.w=Kmeans.sqmwss.w,
mean.y=Kmeans.mean.y, class=Kmeans.class,cluster.method=cluster.method)

}else if(cluster.method=="HACA")
  
  #========================HACA==============================#  
  {
hc <- hclust(dist(w), method=HACA.link)
if (is.null(HACA.cut)){
  HACA.cut<-M
}

HACA.class <- cutree(hc, k=HACA.cut)

HACA.size <- NULL
HACA.mean.w <- NULL
HACA.sqmwss.w <- NULL
HACA.mean.y <- NULL
HACA.wss.w <- NULL

for (i in 1:HACA.cut)
{
        c <- which(HACA.class==i)
        if (length(c)>1)
        {
                mean.w <- apply(w[c,],2,mean)  #the average w within each group:a vector with K elements 
                mean.y <- mean(apply(Y[c,],1,sum))  #the average y of each group:a scalar
                sqmwss.w <- sqrt(mean(apply((w[c,]-matrix(rep(mean.w,length(c)),
                length(c),K, byrow=TRUE))^2,1,sum)))
                wss.w <- sum((w[c,]-matrix(rep(mean.w,length(c)),length(c),
                K, byrow=TRUE))^2)
        }
        if (length(c)==1)
        {
                mean.w <- w[c,]
                mean.y <- sum(Y[c,])
                sqmwss.w <- 0
                wss.w <- 0
        }
        if (length(c)==0)
        {
                mean.w <- c(NA,NA,NA)
                mean.y <- NA
                sqmwss.w <- NA
                wss.w <- NA
        }
        HACA.size <- c(HACA.size, length(c))
        HACA.mean.w <- rbind(HACA.mean.w, mean.w)
        HACA.mean.y <- c(HACA.mean.y, mean.y)
        HACA.sqmwss.w <- c(HACA.sqmwss.w, sqmwss.w)
        HACA.wss.w <- c(HACA.wss.w, wss.w)
}

output <- list(W=w,size=HACA.size, mean.w=HACA.mean.w,
wss.w=HACA.wss.w, sqmwss.w=HACA.sqmwss.w,
mean.y=HACA.mean.y,class=HACA.class,cluster.method=cluster.method)


}
class(output) <- "cd.cluster"
return(output)
}
