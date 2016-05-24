
kmeansN <- function(x, k, variable.names="x",target.name="y", 
                    cluster.name="cluster",n = 100) {
  # initialize
  .id <- NULL; .order_key <- NULL; cluster <- NULL; i <- NULL;
  x$.id <- 1:nrow(x)
  km.base <- ldply(foreach(i=1:n, combine=rbind) %do% {
    x2 <- x
    set.seed(i)
    km <- kmeans(x2[,variable.names],k)
    x2[,cluster.name] <- km$cluster
    m <- aggregate(x2[,target.name],
                   by=list(x2[,cluster.name]),
                   FUN=mean, na.rm=T)
    names(m) <- c(cluster.name, ".order_key")
    m <- arrange(m, .order_key)
    m$new.cluster <- 1:k
    x2 <- merge(x2, m, by=cluster.name)
    data.frame(n=i, .id=x2$.id, cluster=x2$new.cluster, y=x2[,target.name])
  })
  
  res <- ddply(km.base, .(.id), summarize,
               cluster=as.numeric(names(sort(table(cluster), decreasing = T)[1])))
  
  names(res)[2] <- cluster.name
  x <- merge(x, res, by=".id")
  x[,!(names(x) %in% ".id")]
}

kmeansN2 <- function(x, variable.names="x",target.name="y",k.list=3:6,
                     cluster.name="cluster",n = 100){
  stopifnot(class(x)=="data.frame")
  
  k <- NULL
  df.all <- ldply(foreach(k=k.list, combine=rbind) %do% {
    kmn <- kmeansN(x, k, variable.names, target.name, cluster.name, n)
    s <- aggregate(kmn[,target.name],
                   by=list(kmn[,cluster.name]),
                   FUN=sd, na.rm=T)
    kmn$.average.sd <- mean(s$x, na.rm=T)
    kmn
  })
  
  df.all
}

ykmeans <- function(x, variable.names="x",target.name="y",k.list=3:6,
                    cluster.name="cluster",n = 100){
  stopifnot(class(x)=="data.frame")
  
  df.all <- kmeansN2(x, variable.names, target.name, k.list, cluster.name, n)
  df <- df.all[!is.na(df.all$.average.sd) & 
                 df.all$.average.sd==min(df.all$.average.sd,na.rm=T), ]
  df
}
