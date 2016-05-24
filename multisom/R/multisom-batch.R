multisom.batch <-function(data= NULL,xheight, xwidth,
                          topo= c("rectangular", "hexagonal"),
                          min.radius=0.0001,max.radius=0.002,
                          maxit=1000,init=c("random","sample","linear"),
                          radius.type=c("gaussian","bubble","cutgauss","ep"),index="all")
{

  set.seed(10)

  min_nc <- 2
  max_nc <- xheight
  dataa<- data

  indice <- pmatch(index, c("db","dunn","silhouette","ptbiserial","ch","cindex","ratkowsky","mcclain",
                            "gamma","gplus","tau","ccc","scott","marriot","trcovw","tracew","friedman",
                            "rubin","ball","sdbw","dindex","hubert","sv","xie-beni","hartigan","ssi","xu",
                            "rayturi","pbm","banfeld","all"))


  if (xheight != xwidth)
    stop("dimension shoud be equal")

  if(is.null(data))
  {
    stop("\n","data matrix is needed")
  }
  else{
    data <- as.matrix(data)
    numberObsBefore <- dim(data)[1]
    data <- na.omit(data) # returns the object with incomplete cases removed
    nn <- numberObsAfter <- dim(data)[1]
    pp <- dim(data)[2]
    TT <- t(data)%*%data
    sizeEigenTT <- length(eigen(TT)$value)
    eigenValues <- eigen(TT/(nn-1))$value

    # Only for indices using vv : CCC, Scott, marriot, tracecovw, tracew, friedman, rubin
    if((is.na(match(12,indice))==FALSE) || (is.na(match(13,indice))==FALSE) ||(is.na(match(14,indice))==FALSE)||
       (is.na(match(15,indice))==FALSE) || (is.na(match(16,indice))==FALSE)||(is.na(match(17,indice))==FALSE)||
       (is.na(match(18,indice))==FALSE)||(is.na(match(31,indice))==FALSE))
    {
      for (i in 1:sizeEigenTT)
      {
        if (eigenValues[i] < 0) {
          stop("The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.")
        }
      }
      s1 <- sqrt(eigenValues)
      ss <- rep(1,sizeEigenTT)
      for (i in 1:sizeEigenTT)
      {
        if (s1[i]!=0)
          ss[i]=s1[i]
      }
      vv <- prod(ss)
    }
  }

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#                                                                                                                      #
  #                                              Indices                                                                 #
  #                                                                                                                      #
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

  ###############################################
  #                                             #
  #        BALL, CH, HARTIGAN and RATKOWSKY     #
  #                                             #
  ###############################################

  gss <- function(x, cl, d)
  {
    n <- length(cl)
    k <- max(cl)
    centers <- matrix(nrow = k, ncol = ncol(x))
    for (i in 1:k){

      if (ncol(x) == 1){
        centers[i, ] <- mean(x[cl == i, ])
      }
      if (is.null(dim(x[cl == i, ]))){
        bb <- matrix(x[cl == i, ],byrow=FALSE,nrow=1,ncol=ncol(x))
        centers[i, ] <- apply(bb, 2, mean)
      }else {
        centers[i, ] <- apply(x[cl == i, ], 2, mean)
      }

    }
    allmean <- apply(x, 2, mean)
    dmean <- sweep(x, 2, allmean, "-")
    allmeandist <- sum(dmean^2)
    withins <- rep(0, k)
    x.2 <- (x - centers[cl, ])^2
    for (i in 1:k) {
      withins[i] <- sum(x.2[cl == i, ])
    }
    wgss <- sum(withins)
    bgss <- allmeandist - wgss

    results <- list(wgss=wgss,bgss=bgss, centers=centers,gss=withins)
    return(results)
  }

  vargss <- function(x, clsize, varwithins)
  {
    nvar <- dim(x)[2]
    n <- sum(clsize)
    k <- length(clsize)

    varallmean <- rep(0, nvar)
    varallmeandist <- rep(0, nvar)
    varwgss <- rep(0, nvar)
    for (l in 1:nvar) varallmean[l] <- mean(x[, l])
    vardmean <- sweep(x, 2, varallmean, "-")
    for (l in 1:nvar) {
      varallmeandist[l] <- sum((vardmean[, l])^2)
      varwgss[l] <- sum(varwithins[, l])
    }
    varbgss <- varallmeandist - varwgss
    vartss <- varbgss + varwgss
    zvargss <- list(vartss = vartss, varbgss = varbgss)
    return(zvargss)
  }
  varwithinss <- function(x, centers, cluster)
  {
    nrow <- dim(centers)[1]
    nvar <- dim(x)[2]
    varwithins <- matrix(0, nrow, nvar)
    x <- (x - centers[cluster, ])^2
    for (l in 1:nvar) {
      for (k in 1:nrow) {
        varwithins[k, l] <- sum(x[cluster == k, l])
      }
    }
    return(varwithins)
  }



  ##############################
  #                            #
  #           SDbw             #
  #                            #
  ##############################

  centers<-function(cl,x)
  {
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)

    centers <- matrix(nrow = k, ncol = ncol(x))
    for (i in 1:k)
    {
      for (j in 1:ncol(x))
      {
        if (is.na(match(i,cl)))
        {
          centers[i, j] <- 0
        }
        else{
          centers[i, j] <- mean(x[cl ==i, j])
        }
      }
    }

    return(centers)
  }

  Average.scattering <- function (cl, x)
  {
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    centers.matrix <- centers(cl,x)

    cluster.size <- numeric(0)
    variance.clusters <- matrix(0, ncol = ncol(x), nrow = k)
    var <- matrix(0, ncol = ncol(x), nrow = k)

    for (u in 1:k)
      cluster.size[u] <- sum(cl == u)

    for (u in 1:k)
    {
      for (j in 1:ncol(x))
      {
        for(i in 1:n)
        {
          if(cl[i]==u)
            variance.clusters[u,j]<- variance.clusters[u,j]+(x[i, j]-centers.matrix[u,j])^2
        }
      }
    }

    for (u in 1:k)
    {
      for (j in 1:ncol(x))
        variance.clusters[u,j]= variance.clusters[u,j]/ cluster.size[u]
    }

    variance.clusters[is.nan(variance.clusters)] <- 0
    variance.matrix <- numeric(0)
    for(j in 1:ncol(x))
      variance.matrix[j]=var(x[,j])*(n-1)/n


    Somme.variance.clusters<-0
    for (u in 1:k)
      Somme.variance.clusters<-Somme.variance.clusters+sqrt((variance.clusters[u,]%*%(variance.clusters[u,])))


    # Standard deviation
    stdev<-(1/k)*sqrt(Somme.variance.clusters)

    #Average scattering for clusters
    scat<- (1/k)* (Somme.variance.clusters /sqrt(variance.matrix %*% variance.matrix))

    scat <- list(stdev=stdev, centers=centers.matrix, variance.intraclusters= variance.clusters, scatt=scat)
    return(scat)
  }

  density.clusters<-function(cl, x)
  {
    x <- as.matrix(x)
    k <- max(cl)
    n <- length(cl)

    distance <- matrix(0, ncol = 1, nrow = n)
    density <-  matrix(0, ncol = 1, nrow = k)
    centers.matrix<-centers(cl,x)
    stdev<-Average.scattering(cl,x)$stdev
    for(i in 1:n)
    {
      u=1
      while(cl[i] != u )
        u<-u+1
      for (j in 1:ncol(x))
      {
        distance[i]<- distance[i]+(x[i,j]-centers.matrix[u,j])^2
      }
      distance[i]<-sqrt(distance[i])
      if (distance[i] <= stdev)
        density[u]= density[u]+1
    }
    dens<-list(distance=distance, density=density)
    return(dens)

  }


  density.bw<-function(cl, x)
  {
    x <- as.matrix(x)
    k <- max(cl)
    n <- length(cl)
    centers.matrix<-centers(cl,x)
    stdev<-Average.scattering(cl,x)$stdev
    density.bw<- matrix(0, ncol = k, nrow = k)
    u<-1

    for(u in 1:k)
    {
      for(v in 1:k)
      {
        if(v!=u)
        {
          distance<- matrix(0, ncol = 1, nrow = n)
          moy<-(centers.matrix[u,]+centers.matrix[v,])/2
          for(i in 1:n)
          {
            if((cl[i]==u)||(cl[i]==v))
            {
              for (j in 1:ncol(x))
              {
                distance[i]<- distance[i]+(x[i,j]-moy[j])^2
              }
              distance[i]<- sqrt(distance[i])
              if(distance[i]<= stdev)
              {
                density.bw[u,v]<-density.bw[u,v]+1
              }
            }
          }
        }
      }
    }
    density.clust<-density.clusters(cl,x)$density
    S<-0
    for(u in 1:k)
      for(v in 1:k)
      {
        if(max(density.clust[u], density.clust[v])!=0)
          S=S+ (density.bw[u,v]/max(density.clust[u], density.clust[v]))
      }
    density.bw<-S/(k*(k-1))
    return(density.bw)

  }

  ########################################################################
  #                                                                      #
  #                                 DB                                   #
  #                                                                      #
  ########################################################################

  Indice.db <- function (x, cl, d = NULL, centrotypes = "centroids", p = 2, q = 2)
  {

    if (sum(c("centroids") == centrotypes) == 0)
      stop("Wrong centrotypes argument")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d <- as.matrix(d)
      }
      row.names(d) <- row.names(x)
    }
    if (is.null(dim(x))) {
      dim(x) <- c(length(x), 1)
    }
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    dAm <- d
    centers <- matrix(nrow = k, ncol = ncol(x))
    if (centrotypes == "centroids") {
      for (i in 1:k) {
        for (j in 1:ncol(x)) {

          if (is.na(match(i,cl)))
          {
            centers[i, j] <- 0
          }
          else {
            centers[i, j] <- mean(x[cl ==i, j])
          }
        }
      }
    }
    else {
      stop("wrong centrotypes argument")
    }
    S <- rep(0, k)
    for (i in 1:k) {
      ind <- (cl == i)
      if (sum(ind) > 1) {
        centerI <- centers[i, ]
        centerI <- rep(centerI, sum(ind))
        centerI <- matrix(centerI, nrow = sum(ind), ncol = ncol(x),
                          byrow = TRUE)
        S[i] <- mean(sqrt(apply((x[ind, ] - centerI)^2, 1,
                                sum))^q)^(1/q)
      }
      else S[i] <- 0
    }
    M <- as.matrix(dist(centers, p = p))
    R <- array(Inf, c(k, k))
    r = rep(0, k)
    for (i in 1:k) {
      for (j in 1:k) {
        R[i, j] = (S[i] + S[j])/M[i, j]
      }
      r[i] = max(R[i, ][is.finite(R[i, ])])
    }
    DB = mean(r[is.finite(r)])
    resul <- list(DB = DB, r = r, R = R, d = M, S = S, centers = centers)
    resul
  }

  ########################################################################
  #                                                                      #
  #                                DUNN                                  #
  #                                                                      #
  ########################################################################

  Indice.dunn <- function(distance=NULL, clusters, Data=NULL, method="euclidean")
  {
    if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
    if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
    if (class(distance)=="dist") distance <- as.matrix(distance)
    n <- unique(clusters)
    nc <- length(n)
    interClust <- matrix(NA, nc, nc)
    intraClust <- rep(NA, nc)

    for (i in 1:nc) {
      c1 <- which(clusters==n[i])
      for (j in i:nc) {
        if (j==i) intraClust[i] <- max(distance[c1,c1])
        if (j>i) {
          c2 <- which(clusters==n[j])
          interClust[i,j] <- min(distance[c1,c2])
        }
      }
    }
    dunn <- min(interClust,na.rm=TRUE)/max(intraClust)
    return(dunn)
  }

  ########################################################################
  #                                                                      #
  #                         SILHOUETTE                                   #
  #                                                                      #
  ########################################################################

  Indice.silhouette<-function(d,cl,singleObject=0)
  {
    d<-as.matrix(d)
    n <- unique(cl)
    nc <- length(n)
    Si<-0
    for(k in 1:nc)
    {
      if ((sum(cl==n[k]))<=1)
        Sil<-singleObject
      else
      {
        Sil<-0
        for(i in 1:length(cl))
        {
          if(cl[i]==n[k])
          {
            ai<-sum(d[i,cl==n[k]])/(sum(cl==n[k])-1)
            dips<-NULL
            for(j in 1:nc)
              if (cl[i]!=j)
                if(sum(cl==n[j])!=1)
                  dips<-cbind(dips,c( (sum(d[i,cl==n[j]])) /(sum(cl==n[j])) ))
            else
              dips<-cbind(dips,c( (sum(d[i,cl==n[j]]))))
            bi<-min(dips)
            Sil<-Sil+(bi-ai)/max(c(ai,bi))
          }
        }

      }
      Si<-Si+Sil
    }
    Si/length(cl)

  }

  ########################################################################
  #                                                                      #
  #                             POINT-BISERIAL                           #
  #                                                                      #
  ########################################################################

  Indice.ptbiserial <- function (x,md,cl1)
  {

    nn <- dim(x)[1]
    pp <- dim(x)[2]

    md2 <- as.matrix(md)
    m01 <- array(NA, c(nn,nn))
    nbr <- (nn*(nn-1))/2
    pb <- array(0,c(nbr,2))

    m3 <- 1
    for (m1 in 2:nn){
      m12 <- m1-1
      for (m2 in 1:m12){
        if (cl1[m1]==cl1[m2]) m01[m1,m2]<-0
        if (cl1[m1]!=cl1[m2]) m01[m1,m2]<-1
        pb[m3,1] <- m01[m1,m2]
        pb[m3,2] <- md2[m1,m2]
        m3 <- m3+1
      }
    }

    y <- pb[,1]
    x <- pb[,2]

    biserial.cor <- function (x, y, use = c("all.obs", "complete.obs"), level = 1){

      if (!is.numeric(x))
        stop("'x' must be a numeric variable.\n")
      y <- as.factor(y)
      if (length(levs <- levels(y)) > 2)
        stop("'y' must be a dichotomous variable.\n")
      if (length(x) != length(y))
        stop("'x' and 'y' do not have the same length")
      use <- match.arg(use)
      if (use == "complete.obs") {
        cc.ind <- complete.cases(x, y)
        x <- x[cc.ind]
        y <- y[cc.ind]
      }
      ind <- y == levs[level]
      diff.mu <- mean(x[ind]) - mean(x[!ind])
      prob <- mean(ind)
      diff.mu * sqrt(prob * (1 - prob))/sd(x)
    }

    ptbiserial <- biserial.cor(x=pb[,2], y=pb[,1], level = 2)
    return(ptbiserial)
  }

  ########################################################################
  #                                                                      #
  #                              CH                                      #
  #                                                                      #
  ########################################################################

  Indice.ch <- function (x, cl, d = NULL, centrotypes = "centroids")
  {
    nb.cl1 <- table(cl)
    nb1.cl1 <- sum(nb.cl1==1)
    if (nb1.cl1 > 0){
      CH <- NA
    }
    if (sum(c("centroids", "medoids") == centrotypes) == 0)
      stop("Wrong centrotypes argument")
    if ("medoids" == centrotypes && is.null(d))
      stop("For argument centrotypes = 'medoids' d cannot be null")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d <- as.matrix(d)
      }
      row.names(d) <- row.names(x)
    }

    n <- length(cl)
    k <- length(unique(cl))
    CH <- (gss(x, cl, d)$bgss/(k-1))/
      (gss(x, cl, d)$wgss/(n-k))
    return(CH)
  }


  ########################################################################
  #                                                                      #
  #                              C-INDEX                                 #
  #                                                                      #
  ########################################################################

  Indice.cindex <- function (d, cl)
  {

    d <- data.matrix(d)
    DU <- 0
    r <- 0
    v_max <- array(1, max(cl))
    v_min <- array(1, max(cl))
    for (i in 1:max(cl)){
      n <- sum(cl == i)
      if (n > 1){
        t <- d[cl == i, cl == i]
        DU = DU + sum(t)/2
        v_max[i] = max(t)
        if (sum(t == 0) == n)
          v_min[i] <- min(t[t != 0])
        else v_min[i] <- 0
        r <- r + n * (n - 1)/2
      }
    }
    Dmin = min(v_min)
    Dmax = max(v_max)
    if (Dmin == Dmax)
      result <- NA
    else result <- (DU - r * Dmin)/(Dmax * r - Dmin * r)
    result
  }

  ########################################################################
  #                                                                      #
  #                            RATKOWSKY                                 #
  #                                                                      #
  ########################################################################

  Indice.ratkowsky <- function(x, cl, d, centrotypes = "centroids")
  {

    n <- unique(cl)
    qq <- length(n)
    centers <- gss(x, cl, d)$centers
    varwithins <- varwithinss(x, centers, cl)
    zvargss <- vargss(x, cl, varwithins)
    ratio <- mean(sqrt(zvargss$varbgss/zvargss$vartss))
    ratkowsky <- ratio/sqrt(qq)
    return(ratkowsky)
  }

  ########################################################################
  #                                                                      #
  #                           MCCLAIN                                    #
  #                                                                      #
  ########################################################################

  Indice.mcclain <- function (cl1,md)
  {

    m <- unique(cl1)
    cn1 <- length(m)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1){

      cluster.size[u] <- sum(cl1 ==m[u])
      du <- as.dist(dmat[cl1 == m[u], cl1 == m[u]])
      within.dist1 <- c(within.dist1, du)
      for (v in 1:cn1){
        if (v != u) {
          suv <- dmat[cl1 == m[u], cl1 == m[v]]
          if (u < v) {
            separation.matrix[u, v] <- separation.matrix[v,u] <- min(suv)
            between.dist1 <- c(between.dist1, suv)
          }
        }
      }
    }

    nwithin1 <- length(within.dist1)
    nbetween1 <- length(between.dist1)
    meanwithin1 <- mean(within.dist1)
    meanbetween1 <- mean(between.dist1)
    mcclain <- (meanwithin1/nwithin1)/(meanbetween1/nbetween1)
    return(mcclain)
  }


  ########################################################################
  #                                                                      #
  #                              GAMMA                                   #
  #                                                                      #
  ########################################################################

  Indice.gamma <- function (cl1,md)
  {

    n <- unique(cl1)
    cn1 <- length(n)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1){
      x <- n[u]
      cluster.size[u] <- sum(cl1 == x)
      du <- as.dist(dmat[cl1 == x, cl1 == x])
      within.dist1 <- c(within.dist1, du)
      average.distance[u] <- mean(du)
      median.distance[u] <- median(du)
      bv <- numeric(0)
      for (v in 1:cn1){
        if (v != u) {
          x1 <- n[v]
          suv <- dmat[cl1 == x, cl1 == x1]
          bv <- c(bv, suv)
          if (u < v) {
            separation.matrix[u, v] <- separation.matrix[v,u] <- min(suv)
            between.dist1 <- c(between.dist1, suv)
          }
        }
      }
    }

    nwithin1 <- length(within.dist1)

    s.plus <- s.moins <- 0
    for (k in 1: nwithin1){

      s.plus <- s.plus+(colSums(outer(between.dist1,within.dist1[k], ">")))
      s.moins <- s.moins+(colSums(outer(between.dist1,within.dist1[k], "<")))
    }

    Gamma <- (s.plus-s.moins)/(s.plus+s.moins)
    return(Gamma)
  }


  ########################################################################
  #                                                                      #
  #                              GPLUS                                   #
  #                                                                      #
  ########################################################################

  Indice.gplus <- function (cl1,md)
  {

    n <- unique(cl1)
    cn1 <- length(n)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1) {
      x <-  n[u]
      cluster.size[u] <- sum(cl1 == x)
      du <- as.dist(dmat[cl1 == x, cl1 == x])
      within.dist1 <- c(within.dist1, du)
      average.distance[u] <- mean(du)
      median.distance[u] <- median(du)
      bv <- numeric(0)
      for (v in 1:cn1) {
        if (v != u) {
          x1 <- n[v]
          suv <- dmat[cl1 == x, cl1 == x1]
          bv <- c(bv, suv)
          if (u < v) {
            separation.matrix[u, v] <- separation.matrix[v,u] <- min(suv)
            between.dist1 <- c(between.dist1, suv)
          }
        }
      }
    }

    nwithin1 <- length(within.dist1)

    s.plus <- s.moins <- 0
    for (k in 1: nwithin1){

      s.plus <- s.plus+(colSums(outer(between.dist1,within.dist1[k], ">")))
      s.moins <- s.moins+(colSums(outer(between.dist1,within.dist1[k], "<")))
    }

    Gplus <- (2*s.moins)/(n1*(n1-1))
    return(Gplus)
  }

  ########################################################################
  #                                                                      #
  #                              TAU                                     #
  #                                                                      #
  ########################################################################

  Indice.tau <- function (cl1,md)
  {

    n <- unique(cl1)
    cn1 <- length(n)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1) {
      x <-  n[u]
      cluster.size[u] <- sum(cl1 == x)
      du <- as.dist(dmat[cl1 == x, cl1 == x])
      within.dist1 <- c(within.dist1, du)
      average.distance[u] <- mean(du)
      median.distance[u] <- median(du)
      bv <- numeric(0)
      for (v in 1:cn1) {
        if (v != u) {
          x1 <- n[v]
          suv <- dmat[cl1 == x, cl1 == x1]
          bv <- c(bv, suv)
          if (u < v) {
            separation.matrix[u, v] <- separation.matrix[v,u] <- min(suv)
            between.dist1 <- c(between.dist1, suv)
          }
        }
      }
    }

    nwithin1 <- length(within.dist1)
    nbetween1 <- length(between.dist1)

    s.plus <- s.moins <- 0
    for (k in 1: nwithin1){

      s.plus <- s.plus+(colSums(outer(between.dist1,within.dist1[k], ">")))
      s.moins <- s.moins+(colSums(outer(between.dist1,within.dist1[k], "<")))
    }

    t.tau  <- (nwithin1*nbetween1)-(s.plus+s.moins)
    tau <- (s.plus-s.moins)/(((n1*(n1-1)/2-t.tau)*(n1*(n1-1)/2))^(1/2))
    return(tau)
  }


  ########################################################################
  #                                                                      #
  #       CCC,SCOTT,MARRIOT,TRCOVW,TRACEW,FRIEDMAN,RUBIN                 #
  #                                                                      #
  ########################################################################

  Indices.WBT <- function(x,cl,P,s,vv)
  {

    n <- dim(x)[1]
    pp <- dim(x)[2]
    nc <- unique(cl)
    qq <- length(nc)
    z <- matrix(0,ncol=qq,nrow=n)
    clX <- as.matrix(cl)

    for (i in 1:n)
      for (j in 1:qq){
        z[i,j]==0
        if (clX[i,1]== nc[j])
        {z[i,j]=1}
      }

    xbar <- solve(t(z)%*%z)%*%t(z)%*%x
    B <- t(xbar)%*%t(z)%*%z%*%xbar
    W <- P-B
    marriot <- (qq^2)*det(W)
    trcovw <- sum(diag(cov(W)))
    tracew <- sum(diag(W))
    if(det(W)!=0)
      scott<-n*log(sum(diag(P))/sum(diag(W)))
    else {cat("Error: division by zero!")}
    friedman <- sum(diag(solve(W)*B))
    rubin <- sum(diag(P))/sum(diag(W))



    R2 <- 1-sum(diag(W))/sum(diag(P))
    v1 <- 1
    u <- rep(0,pp)
    c <- (vv/(qq))^(1/pp)
    u <- s/c
    k1 <- sum((u>=1)==TRUE)
    p1 <- min(k1,qq-1)
    if (all(p1>0,p1<pp))
    {
      for (i in 1:p1)
        v1 <- v1*s[i]
      c <- (v1/(qq))^(1/p1)
      u <- s/c
      b1 <- sum(1/(n+u[1:p1]))
      b2 <- sum(u[p1+1:pp]^2/(n+u[p1+1:pp]),na.rm=TRUE)
      E_R2 <- 1-((b1+b2)/sum(u^2))*((n-qq)^2/n)*(1+4/n)
      ccc <- log((1-E_R2)/(1-R2))*(sqrt(n*p1/2)/((0.001+E_R2)^1.2))
    }else {
      b1 <- sum(1/(n+u))
      E_R2 <- 1-(b1/sum(u^2))*((n-qq)^2/n)*(1+4/n)
      ccc <- log((1-E_R2)/(1-R2))*(sqrt(n*pp/2)/((0.001+E_R2)^1.2))
    }
    results <- list(ccc=ccc,scott=scott,marriot=marriot,trcovw=trcovw,tracew=tracew,friedman=friedman,rubin=rubin)
    return(results)
  }

  ########################################################################
  #                                                                      #
  #                              BALL                                    #
  #                                                                      #
  ########################################################################

  Indice.ball <- function(x, cl, d = NULL, centrotypes = "centroids")
  {

    wgssB <- gss(x, cl, d)$wgss
    qq <- length(unique(cl))
    ball <- wgssB/qq
    return(ball)
  }

  ########################################################################
  #                                                                      #
  #                             SDBW                                     #
  #                                                                      #
  ########################################################################

  Indice.sdbw<-function(x, cl)
  {

    x <- as.matrix(x)
    Scatt <- Average.scattering(cl,x)$scatt
    Dens.bw <- density.bw(cl,x)
    SDbw <- Scatt+Dens.bw
    return(SDbw)
  }

  ########################################################################
  #                                                                      #
  #                             D INDEX                                  #
  #                                                                      #
  ########################################################################

  Indice.dindex<- function(cl, x)
  {

    x <- as.matrix(x)
    distance<-density.clusters(cl, x)$distance
    n<-length(distance)
    S<-0
    for(i in 1:n)
      S<-S+distance[i]
    inertieIntra<-S/n
    return(inertieIntra)
  }

  ########################################################################
  #                                                                      #
  #                             HUBERT                                   #
  #                                                                      #
  ########################################################################

  Indice.hubert<-function(x, cl)
  {

    k <- max(cl)
    n<-dim(x)[1]
    y <- matrix(0, ncol = dim(x)[2], nrow = n)
    md <- dist(x, method="euclidean")
    P<- as.matrix(md)
    meanP<-mean(P)
    variance.matrix <- numeric(0)

    for(j in 1:n)
      variance.matrix[j]=var(P[,j])*(n-1)/n
    varP<-sqrt(variance.matrix %*% variance.matrix)

    centers.clusters<-centers(cl,x)
    for(i in 1:n){
      for(u in 1:k){
        if(cl[i]==u)
          y[i,]<-centers.clusters[u,]
      }
    }

    Q<- as.matrix(dist(y, method="euclidean"))
    meanQ<-mean(Q)
    for(j in 1:n)
      variance.matrix[j]=var(Q[,j])*(n-1)/n
    varQ<-sqrt(variance.matrix %*% variance.matrix)

    M<-n*(n-1)/2
    S<-0
    n1<-n-1

    for(i in 1:n1){
      j<-i+1
      while(j<=n)
      {
        S<-S+(P[i,j]-meanP)*(Q[i,j]-meanQ)
        j<-j+1
      }

    }
    hubert <-S/(M*varP*varQ)

    return(hubert)
  }

  ########################################################################
  #                                                                      #
  #                                SV                                    #
  #                                                                      #
  ########################################################################

  Indice.sv <- function(cl,data)
  {
    n <- unique(cl)
    nc <- length(n)
    den <- rep(NA, nc)
    num <- rep(NA, nc)
    centers <- centers(cl,data)
    d <- as.matrix(dist(centers, p = 2))

    for (i in 1:nc) {
      c1 <- which(cl == n[i])
      distance <- rep(NA, length(c1))
      for( i1 in 1: length(c1)) {
        distance[i1] <- sqrt(sum((data[c1[i1],] - centers[i, ])^2))
      }
      den[i] <- max(distance)
      dis <- d[i,]
      num[i] <- min(dis[dis != 0])
    }

    sv <- sum(num)/sum(den)
    return(sv)
  }

  ########################################################################
  #                                                                      #
  #                             XIE-BENI                                 #
  #                                                                      #
  ########################################################################

  Indice.XieBeni <- function(x,d, cl)
  {

    distance <- as.matrix(d)
    n <- unique(cl)
    nc <- length(n)
    interClust <- matrix(NA, nc, nc)
    wgssB <- gss(x, cl,d)$wgss
    N <-dim(x)[1]

    for (i in 1:nc) {
      c1 <- which(cl == n[i])
      for (j in i:nc) {
        if (j>i) {
          c2 <- which(cl == n[j])
          interClust[i,j] <- min(distance[c1,c2])
        }
      }
    }

    xiebeni <- (wgssB/N) / min(interClust,na.rm=TRUE)^2
    return(xiebeni)
  }

  ########################################################################
  #                                                                      #
  #                             HARTIGAN                                 #
  #                                                                      #
  ########################################################################

  Indice.hartigan <- function(x,cl,d)
  {
    bgss<- gss(x, cl, d)$bgss
    wgss<- gss(x, cl, d)$wgss
    hart <- log(bgss/wgss)
    return(hart)
  }

  ########################################################################
  #                                                                      #
  #                               SSI                                    #
  #                                                                      #
  ########################################################################

  Indice.ssi <- function (x,cl)
  {

    centers <-centers(cl,x)
    ncl <- dim(centers)[1]
    nvar <- dim(centers)[2]
    n <- sum(cl)

    cmax <- apply(centers, 2, max)
    cmin <- apply(centers, 2, min)
    cord <- apply(centers, 2, order)
    cmaxi <- cord[ncl,]
    cmini <- cord[1,]

    meanmean <- mean(centers)
    absmdif <- abs(apply(centers, 2, mean) - meanmean)
    span <- cmax - cmin
    csizemax <- cl[cmaxi]
    csizemin <- cl[cmini]

    hiest <- nvar
    hiestw <- hiest * max(span) * max(max(csizemax), max(csizemin)) * exp(-min(absmdif))

    sist <- sum(span)/hiest

    sistw <- (span * exp(-absmdif)) %*% sqrt(csizemax*csizemin) / hiestw

    # return(list(ssi=sist, ssiw=sistw))
    return(sistw)

  }

  ########################################################################
  #                                                                      #
  #                             XU                                       #
  #                                                                      #
  ########################################################################

  Indice.xu <- function(x, cl,d)
  {
    n <- sum(cl)
    k <- length(cl)
    d <- dim(x)[2]
    wgss<- gss(x, cl, d)$wgss

    xuindex <- d * log(sqrt(wgss/(d*(n^2)))) + log(k)
    return(xuindex)
  }


  ########################################################################
  #                                                                      #
  #                             RAY-TURI                                 #
  #                                                                      #
  ########################################################################

  Indice.rayturi <- function(x,cl,Dist)
  {
    n <- dim(x)[1]
    centers <- centers(cl,x)
    d<- dist(centers)^2
    rayInter <- min(d[d != 0])*n
    rayIntra <- gss(x,cl,Dist)$wgss
    ray <- rayIntra/rayInter
    return(ray)
  }


  ########################################################################
  #                                                                      #
  #                                 PBM                                  #
  #                                                                      #
  ########################################################################

  Indice.pbm<- function(data,cl)
  {
    n <-dim(data)[1]
    p <-dim(data)[2]
    k <- length(unique(cl))

    centers <-centers(cl,data)
    d <- max(dist(centers, p = 2))
    bar <- apply(data, 2, mean)
    centers1 <-matrix(rep(bar,n),nrow=n,ncol=p,byrow=TRUE)
    ew <- sqrt(sum((data-centers1)^2))
    et <- sqrt(sum((data-centers[cl,])^2))
    pbm <-((et*d)/(k*ew))^2
    return(pbm)
  }

  ########################################################################
  #                                                                      #
  #                          BANFELD-RAFTERY                             #
  #                                                                      #
  ########################################################################

  Indice.banfeld<- function(x, cl,d)
  {
    wgss<- gss(x, cl, d)$gss
    k<- max(cl)
    ban<- rep(0,k)
    for (i in 1:k)
    {
      nk  <-sum(cl == i)
      if (wgss[i] !=0)
        ban[i] <- nk*log(wgss[i]/nk)
    }
    return(sum(ban))
  }


  ######################################
  result <- array(0, c(max_nc-1,30))
  col <- rep(0,max_nc-1)
  nb <- rep(0,max_nc-1)
  n <- max_nc-1
  ne <- list()
  lis <-list()
  lis1 <-list()
  for (i in 1:max_nc){
    col[i-1]<- paste(i,i,sep="x")
  }
  rownames(result)<- rev(col)
  colnames(result) <- c("db","dunn","silhouette","ptbiserial","ch","cindex","ratkowsky","mcclain","gamma","gplus",
                        "tau","ccc","scott","marriot","trcovw","tracew","friedman","rubin","ball","sdbw",
                        "dindex","hubert","sv","xie-beni","hartigan","ssi","xu","rayturi","pbm","banfeld")

  n1<-n2<-n3<-n4<-n5<-n6<-n7<-n8<-n9<-n10<-n11<-n12<-n13<-n14<-n15<-n16<-n17<-
    n18<-n19<-n20<-n21<-n22<-n23<-n24<-n25<-n26<-n27<-n28<-n29<-n30<-1

  ###################################################################################
  repeat {

    res.batch <- BatchSOM(data,grid = somgrid(xheight, xwidth,topo), min.radius ,max.radius,maxit,init,radius.type)
    codes <- res.batch$code
    cl  <- res.batch$classif
    Dist <- dist(data, method="euclidean")

    ##################################### Indice.db###################################
    if(is.na(match(1,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n1,1] <-Indice.db(data, cl, d = NULL, centrotypes = "centroids", p = 2, q = 2)$DB
      lis1[[n1]] <- cl
      n1 <-n1+1
    }

    #################################### Indice.dunn#################################
    if(is.na(match(2,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n2,2]<-Indice.dunn(Dist, cl,Data=NULL, method="euclidean")
      lis1[[n2]] <- cl
      n2 <-n2+1
    }

    ################################## Indice.silhouette############################
    if(is.na(match(3,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n3,3] <-Indice.silhouette(Dist,cl,singleObject=0)
      lis1[[n3]] <- cl
      n3 <-n3+1
    }

    ################################# Indice.ptbiserial#############################
    if(is.na(match(4,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n4,4]<-Indice.ptbiserial(data,Dist,cl)
      lis1[[n4]] <- cl
      n4 <-n4+1
    }

    ################################ Indice.ch #####################################
    if(is.na(match(5,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n5,5]<-Indice.ch(data, cl, d = NULL, centrotypes = "centroids")
      lis1[[n5]] <- cl
      n5 <-n5+1
    }

    ############################### Indice.cindex##################################
    if(is.na(match(6,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n6,6]<-Indice.cindex(Dist, cl)
      lis1[[n6]] <- cl
      n6 <-n6+1
    }

    ############################## Indice.ratkowsky###############################
    if(is.na(match(7,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n7,7]<-Indice.ratkowsky(data,cl,Dist, centrotypes = "centroids")
      lis1[[n7]] <- cl
      n7 <-n7+1
    }

    ############################# Indice.mcclain#################################
    if(is.na(match(8,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n8,8] <-Indice.mcclain(cl,Dist)
      lis1[[n8]] <- cl
      n8 <-n8+1
    }

    ############################ Indice.gamma###################################
    if(is.na(match(9,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n9,9] <-Indice.gamma(cl ,Dist)
      lis1[[n9]] <- cl
      n9 <-n9+1
    }

    ############################ Indice.gplus####################################
    if(is.na(match(10,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n10,10] <-Indice.gplus(cl,Dist)
      lis1[[n10]] <- cl
      n10 <-n10+1
    }

    ########################### Indice.tau######################################
    if(is.na(match(11,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n11,11]<-Indice.tau(cl,Dist)
      lis1[[n11]] <- cl
      n11 <-n11+1
    }

    ########################### Indice.ccc######################################
    if(is.na(match(12,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n12,12] <-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$ccc
      lis1[[n12]] <- cl
      n12 <-n12+1
    }

    ######################## Indice.scott######################################
    if(is.na(match(13,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n13,13]<-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$scott
      lis1[[n13]] <- cl
      n13 <-n13+1
    }

    ####################### Indice.marriot######################################
    if(is.na(match(14,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n14,14]<-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$marriot
      lis1[[n14]] <- cl
      n14 <-n14+1
    }

    ###################### Indice.trcovw######################################
    if(is.na(match(15,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n15,15]<-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$trcovw
      lis1[[n15]] <- cl
      n15 <-n15+1
    }

    ##################### Indice.tracew######################################
    if(is.na(match(16,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n16,16] <-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$tracew
      lis1[[n16]] <- cl
      n16 <-n16+1
    }

    #################### Indice.friedman#####################################
    if(is.na(match(17,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n17,17] <-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$friedman
      lis1[[n17]] <- cl
      n17 <-n17+1
    }

    ################### Indice.rubin##########################################
    if(is.na(match(18,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n18,18]<-Indices.WBT(data,cl,P=TT,s=ss,vv=vv)$rubin
      lis1[[n18]] <- cl
      n18 <-n18+1
    }

    ######################## Indice.ball########################################
    if(is.na(match(19,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n19,19] <-Indice.ball(data, cl, d = NULL, centrotypes = "centroids")
      lis1[[n19]] <- cl
      n19 <-n19+1
    }

    ######################### Indice.sdbw#######################################
    if(is.na(match(20,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n20,20] <-Indice.sdbw(data, cl)
      lis1[[n20]] <- cl
      n20 <-n20+1
    }

    ######################### Indice.dindex#######################################
    if(is.na(match(21,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n21,21]<-Indice.dindex(cl,data)
      lis1[[n21]] <- cl
      n21 <-n21+1
    }

    ######################## Indice.hubert########################################
    if(is.na(match(22,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n22,22] <-Indice.hubert(data, cl)
      lis1[[n22]] <- cl
      n22 <-n22+1
    }

    ########################### Indice.sv#########################################
    if(is.na(match(23,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n23,23]<-Indice.sv(cl,data)
      lis1[[n23]] <- cl
      n23 <-n23+1
    }

    ########################## Indice.xie-beni##########################################
    if(is.na(match(24,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n24,24]<-Indice.XieBeni(data,Dist, cl)
      lis1[[n24]] <- cl
      n24 <-n24+1
    }

    ######################## Indice.Hartigan########################################
    if(is.na(match(25,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n25,25]<-Indice.hartigan(data,cl,Dist)
      lis1[[n25]] <- cl
      n25 <-n25+1
    }

    ########################## Indice.ssi#############################################
    if(is.na(match(26,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n26,26]<-Indice.ssi(data,cl)
      lis1[[n26]] <- cl
      n26 <-n26+1
    }

    ############################ Indice.xu ###########################################
    if(is.na(match(27,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n27,27]<-Indice.xu(data, cl,Dist)
      lis1[[n27]] <- cl
      n27<-n27+1
    }

    ############################ Indice.rayturi########################################
    if(is.na(match(28,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n28,28]<-Indice.rayturi(data,cl)
      lis1[[n28]] <- cl
      n28 <-n28+1
    }

    ############################ Indice.pbm ###########################################
    if(is.na(match(29,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n29,29]<-Indice.pbm(data,cl)
      lis1[[n29]] <- cl
      n29<-n29+1
    }

    ############################ Indice.banfeld ######################################
    if(is.na(match(30,indice))==FALSE||is.na(match(31,indice))==FALSE)
    {
      result[n30,30]<-Indice.banfeld(data, cl,Dist)
      lis1[[n30]] <- cl
      n30 <-n30+1
    }

    data<- codes
    xheight <- xheight - 1
    xwidth <- xwidth - 1

    if (xheight ==1 & xwidth == 1)

      break
  }

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#                                                                                                                      #
  #                                                                                                                      #
  #                                            Number of Clusters                                                        #
  #                                                                                                                      #
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

  for(i in 1:length(lis1))
  {
    lis[[i]]<-unique(lis1[[i]])
  }

  dim <- rep(0,n)
  vide <-a <- a1 <-list()
  j <- j1 <-1
  for (i1 in  max_nc:2)
  {
    dim[j1] <-i1^2
    j1 <-j1+1
    for(i in 1:length(dim))
    {
      for(k in 1:dim[i])
      {
        if(is.na(match(k,lis[[i]]))==TRUE)
          vide[[j]] <-k
        j <-j+1
      }
      a[[i]] <- unlist(vide)
      vide <-list()
      if(length(lis[[i]])== dim[i])
        a[[i]] <-0
    }
  }


  nb[1] <- length(unique(lis[[1]]))
  for(i in 2:n)
  {
    nc <-lis1[[i]]
    c <- length(lis[[i]])
    m <- a[[i-1]]
    m2 <-a[[i]]

    k <-1
    while(k <=length(m))
    {
      n1 <- m[k]
      len <-length(which(lis1[[i]]== nc[n1]))

      if(len ==1){
        c <- c-1
        m2[length(m2)+1]<- nc[n1]
      }
      if (len > 1)
      {
        en <- 0
        m1 <- which(lis1[[i]]==nc[n1])
        for (nn in 1:len)
        {
          if((is.na(match(m1[nn],a[[i-1]]))== FALSE))
          {
            en <- en +1
          }
        }
        if (en == len){
          c <-c-1
          for (nn in 1:len){
            m <- m[-which(m ==m1[nn])]
          }
          m2[length(m2)+1]<- nc[n1]
          k <-k-1
        }
      }
      k <- k+1
    }

    nb[i] <-c
    a[[i]] <- m2
    j <-1
    for(i in 1:length(dim))
    {
      for(k in 1:dim[i])
      {
        if(is.na(match(k,a[[i]]))==TRUE)
          a1[[j]] <-k
        j <-j+1
      }
      v <-unlist(a1)
      ne[[i]] <-v
      a1 <-list()
    }
  }


  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#                                                                                                                      #
  #                                                                                                                      #
  #                                            Best Number of Clusters                                                   #
  #                                                                                                                      #
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

  nc.DB<-indice.DB<-0
  if(is.na(match(1,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    #DB- The minimum value of the index is used to indicate the optimal number of clusters.
    nc.DB <- as.vector(nb[which.min(result[,1])])
    indice.DB <- min(result[,1],na.rm = TRUE)
    best.nc<-nc.DB
  }

  nc.Dunn<-indice.Dunn<-0
  if(is.na(match(2,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    #DUNN - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.Dunn <- as.vector(nb[which.max(result[,2])])
    indice.Dunn <- max(result[,2],na.rm = TRUE)
    best.nc<-nc.Dunn
  }

  nc.Silhouette<-indice.Silhouette<-0
  if(is.na(match(3,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # SILHOUETTE - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.Silhouette <- as.vector(nb[which.max(result[,3])])
    indice.Silhouette <- max(result[,3],na.rm = TRUE)
    best.nc<-nc.Silhouette
  }

  nc.ptbiserial<-indice.ptbiserial<-0
  if(is.na(match(4,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # PTBISERIAL - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.ptbiserial <- as.vector(nb[which.max(result[,4])])
    indice.ptbiserial <- max(result[,4],na.rm = TRUE)
    best.nc<-nc.ptbiserial
  }

  nc.CH<-indice.CH<-0
  if(is.na(match(5,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    #CH- The maximum value of the index is used to indicate the optimal number of clusters.
    nc.CH <- as.vector(nb[which.max(result[,5])])
    indice.CH <- max(result[,5],na.rm = TRUE)
    best.nc<-nc.CH
  }

  nc.cindex<-indice.cindex<-0
  if(is.na(match(6,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # CINDEX - The minimum value of the index is used to indicate the optimal number of clusters.
    nc.cindex <- as.vector(nb[which.min(result[,6])])
    indice.cindex <- min(result[,6],na.rm = TRUE)
    best.nc<-nc.cindex
  }

  nc.Ratkowsky<-indice.Ratkowsky<-0
  if(is.na(match(7,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # TARKOWSKY - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.Ratkowsky <- as.vector(nb[which.max(result[,7])])
    indice.Ratkowsky <- max(result[,7],na.rm = TRUE)
    best.nc<-nc.Ratkowsky
  }

  nc.McClain<-indice.McClain<-0
  if(is.na(match(8,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # MCCLAIN - The minimum value of the index is used to indicate the optimal number of clusters.
    nc.McClain <- as.vector(nb[which.min(result[,8])])
    indice.McClain <- min(result[,8],na.rm = TRUE)
    best.nc<-nc.McClain
  }

  nc.Gamma<-indice.Gamma<-0
  if(is.na(match(9,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # GAMMA - The maximum value of the index is taken to represent the correct number of clusters
    nc.Gamma <- as.vector(nb[which.max(result[,9])])
    indice.Gamma <- max(result[,9],na.rm = TRUE)
    best.nc<-nc.Gamma
  }

  nc.Gplus<-indice.Gplus<-0
  if(is.na(match(10,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # GPLUS - The minimum values of the index are used to determine the optimal number of clusters in the data
    nc.Gplus <- as.vector(nb[which.min(result[,10])])
    indice.Gplus <- min(result[,10],na.rm = TRUE)
    best.nc<-nc.Gplus
  }

  nc.Tau<-indice.Tau<-0
  if(is.na(match(11,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # TAU - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.Tau <- as.vector(nb[which.max(result[,11])])
    indice.Tau <- max(result[,11],na.rm = TRUE)
    best.nc<-nc.Tau
  }

  nc.CCC<-indice.CCC<-0
  if(is.na(match(12,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # CCC - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.CCC <-as.vector(nb[which.max(result[,12])])
    indice.CCC <- max(result[,12],na.rm = TRUE)
    best.nc<-nc.CCC
  }

  #Some indices need to compute difference between hierarchy levels to identify relevant number of clusters
  if((is.na(match(13,indice))==FALSE)||(is.na(match(14,indice))==FALSE)||(is.na(match(15,indice))==FALSE)||
     (is.na(match(16,indice))==FALSE)||(is.na(match(17,indice))==FALSE)||(is.na(match(18,indice))==FALSE)||
     (is.na(match(19,indice))==FALSE)||(is.na(match(21,indice))==FALSE)||(is.na(match(22,indice))==FALSE)||
     (is.na(match(25,indice))==FALSE)||(is.na(match(27,indice))==FALSE)||(is.na(match(31,indice))==FALSE))

  {

    Diff <- array(0, c(max_nc-min_nc+1,12))
    Diff[,1] <- min_nc:max_nc
    for (n in min_nc:max_nc)
    {
      if (n==min_nc)
      {
        Diff[n-min_nc+1,2] <- abs(result[n-min_nc+1,13]-NA)   #Scott
        Diff[n-min_nc+1,3] <- abs(result[n-min_nc+1,14]-NA)   #Marriot
        Diff[n-min_nc+1,4] <- abs(result[n-min_nc+1,15]-NA)   #Trcovw
        Diff[n-min_nc+1,5] <- abs(result[n-min_nc+1,16]-NA)   #Tracew
        Diff[n-min_nc+1,6] <- abs(result[n-min_nc+1,17]-NA)   #Friedman
        Diff[n-min_nc+1,7] <- abs(result[n-min_nc+1,18]-NA)   #Rubin
        Diff[n-min_nc+1,8] <- abs(result[n-min_nc+1,21]-NA)   #D index
        Diff[n-min_nc+1,9] <- abs(result[n-min_nc+1,19]-NA)   #Ball
        Diff[n-min_nc+1,10]<- abs(result[n-min_nc+1,22]-NA)   #Hubert
        Diff[n-min_nc+1,11]<- abs(result[n-min_nc+1,25]-NA)   #hartigan
        Diff[n-min_nc+1,12]<- abs(result[n-min_nc+1,27]-NA)   #xu
      }
      else
      {
        if(n==max_nc)
        {
          Diff[n-min_nc+1,2] <- abs(result[n-min_nc+1,13]-result[n-min_nc,13])
          Diff[n-min_nc+1,3] <- abs(result[n-min_nc+1,14]-NA) # Marriot
          Diff[n-min_nc+1,4] <- abs(result[n-min_nc+1,15]-result[n-min_nc,15]) # trcovw
          Diff[n-min_nc+1,5] <- abs(result[n-min_nc+1,16]-NA) #traceW
          Diff[n-min_nc+1,6] <- abs(result[n-min_nc+1,17]-result[n-min_nc,17])
          Diff[n-min_nc+1,7] <- abs(result[n-min_nc+1,18]-NA) #Rubin
          Diff[n-min_nc+1,8] <- abs(result[n-min_nc+1,21]-NA) # D index
          Diff[n-min_nc+1,9] <- abs(result[n-min_nc+1,19]-result[n-min_nc,19])
          Diff[n-min_nc+1,10]<- abs(result[n-min_nc+1,22]-NA)
          Diff[n-min_nc+1,11]<- abs(result[n-min_nc+1,25]-result[n-min_nc,25])
          Diff[n-min_nc+1,12]<- abs(result[n-min_nc+1,27]-NA)   #xu


        }
        else
        {

          Diff[n-min_nc+1,2] <- abs(result[n-min_nc+1,13]-result[n-min_nc,13])
          Diff[n-min_nc+1,3] <- ((result[n-min_nc+2,14]-result[n-min_nc+1,14])-(result[n-min_nc+1,14]-result[n-min_nc,14]))
          Diff[n-min_nc+1,4] <- abs(result[n-min_nc+1,15]-result[n-min_nc,15])
          Diff[n-min_nc+1,5] <- ((result[n-min_nc+2,16]-result[n-min_nc+1,16])-(result[n-min_nc+1,16]-result[n-min_nc,16]))
          Diff[n-min_nc+1,6] <- abs(result[n-min_nc+1,17]-result[n-min_nc,17])
          Diff[n-min_nc+1,7] <- ((result[n-min_nc+2,18]-result[n-min_nc+1,18])-(result[n-min_nc+1,18]-result[n-min_nc,18]))
          Diff[n-min_nc+1,8] <- ((result[n-min_nc+2,21]-result[n-min_nc+1,21])-(result[n-min_nc+1,21]-result[n-min_nc,21]))
          Diff[n-min_nc+1,9] <- abs(result[n-min_nc+1,19]-result[n-min_nc,19])
          Diff[n-min_nc+1,10]<- abs(result[n-min_nc+1,22]-result[n-min_nc,22])
          Diff[n-min_nc+1,11]<- abs(result[n-min_nc+1,25]-result[n-min_nc,25]) # Hartigan
          Diff[n-min_nc+1,12]<- abs(result[n-min_nc+1,27]-result[n-min_nc,27])   #xu

        }
      }
    }

  }

  nc.Scott<-indice.Scott<-0
  if(is.na(match(13,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # SCOTT - The maximum difference between hierarchy levels was used to suggest the correct number of partitions.
    nc.Scott <-nb[which((Diff[,2])==max((Diff[,2]),na.rm=T))]
    indice.Scott <- max(Diff[,2],na.rm = TRUE)
    best.nc<-nc.Scott
  }

  nc.Marriot<-indice.Marriot<-0
  if(is.na(match(14,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # MARRIOT - The maximum difference between successive levels was used to determine the best partition level.
    nc.Marriot <- nb[which((Diff[,3])==max((Diff[,3]),na.rm=TRUE))]
    indice.Marriot <- max(Diff[,3],na.rm = TRUE)
    best.nc<-nc.Marriot
  }

  nc.TrCovW<-indice.TrCovW<-0
  if(is.na(match(15,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # TRCOVW - To determine the number of clusters in the data, maximum difference scores were used.
    nc.TrCovW <- nb[which((Diff[,4])==max((Diff[,4]),na.rm=T))]
    indice.TrCovW <- max(Diff[,4],na.rm = TRUE)
    best.nc<-nc.TrCovW
  }


  nc.TraceW<-indice.TraceW<-0
  if(is.na(match(16,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # TRACEW - To determine the number of clusters in the data, maximum difference scores were used.
    nc.TraceW <- nb[which((Diff[,5])==max((Diff[,5]),na.rm=T))]
    indice.TraceW <- max(Diff[,5],na.rm = TRUE)
    best.nc<-nc.TraceW
  }

  nc.Friedman<-indice.Friedman<-0
  if(is.na(match(17,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # FRIEDMAN - The maximum difference in values of trace W-1B criterion was used to indicate the optimal number of clusters.
    nc.Friedman <- nb[which((Diff[,6])==max((Diff[,6]),na.rm=T))]
    indice.Friedman <- max(Diff[,6],na.rm = TRUE)
    best.nc<-nc.Friedman
  }

  nc.Rubin<-indice.Rubin<-0
  if(is.na(match(18,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # RUBIN - The difference between levels was used.
    nc.Rubin <- nb[which((Diff[,7])==min((Diff[,7]),na.rm=T))]
    indice.Rubin <- min(Diff[,7],na.rm = TRUE)
    best.nc<-nc.Rubin
  }

  nc.Ball<-indice.Ball<-0
  if(is.na(match(19,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # BALL - The largest difference between levels was used to indicate the optimal solution.
    nc.Ball <- nb[which((Diff[,9])==max((Diff[,9]),na.rm=T))]
    indice.Ball <- max(Diff[,9],na.rm = TRUE)
    best.nc<-nc.Ball
  }


  nc.Hartigan<-indice.Hartigan<-0
  if(is.na(match(25,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # HARTIGAN - The maximum differences between hierarchy levels were taken as indicating the correct number of clusters in the data.
    nc.Hartigan <- nb[which((Diff[,11])==max((Diff[,11]),na.rm=T))]
    indice.Hartigan <- max(Diff[,11],na.rm = TRUE)
    best.nc<-nc.Hartigan
  }

  nc.xu<-indice.xu<-0
  if(is.na(match(27,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    #XU - The maximum value of the second differences is taken as the proposed number of clusters
    nc.xu <- nb[which((Diff[,12])==max((Diff[,12]),na.rm=T))]
    indice.xu<- max(Diff[,12],na.rm = TRUE)
    best.nc<-indice.xu
  }

  nc.SDbw<-indice.SDbw<-0
  if(is.na(match(20,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # SDBW - The minimum value of the index is used to indicate the optimal number of clusters.
    nc.SDbw <- as.vector(nb[which.min(result[,20])])
    indice.SDbw<- min(result[,20],na.rm = TRUE)
    best.nc<-nc.SDbw
  }

  nc.sv <-indice.sv<-0
  if(is.na(match(23,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # SV - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.sv <- as.vector(nb[which.max(result[,23])])
    indice.sv<- max(result[,23],na.rm = TRUE)
    best.nc<-nc.sv

  }

  nc.xiebeni <-indice.xiebeni<-0
  if(is.na(match(24,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    #XIE-BENI - The minimum value of the index is used to indicate the optimal number of clusters.
    nc.xiebeni <- as.vector(nb[which.min(result[,24])])
    indice.xiebeni<- min(result[,24],na.rm = TRUE)
    best.nc<-nc.xiebeni

  }

  nc.ssi<-indice.ssi<-0
  if(is.na(match(26,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # SSI - The maximum value of the index is used to indicate the optimal number of clusters.
    nc.ssi <- as.vector(nb[which.max(result[,26])])
    indice.ssi<- max(result[,26],na.rm = TRUE)
    best.nc<-nc.ssi
  }

  nc.rayturi<-indice.rayturi<-0
  if(is.na(match(28,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # RAY-TURI - The minimum  value of the index is used to indicate the optimal number of clusters.
    nc.rayturi <- as.vector(nb[which.min(result[,28])])
    indice.rayturi<- min(result[,28],na.rm = TRUE)
    best.nc<-nc.rayturi
  }


  nc.Pbm<-indice.Pbm<-0
  if(is.na(match(29,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # PBM -The maximum value of the index is used to indicate the optimal number of clusters
    nc.Pbm <-as.vector(nb[which.max(result[,29])])
    indice.Pbm <- max(result[,29],na.rm = TRUE)
    best.nc<-nc.Pbm
  }

  nc.Ban<-indice.Ban<-0
  if(is.na(match(30,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    # BANFELD - The minimum value of the index is used to indicate the optimal number of clusters.
    nc.Ban <-as.vector(nb[which.min(result[,30])])
    indice.Ban <- min(result[,30],na.rm = TRUE)
    best.nc<-nc.Ban
  }

  nc.Hubert<-indice.Hubert<-0
  if(is.na(match(22,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {
    nc.Hubert  <- 0.00
    indice.Hubert  <- 0.00
    options(graphics.record = TRUE)
    plot(nb,result[,22], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert Statistic values")))
    plot(Diff[,1],Diff[,10], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert statistic second differences")))
    cat(paste ("*** : The Hubert index is a graphical method of determining the number of clusters.
               In the plot of Hubert index, we seek a significant knee that corresponds to a
               significant increase of the value of the measure i.e the significant peak in Hubert
               index second differences plot.", "\n", "\n"))
  }

  nc.Dindex<-indice.Dindex<-0
  if(is.na(match(21,indice))==FALSE||is.na(match(31,indice))==FALSE)
  {

    nc.Dindex <- 0.00
    indice.Dindex<- 0.00
    plot(nb,result[,21], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Dindex Values")))
    plot(Diff[,1],Diff[,8], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Second differences Dindex Values")))
    cat(paste ("*** : The D index is a graphical method of determining the number of clusters.
               In the plot of D index, we seek a significant knee (the significant peak in Dindex
               second differences plot) that corresponds to a significant increase of the value of
               the measure.", "\n", "\n"))
  }


  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#                                                                                                                      #
  #                                                                                                                      #
  #                                            Displaying results                                                        #
  #                                                                                                                      #
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

  if((is.na(match(1,indice))==FALSE) || (is.na(match(2,indice))==FALSE) ||(is.na(match(3,indice))==FALSE)||
     (is.na(match(4,indice))==FALSE) || (is.na(match(5,indice))==FALSE)||(is.na(match(6,indice))==FALSE)||
     (is.na(match(7,indice))==FALSE) ||(is.na(match(8,indice))==FALSE)||(is.na(match(9,indice))==FALSE)||
     (is.na(match(10,indice))==FALSE)||(is.na(match(11,indice))==FALSE)||(is.na(match(12,indice))==FALSE)||
     (is.na(match(13,indice))==FALSE)||(is.na(match(14,indice))==FALSE)||(is.na(match(15,indice))==FALSE)||
     (is.na(match(16,indice))==FALSE)||(is.na(match(17,indice))==FALSE)||(is.na(match(18,indice))==FALSE)||
     (is.na(match(19,indice))==FALSE) ||(is.na(match(20,indice))==FALSE)||(is.na(match(23,indice))==FALSE)||
     (is.na(match(24,indice))==FALSE)||(is.na(match(25,indice))==FALSE)||(is.na(match(26,indice))==FALSE)||
     (is.na(match(27,indice))==FALSE)||(is.na(match(28,indice))==FALSE)||(is.na(match(29,indice))==FALSE)||
     (is.na(match(30,indice))==FALSE)||(is.na(match(31,indice))==FALSE) )
  {

    results <- c(nc.DB,indice.DB,nc.Dunn,indice.Dunn,nc.Silhouette,indice.Silhouette,nc.ptbiserial,
                 indice.ptbiserial,nc.CH,indice.CH,nc.cindex,indice.cindex,nc.Ratkowsky,indice.Ratkowsky,
                 nc.McClain,indice.McClain,nc.Gamma,indice.Gamma,nc.Gplus,indice.Gplus,nc.Tau,indice.Tau,
                 nc.CCC,indice.CCC,nc.Scott,indice.Scott,nc.Marriot,indice.Marriot,nc.TrCovW,indice.TrCovW,
                 nc.TraceW,indice.TraceW,nc.Friedman,indice.Friedman,nc.Rubin,indice.Rubin,nc.Ball,indice.Ball,
                 nc.SDbw,indice.SDbw,nc.Dindex,indice.Dindex,nc.Hubert,indice.Hubert,nc.sv,indice.sv,nc.xiebeni,indice.xiebeni,
                 nc.Hartigan,indice.Hartigan,nc.ssi,indice.ssi,nc.xu,indice.xu,nc.rayturi,indice.rayturi,nc.Pbm,indice.Pbm,nc.Ban,indice.Ban)


    result <-round(result, digits=4)
    results1<-matrix(c(results),nrow=2,ncol=30)
    resultats <- matrix(c(results),nrow=2,ncol=30,dimnames = list(c("Number_clusters","Value_Index"),
                                                                  c("db","dunn","silhouette","ptbiserial","ch","cindex","ratkowsky","mcclain","gamma",
                                                                    "gplus","tau","ccc","scott","marriot","trcovw","tracew","friedman","rubin","ball",
                                                                    "sdbw","dindex","hubert","sv","xie-beni","hartigan","ssi","xu","rayturi","pbm","banfeld")))


    resultats<-round(resultats, digits=3)

  }

  ######################## Summary results #####################################
  if (is.na(match(31,indice))==FALSE)
  {
    cat("*******************************************************************", "\n")
    cat("* Among all indices:                                               ", "\n")
    x <- unique(results1[1,])
    x <- sort(x)
    x<- x[-1]
    c=0
    for(i in 1:length(x))
    {
      vect<-which(x[i]== results1[1,])
      if(length(vect)>0)
        cat("*",length(vect), "proposed", x[i],"as the best number of clusters", "\n")

      if(c<length(vect))
      {
        j=x[i]
        c<-length(vect)
      }
    }

    pos<-which(nb==j)
    if(length(pos)>1){
      pos <-pos[1]
    }
    niveau<-max_nc-pos+1
    cat("\n","                  ***** Conclusion *****                           ", "\n", "\n")
    cat("* According to the majority rule, the best number of clusters is ",j , "\n")
    cat("and the best layer is ",paste(niveau,niveau,sep="x"), "\n", "\n")
    cat("*******************************************************************", "\n")
    jj <-j
    cl <-lis1[[1]]
    par <- nrow(dataa)
    partition <- rep(0,par)
    x <-list()
    b <- list()
    num <- ne[[pos]]
    for(i in 1:nb[[pos]])
    {
      y <- which(lis1[[pos]]==num[i])
      if ((pos-1 ==1) || (pos-1 ==0)){
        b[[i]] <- y
      }

      if(pos-1 > 1)
      {
        for(j in (pos-1):2)
        {
          x <-list()
          l <- 1
          for(k in 1:length(y))
          {
            x[[l]] <- which(lis1[[j]]==y[k])
            l <- l+1
          }
          y <- unlist(x)
        }
        b[[i]] <- y
      }

      for(j in 1:length(y))
      {
        for(k in 1:par)
        {
          if (cl[k]==y[j])
            partition[k]<- i

        }
      }
    }

    ###############################################################
    grid = somgrid(max_nc,max_nc,topo)
    plot(grid)
    symbols(grid$pts[, 1], grid$pts[, 2],circles = rep(0.5, nrow(grid$pts)),
            inches = FALSE, add = TRUE, bg = "white")
    main <- "Mapping plot"
    title.y <- max(grid$pts[,2]) + 1.2
    text(mean(range(grid$pts[,1])),title.y,main)


    if(pos-1 ==0){
      for (i in 1:length(num))
      {
        n <-num[i]
        lc <- length(which(cl==n))
        b <-paste(i,lc,sep=":")
        text(grid$pts[n,1],grid$pts[n, 2],b)
      }
    }
    else{
      for(i in 1:jj)
      {
        x <-b[[i]]
        lc<- c()
        n <-1
        for(j1 in 1:length(x))
        {
          if (is.na(match(x[j1],cl))==FALSE)
          {
            lc[n]<- length(which(cl==x[j1]))
            b1 <-paste(i,lc[n],sep=":")
            text(grid$pts[x[j1],1],grid$pts[x[j1], 2],b1)
            n <-n+1
          }

        }
      }


    }
  }


  if (is.na(match(31,indice))==FALSE)
  {
    results.final <- list(All.index.by.layer=result,Best.nc=resultats,Best.partition=partition)
  }

  if((is.na(match(22,indice))==FALSE) || (is.na(match(21,indice))==FALSE))

    results.final <- list(All.index.by.layer=result[,indice])


  if((is.na(match(1,indice))==FALSE)||(is.na(match(2,indice))==FALSE) ||(is.na(match(3,indice))==FALSE)||
     (is.na(match(4,indice))==FALSE)||(is.na(match(5,indice))==FALSE)||(is.na(match(6,indice))==FALSE)||
     (is.na(match(7,indice))==FALSE)||(is.na(match(8,indice))==FALSE)||(is.na(match(9,indice))==FALSE)||
     (is.na(match(10,indice))==FALSE)||(is.na(match(11,indice))==FALSE)||(is.na(match(12,indice))==FALSE)||
     (is.na(match(13,indice))==FALSE)||(is.na(match(14,indice))==FALSE)||(is.na(match(15,indice))==FALSE)||
     (is.na(match(16,indice))==FALSE)||(is.na(match(17,indice))==FALSE)||(is.na(match(18,indice))==FALSE)||
     (is.na(match(19,indice))==FALSE)||(is.na(match(20,indice))==FALSE)||(is.na(match(23,indice))==FALSE)||
     (is.na(match(24,indice))==FALSE)||(is.na(match(25,indice))==FALSE)||(is.na(match(26,indice))==FALSE)||
     (is.na(match(27,indice))==FALSE)||(is.na(match(28,indice))==FALSE)||(is.na(match(29,indice))==FALSE)||
     (is.na(match(30,indice))==FALSE))


    results.final <- list(All.index.by.layer=result[,indice],Best.nc=resultats[,indice])


  return(results.final)

  }



