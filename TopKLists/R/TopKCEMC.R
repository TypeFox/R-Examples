###New set of functions for generate combined ranked list using CEMC
###taking into account difference spaces of input ranked lists
###Last modified: 04/28/10

`CEMC` <-
function(input,space=NULL,k=NULL,dm="k",kp=0.5,N=NULL, N1=NULL,rho=0.1,
         e1=0.1,e2=1,w=0.5,b=0,init.m="p",init.w=0, d.w=NULL,input.par=NULL,
         extra=0){  

  ##input: a list of several top K lists, may have different length
  ##space: underlying spaces for the lists
  ##k: desired length of combined list
  ##dm: distance measure, "s" for Spearman, "k" for kendall(p=0)
  ##kp: partial distance used in Kendall's tau distance measure
  ##N: number of samples generated in each iterate
  ##N1: number of samples retained after each iterate 
  ##rho: proportion of samples used to estimate new probability matrix
  ##e1: stop criteria w.r.t. p
  ##e2: stop criteria w.r.t. y
  ##w: weight of the new probabilities estimate each iterate
  ##b: parameter used in .blur function
  ##init.m: initialization method, see the function init.p for details
  ##init.w: probability matrix initialization:
  ##      (1-init.w) * uniform +  init.w * estimated from input lists
  ##d.w: weights for distances from different input lists
  ##input.par: input parameters in a dataframe

  if (missing(input))
    stop("You need to input the top-k lists to be aggregated")
  if (!is.null(input.par)) {
    for (p.n in names(input.par))
      assign(p.n,input.par[1,p.n])
  }
  
  time.start <- proc.time()

  topK.input <- input
  space.input <- space
  
  a <- length(input) # number of input top K lists
  k.a <- numeric(a) # lengths of input lists
  for (i in 1:a)
    k.a[i] <- length(input[[i]])
  item <- sort(unique(unlist(input))) # all items in the input top K lists
  n <- length(item)
  ##added to allow for k to be NULL
  if (is.null(k)) 
    k <- n
  else if (k>n) 
    k <- n
  if (extra > 0) {
    n <- n + extra
    k <- k + extra
  }
  ##recode ranks: ranked=1:k.a, unranked but in space=k.a+1, not in space=0 
  rank.a <- matrix(0,nrow=n,ncol=a)
  if (is.null(space)) { 
    for (i in 1:a) {
      input[[i]] <- match(input[[i]],item) # recode items from 1 to n
      rank.a[,i] <- 1+k.a[i]
      rank.a[input[[i]],i] <- 1:k.a[i]
    }
  }
  else {
    for (i in 1:a) {
      input[[i]] <- match(input[[i]],item) # recode items from 1 to n
      space[[i]] <- match(space[[i]],item) 
      rank.a[space[[i]],i] <- 1+k.a[i]
      rank.a[input[[i]],i] <- 1:k.a[i]
    }
  }
  
  ##set default values for N and N1
  if (is.null(N))
    N <- 50 * n

  if (is.null(N1))
    N1 <- round(0.1 * N)

  p <- init.p(input,n,k,init.m,init.w)
  p2 <- p
  y <- 0
  samp2 <- TopKSample.c(.blur(p,b),N1)[1:k,]
  iter <- 0
  y.count <- 0 # counter used in stopping rule
  
  if (is.null(d.w))
    d.w <- rep(1,a) # equal weight as default
  
  repeat {
    iter <- iter + 1
    samp <- cbind(samp2,TopKSample.c(.blur(p,b),N-N1)[1:k,])

 ##   if (iter == 10)
 ##     samp[,N] <- new.c
    
    rank.b <- matrix(k+1,nrow=n,ncol=N)
    rank.b[cbind(as.vector(samp),rep(1:N,each=k))] <- 1:k
    
    dist <- numeric(N)
    for (i in 1:a) {
      if (dm=="s")
        dist <- dist + Spearman(rank.a[,i],rank.b,k.a[i],k,n) * d.w[i]
      else if (dm=="k") 
        dist <- dist + Kendall2Lists.c(rank.a[,i],rank.b,k.a[i],k,n,kp) * d.w[i]
      else stop("Invalid distance measure")
    }
    y2 <- sort(dist)[round(N*rho)]
    samp2 <- samp[,order(dist)[1:N1]]
    samp3 <- samp[,dist<=y2]
    n.samp3 <- dim(samp3)[2]
#    print(samp2[,1])
#    print(min(dist))
    for (i in 1:k)
      for (j in 1:n)
        p2[j,i] <- sum(samp3[i,]==j)/n.samp3
    if (k < n) 
      p2[,k+1] <- 1 - apply(p2[,1:k],1,sum)
    if (abs(y-y2) < e2) y.count <- y.count + 1
    else y.count <- 0
    y <- y2
    p <- p*(1-w)+p2*w
    if (y.count >= 5) break
  }

  result <- samp[,order(dist)[1]]
  rank.result <- rep(k+1,n)
  rank.result[result] <- 1:k
  if (k < n) 
    dimnames(p) <- list(c(item,rep(0,extra)),1:(k+1))
  else
    dimnames(p) <- list(c(item,rep(0,extra)),1:k)

  diff.p <- mean(abs(p[result,1:k] - diag(k)))
  #difference between sorted p and the identity matrix

  uc <- mean(1-p[cbind(result,1:k)])    
  #mean uncertainty for the final list
  
  result.ori <- item[result[result <= n-extra]] #result in original coding
  
  dist.s <- 0
  dist.k <- 0
  for (i in 1:a) {
    dist.s <- dist.s + Spearman(rank.a[,i],rank.result,k.a[i],k,n) * d.w[i]
    dist.k <- dist.k + Kendall2Lists.c(rank.a[,i],rank.result,k.a[i],k,n,kp) *d.w[i]
  }
  
  time.end <- proc.time()
  time.use <- sum((time.end-time.start)[-3])
  
  list(TopK=result.ori,ProbMatrix=p,
       #output.other=data.frame(diff.p=diff.p, uc=uc,dist.s=dist.s,
       #dist.k=dist.k, iter=iter, time=time.use),
       #input.list=topK.input,input.space=space.input,d.w=d.w,
       input.par=data.frame(k=k,dm=I(dm),kp=kp,N=N,N1=N1,extra=extra,
       rho=rho,e1=e1,e2=e2,w=w,b=b,init.m=I(init.m),init.w=init.w))
}

`TopKSample` <-
function(p,N) {
  #Sampler to generate N top K lists according to p
  #p: matrix of dimension n*k, n is the number of items
  #N: the number of samples

  k <- dim(p)[2]
  n <- dim(p)[1]
  item <- 1:n
  samp <- matrix(0,nrow=k,ncol=N)

  i <- 1
  fail <- FALSE
  repeat {
    samp[1,i] <- which(rmultinom(1,1,p[,1])==1)
    for (j in 2:k) {
      remain <- item[-samp[1:(j-1),i]]
      if (sum(p[remain,j]) == 0) {
        fail <- TRUE
        break
      }
      samp[j,i] <- remain[which(rmultinom(1,1,p[remain,j])==1)]
    }
    if (!fail) i <- i + 1
    if (i == N+1) break
  }
  samp
}

`TopKSample.c` <-
function(p,N) {

  #this version calls a C subroutine for faster sampling

  k <- dim(p)[2]
  n <- dim(p)[1]
  samp <- matrix(0,nrow=k,ncol=N)
  seed <- round(runif(1,10000,3000000))
  
  samp[] <- .C("topksamplec",as.double(p),as.integer(k),as.integer(n),
               as.integer(N),samp=as.integer(samp),as.integer(seed))$samp

  samp

}

`.blur` <-
function(p,b) {

  #p: vector of probabilites
  #b: exchange proportions, see below

  n <- dim(p)[1]
  k <- dim(p)[2]
  p2 <- matrix(0,nrow=n,ncol=k)

  p2[,1] <- p[,1]*(1-b) + p[,2]*b
  
  for (i in 2:(k-1))
    p2[,i] <- p[,i-1]*b+p[,i]*(1-2*b)+p[,i+1]*b

  p2[,k] <- p[,k]*(1-b) + p[,k-1]*b
  p2
}

`.cc.rank` <-
function(input.list) {

  #calculate composite rank
  #input.list: a list of topk lists, may have different length

  n.list <- length(input.list) #number of lists
  nn.list <- unlist(lapply(input.list, length)) #lengths of individual lists

  item <- unique(unlist(input.list))
  n.item <- length(item)
  
  score <- matrix(0,nrow=n.item, ncol=n.list)

  for (i in 1:n.list) {
    list.code <- match(input.list[[i]],item) #coded list
    score[,i] <- nn.list[i]+1 #scores for items not in the list
    score[list.code,i] <- 1:nn.list[i] #scores for item in the list
    score[,i] <- score[,i]/nn.list[i] #weighted score
  }

  c.score <- apply(score,1,sum) #composite scores
  item[order(c.score)] #order w.r.t the composite scores

}

`init.p` <-
function(topK,n,k,init.m="p",init.w=0) {

  #topK: a list of input lists, with items coded from 1 to n
  #n: total number of items
  #k: length of target list
  #init.m: initialization method, "p" point mass, "s" smooth
  #        "cp" point mass using composite ranks
  #        "cs" smooth using composite ranks
  #init.w: initialization weight
 
  if (k < n) {
    p.u <- matrix(1/n,nrow=n,ncol=k+1)
    p.u[,k+1] <- 1-k/n
  }
  else {
    p.u <- matrix(1/n,nrow=n,ncol=k)
  }
  
  a <- length(topK)
  
  if (init.m %in% c("p","s")) {
    p.e <- matrix(0,nrow=n,ncol=n)
    for (i in 1:a) {
      k.a <- length(topK[[i]])
      p.e[cbind(topK[[i]],1:k.a)] <- 1/a + p.e[cbind(topK[[i]],1:k.a)]
      if (k.a < n) {
	p.e[(1:n)[-match(topK[[i]],1:n)],(k.a+1):n] <-
	  1/((n-k.a)*a) + p.e[(1:n)[-match(topK[[i]],1:n)],(k.a+1):n]
      }
    }
#    if ((k+1)<n)
#      p.e <- cbind(p.e[,1:k],apply(p.e[,(k+1):n],1,sum))
    }
  else if (init.m %in% c("cp","cs")) {
    p.e <- matrix(0,nrow=n,ncol=n)
    c.list <- .cc.rank(topK)
    p.e[cbind(c.list,1:n)] <- 1
#    p.e[,(k+1):n] <- (1-apply(p.e[,1:k],1,sum))/(n-k)
    }
  else stop("Invalid initialization method")
  
  if (init.m %in% c("s","cs")) {
      dist <- 0
      for (i in 1:(a-1)) {
        for (j in (i+1):a) {
          dist <- dist + Spearman(topK[[i]],topK[[j]],n)
        }
      }
      dist.avg <- dist/(a*(a-1)/2)/n #average rank difference
      normal.c <- sum(dnorm(1:n,mean=(n+1)/2,sd=dist.avg/2))
                                     #normalizing constant
      weight <- matrix(0,nrow=n,ncol=n)
      for (i in 1:n) {
        weight[,i] <- dnorm(1:n,mean=i,sd=dist.avg/2)/normal.c
        weight[i,i] <- weight[i,i]+(1-sum(weight[,i])) #make weights sum to 1
      }
      p.e <- p.e %*% weight
     }
  if (k < n) 
    p.e <- cbind(p.e[,1:k],1-apply(p.e[,1:k],1,sum))
  p.u*(1-init.w) + p.e*init.w
}

`Kendall2Lists` <-
function(rank.a,rank.b,k.a,k.b,n,p=0) {
 
  ##Kendall's tau distance between top K lists
  ##n: total number of objects, numbered from 1 to n
  ##rank.a: a single top k list
  ##rank.b: a vector of matrix form of top k list(s) to be compared to list a
  ##k.a: value of k for rank.a
  ##k.b: value of k for rank.b
  ##p: distance added for tied pair (potential problem when p != 0)
 
  if (is.vector(rank.b)) {
    n.b <- 1
  }
  else {
    n.b <- ncol(rank.b) # number of lists in b
  }
 
  dist <- numeric(n.b)
  d.a <- rank.a - rep(rank.a,each=n) # compare all items to each other
  for (i in 1:n.b) {
     if (n.b == 1) d.b <- rank.b - rep(rank.b,each=n)
    else d.b <- rank.b[,i] - rep(rank.b[,i],each=n)
    count <- table(c(sign(d.a)*sign(d.b),-1,0,1)) #make sure all 3 values exist
    dist[i] <- (count["-1"]-1)*1/2 + (count["1"]-1)*0/2 + (count["0"]-n-1)*p/2
  }  
  dist
}

## `kendall.c` <-
## function(rank.a,rank.b,k.a,k.b,n,p=0) {

##   #version that calls a C subroutine to speed up calculation
##   #Kendall's tau distance between top K lists
##   ##rank.a: a single top k list
##   ##rank.b: a vector of matrix form of top k list(s) to be compared to list a
##   ##k.a: value of k for rank.a
##   ##k.b: value of k for rank.b
##   ##p: distance added for tied pair (potential problem when p != 0)
##   ##n: total number of objects, numbered from 1 to n
  
##   if (is.vector(rank.b)) {
##     n.b <- 1
##   }
##   else {
##     n.b <- ncol(rank.b) # number of lists in b
##   }

##   dist <- numeric(n.b)
##   .C("kendallc",as.integer(rank.a),as.integer(rank.b),as.integer(n),
##      as.integer(n.b),as.double(p),dist=as.double(dist))$dist
## }

`Kendall2Lists.c` <-
function(rank.a,rank.b,k.a,k.b,n,p=0) {

  ##Kendall's tau distance between top K lists with different underlying spaces
  ##rank.a: a single top k list
  ##rank.b: a vector of matrix form of top k list(s) to be compared to list a
  ##k.a: value of k for rank.a
  ##k.b: value of k for rank.b
  ##p: distance added for tied pair (potential problem when p != 0)
  ##n: total number of objects, numbered from 1 to n
  
  if (is.vector(rank.b)) 
    n.b <- 1
  else
    n.b <- ncol(rank.b) # number of lists in b

  dist <- numeric(n.b)
  .C("kendall2c",as.integer(rank.a),as.integer(rank.b),as.integer(n),
     as.integer(n.b),as.integer(k.a),as.integer(k.b),
     as.double(p),dist=as.double(dist))$dist
}

`Spearman` <-
function(rank.a,rank.b,k.a,k.b,n) {
  
  ##Spearman's footrule distance between two top K lists
  ##rank.a: a single top k list
  ##rank.b: a vector of matrix form of top k list(s) to be compared to list a
  ##k.a: value of k for rank.a
  ##k.b: value of k for rank.b
  ##n: total number of objects, numbered from 1 to n
  
  if (is.vector(rank.b)) {
    n.b <- 1
  }
  else {
    n.b <- ncol(rank.b) # number of lists in b
  }
 
  if (n.b == 1)
    d <- sum(abs(rank.a-rank.b)*((rank.a<=k.a) | (rank.b<=k.b)))
    # for each list in b, only calculate the distance using the items
    # exist in list a and/or the list in b
  else
    d <- apply(abs(rank.a-rank.b)*((rank.a<=k.a) | (rank.b<=k.b)),2,sum)
  d
}

