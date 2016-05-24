obliclus <-
function (A, cluster=NULL, normalize=TRUE, N.random=10, fixed=FALSE,
                      N.cluster=NULL, ini.cluster=NULL, maxit=1000, alpha=1,
                      method="oblimin", geomin.par=0.01, oblimin.index=0,
                      fit.cr=1e-04)
{
  if (!is.matrix(A))
    stop("'A' must be matrix form.")

  if (is.null(N.cluster))
    N.cluster <- dim(A)[2]


  N.var <- dim(A)[1]; N.fac <- dim(A)[2]
  info <- list("N.var"=N.var, "N.fac"=N.fac, "N.random"=N.random, "fixed"=fixed,
               "N.cluster"=N.cluster, "alpha"=alpha, "method"=method,
               "geomin.par"=geomin.par, "oblimin.index"=oblimin.index,
               "fit.cr"=fit.cr)

  ##Kaiser's Normalization
  if (normalize) {
    H.norm <- diag(apply(A^2, 1, sum))^(1/2)
    A <- solve(H.norm) %*% A
  }

  ##Different N.random initial starts of cluster parameter
  value.criteria <- 0
  for (r in 1:info$N.random) {
    if (!fixed) {
      cluster <- random.cluster.C(c(N.var, N.cluster, 1, 1))
      while (any(cluster == 0))
        cluster <- random.cluster.C(c(N.var, N.cluster, 1, 1))
    }
    if(!is.null(ini.cluster)){
      cluster <- ini.cluster
    }

    conv <- FALSE

    ##Use GP algorithm and k-means algorithm
    for(i in 1:maxit){
      if(i != 1)
        criteria.old <- criteria.new

      if (i == 1) ##Initialize T
        T.n <- qr.Q(qr(matrix(rnorm(N.fac * N.fac), N.fac)))

      if (fixed)
        T <- T.n

      ##Update T
      if (!fixed) {
        for (r.j in 1:10) {
          if (r.j == 1) {
            T.n.j <- T.n
          }
          else {
            T.n.j <- qr.Q(qr(matrix(rnorm(N.fac * N.fac), N.fac)))
          }

          for (j in 1:maxit) {
            T <- T.n.j

            if (j != 1)
              criteria.old.j <- criteria.new.j

            ##Gradient of criteria at T
            G <- analy.deriv(A, T, cluster, info)

            ##Calculate optimal alpha
            alpha <- alpha.optim(A, T, G, cluster, info)

            ##Update T.j
            X <- T - alpha * G
            T.n.j <- X %*% sqrt(solve(diag(diag(t(X) %*% X))))
            criteria.new.j <- criteria(A, T.n.j, cluster, info)

            if (j != 1)
              if (criteria.new.j > criteria.old.j)
                break

            ##Check convergence
            if (j == 1) {
              criteria.old.j <- criteria.new.j
              s <- 1
            }
            else {
              s <- abs(criteria.old.j - criteria.new.j)
            }
            if (s < info$fit.cr) {
              conv <- TRUE
              T.j <- T.n.j
              break
            }
          }

          if (r.j == 1) {
            value.criteria.j <- criteria(A, T.j, cluster, info)
            T.out.j <- T.j
            conv.out.j <- conv
          }
          else {
            if (!is.numeric(alpha))
              next
            value.criteria.new.j <- criteria(A, T.j, cluster, info)
            if (value.criteria.j > value.criteria.new.j) {
              T.out.j <- T.j
              value.criteria.j <- value.criteria.new.j
              conv.out.j <- conv
            }
          }
        }
      }
      else { ##Case where cluster structure is fixed
        ##Gradient of criteria at T
        G <- analy.deriv(A, T, cluster, info)

        ##Calculate optimal alpha
        alpha <- alpha.optim(A, T, G, cluster, info)

        ##Update T.j
        X <- T - alpha * G
        T.out.j <- X %*% sqrt(solve(diag(diag(t(X) %*% X))))
        criteria.new.j <- criteria(A, T.out.j, cluster, info)
      }

      T.n <- T.out.j
      criteria.new <- criteria(A, T.n, cluster, info)

      ##Update cluster parameter
      if(!fixed){
        cluster.tmp <- cluster
        B <- A %*% solve(t(T.n))
        cluster <- kmeans(B^2, info$N.cluster, nstart=100, iter.max=100, algorithm="Lloyd")$cluster
        tmp <- criteria(A, T.n, cluster, info)

        if(tmp > criteria.new){
          cluster <- cluster.tmp
        }

        if(!is.numeric(cluster)){
          T <- T.n
          conv <- TRUE
          break
        }
      }

      criteria.new <- criteria(A, T.n, cluster, info)


      ##Check convergence
      if(i == 1){
        criteria.old <- criteria.new
        s <- 1
      }
      else {
        s <- abs(criteria.old - criteria.new)
      }
      if (s < info$fit.cr) {
        conv <- TRUE
        T <- T.n
        break
      }
    }

    if (r == 1) {
      value.criteria <- criteria(A, T, cluster, info)
      T.out <- T
      cluster.out <- cluster
      num.out <- 1
      conv.out <- conv
    }
    else {
      if (!is.numeric(alpha))
        next
      value.criteria.new <- criteria(A, T, cluster, info)
      if (value.criteria > value.criteria.new) {
        T.out <- T
        cluster.out <- cluster
        num.out <- r
        value.criteria <- value.criteria.new
        conv.out <- conv
      }
    }
  }

  T <- T.out
  A <- A %*% solve(t(T))
  A.out <- A
  E <- t(T) %*% T
  clust <- cluster.out

  ##Denormalization
  if (normalize) {
    A.out <- H.norm %*% A.out
  }

  return(list("A"=A.out, "T"=T, "E"=E, "conv"=conv.out, "cluster"=clust, "cr"=value.criteria))
}
