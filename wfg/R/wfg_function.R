
wfg <- function(net, attr=NULL, under.sample=FALSE, prioritize=FALSE)
{
  result <- list()

  n <- dim(get.adjacency(net))[1]
  n_edge <- sum(get.adjacency(net))/2

  A <- get.adjacency(net)
  A <- as.matrix(A)
  eff <- matrix(1, nrow=n, ncol=n)

  if (!is.null(attributes))
  {
    # vector of edges
    y <- A[upper.tri(A)]

    df <- data.frame(y)

    n_att <- dim(attr)[2]
    model <- list()

    for (a in 1:n_att)
    {
      tp <- typeof(attr[,a])
      x <- c()
      from.mat <- replicate(n, attr[,a])
      to.mat <- t(from.mat)
      fromv <- from.mat[upper.tri(from.mat)]
      tov <- to.mat[upper.tri(to.mat)]

      if (tp=="double") # numeric
      {
        x <- abs(fromv-tov)
      }
      else # otherwise
      {
        x <- 1-(fromv==tov)*1
      }
      df <- data.frame(df,x)
      df.log <- df
    }

    if (under.sample)
    {
      # sample from data for logistric regression
      df.y1 <- df[df$y==1,]
      df.y0 <- df[df$y==0,]

      # number of sampling
      n.resamp <- n_edge
      resamp0 <- sample(seq(1:  dim(df.y0)[1]  ), n.resamp, replace = FALSE)
      df.y0.resample <- df.y0[resamp0,]

      df.log <- rbind(df.y1, df.y0.resample)
    }

    model <- glm(y ~ ., data = df.log, family = "binomial")

    beta <- model$coefficients
    predictor <- as.matrix(df[,2:(n_att+1)])
    lp <- predictor %*% beta[2:(n_att+1)]
    elp <- exp(lp)

    eff[upper.tri(eff)] <- elp
    eff <- eff+t(eff)
  }
  A.weighted <- A*eff
  diag(A.weighted) <- 0
  net2 <- graph.adjacency(A.weighted, mode="undirected", weighted=TRUE)
  cmt2 <- fastgreedy.community(net2)
  memb <- membership(cmt2)

  if (!is.null(attributes))
  {
    names(beta) <- c("Intercept", colnames(attr))
    result$beta <- beta
  }

  if (prioritize){

    if (!is.null(attributes))
    {
      from.mat <- replicate(n,memb)
      to.mat <- t(from.mat)

      fromv <- from.mat[upper.tri(from.mat)]
      tov <- to.mat[upper.tri(to.mat)]

      betas0 <- c()
      betas1 <- c()

      memb.label <- unique(memb)

      beta.matrix <- matrix(0, nrow=length(memb.label), ncol=n_att)

      for (k in 1:n_att)
      {
        d <- dist(attr[,k], method = "manhattan")
        d.mat <- as.matrix(d)
        x.pri <- d.mat[upper.tri(d.mat)]
        string <- data.frame(x.pri, fromv, tov)
        i <- 1
        for (label in memb.label)
        {
          within <- string[string$fromv==label & string$tov==label, ]
          reachout <- string[string$fromv!=label | string$tov!=label, ]

          n_within <- dim(within)[1]
          index <- sample(1:dim(reachout)[1], n_within)
          y.pri <- c(rep(1,n_within),rep(0,length(index)))

          df.intermediate <- rbind(within, reachout[index,])
          df.log <- cbind(df.intermediate,y.pri)

          model <- glm(y.pri~x.pri, data = df.log, family = "binomial")
          beta.pri <- model$coefficients

          betas0[i] <- beta.pri[1]
          betas1[i] <- beta.pri[2]

          i <- i+1
        }

        beta.matrix[,k] <- betas1
      }
    }
    row.names(beta.matrix) <- memb.label
    colnames(beta.matrix) <- colnames(attr)
    result$beta.matrix <- beta.matrix
  }

  result$memb <- memb
  return(result)
}
