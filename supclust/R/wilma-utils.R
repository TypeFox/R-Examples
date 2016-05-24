## The function implementing the cleaning stage for Wilma
back.search <- function(genes, x, y, verbose = FALSE)
{
    reduced <- TRUE
    while (reduced && (N <- length(genes)) >= 2)
    {
	reduced <- FALSE
	if (score(mx <- rowMeans(x[, genes]), y) == 0)
	{
            margins <- sapply(seq(along=genes), function(i)
                              margin(rowMeans(x[, genes[-i], drop=FALSE]), y))
	    if (max(margins) > margin(mx, y)) {
                imax <- which.max(margins)
		if(verbose)
		    cat("Gen: ", genes[imax],
			"Score:	  0",
			"Margin: ", round(max(margins),3),"\n")
		genes	   <- genes[-imax]
		reduced	   <- TRUE
	    }
	}

	if ((sc <- score(mx <- rowMeans(x[, genes]), y)) > 0)
	{
	    old.score <- sc
            scores	<- numeric(N)
	    for (i in 1:length(genes))
		scores[i]  <- score(rowMeans(x[, genes[-i], drop=FALSE]), y)

            minsc <- min(scores)
	    indices <- which(scores == minsc)# all of them
            mgs <- sapply(indices, function(ii)
                          margin(rowMeans(x[, genes[- ii], drop=FALSE]), y))

	    if (minsc < old.score || max(mgs) > margin(mx, y)) {
                imax <- which.max(mgs)
		if(verbose)
		    cat("Gen: ", indices[imax],
			"Score: ", minsc,
			"Margin: ", round(mgs[imax],3),"\n")
		genes	   <- genes[-indices[imax]]
		reduced	   <- TRUE
              }
          }
      }
    genes
  }


## The margin function for Wilma
margin <- function(x, resp)
  {
    .C(R_margin,
       as.double(x[order(resp)]),
       as.integer(sum(resp==0)),
       as.integer(sum(resp==1)),
       re = double(1))$re
  }


## The score function for Wilma
score <- function(x, resp)
  {
    .C(R_score,
       as.double(x),
       as.integer(resp),
       as.integer(length(x)),
       re = double(1))$re
  }


## The function implementing the sign-flip for Wilma -- using score()
sign.flip  <- function(x, y)
  {
    y <- as.integer(y)
    O <- as.integer(0)
    if(any(y < O | y > 1L)) stop("'y' must be 0 or 1!")
    n <- length(y)
    scores <- apply(x, 2, score, y)
    n0 <- sum(y == O)
    middle <- n0 * (n - n0) / 2
    lrg  <- scores > middle
    x[,lrg] <- -x[,lrg]
    list(flipped.matrix = x, signs = ((-2)*lrg)+1)
  }


## Classification with Wilma's clusters - nearest neighbor rule
nnr <- function(xlearn, xtest, ylearn)
{
    as.numeric(knn(xlearn, xtest, ylearn))-1
}


## Classification with Wilma's clusters - diagonal linear discriminant analysis
dlda <- function (xlearn, xtest, ylearn)
  {
    ## Definition of variables
    n      <- nrow(xlearn)
    p      <- ncol(xlearn)
    nk     <- rep(0, max(ylearn) - min(ylearn) + 1)
    K      <- length(nk)
    m      <- matrix(0, K, p)
    v      <- matrix(0, K, p)
    disc   <- matrix(0, nrow(xtest), K)

    ## Computing mean and variances
    for (k in (1:K))
      {
        which   <- (ylearn == k + min(ylearn) - 1)
        nk[k]   <- sum(which)
        m[k, ]  <- apply(xlearn[which, , drop = FALSE], 2, mean)
        v[k,]   <- apply(xlearn[which, , drop = FALSE], 2, var)
      }

    ## Computing the pooled variance
    vp <- apply(v, 2, function(z) sum((nk - 1) * z)/(n - sum(nk!=0)))

    ## Computing the discriminant function
    for (k in (1:K))
      {
        disc[,k] <- apply(xtest, 1, function(z) sum((z-m[k,])^2*(1/vp)))
      }

    ## Prediction and output
    pred <- apply(disc, 1, function(z) (min(ylearn):max(ylearn))[order(z)[1]])
    pred
   }


## Classification with Wilma's clusters - logistic regression
logreg <- function(xlearn, xtest, ylearn)
  {
    xvalues    <- xlearn
    op         <- options(warn = -1)
    model      <- glm(ylearn~., data = data.frame(xvalues), family = binomial)
    xvalues    <- xtest
    predic     <- predict(model,new = data.frame(xvalues), type = "response")
    options(op)
    as.numeric(predic >= 0.5)
  }


## Classification with Wilma's clusters - aggregated trees
aggtrees <- function(xlearn, xtest, ylearn)
  {
    noc     <- ncol(xtest)
    predic  <- matrix(0, nrow(xtest), noc)
    for (j in 1:noc)
      {
        xvalues    <- xlearn
        model      <- rpart(ylearn ~ ., data = data.frame(xvalues))
        xvalues    <- xtest
        predic[,j] <- predict(model, newdata = data.frame(xvalues))
      }
    final <- rowSums(predic)
    for (j in 1:nrow(xtest))
      {
        if (final[j] == (noc/2)) final[j] <- final[j]+(predic[j,1]-0.5)
      }
    as.numeric(final > (round((noc-0.1)/2)))
  }

