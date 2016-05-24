
# ===================================================
# Creating a SMOTE training sample for classification problems
# L. Torgo, Feb 2010, Nov 2014
# ---------------------------------------------------
smote <- function(form,data, perc.over=2,k=5,perc.under=2)
  
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # per.over is the number of new cases (smoted cases) generated
  #          for each rare case. If perc.over < 1 a single case
  #          is generated uniquely for a randomly selected perc.over
  #          of the rare cases
  # k is the number of neighbours to consider as the pool from where
  #   the new examples are generated
  # perc.under is the number of "normal" cases that are randomly
  #            selected for each smoted case

{

  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  minCl <- names(which.min(table(data[,tgt])))
  
  # get the cases of the minority class
  minExs <- which(data[,tgt] == minCl)

  # generate synthetic cases from these minExs
  if (tgt < ncol(data)) {
      orig.order <- colnames(data)
      cols <- 1:ncol(data)
      cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
      data <-  data[,cols]
  }
  newExs <- smote.exs(data[minExs,],ncol(data),perc.over,k)
  if (tgt < ncol(data)) {
      newExs <- newExs[,cols]
      data <- data[,cols]
  }
  
  # get the undersample of the "majority class" examples
  selMaj <- sample((1:NROW(data))[-minExs],
                   as.integer(perc.under*nrow(newExs)),
                   replace=T)

  # the final data set (the undersample+the rare cases+the smoted exs)
  newdataset <- rbind(data[selMaj,],data[minExs,],newExs)

  if (tgt < ncol(data)) newdataset <- newdataset[,orig.order]
  newdataset

}



# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
# L. Torgo, Feb 2010
# ---------------------------------------------------
smote.exs <- function(data,tgt,N,k)
  # INPUTS:
  # data are the rare cases (the minority "class" cases)
  # tgt is the name of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours to use for the generation
  # OUTPUTS:
  # The result of the function is a N*T set of generated
  # examples with rare values on the target
{
  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2]-1)
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]
  
  if (N < 1) { # only a percentage of the T cases will be SMOTEd
    nT <- NROW(T)
    idx <- sample(1:nT,as.integer(N*nT))
    T <- T[idx,]
    N <- 1
  }

  p <- dim(T)[2]
  nT <- dim(T)[1]

  ranges <- apply(T,2,max)-apply(T,2,min)
  
  nexs <-  as.integer(N) # this is the number of artificial exs generated
                                        # for each member of T
  new <- matrix(nrow=nexs*nT,ncol=p)    # the new cases

  for(i in 1:nT) {

    # the k NNs of case T[i,]
    xd <- scale(T,T[i,],ranges)
    for(a in nomatr) xd[,a] <- xd[,a]==0
    dd <- drop(xd^2 %*% rep(1, ncol(xd)))
    kNNs <- order(dd)[2:(k+1)]

    for(n in 1:nexs) {
      # select randomly one of the k NNs
      neig <- sample(1:k,1)

      ex <- vector(length=ncol(T))

      # the attribute values of the generated case
      difs <- T[kNNs[neig],]-T[i,]
      new[(i-1)*nexs+n,] <- T[i,]+runif(1)*difs
      for(a in nomatr)
        new[(i-1)*nexs+n,a] <- c(T[kNNs[neig],a],T[i,a])[1+round(runif(1),0)]

    }
  }
  newCases <- data.frame(new)
  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))

  newCases[,tgt] <- factor(rep(data[1,tgt],nrow(newCases)),levels=levels(data[,tgt]))
  colnames(newCases) <- colnames(data)
  newCases
}

 
