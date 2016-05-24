#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

colAUC = function (X, y, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
{   
  # input:
  #   X - 2D matrix (of feature columns and samples rows)
  #   y - 1D array identifying which class each sample belongs to (csource(lass numbers start from 1)
  
  #=======================================
  # make sure inputs are in correct format
  #=======================================
  y    = as.factor(y)
  X    = as.matrix(X)               # make sure inputs are in correct format
  alg  = match.arg(alg)
  if (nrow(X)==1) X = t(X)
  if (plotROC) alg = "ROC"
  
  #=======================================
  # Prepare for calculations & error-check
  #=======================================
  nR   = nrow(X)                    # number of samples
  nC   = ncol(X)                    # get dimentions of the data set
  nY   = table(y)                   # number of elements per label
  uL   = as.factor(rownames(nY))    # find all the classes among the labels 
  nL   = length(nY)                 # number of unique classes
  if (nL<=1) 
    stop("colAUC: List of labels 'y' have to contain at least 2 class labels.")
  if (!is.numeric(X)) stop("colAUC: 'X' must be numeric")
  if (nR!=length(y))  stop("colAUC: length(y) and nrow(X) must be the same")
  L    = matrix(rep(uL,each=nR),nR,nL)  # store vector L as row vector and copy it into nR rows
  per  = combs(1:nL,2)                  # find all possible pairs of L columns
  nP   = nrow(per)                      # how many possible pairs were found?
  Auc  = matrix(0.5,nP,nC)              # initialize array to store results
  rownames(Auc) = paste(uL[per[,1]]," vs. ",uL[per[,2]], sep="")
  colnames(Auc) = colnames(X)
  
  #=======================================
  # prepare the plot, if needed
  #=======================================
  if (plotROC) {                        # initialize the plot
    plot(c(0,1), c(0,1), type='n', xaxs="i", yaxs="i",
     xlab="probability of false alarm", sub="(1-Specificity)", 
     ylab="probability of detection\n(Sensitivity)")
    title("ROC Curves")
    abline(h=0:10/10, v=0:10/10, col = "lightgray") # grid on
    if (nC*nP<20) { # if too many curves than skip the labels
      S = colnames(Auc)
      if (is.null(S)) S=paste('col',1:nC);
      if (nP>1) S = paste(rep(S,each=nP), "[", rownames(Auc), "]")
      legend("bottomright", S,  col=1:(nC*nP),  lty=1, lwd=1, pch=20, 
           merge=TRUE, inset=0.01, bg="white")
    }
    nClr = 1
  }
  
  #=============================================
  # Calculate AUC by using Wilcox test algorithm
  #=============================================
  if(alg=='Wilcoxon') {
    idxL = vector(mode="list", length=nL) 
    for (i in 1:nL) idxL[[i]] = which(y==uL[i])
    for (j in 1:nC) {                   # for each column representing a feature
      for (i in 1:nP) {                 # go through all permutations of columns in d
        c1 = per[i,1]                   # and identify 2 classes to be compared
        c2 = per[i,2]
        n1 = as.numeric(nY[c1])        
        n2 = as.numeric(nY[c2])
        if (n1>0 & n2>0) {
          r = rank(c(X[idxL[[c1]],j], X[idxL[[c2]],j]))
          Auc[i,j] = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
        }
      } # end of 'for i' loop
    } # end of 'for j' loop
  } 
  
  #==============================================
  # Calculate AUC by using integrating ROC curves
  #==============================================
  if(alg=='ROC') {                        # use 'ROC' algorithm
    for (j in 1:nC) {                     # for each column representing a feature
      x = sort(X[, j], index=TRUE)        # sort all columns and store each one in x[[1]]. x[[2]] stores original positions
      nunq = which(diff(x$x)==0)          # find non-unique A's in column j (if vector is [1 1] nunq=1
      nTies = length(nunq)                # number of non-unique values
      if (nTies<nR-1) {                   # make sure all numbers in A column are not the same
        idx = y[x$ix]                     # reorder label vector in the same order as b, or associate label with each number in b
        # assign column for each label (class) and for each point add 1 in the column corresponding to its class
        d = ( matrix(rep(idx,nL),nR,nL) == L ) 
        for (i in 1:nL) d[,i] = cumsum(d[,i])  # cumulative sum of d columns
        if (nTies) d = d[-nunq, ]         # remove non unique rows if any
        d = rbind( matrix(0,1,nL), d )    # append row of zeros at the beggining
        nD = nrow(d)
        # assume that col#1 ploted on x axis is correct clasification and col#2 (y) is false find AUC
        for (i in 1:nP) {                 # go through all permutations of columns in d
          c1 = per[i,1]                   # and identify 2 classes to be compared
          c2 = per[i,2]
          n  = d[nD,c1]*d[nD,c2]          # normalize area to 1 at the maximum
          if (n>0) Auc[i,j] = trapz(d[,c1], d[,c2])/n  # Trapezoidal numerical integration
          if (plotROC) {                  # plot the results
            xx = if(n>0) d[,c1]/d[nD,c1] else c(0,1)         
            yy = if(n>0) d[,c2]/d[nD,c2] else c(0,1)
            if (2*Auc[i,j]<1) { xx=1-xx; yy=1-yy; } # if auc<0.5 than mirror it to the other side of 0.5
            lines(xx, yy, col=nClr, type='o', pch=20)
            nClr = nClr+1                 # next color
          }
        } # end of 'for i' loop
      }
    } # end of 'for j' loop
  }
  Auc = pmax(Auc, 1-Auc) # if any auc<0.5 than mirror it to the other side of 0.5
  return (Auc)
}
