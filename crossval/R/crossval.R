### crossval.R  (2014-03-29)
###
###    Generic Function for Cross Valdidation 
###
### Copyright 2009-14  Korbinian Strimmer
###
###
### This file is part of the `crossval' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# general routines for conducting cross validation


#  predfun(Xtrain, Ytrain, Xtest, Ytest, ...)
#  must be given

# B times repeated K-fold cross-validation 
crossval = function(predfun,
                       X, Y, 
                       K=10,  # number of folds
                       B=20,  # number of repetitions
                       verbose=TRUE,
                       ...    # optional arguments for predfun
)
{
  ygrouped = group.samples(Y)
 
  # makes no sense to have more folds than entries in largest group
  groupsize = sapply(ygrouped, length)
  nfolds = min(K, max(groupsize))
  if (verbose) cat("Number of folds:", nfolds, "\n")


  allfolds = B*nfolds
  if (verbose) cat("Total number of CV fits:", allfolds, "\n")
 
  stat.cv = NULL
 
  i = 1
  for (b in 1:B)
  {
    if (verbose) cat("\nRound #", b, "of", B, "\n")

    folds = get.folds(ygrouped, K=nfolds)

    for (f in 1:nfolds)
    {
      if (verbose) cat("CV Fit #", i, "of", allfolds, "\n")

      #### prepare  test and training data set ####

      test.idx = folds[[f]]
      train.x = X[-test.idx, , drop=FALSE]                     
      train.y = Y[-test.idx]
      test.x = X[test.idx, , drop=FALSE] 
      test.y = Y[test.idx]

      ### learn predictor and compute test error ####
      stat.new = predfun(train.x, train.y, test.x, test.y, ...)
      stat.cv = rbind(stat.cv, stat.new)

      rownames(stat.cv)[i] = paste0("B",b,".F",f)

      i = i+1
    }   
  }

  stat = apply(stat.cv, 2, mean)
  stat.se = apply(stat.cv, 2, sd) / sqrt(allfolds)

  return(list(stat.cv=stat.cv, stat=stat, stat.se=stat.se))
}


########################################################################
# private functions
########################################################################

# return a list with samples arranged by group
group.samples = function(y)
{
  # split samples into groups
  if(is.factor(y))
  {
    ygrouped = split(seq(y), y)
  }
  else
  {
    ygrouped = list(all=seq(length(y)))
  }

  return( ygrouped )
}

# divide samples into sets of similar size with 
# evenly distributed samples (balanced per group) 
get.folds = function(ygrouped, K)
{
  # makes no sense to have more folds than entries in largest group
  groupsize = sapply(ygrouped, length)
  nfolds = min(K, max(groupsize))
  if (K != nfolds) cat("Number of folds:", nfolds, "\n")
  
  # assign the members of each group evenly to the folds
  m = NULL
  for (i in 1:length(ygrouped) )
  {
     a = ceiling(groupsize[i]/nfolds)
     ridx = sample.int(groupsize[i], groupsize[i])
     v = c( rep(NA, nfolds*a-groupsize[i]), ygrouped[[i]][ridx] ) # pad with NAs
     ridx = sample.int(nfolds, nfolds) # reshuffle column containing the NAs
     v[1:nfolds] = v[ridx]
     m = c(m,v)
  }
  m = matrix(m, nrow=nfolds) # note that all NAs of a group are all in one column

  folds =  vector("list", nfolds)
  for(j in 1:nfolds)
  {
    keep = !is.na(m[j , ])
    folds[[j]] = m[j, keep]
  }

  return( folds )
}




