### Detect and test most linear or most gap-linear subsequences.
### Last modified: 09-28-2007.
### Author: Yanming Di.

## See Di and Michael (2007) for the detailed description of the problem.

## Let Tn = (T1, ..., Tn) be the arrival times of a homogeneous
## Poisson process. We are interested in testing whether Tn possesses
## a sub-sequence that has spacings that are too regular. Two measures of
## regularity are considered, almost-linearity and almost gap-linearity.

## For a sequence Tn, almost linear subsequences of length k can be described
## as follows.  We look at all sub-sequences of length k + 1 of the sequence
## Tn (which have k spacings). For each of the sub-sequence, we look at the
## normalized spacings of it. If these spacings are too close to (1/k, ...,
## 1/k), we will suspect that the subsequence is too "linear".  If we use L1
## norm to measure the distance between the spacings to (1/k, ..., 1/k), it is
## equivalent to look at the (minus) minimum of the spacings. Lastly, to
## incorporate multiple testing effect, we take the maximum of these minimums.

## Almost gap-linear subsequences can be described as follows. Let \sigma be a
## subsequence of (0, ..., n). \sigma can be viewed as the indices of a
## subsequence of Tn. We still look at the spacings of this subsequence of Tn,
## but now we compare the spacings to the standardized spacings of the
## indicies \sigma. For example, when
##     Tn = (11, 14, 15, 20),
## and
##     \sigma = (0, 1, 3)
## we wil compare the standardized spacings of the subsequence (11, 14, 20),
## to the standardized spacings of the corresponding indicies (0, 1, 3).

## Recursive algorithms are used to compute the test statistics.  For each
## pair of (N, k), it's an
##   O(m N^3 k)
## algorithm for the linearity test and
##   O(m N^4 k)
## algorithm for the gap-linearity test. Here
##   N: N+1 is the sequence length;
##   k: k+1 is subsequence length;
##   m: simulation sample size.

## For the linearity and gap-linearity tests described above, we provide the
## following functions.

## test.lin.t(Tn, k): find the most linear length (k+1) subsequence of Tn and
##   compute the test statistic for testing linearity for this subsequence.
##
## test.lin.p(t, n, k): compute the linearity p-value corresponding a test
##   statistic t, which is computed for a length k+1 subsequence of a length
##   n+1 sequence. This p-value is computed by comparing the test statistic t
##   to a precomputed quantile table.
##
## test.lin(Tn): compute the test statistics for Tn and for all k, and
##   compute the corresponding p-values (a seperate test statistic need to be
##   calculated for each k).  This might be the only function need be called
##   directly by an end user.
## 
## test.gaplin.t(Tn, k):
## test.gaplin.p(t, n, k):
## test.gaplin(Tn):
## 
## are the corresponding functions for testing gap-linearity.

## load("tables/q.testlin.rda");
## load("tables/q.testgaplin.rda");

## Oct. 31, 2005

## Now think about the test under interpretation B. For a sequence x
## of length N + 1, we want to find the most "linear" length k + 1
## subsequence of it.  But now, for a subsequence sigma of (1, ...,
## N+1), we compare the standardized spaceings of x[sigma[]] to the
## standardized spacings of sigma[] instead of to (1/k, ..., 1/k).

## Things get more complicated, since some good monotonicity property
## is lost. Things still can be done. But it now takes 5 loops and we
## will have a O(N^5) algorithm. 

## Consider the subsequence having x[s] and x[e] as end points.  We
## focus on this subsequence, and subsubsequences of it. The test
## statistic we need to find is equivalent to max(min(spacings)) of
## the sequence
##   y[] = standardized x[s:e] - 0:N /N,
## where N+1 is the length of the x[s:e].
## We know that y[1] = y[N+1] = 0. Some y[i]'s are negative. 


## Let mmspse[s, e0, k] be the max(min(spacings)) of all length k +1
## subsequence of y[].  Then we have the following recusive relation
## mmspse[s, e0, k+1] = 
##   max_j(min(mmsspse[s, e0-j, k], d[e0-j, e0]),
## where d[a, b] = y[b] - y[a]. d[a,b] can be negative. We need
##   mmsp[s, e, k+1] = mmspse[s, e, k+1].
## Note here, we need the extra loop for e0.

## The above recursion can help us reduce the order of the algorithm
## to O(N^5).  However, it will only find the test statistic and will
## not find the sequence that scores the highest. Modifications can be
## made to achieve that goal, but it may need a lot of memory.

## Extra savings of time comes from the fact that we can find the test
## statistics correspoind to different k (from 2 to N), in one run. We
## even can do all (N, k) pairs in one run, if our goal is to contruct
## the quantile tables for the test statistics. On the other hand,
## even if we just need the statistic of a particular (N, k), we can't
## save too much time.

## The chanllenge is to arrange the loops in a smart order. Another
## challenge is to be able to do the computing incrementally. For
## example, when we calculate the quantile table for N = 100,
## essentially, we will compute all the tables for N = 1, ..., 99. Can
## we save some essential data, so that when later we want to
## calculate N = 200, we can use those to start from 101 in stead of
## 1.

## We need 5 loops to archieve our goal. A loop for the starting point
## s of the subsequence, a loop for the end pioint e0, for each e0, we
## need another loop for e, a loop for the mid point, we call that
## midpoint, and a loop for k. The midpoint divides a subsequence of
## k+1 into a subsequence of length k and subsequence of length 2. The
## ranges for the loops are
##   e = 3:(N+1);
##   s = 1:(e-2);
##   k = 2:(e-s);
##   e0 = (s+k):e;
##   mid=(s+k-1):(e0-1).
## So it is a O(N^5) algorithm and test statistics for ALL k are
## computed.
## Note, we have
##  1 <= s < s+k <= e0 <= e <= N+1,
## and
##  2 <= k <= N.

## Storage, we need
##   d[i,j]: N by N
##   mssp[k, s, e]:
##   mmspse[k, s, e0]:
##         N x (N+1) x (N+1)
##   tB[N,k]: N by N
recursive.B = function (x, indices=FALSE) {
  N0.plus.1 = length(x);
  N0 = N0.plus.1 - 1;

  mmsp = vector("numeric", N0*N0.plus.1*N0.plus.1);
  mmsp[] = -1;
  dim(mmsp) = c(N0, N0.plus.1, N0.plus.1);
  mmspse = vector("numeric", N0*N0.plus.1);
  dim(mmspse) = c(N0, N0.plus.1);
  d = vector("numeric", N0.plus.1 * N0.plus.1);
  dim(d) = c(N0.plus.1, N0.plus.1);
  tB = vector("numeric", N0 * N0);
  dim(tB) = c(N0, N0);
  J = vector("numeric", N0.plus.1);
  J[] = 1;
  dim(J) = c(N0.plus.1, 1);
  for (e in 3:(N0.plus.1)) {
    for (s in 1:(e-2)) {
      N = e - s;
      N.plus.1 = N + 1;
      y = (x[s:e] - x[s])/(x[e] - x[s]) - (0:N)/N;  
      dim(y) = c(1, N.plus.1);
      tmp0 = J[s:e] %*% y;
      d[s:e, s:e] = tmp0 - t(tmp0); 
      mmspse[]  = -1;
      mmspse[1, ] = d[s, ]; 
      for (k in 2:(e-s)) {
        k1 = k - 1;
        for (e0 in (s+k):e) {
          tmp = -1;
          for (mid in (s+k1):(e0-1)) {
            tmp1 = min(mmspse[k1, mid], d[mid, e0]);
            if (tmp < tmp1) tmp = tmp1;
          } # mid loop
          mmspse[k, e0] = tmp;
        } # e0 loop
      }# k loop
      mmsp[, s, e] = mmspse[, e];
    }# s loop
    for (k0 in 2:(e-1)) { tB[e-1,k0] = max(mmsp[k0,,]) + 1/k0;}
  }# e loop
  if (indices) {
    sigma = array(0, c(N0, 2));
    for (k in 2:N0) {
      tmp = which.max(mmsp[k,,]);
      sigma[k, 2] = (tmp - 1) %/% N0.plus.1  + 1;
      sigma[k, 1] = (tmp - 1) %% N0.plus.1 + 1; 
    }
    return(list(tB=tB, sigma = sigma));
  } else {
    return(tB);
  }
}

most.gaplin.sub  = function(x, k, T) {
  ## Search among all length k+1 subsequence of x the most "linear" one.
  ## Subsequences searched here have the same end points as x.  Here, we
  ## consider "linear" under interpretation B. For the subsequence
  ##   x[i0], ..., x[ik]
  ## we compare the standardized spacings of it to the standardized
  ## spacings of
  ##   i0, ..., ik
  ## Only those with min(spacings + 1/k - e_sigma) > T are considered, 
  ## where
  ##   e_sigma[j] = (i[j] - i[0]) / (i[k] - i[0]); 
  ##
  ## Input:
  ##   x=(x0, x1, ..., xn) vector of length n+1.
  ##   k: see above. 
  ##   T: current best max(min(spacings + 1/k - e_sigma)). 
  ## 
  ## Output:
  ##   improved: whether we've found a better deal! True if a
  ## subsequence with min(spacings) > T has been found.
  ##   sub: the highest scored subsequence.
  ##     T: the new max(min(spacings)).

  ## Note: In Michael's note. x = (x0, ..., xn), but the indices in R
  ## starting with 1.  So it might be confusing sometimes.

  ## Here, we just borrow some code from "recusive.B".

  N.plus.1 = length(x);
  N = N.plus.1 - 1;

  J = vector("numeric", N.plus.1);
  J[] = 1;
  dim(J) = c(N.plus.1, 1);

  mmspse = vector("numeric", N*N.plus.1);
  dim(mmspse) = c(N, N.plus.1);

  sigma = mmspse;
  sigma[1,] = 1;

  y = (x[] - x[1])/(x[N.plus.1] - x[1]) - (0:N)/N;  
  dim(y) = c(1, N.plus.1);
  tmp0 = J %*% y;
  d = tmp0 - t(tmp0); 

  mmspse[]  = -1;
  mmspse[1, ] = d[1, ]; 

  
  for (k0 in 2:k) {
    k1 = k0 - 1;
    for (e0 in (1+k0):N.plus.1) {
      tmp = -1;
      for (mid in k0:(e0-1)) {
        tmp1 = min(mmspse[k1, mid], d[mid, e0]);
        if (tmp < tmp1) {
          tmp = tmp1;
          sigma[k0, e0] = mid;
        }
      } # mid loop
      mmspse[k0, e0] = tmp;
    } # e0 loop
  }# k loop
  t = mmspse[k, N.plus.1] + 1/k;

  sigma0 = vector("numeric", k+1);
  sigma0[k+1] = N.plus.1;
  for (k0 in k:1) {
    sigma0[k0] = sigma[k0, sigma0[k0+1]];
  }
  list(improved=(t>T), t=t, sub=x[sigma0], sigma = sigma0);
}

most.linear.sub = function(x, k, t) {
## Search among all length k+1 subsequence of x the most "linear" one.
## The subsequences searched here have the same end points as x.
## Only those with min(spacings) > T are considered. 

## Input:
##   x=(x0, x1, ..., xn) vector of length n+1.
##   k: see above. 
##   T: current best max(min(spacings)). 
## 
## Output:
##   improved: whether we've found a better deal! True if a
## subsequence with min(spacings) > T has been found.
##   sub: the highest scored subsequence.
##     T: the new max(min(spacings)).

## Note: In Michael's note. x = (x0, ..., xn), but the indices in R starting with 1. 
## So it might be confusing sometimes.

## The key idea here is we only want to consider subsequences with
## min(spacings) > T and we elmininate those unqualified ones as early
## as possible. So only Nk key operations are needed and most of them
## are comparisons and "+".

  N.plus.1 = length(x);
  N = N.plus.1 - 1;
  ## assert(k <= N) # k = N is OK now.
  ## x should be sorted!
  y =  x[] / (x[N.plus.1] - x[1]);
  
  ## sigma = array(0, k+1); # array is a slow thing to do in R
  sigma = sigma0 = 1:(k+1);
  ## sigma[1] = 1;
  sigma[k+1] = N.plus.1;
  ## The end points of the subsequences are fixed and are the same as x's.

  ## Uncomment the following for debugging
  ## plot(cbind(x, 0));
  ## yy = 0;

  eps = .Machine$double.eps;
  max.min.sp = t;
  threshold = max.min.sp + eps;

  ## We will only consider subsequences with min(spacings) greater
  ## than the threshold value. Once we find a more linear subsquence,
  ## we will update the thresshold value accordingly. E.g, if
  ##   sigma = (s0, s1, ..., si, ..., sk)
  ## are the indices to the currently best subsequence and the spacing i
  ##   (y[si-1], y[si])
  ## is the min spacing of it. Then during our next update. we can use this
  ## min spacing as the new threshold.  Also, si will be the first
  ## index, the algorithm needs to adjust.
  
  ind = 1;
  ind.sigma = 1;
  starting.from = 1;
  repeat{
    for (i in ind.sigma:k) {
      ind = ind + 1;
      while ((y[ind] - y[sigma[i]] < threshold) && (ind <= N.plus.1))
        { ind = ind +1;}
      if (ind > N.plus.1) {
        break;
      }
      sigma[i+1] = ind;
    }
    if (ind > N.plus.1) {
      break
    }
    ## if the we get here, we've got a better deal!
    sigma[k+1]=N.plus.1; # sigma[k+1] might be changed in the for loop
    sigma0 = sigma; # x[sigma] is a more linear "subsequence" 
    ## max.min.sp = min(diff(y[sigma]));
    spacings = diff(y[sigma]);
    ind.sigma = which.min(spacings);
    max.min.sp = spacings[ind.sigma];
    threshold = max.min.sp + eps;
    ind = sigma[ind.sigma + 1] - 1;
    ## debug
    ## points(cbind(x[sigma], yy), col="red");
    ## yy=yy + 0.001 
  }
  ## debug
  ## print(yy); 
  list(improved=(max.min.sp > t), t=max.min.sp, sub=x[sigma0]);
}

test.lin.t = function(Tn, k) {
  ## Calculate the test statistic for the sequence Tn for subsequences
  ## of length k+1 (which have k spacings). See (33) of Michael's note
  ## for details.

  ## For the sequence Tn, the test statistic can be calculated as
  ## follows. We look at all sub-sequences of length k+1 of Tn. For
  ## each such sbu-sequence, we look at the normalized spacings of
  ## it. If this spacings are too close to (1/k, ..., 1/k) under |.|
  ## norm. We will suspect that it is too "regular" or "linear". It is
  ## equivalent to look at the (minus) minimum spacing. We then find
  ## the maximum of all the minimums.

  ## Arguments:
  ## Input:
  ##   Tn: the sequence to be tested
  ##   k: k + 1 is the length of subquences
  ## Output is a list contains:
  ##   t: the test statistic
  ##   sub: (one of) the highest scored subsequence, the most "linear"
  ##     subsequence

  Tn = sort(Tn);
  cur = most.linear.sub(Tn, k, 0);
  ## Find one subsequence first, then we'll only consider sequences
  ## more linear than this one. The end points of the subsequences
  ## here are the same as Tn's
  N.plus.1 = length(Tn);
  N = N.plus.1 - 1;
  
  ## Search subsequences with other possible end points.
  for (i in 1:(N-k+1))
    for (j in (i+k):(N.plus.1)) {
      try = most.linear.sub(Tn[i:j], k, cur$t);
      if (try$improved) {
         cur = try;
      }
   };

  list(t = cur$t, sub = cur$sub);  
}

test.gaplin.t = function(Tn, k) {
  ## Calculate the test statistic for the sequence Tn for subsequences
  ## of length k+1 (which have k spacings). See (33) of Michael's note
  ## for details.

  ## For the sequence Tn, the test statistic can be calculated as
  ## follows. We look at all sub-sequences of length k+1 of Tn. For
  ## each such sbu-sequence, we look at the normalized spacings of
  ## it. If this spacings are too close to (1/k, ..., 1/k) under |.|
  ## norm. We will suspect that it is too "regular" or "linear". It is
  ## equivalent to look at the (minus) minimum spacing. We then find
  ## the maximum of all the minimums.

  ## Arguments:
  ## Input:
  ##   Tn: the sequence to be tested
  ##   k: k + 1 is the length of subquences
  ## Output is a list contains:
  ##   t: the test statistic
  ##   sub: (one of) the highest scored subsequence, the most "linear"
  ##     subsequence

  Tn = sort(Tn);
  cur = most.gaplin.sub(Tn, k, 0);
  ## Find one subsequence first, then we'll only consider sequences
  ## more linear than this one. The end points of the subsequences
  ## here are the same as Tn's
  N.plus.1 = length(Tn);
  N = N.plus.1 - 1;
  
  ## Search subsequences with other possible end points.
  for (i in 1:(N-k+1))
    for (j in (i+k):(N.plus.1)) {
      try = most.gaplin.sub(Tn[i:j], k, cur$t);
      if (try$improved) {
         cur = try;
      }
   };

  list(t = cur$t, sub = cur$sub);  
}

test.lin.p = function(t, n, k) {
##   Calculate the p-value corresponding to test statistic t. A table
##   need to be constructed (using test.lin.table) before calling this
##   function.

##  table.name = sprintf('tables/test.A.table.n%d', as.integer(n));
##   q = as.matrix(read.table(table.name, sep=","));
  q = get(sprintf("q.testlin.n%d", as.integer(n)));
  i = sum(q[k,] < t);
  p = 1 - q[1,i];

}

test.gaplin.p = function(t, n, k) {
##   Calculate the p-value corresponding to test statistic t. A table
##   need to be constructed (using test.gaplin.table) before calling this
##   function.
  ## table.name = sprintf("tables/test.B.table.n%d", as.integer(n));
  ## q = as.matrix(read.table(table.name, sep=","));
  q = get(sprintf("q.testgaplin.n%d", as.integer(n)));
  i = sum(q[k,] < t);
  p = 1 - q[1,i];
}

test.lin = function(Tn) {
  ## Calculate the test for Tn for all k and give the p-values. (A
  ## seperate test statistic need to be calculated for each k.) This
  ## might be the only function need be called directly by the user.
  n = length(Tn) - 1;
  cat("Tn=(");
  cat(Tn, sep=",");
  cat(")\n");
  
  cat(sprintf('%2s  %5s  %8s  %s\n', 'k', 't_lin', 'p-value',
  'highest scored subsequence'));
  # do we need k=n?
  for (k in 2:n) {
    tea = test.lin.t(Tn, k);
    p = test.lin.p(tea$t, n, k);
    cat(sprintf('%2d  %5.3f  %8.3f   (', as.integer(k), tea$t, p));
    cat(tea$sub, sep=",");
    cat(")\n");
  }

  invisible(NULL);
}

test.gaplin = function(Tn) {
  ## Calculate the test for Tn for all k and give the p-values. (A
  ## seperate test statistic need to be calculated for each k.) This
  ## might be the only function need be called directly by the user.
  Tn = sort(Tn);
  N = length(Tn) - 1;
  cat("Tn=(");
  cat(Tn, sep=",");
  cat(")\n");

  tea = recursive.B(Tn, TRUE);
  tB = tea$tB[N,];

  cat(sprintf('%2s  %5s  %8s  %s\n', 'k', 't_gaplin', 'p-value',
  'highest scored subsequence'));
  for (k in 2:N) {
    ## tea = test.gaplin.t(Tn, k);
    ## p = test.gaplin.p(tea$t, N, k);
    ## cat(sprintf('%2d  %5.3f  %8.3f   (', as.integer(k), tea$t, p));
    ## cat(tea$sub, sep=",");

    p = test.gaplin.p(tB[k], N, k);
    sub = most.gaplin.sub(Tn[tea$sigma[k,1]:tea$sigma[k,2]], k, 0)$sub;
    cat(sprintf('%2d  %5.3f  %8.3f   (', as.integer(k), tB[k], p));
    cat(sub, sep=",");
    cat(")\n");
  }
  invisible(NULL);
}
