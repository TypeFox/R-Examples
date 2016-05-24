# author: Kobi Perl
# Based on the following thesis:
#   Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences. Haifa. 
#   Retrieved from http://bioinfo.cs.technion.ac.il/people/zohar/thesis/eran.pdf

# mHG definition:
#   mHG(lambdas) = min over 1 <= n <= N of HGT (b_n(lambdas); N, B, n)
# Where HGT is the hypergeometric tail:
#   HGT(b; N, B, n) = Probability(X >= b)
# And:
#   b_n = sum over 1 <= i <= n of lambdas[i]
# Fields:
#   mHG - the statistic itself
#   n - the index for which it was obtained
#   b (short for b_n) - sum over 1 <= i <= n of lambdas[i]
setClass("mHG.statistic.info", slots = c(mHG = "numeric", n = "numeric", b = "numeric"))

# We define EPSILON to account for small changes in the calculation of p-value
# between the function calculating the statistic and this function calculating the p-values
# Specifically, if (statistic + EPSILON) is higher than the hypergeometric tail associated by a cell
# in the W*B path matrix, used in the p-value calculation, then the cell is not in the "R region"
# We warn ifthe mHG statistic gets below -EPSILON
EPSILON <- 1e-10 

mHG.test <- function(lambdas, n_max = length(lambdas)) {
  # Performs a minimum-hypergeometric test.
  # Test is based on the following thesis:
  #   Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences. Haifa. 
  #   Retrieved from http://bioinfo.cs.technion.ac.il/people/zohar/thesis/eran.pdf
  # 
  # The null-hypothesis is that the 1s in the lambda list are randomly and uniformly 
  # distributed in the lambdas list. The alternative hypothesis is that the 1s tend
  # to appeard in the top of the list. As the designation of "top" is not a clear-cut
  # multiple hypergeometric tests are performed, with increasing length of lambdas being
  # considered to be in the "top". The statistic is the minimal p-value obtained in those
  # tests. A p-value is calculated based on the statistics.
  #
  # Input:
  #   lambdas - {0,1}^N, sorted from top to bottom.
  #   n_max - the algorithm will only consider the first n_max partitions.
  # Output:
  #   "htest" S3 class - with the following fields:
  #     statistic - The mHG statistic. Definition:
  #       mHG(lambdas) = min over 1 <= n < n_max of HGT (b_n(lambdas); N, B, n)
  #       Where HGT is the hypergeometric tail:
  #         HGT(b; N, B, n) = Probability(X >= b)
  #       And:
  #         b_n = sum over 1 <= i <= n of lambdas[i]
  #     n - the index for which the statistic was obtained
  #     b (short for b_n) - sum over 1 <= i <= n of lambdas[i]
  #     parameters:
  #       N - total number of white and black balls.
  #       B - number of black balls.
  #       n_max
  #     p.value
  N <- length(lambdas)
  B <- sum(lambdas)
  W <- length(lambdas) - B
  
  mHG.statistic.info <- mHG.statistic.calc(lambdas, n_max) # The uncorrected for MHT p-value
  p <- mHG.pval.calc(mHG.statistic.info@mHG, N, B, n_max)
  result <- list()
  class(result) <- "htest"
  result$statistic <- c("mHG" = mHG.statistic.info@mHG)
  result$parameters <- c("N" = N, "B" = B, "n_max" = n_max)
  result$p.value <- p

  result$n <- mHG.statistic.info@n # Not an official field of htest
  result$b <- mHG.statistic.info@b # Not an official field of htest
  
  return(result)
}

mHG.statistic.calc <- function(lambdas, n_max = length(lambdas)) {
  # Calculates the mHG statistic.
  # mHG definition:
  #   mHG(lambdas) = min over 1 <= n < N of HGT (b_n(lambdas); N, B, n)
  # Where HGT is the hypergeometric tail:
  #   HGT(b; N, B, n) = Probability(X >= b)
  # And:
  #   b_n = sum over 1 <= i <= n of lambdas[i]
  # In R, HGT can be obtained using:
  #   HGT(b; N, B, n) = phyper((b-1), B, N - B, n, lower.tail = F)
  #
  # Input:
  #   lambdas - sorted and labeled {0,1}^N.
  #   n_max - the algorithm will only consider the first n_max partitions.
  # Output: mHG.statistic
  # 
  # Statistic is defined in the following thesis:
  #   Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences. Haifa. 
  #   Retrieved from http://bioinfo.cs.technion.ac.il/people/zohar/thesis/eran.pdf
  # 
  # If several n gives the same mHG, then the lowest one is taken.
  
  # Input check
  stopifnot(n_max > 0)
  stopifnot(n_max <= length(lambdas))
  stopifnot(length(lambdas) > 0)
  stopifnot(all(lambdas == 1 | lambdas == 0))
  
  N <- length(lambdas)
  B <- sum(lambdas)
  W <- N - B
  
  mHG <- 1
  mHG.n <- 0
  mHG.b <- 0
  m <- 0 # Last time we saw a one
  HG_row <- numeric(B + 1) # The first B + 1 hypergeometric probabilities, HG[i] = Prob(X == (i - 1))
                       # updated for the current number of tries n.
  HG_row[1] <- 1 # For n = 0, b = 0
  b <- 0
  n <- 0
  while (n < n_max) {  # iterating on different N to find minimal HGT
    n <- n + 1
    b <- b + lambdas[n]
    
    if (lambdas[n] == 1) { # Only then HGT can decrease (see p. 19 in thesis)
      HG_row <- HG_row_n.calc(HG_row, m, n, b, N, B)
      m <- n
      
      HGT <- 1 - sum(HG_row[1:b]) # P(X >= b) = 1 - P(X < b)
      # statistic
      if (HGT < mHG) {
        mHG <- HGT
        mHG.n <- n
        mHG.b <- b
      }      
    }    
  }
  return(new("mHG.statistic.info", mHG = mHG, n = mHG.n, b = mHG.b))
}

mHG.pval.calc <- function(p, N, B, n_max = N) {
  # Calculates the p-value associated with the mHG statistic.
  # Guidelines for the calculation are to be found in:
  #   Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences. Haifa. 
  #   Retrieved from http://bioinfo.cs.technion.ac.il/people/zohar/thesis/eran.pdf
  # (pages 11-12, 19-20)
  # Input:
  #   p - the mHG statistic. Marked as p, as it represenets an "uncorrected" p-value.
  #   N - total number of white and black balls (according to the hypergeometric problem definition).
  #   B - number of black balls.
  #   n_max - the algorithm will calculate the p-value under the null hypothesis that only the 
  #           first n_max partitions are taken into account in determining in minimum.
  # Output: p-value.
  
  # Input check
  stopifnot(n_max > 0)
  stopifnot(n_max <= N)
  stopifnot(N >= B)
  stopifnot(p <= 1)
  
  if (p < -EPSILON) {
	  warning("p-value calculation will be highly inaccurate due to an extremely small mHG statistic")
  }
  
  # p - the statistic.
  # N\B - number of all \ black balls.
  W <- N - B
  
  R_separation_line <- R_separation_line.calc(p, N, B, n_max)
  
  pi_r <- pi_r.calc(N, B, R_separation_line)
  p_corrected <- 1 - pi_r[W+2,B+2]
  return(p_corrected)
}

R_separation_line.calc <- function(p, N, B, n_max) {
  # Determine R separation line - This is the highest (broken) line crossing the B*W matrix horizontally, that underneath it all
  # the associated p-values are higher than p, or w + b > n_max.
  #
  # (This is a bit different from the original definition to make the calculation more efficient)
  #
  # Input:
  #   p - the mHG statistic. Marked as p, as it represenets an "uncorrected" p-value.
  #   N - total number of white and black balls (according to the hypergeometric problem definition).
  #   B - number of black balls.
  #   n_max - Part of the constraint on the line, the null hypothesis is calculated under
  #           the assumption that the first n_max partitions are taken into account in determining the minimum.
  # Output:
  #   R_separation_line - represented as a vector size B + 1, index b + 1 containing 
  #                       the first (high enough) w to the right of the R separation line (or W + 1 if no such w exist).
  # See:
  #   Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences. Haifa. 
  #   Retrieved from http://bioinfo.cs.technion.ac.il/people/zohar/thesis/eran.pdf
  #   (pages 11-12)
  W <- N - B
  R_separation_line <-  rep(W + 1, times = B + 1)
  
  HG_row <- numeric(B) # First B in HG_row
  HG_row[1] <- 1 # For n = 0, b = 0
  b <- 0
  w <- 0
  HGT <- 1 # Initial HGT
  
  # We are tracing the R line - increasing b until we get to a cell where the associated p-values are smaller
  # than p, and then increasing w until we get to a cell where the associated p-values are bigger than p
  should_inc_w <- (HGT <= (p + EPSILON)) && (w < W) && (b <= (n_max - w))
  should_inc_b <- (HGT > (p + EPSILON)) && (b < B) && (b <= (n_max - w))
  
  while (should_inc_w || should_inc_b) {
    while (should_inc_b) { # Increase b until we get to the R line (or going outside the n_max zone)
      R_separation_line[b+1] <- w
      b <- b + 1
      HG_row[b+1] <- HG_row[b] * d_ratio(b + w, b, N, B)
      HG_row[1:b] <- HG_row[1:b] * v_ratio(b + w, 0:(b-1), N, B)
      HGT <- 1 - sum(HG_row[1:b]) # P(X >= b) = 1 - P(X < b)
      should_inc_b <- (HGT > (p + EPSILON)) && (b < B) && (b <= (n_max - w))
    }
    if (b > (n_max - w)) {
      # We can stop immediately and we do not need to calculate HG_row anymore
      R_separation_line[(b+1):(B+1)] <- w
      should_inc_w <- F
    } else {
      should_inc_w <- (HGT <= (p + EPSILON)) && (w < W)
      while (should_inc_w) { # Increase w until we get outside the R line (or going outside the n_max zone)
        w <- w + 1
        HG_row[1:(b+1)] <- HG_row[1:(b+1)] * v_ratio(b + w, 0:b, N, B)
        HGT <- 1 - sum(HG_row[1:b]) # P(X >= b) = 1 - P(X < b)
        should_inc_w <- (HGT <= (p + EPSILON)) && (w < W) && (b <= (n_max - w))
      }
      if (b > (n_max - w)) {
        # We can stop immediately and we do not need to calculate HG_row anymore
        R_separation_line[(b+1):(B+1)] <- w
        should_inc_b <- F        
      } else {
        should_inc_b <- (HGT > (p + EPSILON)) && (b < B) && (b <= (n_max - w))
      }
    }
  }
  if (HGT > (p + EPSILON)) # Last one
    R_separation_line[b+1] <- w
  return(R_separation_line)
}

pi_r.calc <- function(N, B, R_separation_line) {
  # Consider an urn with N balls, B of which are black and W white. pi_r stores 
  # The probability of drawing w white and b black balls in n draws (n = w + b)
  # with the constraint of P(w,b) = 0 if (w, b) is on or above separation line.
  # Row 1 of the matrix represents w = -1, Col 1 represents b = -1.
  #
  # Input:
  #   N - total number of white and black balls (according to the hypergeometric problem definition).
  #   B - number of black balls.
  #   R_separation_line - represented as a vector size B + 1, index b + 1 containing 
  #                       the first (high enough) w to the right of the R separation line.
  # See:
  #   Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences. Haifa. 
  #   Retrieved from http://bioinfo.cs.technion.ac.il/people/zohar/thesis/eran.pdf
  #   (pages 20)
  W <- N - B
  
  pi_r <- matrix(data = 0, nrow = W + 2, ncol = B + 2)
  pi_r[1,] <- 0 # NOTE: Different from the thesis (see page 20 last paragraph),
                # should be 1 according to that paragraph, but this seems wrong.
  pi_r[,1] <- 0
  
  for (b in 0:B) {
    w <-  R_separation_line[b+1]
    while (w < (W + 1)) {
      if ((w == 0) && (b == 0)) {
        pi_r[2,2] <- 1 # Note, this cell will be 0 if it's left to the R separation line (should not occure) 
      } else {
        # Apply the recursion rule:
        # P(w,b) = P((w,b)|(w-1,b))*P(w-1,b)+P((w,b)|(w,b-1))*P(w,b-1)
        pi_r[w + 2, b + 2] <- (W - w + 1) / (B + W - b - w + 1) * pi_r[w + 1, b + 2] +
          (B - b + 1)/ (B + W - b - w + 1) * pi_r[w + 2, b + 1]
      }
      w <- w + 1
    }
  }
  return(pi_r)
}

HG_row_n.calc <- function(HG_row_m, m, n, b_n, N, B) {
  # Calculate HG row n. This row contains the first (b_n  + 1)
  # hypergeometric probabilities, HG[i] = Prob(X == (i - 1)), for number of tries n.
  # Does so given an updated HG row m (m < n), which contains the first (b_n)
  # hypergeometric probabilities.
  #
  # Input:
  #   HG_row_m - updated HG row m (m < n), which contains the first (b_n)
  #              hypergeometric probabilities.
  #   m - the number of tries (m < n) for which the HG_row_m fits.
  #   n - the number of tries (n > m) for which we want to calculate the HG row
  #   b_n - The maximal b for which we need to calculate the hypergeometric probabilities.
  #   N - total number of white and black balls (according to the hypergeometric problem definition).
  #   B - number of black balls.
  
  # The function directs the calculation to an iteration solution (with the cost of B(n-m)) 
  # or a recursive solution (with the cost B * log(B)). This multiplier helps to determine
  # when to use the recursion solution - it is not a theoretical result, but an empirical one.
  RECURSION_OVERHEAD_MULTIPLIER <- 20
  
  HG_row_n.calc.func <- NULL
  if ((n - m) <= (RECURSION_OVERHEAD_MULTIPLIER * log2(b_n))) {
	  HG_row_n.calc.func <- HG_row_n.calc.iter
  } else {
  	HG_row_n.calc.func <- HG_row_n.calc.recur
  }
  return(HG_row_n.calc.func(HG_row_m, m, n, b_n, N, B))
}


HG_row_n.calc.recur <- function(HG_row_m, m, n, b_n, N, B) {
  # Calculate HG row n recursively.
  # See function documentation for "HG_row_n.calc", to gain insight on input and outputs. 
  
  # NOTE: This implementation is my interpretation of a very unclear statement in Eden's thesis that a row can be
  # calculated in O(B) without considering previous rows. I did this in O(B * log(B)).
  
  HG_row_n.calc.recur.inner <- function(b_n_start, b_n_end, m_start)   {
    # The code works directly on HG_row_m - updating it recursively, filling the right
    # values from right to left. Filling this row in a directed manner allow us to update a cell
    # in a recursive manner according to the one left to it, without being concerned that the cell
    # to the left was already updated to a "too big" number of tries.
    #
    # The code work on the subtree induced by the limits:
    # Rows (m_start, n) and columns (b_n_start, b_n_end),
    # updating the subtree in a recursive manner until row_n.
    # 
    # The code assumes that:
    #   HG_row_m[b_n_start] has the correct hypergeometric probability value for number of tries (row) 
    #   m_start.
      
    # Stop condition
    if ((n - m_start) == 0) {
      return()
    } else {
      # split the tree to two subtrees to be evaluated separately.
      # R_tree Will be used to calculate the HG_row_n entries corresponding to (b_n_start:b_n_start + r_split),
      # and L_tree will be used to calculate the rest.
      r_split <- floor((b_n_end - b_n_start + 1) / 2)
      l_split <- (b_n_end - b_n_start + 1) - r_split
      
      # If m is above the root of the tree we are working on - then some rows were already
      # calculated and we can take this to our advantage.
      rows_already_calc <- max(m - m_start, 0)
      
      # Go diagonally (increasing both b_n_start and m) until we get to the root of the right tree.
      i <- rows_already_calc
      while (i < l_split) { # Needs to occur sequentially
        i <- i + 1
        HG_row_m[b_n_start + i + 1] <<- HG_row_m[b_n_start + i] * d_ratio(m_start + i, b_n_start + i, N, B)
      }
      # Calculate the right tree.
      HG_row_n.calc.recur.inner(b_n_start + l_split, b_n_end, m_start + l_split)
      
      # Go upwards (increasing only m) until we get to the root of the left tree.      
      i <- rows_already_calc
      while (i < (n + 1 - l_split - m_start)) { # Needs to occur sequentially
        i <- i + 1
        HG_row_m[b_n_start + 1] <<- HG_row_m[b_n_start + 1] * v_ratio(m_start + i, b_n_start, N, B)
      }
      # Calculate the left tree.
      HG_row_n.calc.recur.inner(b_n_start, b_n_start + l_split - 1, n + 1 - l_split)
    }
  }

  # HG_row_n is calculated from a tree, for which the root is located in
  # b_n rows beneath - we will mark this as m_start.
  # If row m is "beneath" this root, We "go up" and calculate it.
  # If we are "above" this root, we will ignore this for now, and use the fact that
  # we have some rows calculated above m_start, in the recursion.
  # NOTE: We can technically initialize HG_row_m[2:] to 0, but knowing that we will
  # only use the cell in index 1, we do not bother.
  m_start <- n - b_n
  while (m < m_start) {
    m <- m + 1
    HG_row_m[1] <- HG_row_m[1] * v_ratio(m, 0, N, B)
  }
  HG_row_n.calc.recur.inner(0, b_n, m_start) # NOTE: This code works on HG_row_m directly.
  return(HG_row_m)
}

HG_row_n.calc.iter <- function(HG_row_m, m, n, b_n, N, B) {
  # Calculate HG row n recursively.
  # See function documentation for "HG_row_n.calc", to gain insight on input and outputs. 
  
  # NOTE: The code works directly on HG_row_m, m - increasing m until it becomes n.

  # Go upwards (increasing only m) until we get to row n-1.    
  b_to_update <- 0:(b_n - 1)
  while (m < (n - 1)) {
    m <- m + 1
    HG_row_m[b_to_update+1] <- HG_row_m[b_to_update+1] * v_ratio(m, b_to_update, N, B)    
  }  
  m <- m + 1
  # Last row to go - first update b_n from the diagonal, then the rest vertically
  HG_row_m[b_n+1] <- HG_row_m[b_n] * d_ratio(m, b_n, N, B)
  HG_row_m[b_to_update+1] <- HG_row_m[b_to_update+1] * v_ratio(m, b_to_update, N, B)
  
  return(HG_row_m)    
}

d_ratio <- function(n, b, N, B) {
  # The ratio between HG(n,b,B,N) and HG(n-1,b-1,B,N)
  # See page 19 in Eden's theis.
  return(n * (B - (b-1)) / (b * (N - (n-1))))
}
v_ratio <- function(n, b, N, B) {
  # The ratio between HG(n,b,B,N) and HG(n-1,b,B,N)
  # See page 19 in Eden's theis
  return((n*(N-n-B+b+1)) / ((n - b) * (N - n + 1)))
}
