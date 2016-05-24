#==============================================================================
#' @importFrom stats sd
InitExprDis = function(x, mut=NULL) {
  # Initiate Gaussian Mixute Model for MAP estimation
  #
  # Args:
  #   x: the mRNA expression vector of gene with patients as the names
  #   mut: the mutation vector for gene with patients as the names
  # Return:
  #   model: the intial model
  # 
  # Date: 
  #   updated date: March 3, 2012
  #   Revised: February 15, 2015
  #
  # Author: Jiarui Ding <jiaruid@cs.ubc.ca>
  #   Department of Computer Science, UBC
  #   Department of Molecular Oncology, BC Cancer Agency 
  #
  
  x = x[!is.na(x)]
  N = length(x)
  if (N < 20)
    stop("The number of patients is too small to robustly estimate 
          the gene expression distributions")
    
  # Take the min and max as outliers
  id.min = which.min(x)
  id.max = which.max(x)
  x.min  = x[id.min] 
  x.max  = x[id.max]
  id.out = c(id.min, id.max)
  x.rm   = x[-id.out]
  
  mut = mut[intersect(names(mut), names(x.rm))]
  len = sum(mut)
  
  # Maybe only part of the patients with mutation data and 
  # make sure not all the patients have this mutation
  if (len > 0 & len < length(mut) - 2) {
    patient  = intersect(names(x.rm), names(which(mut>0)))
    x.mut    = x.rm[patient]
    patient.nonmut = setdiff(names(mut), patient)  
    x.nonmut = x.rm[patient.nonmut]
    
    # 
    mu.nonmut = mean(x.nonmut)
    sd.nonmut = sd(x.nonmut)
    
    x.mut.up   = x.mut[x.mut >= mu.nonmut + sd.nonmut*1.0]
    x.mut.down = x.mut[x.mut <= mu.nonmut - sd.nonmut*1.0]
    
    if (length(x.mut.up) == 0) {
      x.up = x.max 
    } else {
      x.up = c(x.mut.up, x.max)
    }

    if (length(x.mut.down) == 0) {
      x.down = x.min
    } else {
      x.down = c(x.mut.down, x.min)
    }
    x.neutral = x.nonmut
  } else {
    # Find the outliers based on boxplot rule
    box.out = DetectBoxplotOutlier(x.rm)   
    outlier = box.out$outlier.value
    
    if (length(outlier) > length(x.rm) / 2) {
      x.up      = x.max
      x.down    = x.min
      x.neutral = x.rm
    } else {
      id.up     = outlier > box.out$cu
      x.up      = c(outlier[id.up], x.max)
      x.down    = c(outlier[!id.up], x.min)
      x.neutral = x.rm[box.out$keep.id]
    }
  }
  mu = c(mean(x.down), mean(x.neutral), mean(x.up))
  
  num.up      = length(x.up)
  num.down    = length(x.down)
  num.neutral = length(x.neutral) 
  lambda = c(num.down, num.neutral, num.up) /
           (num.down + num.neutral + num.up)
  
  sigma.all     = sd(x)
  sigma.up      = sd(x.up)
  sigma.down    = sd(x.down)
  sigma.neutral = sd(x.neutral)
  
  if (sigma.all == 0) {
    stop("The gene is not expressed in all patients or has the same 
          expression values in all patients")
  }
  
  if (mu[2] - mu[1] < sigma.all) {
    lambda[1] = 0.01
  }
  if (mu[3] - mu[2] < sigma.all) {
    lambda[3] = 0.01
  }
  lambda = lambda / sum(lambda)
  
  
  min.reads = 5
  if (num.neutral < min.reads || sigma.neutral==0) 
    sigma.neutral = sigma.all
  
  if (num.up < min.reads) 
    sigma.up = sigma.neutral
  
  if (num.down < min.reads) 
    sigma.down = sigma.neutral
  
  sigma = c(sigma.down, sigma.neutral, sigma.up)
  
  # The model
  model = list()
  model$lambda = lambda
  model$mu     = mu
  model$sigma  = sigma
  
  return(model)
}


#==============================================================================
# Detecting outliers based on the
# boxplot rule, the ideal fourths are used to estimate the quartiles
#
# Ref: Frigge, et al. "Some implementations of the Boxplot"
#
DetectBoxplotOutlier = function(x, k=1.5) {
  x = x[!is.na(x)]
  N = length(x)
  
  ideal.fourth = ComputeIdealFourth(x)
  range = k * (ideal.fourth$qu - ideal.fourth$ql)
  cl = ideal.fourth$ql - range
  cu = ideal.fourth$qu + range
  
  outlier.id    = which(x < cl | x > cu)
  keep.id       = setdiff(seq(1:N), outlier.id)
  outlier.value = x[outlier.id]
  
  list(outlier.value=outlier.value, 
       outlier.id=outlier.id, 
       keep.id=keep.id, 
       cl=cl, 
       cu=cu)
}

#==============================================================================
# Compute the ideal (machine) fourths of x
#
# Frigge. "Resistant outlier rules and the non-Gaussian case"
# Frigge, et al. "Some implementations of the Boxplot"
#
ComputeIdealFourth = function(x, na.rm=FALSE) {
  if (na.rm == TRUE) {
    x = x[!is.na(x)]
  }
  x = sort(x)
  N = length(x)
  
  j  = floor(N / 4 + 5 / 12)
  g  = N / 4 - j + 5 / 12
  ql = (1 - g) * x[j] + g * x[j + 1]
  
  k  = N - j + 1
  qu = (1 - g) * x[k] + g * x[k - 1]
  
  list(ql=ql, qu=qu)
}
