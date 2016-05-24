my_ascore_method <-
function(U, V)
{
  # Internal function for ASCORE method
  # get inverse of V
  V.eigen = eigen(V)
  vals = V.eigen$values
  invals = ifelse(abs(vals) > 1e-9, 1/vals, 1e-9)
  M = t(V.eigen$vectors)
  V.sqr = t(M) %*% diag(sqrt(invals)) %*% M
  
  # standardize U
  Ustd = V.sqr %*% U
  Ustd2 = Ustd^2
  # Ustd in decreasing order
  Uord = Ustd[order(abs(Ustd), decreasing=TRUE)]
  Uord2 = Uord^2
  
  ## get scores
  k = length(U) 
  score <- score.ords <- rep(0, length(U))
  for (j in 1:k)
  {
    score[j] = sum(Ustd2[1:j] - 1) / sqrt(2*j)
    score.ords[j] = sum(Uord2[1:j] - 1) / sqrt(2*j)
  }
  S1 = max(score)
  S2 = max(score.ords)
  c(S1, S2)
}

my_assu_method <-
function(U, V)
{
  # Internal function for ASSU method
  # get inverse of V
  V.eigen = eigen(V)
  vals = V.eigen$values
  invals = ifelse(abs(vals) > 1e-9, 1/vals, 1e-9)
  M = t(V.eigen$vectors)
  V.inv = t(M) %*% diag(invals) %*% M
  V.sqr = t(M) %*% diag(sqrt(invals)) %*% M
  
  # standardize U
  Ustd = V.sqr %*% U
  # Ustd in decreasing order
  Uord = Ustd[order(abs(Ustd), decreasing=TRUE)]
  # decreasing order of U
  Ord = order(U^2, decreasing=TRUE)
  # ord cov matrix
  V.ord = V[Ord, Ord]
  
  ## get scores
  k = length(U) 
  scores <- scores.ord <- rep(0, k)
  pvals <- pvals.ord <- rep(0, k)
  for (j in 1:k)
  {
    aux = my_ssu_method(U[1:j], V[1:j,1:j])
    scores[j] = aux[1]
    pvals[j] = aux[2]
    aux.ord = my_ssu_method(Uord[1:j], V.ord[1:j,1:j])
    scores.ord[j] = aux.ord[1]
    pvals.ord[j] = aux.ord[2]
  }
  p1 = min(pvals)
  p2 = min(pvals.ord)
  S1 = scores[pvals==p1]
  S2 = scores.ord[pvals.ord==p2] 
  c(S1, p1, S2, p2)
}

my_assuw_method <-
function(U, V)
{
  # Internal function for ASSUW method
  # decreasing order of U
  Ordw = order(U^2/diag(V), decreasing=TRUE)
  # Ustd in decreasing order
  Uordw = U[Ordw]
  # ord cov matrix
  V.ordw = V[Ordw, Ordw]
  
  ## get scores
  k = length(U) 
  scores <- scores.ord <- rep(0, k)
  pvals <- pvals.ord <- rep(0, k)
  for (j in 1:k)
  {
    aux = my_ssuw_method(U[1:j], V[1:j,1:j])
    scores[j] = aux[1]
    pvals[j] = aux[2]
    aux.ord = my_ssuw_method(Uordw[1:j], V.ordw[1:j,1:j])
    scores.ord[j] = aux.ord[1]
    pvals.ord[j] = aux.ord[2]
  }
  p1 = min(pvals)
  p2 = min(pvals.ord)
  S1 = scores[pvals==p1]
  S2 = scores.ord[pvals.ord==p2] 
  c(S1, p1, S2, p2)
}

my_asum_method <-
function(U, V)
{
  # Internal function for ASUM method
  # decreasing order of U
  Ords = order(U/sqrt(diag(V)), decreasing=TRUE)
  # Ustd in decreasing order
  Uords = U[Ords]
  # ord cov matrix
  V.ords = V[Ords, Ords]
  
  ## get scores
  k = length(U) 
  scores <- scores.ord <- rep(0, k)
  pvals <- pvals.ord <- rep(0, k)
  for (j in 1:k)
  {
    aux = my_sum_method(U[1:j], V[1:j,1:j])
    scores[j] = aux[1]
    pvals[j] = aux[2]
    aux.ord = my_sum_method(Uords[1:j], V.ords[1:j,1:j])
    scores.ord[j] = aux.ord[1]
    pvals.ord[j] = aux.ord[2]
  }
  p1 = min(pvals)
  p2 = min(pvals.ord)
  S1 = scores[pvals==p1]
  S2 = scores.ord[pvals.ord==p2] 
  c(S1, p1, S2, p2)
}

my_bst_method <-
function(U, V)
{
  # Internal function for BST method
  # score statistic
  if (is.null(dim(V))) 
  {
    score = 0.5 * (sum(U*U) - V)
    if (is.na(score) || is.infinite(score) || is.nan(score))
      score = 0
  } else {
    # goeman statistic (unknown distribution)
    score = 0.5 * (sum(U*U) - sum(diag(V)))
  }
  score
}

my_calpha_method <-
function(casecon, gen)
{
  # Internal function for BST method
  nA = sum(casecon)
  nU = sum(casecon==0)
  p0 = nA / (nA + nU)
  
  m = ncol(gen)
  # copies of the i-th variant type
  n = apply(gen, 2, function(x) sum(x>0, na.rm=TRUE))
  # copies of the i-th variant type in the cases
  g = apply(gen[casecon==1,], 2, function(x) sum(x>0, na.rm=TRUE))
  # Test statistic 
  Talpha = sum((g - n*p0)^2 - (n * p0 * (1-p0)))
  Talpha
}

my_cmat_method <-
function(casecon, gen, weights)
{
  # Internal function for CMAT method
  nA = sum(casecon)
  nU = length(casecon) - nA
  
  if (is.vector(gen))
  {
    # weighted minor-allele counts in cases 'A' and controls 'U'
    mA = sum(gen[casecon==1] * weights, na.rm=TRUE)
    mU = sum(gen[casecon==0] * weights, na.rm=TRUE)
    # weighted major-allele counts in cases 'A' and controls 'U'
    MA = sum((2 - gen[casecon==1]) * weights, na.rm=TRUE)
    MU = sum((2 - gen[casecon==0]) * weights, na.rm=TRUE)    	
  } else {
    # matrix of weights
    W = diag(weights)	
    # weighted minor-allele counts in cases 'A' and controls 'U'
    mA = sum(gen[casecon==1,] %*% W, na.rm=TRUE)
    mU = sum(gen[casecon==0,] %*% W, na.rm=TRUE)
    # weighted major-allele counts in cases 'A' and controls 'U'
    MA = sum((2 - gen[casecon==1,]) %*% W, na.rm=TRUE)
    MU = sum((2 - gen[casecon==0,]) %*% W, na.rm=TRUE)    	
  }
  # CMAT statistic	
  cmat1 = (nA + nU) / (2 * nA * nU * sum(weights))
  cmat2 = ((mA * MU - mU * MA)^2) / ((mA + mU) * (MA + MU))
  cmat.stat = cmat1 * cmat2
  cmat.stat
}

my_cmc_method <-
function(casecon, X.new)
{
  # Internal function for CMAT method
  ## number of individuals N, cases nA, controls nU
  N = nrow(X.new)
  nA = sum(casecon)
  nU = N - nA
  ## matrix of genotypes in cases
  Xx = X.new[casecon==1,]  
  ## matrix of genotypes in controls  
  Yy = X.new[casecon==0,] 
  ## get means
  Xx.mean = colMeans(Xx, na.rm=TRUE)
  Yy.mean = colMeans(Yy, na.rm=TRUE)
  ## center matrices Xx and Yy
  Dx = sweep(Xx, 2, Xx.mean)
  Dy = sweep(Yy, 2, Yy.mean)
  
  ## pooled covariance matrix
  if (sum(complete.cases(X.new)) == N)  # no missing values
  {  
    COV = (t(Dx) %*% Dx + t(Dy) %*% Dy) / (N-2)
  } else {  # with missing values
    ## covariance matrix of cases
    tDx = t(Dx)
    Sx = matrix(0, ncol(X.new), ncol(X.new))
    for (i in 1:nrow(tDx))
    {
      for (j in i:ncol(Dx))
      {
        Sx[i,j] = sum(tDx[i,] * Dx[,j], na.rm=TRUE)
      }
    }
    sx.diag = diag(Sx)
    Sx = Sx + t(Sx)
    diag(Sx) = sx.diag
    ## covariance matrix of controls
    tDy = t(Dy)
    Sy = matrix(0, ncol(X.new), ncol(X.new))
    for (i in 1:nrow(tDy))
    {
      for (j in i:ncol(Dy))
      {
        Sy[i,j] = sum(tDy[i,] * Dy[,j], na.rm=TRUE)
      }
    }
    sy.diag = diag(Sy)
    Sy = Sy + t(Sy)
    diag(Sy) = sy.diag
    ## pooled covariance matrix
    COV = (1/(N-2)) * (Sx + Sy)	
  }
  
  ## general inverse
  if (nrow(COV) == 1) # only one variant
  { 
    if (COV < 1e-8) COV = 1e-8
    COV.inv = 1 / COV
  } else {
    COV.eigen = eigen(COV)
    eig.vals = COV.eigen$values  
    inv.vals = ifelse(abs(eig.vals) <= 1e-8, 0, 1/eig.vals)
    EV = solve(COV.eigen$vectors)
    COV.inv = t(EV) %*% diag(inv.vals) %*% EV
  }	
  
  ## Hotellings T2 statistic
  stat = t(Xx.mean - Yy.mean) %*% COV.inv %*% (Xx.mean - Yy.mean) * nA * nU / N
  as.numeric(stat)
}

my_gdbr_fstat <- 
function(casecon, G)
{
  # Internal function for GDBR method
  # center y
  y.new = casecon - mean(casecon)
  I = diag(1, length(casecon))
  # get projection 'hat' matrix
  H = y.new %*% solve((t(y.new) %*% y.new)) %*% t(y.new)
  # calculate F statistic
  Fstat.num = sum(diag(H %*% G %*% H))
  Fstat.denom = sum(diag((I - H) %*% G %*% t(I - H)))
  Fstat = Fstat.num / Fstat.denom	
  Fstat
}

my_getUV <-
function(y, X)
{
  # Internal function for getUV
  # center phenotype y
  y.new = y - mean(y)
  # get score vector
  U = colSums(y.new * X, na.rm=TRUE)
  # get covariance matrix
  X.new = scale(X, scale=FALSE)
  if (sum(complete.cases(X)) != length(y))  # missing data
  {
    tX.new = t(X.new)
    Sx = matrix(0, ncol(X), ncol(X))
    for (i in 1:nrow(tX.new))
    {
      for (j in i:ncol(X.new))
      {
        Sx[i,j] = sum(tX.new[i,] * X.new[,j], na.rm=TRUE)
      }
    }
    sx.diag = diag(Sx)
    Sx = Sx + t(Sx)
    diag(Sx) = sx.diag
    V = mean(y) * (1 - mean(y)) * Sx
  } else {  # no missing data
    V = mean(y) * (1 - mean(y)) * (t(X.new) %*% X.new)
  }
  # results
  res.uv = list(U=U, V=V)
  return(res.uv)
}

my_orwss_method <-
function(casecon, gen, cs)
{	
  # Internal function for ORWSS method
  ## num of variants
  p = ncol(gen)
  ## calculate amended log(OR)
  gama = rep(0, p)
  for (j in 1:p)
  {
    tj = table(casecon, gen[,j]) + 0.5	
    # log odds ratio
    gama[j] = log((tj[1,1]*tj[2,2]) / (tj[1,2]*tj[2,1]))
  }
  w = gama
  if (!is.null(cs))
    w[abs(gama - mean(gama)) <= cs * sd(gama)] = 0
  ## calculate genetic score	
  score = rowSums(gen %*% diag(w), na.rm=TRUE) 
  rank.score = order(score)
  
  ## sum of ranks of cases
  x = sum(rank.score[casecon==1])
  x
}

my_rarecov_method <- 
function(casecon, gen, dif)
{
  # Internal function for RARECOV method
  # rare cover algorithm
  M = ncol(gen)
  selected = NULL
  temp.stats = rep(0, M)
  temp.pvals = rep(0, M)
  
  ## get the first variant
  for (j in 1:M)
  {
    temp = chisq.test(table(gen[,j], casecon))
    temp.stats[j] = temp$statistic
    temp.pvals[j] = temp$p.value
  }
  win.j = which(temp.stats == max(temp.stats))
  selected = c(selected, win.j[1])
  Xcorr = temp.stats[win.j[1]]
  Xpval = temp.pvals[win.j[1]]
  rest = setdiff(1:M, selected)
  
  ## get rest of variants
  Xcorr.dif = 1
  while (Xcorr.dif > dif)
  {   
    temp.stats = rep(0, M) 
    temp.pvals = rep(0, M)
    for (j in 1:length(rest))
    {
      gen.new = rowSums(gen[,c(selected,rest[j])], na.rm=TRUE)
      gen.new[gen.new != 0] = 1
      temp = chisq.test(table(gen.new, casecon))
      temp.stats[rest[j]] = temp$statistic
      temp.pvals[rest[j]] = temp$p.value
    }
    win.j = which(temp.stats == max(temp.stats))
    Xcorr.new = temp.stats[win.j[1]]
    Xpval.new = temp.pvals[win.j[1]]
    Xcorr.dif = Xcorr.new - Xcorr
    if (Xcorr.dif > dif)
    {
      selected = c(selected, win.j)
      Xcorr = Xcorr.new
      Xpval = Xpval.new
      rest = setdiff(1:M, selected)
    }
  }
  list(stat=Xcorr, pval=Xpval, sel=selected)
}

my_rbt_method <-
function(casecon, gen)
{
  # Internal function for RBT method
  ## num of counts of rare allele in cases A and controls U
  fA = colSums(gen[casecon==1,], na.rm=TRUE)
  fU = colSums(gen[casecon==0,], na.rm=TRUE)
  
  ## RBT S+ (enrichment of mutations in cases)
  ## get those A larger than U, and sort them
  AltU = fA > fU  # get those A larger than U
  if (sum(AltU) > 0)
  {
    kA.plus = sort(fA[AltU])
    kU.plus = fU[AltU][order(fA[AltU])]
    ## how many unique A > U
    uniqA.plus = unique(sort(fA[AltU]))
    uniqU.plus = unique(sort(fU[AltU]))
    if (length(uniqA.plus) == 1)
    {
      fua = (uniqU.plus + uniqA.plus) / 2
      pU = ppois(uniqU.plus, fua)
      pA = 1 - ppois(uniqA.plus - 1, fua) 
      PUA.plus = -log(pU * pA)
      NUA.plus = 1
    } else {
      ## matrix of observed counts
      ncolsA.plus = length(uniqA.plus)
      nrowsU.plus = length(uniqU.plus)		
      NUA.plus = matrix(0, nrowsU.plus, ncolsA.plus)
      
      for (j in 1:ncolsA.plus) {
        for (i in kU.plus[kA.plus==uniqA.plus[j]]) {
          NUA.plus[which(uniqU.plus==i),j] = sum(kU.plus[kA.plus==uniqA.plus[j]] == i)
        }
      }
      # dimnames(NUA.plus) = list(uniqU.plus, uniqA.plus)
      
      # matrix of poisson probabilities
      PUA.plus = matrix(0, nrowsU.plus, ncolsA.plus)
      for (i in 1:nrowsU.plus) {
        for (j in 1:ncolsA.plus) {
          fua = (uniqU.plus[i] + uniqA.plus[j]) / 2
          pU = ppois(uniqU.plus[i], fua)
          pA = 1 - ppois(uniqA.plus[j]-1, fua) 
          PUA.plus[i,j] = -log(pU * pA)
        }
      }
      # dimnames(PUA.plus) = list(uniqU.plus, uniqA.plus)
    }
  } else { # sum(AltU) <= 0
    NUA.plus = 0
    PUA.plus = 0
  } 
  
  ## RBT S- (enrichment of mutations in controls)
  ## get those A smaller than U, and sort them
  AstU = fA < fU  # get those A smaller than U
  if (sum(AstU) > 0)
  {
    kA.minus = sort(fA[AstU])
    kU.minus = fU[AstU][order(fA[AstU])]
    ## how many unique A > U
    uniqA.minus = unique(sort(fA[AstU]))
    uniqU.minus = unique(sort(fU[AstU]))
    if (length(uniqA.minus) == 1)
    {
      fua = (uniqU.minus + uniqA.minus) / 2
      pU = ppois(uniqU.minus, fua)
      pA = 1 - ppois(uniqA.minus - 1, fua) 
      PUA.minus = -log(pU * pA)
      NUA.minus = 1
    } else {				
      ## matrix of observed counts
      ncolsA.minus = length(uniqA.minus)
      nrowsU.minus = length(uniqU.minus)
      NUA.minus = matrix(0, nrowsU.minus, ncolsA.minus)
      
      for (j in 1:ncolsA.minus) {
        for (i in kU.minus[kA.minus==uniqA.minus[j]]) {
          NUA.minus[which(uniqU.minus==i),j] = sum(kU.minus[kA.minus==uniqA.minus[j]] == i)
        }
      }
      # dimnames(NUA.minus) = list(uniqU.minus, uniqA.minus)
      
      # matrix of poisson probabilities
      PUA.minus = matrix(0, nrowsU.minus, ncolsA.minus)
      for (i in 1:nrowsU.minus) {
        for (j in 1:ncolsA.minus) {
          fua = (uniqU.minus[i] + uniqA.minus[j]) / 2
          pU = ppois(uniqU.minus[i], fua)
          pA = 1 - ppois(uniqA.minus[j]-1, fua) 
          PUA.minus[i,j] = -log(pU * pA)
        }
      }
      # dimnames(PUA.minus) = list(uniqU.minus, uniqA.minus)
    }
  } else {
    NUA.minus = 0
    PUA.minus = 0
  }
  
  ## RBT S+ (enrichment of mutations in cases)
  rbt.plus = sum(NUA.plus * PUA.plus)
  ## RBT S- (enrichment of mutations in controls)
  rbt.minus = sum(NUA.minus * PUA.minus)
  ## RBT statistic	
  stat = max(rbt.plus, rbt.minus)
  stat
}



my_rvt1_method <- 
function(casecon, gen, ymean)
{
  # Internal function for RVT1 method
  ## casecon is already centered
  ## gen is a matrix with recoded RVs
  
  ## get proportion of rare variants
  if (is.vector(gen)) {
    x.prop = mean(gen, na.rm=TRUE)
  } else {
    x.prop = rowMeans(gen, na.rm=TRUE)
  }
  # get score vector U
  U = sum(casecon * x.prop)
  ## V
  xv = sum((x.prop - mean(x.prop))^2)
  V = ymean * (1 - ymean) * xv
  ## get score
  score = sum(U^2 / V)
  score
}

my_rvt2_method <- 
function(casecon, gen, ymean)
{
  # Internal function for RVT2 method
  ## casecon is already centered
  ## gen is a matrix with recoded RVs
  
  ## create indicator dummy
  if (is.vector(gen)) {
    x.ind = gen
  } else {
    x.ind = rowSums(gen, na.rm=TRUE)
  }
  x.ind[x.ind != 0] = 1
  # get score vector U
  U = sum(casecon * x.ind, na.rm=TRUE)
  ## V
  xv = sum((x.ind - mean(x.ind))^2)
  V = ymean * (1 - ymean) * xv
  ## get score
  score = sum(U^2 / V)
  if (is.na(score) || is.infinite(score) || is.nan(score))
    score = 0
  return(score)
}

my_rwas_method <- function(casecon, gen, weights)
{
  # Internal function for RWAS method
  # number of individuals and variants
  N.pop = nrow(gen)
  
  if (is.vector(gen))
  {
    p.pop = mean(gen, na.rm=TRUE) / 2
    N.cases = sum(gen[casecon==1], na.rm=TRUE)
    N.controls = sum(gen[casecon==0], na.rm=TRUE)
    p.cases = mean(gen[casecon==1], na.rm=TRUE) / 2
    p.controls = mean(gen[casecon==0], na.rm=TRUE) / 2
  }
  if (is.matrix(gen))
  {
    # calculation of maf 
    p.pop = colSums(gen) / (2*N.pop)
    if (all(!is.na(p.pop))) # no missing values
    {
      N.cases = sum(casecon)
      N.controls = sum(casecon==0)
      p.cases = colSums(gen[casecon==1,]) / (2 * N.cases)
      p.controls = colSums(gen[casecon==0,]) / (2 * N.controls)
    } else { # with missing data
      N.pop = apply(gen, 2, function(x) sum(!is.na(x)))
      N.cases = apply(gen[casecon==1,], 2, function(x) sum(!is.na(x)))
      N.controls = apply(gen[casecon==0,], 2, function(x) sum(!is.na(x)))
      p.pop = colSums(gen, na.rm=TRUE) / (2 * N.pop)
      p.cases = colSums(gen[casecon==1,], na.rm=TRUE) / (2 * N.cases)
      p.controls = colSums(gen[casecon==0,], na.rm=TRUE) / (2 * N.controls)
    }
  }
  
  # average minor allele frequencies (cases and controls)
  p.total = p.cases*(N.cases/N.pop) + p.controls*(N.controls/N.pop) 
  #p.total = (p.cases + p.controls) / 2 
  
  # calculation of z-scores and weights
  #z.scores = (p.cases - p.controls) / (sqrt(2/N)*sqrt(p.total*(1-p.total)))
  z.num = p.cases - p.controls
  z.denom = sqrt(N.pop/(2 * N.cases * N.controls)) * sqrt(p.total * (1-p.total))
  z.scores = z.num / z.denom
  w = sqrt((1 - p.pop) / p.pop)
  
  # calculation of RWAS
  rwas = sum(weights * w * z.scores) / sqrt(sum(w^2))
  rwas
}

my_score_method <-
function(U, V)
{
  # score statistic and p-value
  if (is.null(dim(V))) 
  {
    score = sum(U^2 / V)
    if (is.na(score) || is.infinite(score) || is.nan(score))
      score = 0
    pval = 1 - pchisq(score, 1)
  } else {
    # get inverse of V
    V.eigen = eigen(V)
    vals = V.eigen$values
    V.rank = sum(abs(vals) > 1e-9)
    invals = ifelse(abs(vals) > 1e-9, 1/vals, 0)
    M = solve(V.eigen$vectors)
    V.inv = t(M) %*% diag(invals) %*% M
    score = t(U) %*% V.inv %*% U
    pval = 1 - pchisq(score, df=V.rank)
  }
  res = c(score, pval)
  res
}

my_ssu_method <-
function(U, V)
{
  # score statistic and p-value
  if (is.null(dim(V)))
  {
    score = sum(U^2 / V)
    if (is.na(score) || is.infinite(score) || is.nan(score))
      score = 0
    pval = 1 - pchisq(score, 1)
  } else {
    if (all(abs(U) < 1e-20)) {
      score = 0
      pval = 1
    } else {
      # get diag of V
      diagV = diag(V)
      diagV = ifelse(diagV > 1e-9, diagV, 1e-9)
      score = as.numeric(t(U) %*% U)
      # distrib of score 
      V.eigen = eigen(V)
      vals = V.eigen$values
      a = sum(vals^3) / sum(vals^2)
      b = sum(vals) - (sum(vals*vals)^2 / sum(vals^3))
      d = (sum(vals^2)^3) / (sum(vals^3)^2)
      pval = 1 - pchisq((score - b)/a, df=d)
    }
  }
  res = c(score, pval)
  res   
}

my_ssuw_method <-
function(U, V)
{
  # score statistic and p-value
  if (is.null(dim(V)))
  {
    score = sum(U^2 / V)
    if (is.na(score) || is.infinite(score) || is.nan(score))
      score = 0
    pval = 1 - pchisq(score, 1)
  } else {
    if (all(abs(U) < 1e-20)) {
      score = 0
      pval = 1
    } else {
      # get diag of V
      diagV = diag(V)
      diagV = ifelse(diagV > 1e-9, diagV, 1e-9)
      score = sum(U^2 / diagV)
      # distrib of score 
      V.eigen = eigen(V %*% diag(1/diagV))
      vals = V.eigen$values
      a = sum(vals^3) / sum(vals^2)
      b = sum(vals) - (sum(vals*vals)^2 / sum(vals^3))
      d = (sum(vals^2)^3) / (sum(vals^3)^2)
      pval = 1 - pchisq((score - b)/a, df=d)
    }
  }
  res = c(score, pval)
  res   
}

my_sum_method <-
function(U, V)
{
  # score statistic and p-value
  if (abs(sum(U)) < 1e-20) {
    score = 0
    pval = 1
  } else {
    ones = rep(1, length(U))
    score = sum(U) / sqrt(sum(V))
    pval = 1 - pchisq(score^2, 1)
  }
  res = c(score, pval)
  res
}

my_uminp_method <-
function(U, V)
{
  ## Requires library "mvtnorm"!!!
  # score statistic and p-value
  if (is.null(dim(V)))
  {
    score = sum(U^2 / V)
    if (is.na(score) || is.infinite(score) || is.nan(score))
      score = 0
    pval = 1 - pchisq(score, 1)
  } else {
    score = abs(U) / (sqrt(diag(V)) + 1e-20)
    M = length(U)
    Vnew = matrix(0, M, M)
    for (i in 1:M)
    {
      for (j in 1:M)
      {
        if(abs(V[i,j] > 1e-20))
          Vnew[i,j] = V[i,j] / sqrt(V[i,i] * V[j,j])
        else Vnew[i,j] = 1e-20
      }
    }
    nd = nrow(Vnew)
    q = max(score)
    pval = as.numeric(1 - pmvnorm(lower=c(rep(-q,nd)), upper=c(rep(q,nd)), mean=c(rep(0,nd)), sigma=Vnew))
  }
  res = c(score, pval)
  names(res) = NULL
  res
}

my_uni_score <- 
function(y, x)
{
  # univariate score statistic 
  u = sum((y - mean(y)) * (x - mean(x)))
  v = mean(y) * (1 - mean(y)) * sum((x - mean(x))^2)
  if (abs(sum(u)) < 1e-20) {
    score = 0
  } else {
    score = (u^2) / v
  }
  if (is.na(score) || is.infinite(score) || is.nan(score))
    score(score)
  return(score)
}

my_vt_method <-
function(casecon, gen, mafs, h.maf)
{
  # mafs: minor allele frequencies
  # h.maf: unique mafs
  z.scores = rep(0, length(h.maf)-1)
  y.new = casecon - mean(casecon)
  for (i in 1:(length(h.maf)-1))
  {
    z.num = sum(gen[,mafs<h.maf[i+1]] * y.new, na.rm=TRUE)
    z.denom = sqrt(sum((gen[,mafs<h.maf[i+1]])^2, na.rm=TRUE))
    z.scores[i] = z.num / z.denom
  }
  stat = max(z.scores)
  stat
}

my_weights_wss <- 
function(casecon, gen)
{
  controls = !casecon
  n.loc = ncol(gen)
  ## calculate weights
  w <- rep(0, n.loc)
  for (j in 1:n.loc)
  {
    nNA = is.na(gen[,j])
    mU = sum(gen[controls,j], na.rm=TRUE)
    nU = sum(controls[!nNA])		
    q = (mU+1) / (2*nU+2)
    n = sum(!nNA)
    w[j] = 1 / sqrt(n * q * (1-q))
  }
  w
}

my_wss_method <-
function(casecon, gen)
{	
  controls = !casecon
  n.loc = ncol(gen)
  
  ## calculate weights
  w <- rep(0, n.loc)
  for (j in 1:n.loc)
  {
    nNA = is.na(gen[,j])
    mU = sum(gen[controls,j], na.rm=TRUE)
    nU = sum(controls[!nNA])		
    q = (mU+1) / (2*nU+2)
    n = sum(!nNA)
    w[j] = sqrt(n * q * (1-q))
  }
  
  ## calculate genetic score	
  score = rowSums(gen %*% diag(1/w), na.rm=TRUE) 
  # rank.score = order(score)
  rank.score = rank(score)
  
  ## sum of ranks of cases
  x = sum(rank.score[casecon==1])
  return(x)
}

my_wst_method <-
function(U, V)
{
  ## weights
  k = 1:length(U)
  w = (1 / (k + 1))^2
  ## statistic
  tw.num = as.numeric(t(w) %*% U)
  tw.denom = as.numeric(t(w) %*% V %*% w)
  stat = tw.num / tw.denom
  stat
}

my_weights_wss <- 
function(casecon, gen)
{
  # internal function for CARV
  # weights a la madsen & browning
  controls = !casecon
  n.loc = ncol(gen)
  ## calculate weights
  w <- rep(0, n.loc)
  for (j in 1:n.loc)
  {
    nNA = is.na(gen[,j])
    mU = sum(gen[controls,j], na.rm=TRUE)
    nU = sum(controls[!nNA])		
    q = (mU+1) / (2*nU+2)
    n = sum(!nNA)
    w[j] = 1 / sqrt(n * q * (1-q))
  }
  w
}


my_score_carv <- 
function(casecon, X.cen, w)
{
  # internal function for CARV
  if (length(w) == 1) 
  {
    U = casecon * (X.cen * w)
    ## statistic
    sco.num = (sum(U, na.rm=TRUE))^2
    sco.denom = sum(U^2, na.rm=TRUE)
    score = sco.num / sco.denom 
  }
  if (length(w) >  1)
  {
    # get score vector U
    U = casecon * (X.cen  %*% diag(w))
    ## statistic
    sco.num = (sum(U, na.rm=TRUE))^2
    sco.denom = sum((rowSums(U, na.rm=TRUE))^2)
    score = sco.num / sco.denom 
  }   
  score	
}


my_hard_approach <- 
function(y, X, w, MAFs, maf)
{
  # internal function for CARV
  # y and X are already centered
  # w contains ak * sk
  # MAFs contains minor allele frequencies
  # maf contains hard maf
  
  # get only rare variants below maf
  w = w[MAFs < maf]
  X.new = X[, MAFs < maf] 	
  stat = my_score_carv(y, X.new, w)
  stat
}


my_variable_approach <- 
function(y, X, w, MAFs, sort.maf)
{
  # internal function for CARV
  # y and X are already centered
  # w contains ak * sk
  # MAFs contains minor allele frequencies
  # sort.maf contains sorted unique MAFs
  
  if (length(sort.maf) > 1)
    sort.maf = sort.maf[-1]
  scores = rep(0, length(sort.maf))
  for (i in 1:length(sort.maf))
  {
    # get only rare variants
    w.new = w[MAFs < sort.maf[i]]
    X.new = X[,MAFs < sort.maf[i]]
    # get score
    scores[i] = my_score_carv(y, X.new, w.new)
  }	
  max(scores)[1]
}


my_stepup_approach <- 
function(y, X, w)
{
  # internal function for CARV
  # variable selection by step-up 
  # y and X are already centered
  # w contains ak * sk
  
  # initial parameters
  M = ncol(X)
  
  # compute univariate test statistic for each variant
  temp.scores = rep(0, M)
  for (j in 1:M) {
    temp.scores[j] = my_score_carv(y, X[,j], w[j])
  }
  # who is the best
  best.score = max(temp.scores)[1]
  selected = which(temp.scores == best.score)[1]
  # leftover variants
  rest = setdiff(1:M, selected)
  
  # build on the model with leftover variants
  decision = 1
  while (decision > 0)
  {
    temp.scores = rep(0, length(rest))
    for (j in 1:length(rest))
    {
      temp = c(selected, rest[j])
      temp.scores[j] = my_score_carv(y, X[,temp], w[temp])
    }
    candidate.score = max(temp.scores)[1]
    candidate = which(temp.scores == candidate.score)[1]
    # compare with previous statistic
    decision = candidate.score - best.score
    if (decision > 0) 
    {
      best.score = candidate.score
      selected = c(selected, candidate)
      rest = setdiff(1:M, selected)
    } #else break 
  }
  best.score
}

