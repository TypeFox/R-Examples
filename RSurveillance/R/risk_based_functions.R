############################################################################
# risk-based surveillance
###########################################################################
# adj.risk
# epi.calc
# sep.rb.bin
# sep.rb.hypergeo
# sep.rb.bin.varse
# sep.rb.hypergeo.varse
# sep.rb2.bin
# sep.rb2.hypergeo
# sse.rb.2stage
# sse.combined
# n.rb
# n.rb.varse

# include freedom functions

##' Adjusted risk
##' @description Calculates adjusted risk for given 
##'   relative risk and population proportions. This is an intermediate calculation
##'   in the calculation of effective probability of infection for risk-based 
##'   surveillance activities
##' @param rr relative risk values (vector of values corresponding to the number of risk strata)
##' @param ppr population proportions corresponding to 
##'   rr values (vector of equal length to rr)
##' @return vector of adjusted risk values (in order corresponding to rr)
##' @keywords methods
##' @export
##' @examples 
##' # examples for adj.risk
##' adj.risk(c(5, 1), c(0.1, 0.9))
##' adj.risk(c(5, 3, 1), c(0.1, 0.1, 0.8))
adj.risk<- function(rr, ppr) {
  sum.prod<- sum(rr*ppr)
  ar<- rr/sum.prod
  return(ar)
}


##' Effective probability of infection (EPI)
##' @description Calculates effective probability of infection (adjusted design prevalence)
##'   for each risk group for risk-based surveillance activities
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector of values corresponding to 
##'   the number of risk strata)
##' @param ppr population proportions corresponding to rr values 
##'   (vector of equal length to rr)
##' @return list of 2 elements, a vector of EPI values and a vector of corresponding
##'   adjusted risks (in corresponding order to rr)
##' @keywords methods
##' @export
##' @examples 
##' # examples for epi.calc
##' epi.calc(0.1, c(5, 1), c(0.1, 0.9))
##' epi.calc(0.02, c(5, 3, 1), c(0.1, 0.1, 0.8))
epi.calc<- function(pstar, rr, ppr) {
  ar<- adj.risk(rr, ppr)
  epi<- pstar*ar
  return(list(epi=epi, adj.risk=ar))
}


##' Binomial risk-based population sensitivity
##' @description Calculates risk-based population sensitivity with a 
##'   single risk factor, using binomial method (assumes a large population),
##'   allows for unit sensitivity to vary among risk strata
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector of values corresponding to the number of risk strata)
##' @param ppr population proportions corresponding to rr values 
##'   (vector of equal length to rr)
##' @param n sample size per risk category (vector same length as 
##'   rr and ppr)
##' @param se unit sensitivity, can vary among risk strata (fixed value or 
##'   vector same length as rr, ppr, n)
##' @return list of 3 elements, a scalar of population-level sensitivity
##'   a vector of EPI values and a vector of corresponding adjusted risks
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.rb.bin
##' sep.rb.bin(0.1, c(5, 3, 1), c(0.1, 0.1, 0.8), c(5, 5, 5), 0.9)
##' sep.rb.bin(0.1, c(5, 1), c(0.1, 0.9), c(10, 5), c(0.95, 0.9))
##' sep.rb.bin(0.1, c(5, 1), c(0.1, 0.9), c(10, 5), c(0.9, 0.9))
##' sep.rb.bin(0.01, c(5, 1), c(0.1, 0.9), c(90, 50), c(0.9, 0.9))
sep.rb.bin<- function(pstar, rr, ppr, n, se) {
  epi<- epi.calc(pstar, rr, ppr)
  p.all.neg<- (1 - se*epi[[1]])^n
  sep<- 1 - prod(p.all.neg)
  return(list(sep=sep, epi=epi[[1]], adj.risk=epi[[2]]))
}


##' Hypergeometric risk-based population sensitivity
##' @description Calculates risk-based population sensitivity with a 
##'   single risk factor, using the hypergeometric method 
##'   (assuming a finite and known population size),
##'   allows for unit sensitivity to vary among risk strata
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector of values corresponding 
##'   to the number of risk strata)
##' @param n sample size per risk category (vector same length as 
##'   rr and ppr)
##' @param N Population size per risk category (vector same length 
##'   as rr and ppr)
##' @param se unit sensitivity, can vary among risk strata (fixed value or a vector the same 
##'   length as rr, ppr, n)
##' @return list of 3 elements, a scalar of population-level sensitivity
##'   a vector of EPI values and a vector of corresponding adjusted risks
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.rb.bin
##' sep.rb.hypergeo(0.1, c(5, 3, 1), c(10, 10, 80), c(5, 5, 5), 0.9)
##' sep.rb.hypergeo(0.1, c(5, 1), c(15, 140), c(10, 5), c(0.95, 0.9))
##' sep.rb.hypergeo(0.1, c(5, 1), c(23, 180), c(10, 5), c(0.9, 0.9))
##' sep.rb.hypergeo(0.01, c(5, 1), c(100, 900), c(90, 50), c(0.9, 0.9))
sep.rb.hypergeo<- function(pstar, rr, N, n, se) {
  ppr<- N/sum(N)
  epi<- epi.calc(pstar, rr, ppr)
  p.all.neg<- (1 - se*n/N)^(epi[[1]]*N)
  sep<- 1 - prod(p.all.neg)
  return(list(sep=sep, epi=epi[[1]], adj.risk=epi[[2]]))
}


##' Binomial risk-based population sensitivity for varying unit sensitivity
##' @description Calculates population sensitivity for a single risk factor 
##'   and varying unit sensitivity using binomial method (assumes large population)
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector of values corresponding 
##'   to the number of risk strata)
##' @param ppr population proportions corresponding to rr values 
##'   (vector of equal length to rr)
##' @param df dataframe of values for each combination of risk stratum and 
##'   sensitivity level, 
##'   col 1 = risk group index, col 2 = unit Se, col 3 = n 
##'   (sample size for that risk group and unit sensitivity)
##' @return list of 3 elements, a scalar of population-level sensitivity
##'   a vector of EPI values and a vector of corresponding adjusted risks
##' @keywords methods
##' @export
##' @examples
##' # examples for sep.rb.bin.varse
##' rg<- c(1, 1, 2, 2)
##' se<- c(0.92, 0.85, 0.92, 0.85)
##' n<- c(80, 30, 20, 30)
##' df<- data.frame(rg, se, n)
##' sep.rb.bin.varse(0.01, c(5, 1), c(0.1, 0.9), df)
##' 
##' rg<- c(1, 1, 2, 2)
##' se<- c(0.95, 0.8, 0.95, 0.8)
##' n<- c(20, 10, 10, 5)
##' df<- data.frame(rg, se, n)
##' sep.rb.bin.varse(0.05, c(3, 1), c(0.2, 0.8), df)
##' 
##' rg<- c(rep(1, 30), rep(2, 15))
##' se<- c(rep(0.95, 20), rep(0.8, 10), rep(0.95, 10), rep(0.8, 5))
##' n<- rep(1, 45)
##' df<- data.frame(rg, se, n)
##' sep.rb.bin.varse(0.02, c(3, 1), c(0.2, 0.8), df)
##' 
##' rg<- c(1, 2, 3, 1, 2, 3)
##' se<- c(0.95, 0.95, 0.95, 0.8, 0.8, 0.8)
##' n<- c(20, 10, 10, 30, 5, 5)
##' df<- data.frame(rg, se, n)
##' sep.rb.bin.varse(0.01, c(5, 3, 1), c(0.1, 0.3, 0.6), df)
sep.rb.bin.varse<- function(pstar, rr, ppr, df) {
  epi<- epi.calc(pstar, rr, ppr)
  p.all.neg<- (1 - df[,2]*epi[[1]][df[,1]])^df[3]
  sep<- 1 - prod(p.all.neg)
  return(list(sep=sep, epi=epi[[1]], adj.risk=epi[[2]]))
}


##' Hypergeometric risk-based population sensitivity for varying unit sensitivity
##' @description Calculates population sensitivity for a single risk factor 
##'   and varying unit sensitivity using hypergeometric approximation method 
##'   (assumes known population size)
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector of values corresponding 
##'   to the number of risk strata)
##' @param N vector of population size for each risk group, corresponding to rr values 
##' (vector of equal length to rr)
##' @param df dataframe of values for each combination of risk stratum and
##'   sensitivity level, 
##'   col 1 = risk group index, col 2 = unit Se, col 3 = n 
##'   (sample size for risk group and unit sensitivity)
##' @return list of 5 elements, a scalar of population-level sensitivity
##'   a vector of EPI values, a vector of corresponding Adjusted risks
##'   a vector of sample sizes (n) per risk group and a vector of 
##'   mean unit sensitivities per risk group
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.rb.hypergeo.varse
##' rg<- c(1, 1, 2, 2)
##' se<- c(0.92, 0.85, 0.92, 0.85)
##' n<- c(80, 30, 20, 30)
##' df<- data.frame(rg, se, n)
##' sep.rb.hypergeo.varse(0.01, c(5, 1), c(200, 1800), df)
##' 
##' rg<- c(1, 1, 2, 2)
##' se<- c(0.95, 0.8, 0.95, 0.8)
##' n<- c(20, 10, 10, 5)
##' df<- data.frame(rg, se, n)
##' sep.rb.hypergeo.varse(0.05, c(3, 1), c(100, 400), df)
##' 
##' rg<- c(rep(1, 30), rep(2, 15))
##' se<- c(rep(0.95, 20), rep(0.8, 10), rep(0.95, 10), rep(0.8, 5))
##' n<- rep(1, 45)
##' df<- data.frame(rg, se, n)
##' sep.rb.hypergeo.varse(0.02, c(3, 1), c(100, 400), df)
##' 
##' rg<- c(1, 2, 3, 1, 2, 3)
##' se<- c(0.95, 0.95, 0.95, 0.8, 0.8, 0.8)
##' n<- c(20, 10, 10, 30, 5, 5)
##' df<- data.frame(rg, se, n)
##' sep.rb.hypergeo.varse(0.01, c(5, 3, 1), c(100, 300, 600), df)
sep.rb.hypergeo.varse<- function(pstar, rr, N, df) {
  ppr<- N/sum(N)
  epi<- epi.calc(pstar, rr, ppr)
  n<- numeric(length(rr))
  se<- n
  for (r in 1:length(rr)) {
    n[r]<- sum(df[df[,1] == r, 3])
    se[r]<- mean(df[df[,1] == r, 2])
  }
  p.all.neg<- (1-se*n/N)^(epi[[1]]*N)
  sep<- 1 - prod(p.all.neg)
  return(list(sep=sep, epi=epi[[1]], adj.risk=epi[[2]], n=n, se=se))
}


##' Binomial risk-based population sensitivity for 2 risk factors
##' @description Calculates risk-based population sensitivity for 
##'   two risk factors, using binomial method (assumes a large population)
##' @param pstar design prevalence (scalar)
##' @param rr1 relative risks for first level risk factor (vector of values corresponding 
##'   to the number of risk strata)
##' @param rr2 relative risks for second level risk factor, 
##'   matrix, rows = levels of rr1, cols = levels of rr2
##' @param ppr1 population proportions for first level risk factor (vector of
##'   same length as rr1)
##' @param ppr2 population proportions for second level 
##'   risk factor, matrix, rows = levels of rr1, cols = levels of rr2
##' @param n matrix of number tested for each risk group 
##'   (rows = levels of rr1, cols = levels of rr2)
##' @param se test unit sensitivity (scalar)
##' @return list of 4 elements, a scalar of population-level sensitivity
##'   a matrix of EPI values, a vector of corresponding Adjusted risks for
##'   the first risk factor and a matrix of adjusted risks for the second 
##'   risk factor
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.rb2.binom
##' pstar<- 0.01
##' rr1<- c(3, 1)
##' ppr1<- c(0.2, 0.8)
##' rr2<- rbind(c(4,1), c(4,1))
##' ppr2<- rbind(c(0.1, 0.9), c(0.3, 0.7))
##' se<- 0.8
##' n<- rbind(c(50, 20), c(20, 10))
##' sep.rb2.binom(pstar, rr1, ppr1, rr2, ppr2, n, se)
sep.rb2.binom<- function(pstar, rr1, ppr1, rr2, ppr2, n, se) {
  ar1<- adj.risk(rr1, ppr1)
  ar2<- array(0, dim = dim(rr2))
  rownames(ar2)<- paste("RR1",1:length(rr1), sep = "=")
  colnames(ar2)<- paste("RR2",1:ncol(rr2), sep = "=")
  epi<- ar2
  p.neg<- ar2
  if (length(se) == 1) se<- array(se, dim = dim(rr2))
  for (i in 1:length(rr1)) {
    ar2[i,]<- adj.risk(rr2[i,], ppr2[i,])
    epi[i,]<- ar1[i]*ar2[i,]*pstar
    p.neg[i,]<- (1 - epi[i,]*se[i,])^n[i,]
  }
  sep<- 1 - prod(p.neg)
  return(list(sep=sep, epi=epi, ar1=ar1, ar2=ar2))
}



##' Hypergeometric risk-based population sensitivity for 2 risk factors
##' @description Calculates risk-based population sensitivity for 
##'   two risk factors, using hypergeometric approximation method 
##'   (assumes a known population size)
##' @param pstar design prevalence (scalar)
##' @param rr1 relative risks for first level risk factor (vector of values corresponding 
##'   to the number of risk strata) 
##' @param rr2 relative risks for second level risk factor, 
##'   matrix, rows = levels of rr1, cols = levels of rr2
##' @param N matrix of population size for each risk group 
##'   (rows = levels of rr1, cols = levels of rr2)
##' @param n matrix of number tested (sample size) for each risk group 
##'   (rows = levels of rr1, cols = levels of rr2)
##' @param se test unit sensitivity (scalar)
##' @return list of 6 elements, a scalar of population-level sensitivity
##'   a matrix of EPI values, a vector of corresponding Adjusted risks for
##'   the first risk factor and a matrix of adjusted risks for the second risk factor,
##'   a vector of population proportions for the first risk factor 
##'   and a matrix of population proportions for the second risk factor
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.rb2.hypergeo
##' pstar<- 0.01
##' rr1<- c(3, 1)
##' rr2<- rbind(c(4,1), c(4,1))
##' N<- rbind(c(100, 500), c(300, 1000))
##' n<- rbind(c(50, 20), c(20, 10))
##' se<- 0.8
##' sep.rb2.hypergeo(pstar, rr1, rr2, N, n, se)
sep.rb2.hypergeo<- function(pstar, rr1, rr2, N, n, se) {
  ppr1<- rowSums(N)/sum(N)
  ppr2<- array(0, dim = dim(rr2))
  rownames(ppr2)<- paste("RR1",1:length(rr1), sep = "=")
  colnames(ppr2)<- paste("RR2",1:ncol(rr2), sep = "=")
  ar1<- adj.risk(rr1, ppr1)
  ar2<- array(0, dim = dim(rr2))
  rownames(ar2)<- rownames(ppr2)
  colnames(ar2)<- colnames(ppr2)
  epi<- ar2
  p.neg<- ar2
  if (length(se) == 1) se<- array(se, dim = dim(rr2))
  for (i in 1:length(rr1)) {
    ppr2[i,]<- N[i,]/sum(N[i,])
    ar2[i,]<- adj.risk(rr2[i,], ppr2[i,])
    epi[i,]<- ar1[i]*ar2[i,]*pstar
    p.neg[i,]<- (1 - se[i,]*n[i,]/N[i,])^(epi[i,]*N[i,])
  }
  sep<- 1 - prod(p.neg)
  return(list(sep=sep, epi=epi, 
              ar1=ar1, ar2=ar2, 
              ppr1=ppr1, ppr2=ppr2))
}


##' Two-stage risk-based system sensitivity
##' @description Calculates system sensitivity for 2 stage risk-based 
##'   sampling, llowing for a single risk factor at each stage and
##'   using either binomial or hypergeometric approxiation
##' @param C Population size (number of clusters), NA = unknown (default)
##' @param pstar.c cluster level design prevalence (scalar)
##' @param pstar.u unit level design prevalence (scalar)
##' @param rr.c cluster level relative risks (vector with length 
##'   corresponding to the number of risk strata), 
##'   use rr.c = c(1,1) if risk factor does not apply  
##' @param rr.u unit level relative risks (vector with length 
##'   corresponding to the number of risk strata), 
##'   use rr.u = c(1,1) if risk factor does not apply  
##' @param ppr.c cluster level population proportions for risk 
##'   categories (vector), NA if no cluster level risk factor
##' @param N population size per risk group for each cluster, 
##'   NA or matrix of N for each risk group 
##'   for each cluster, N=NA means cluster sizes not provided
##' @param rg vector of cluster level risk group (index) for each cluster
##' @param n sample size per risk group for each cluster sampled,
##'   matrix, 1 row for each cluster, columns = unit level risk groups
##' @param ppr.u unit level population proportions for each risk group (optional) 
##'   matrix, 1 row for each cluster, columns = unit level risk groups, 
##'   not required if N is provided
##' @param se unit sensitivity for each cluster, scalar or 
##'   vector of values for each cluster, equal in length to n
##' @return list of 2 elements, a scalar of population-level (surveillance system) 
##'   sensitivity and a vector of cluster-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples for sse.rb.2stage
##' pstar.c<- 0.02
##' pstar.u<- 0.1
##' rr.c<- c(5, 1)
##' ppr.c<- c(0.1, 0.9)
##' rr.u<- c(3, 1)
##' se<- 0.9
##' n<- cbind(rep(10, 50), rep(5, 50))    
##' rg<- c(rep(1, 30), rep(2, 20))
##' ppr.u<- cbind(rep(0.2, 50), rep(0.8, 50))
##' N<- cbind(rep(30, 50), rep(120, 50))
##' C<- 500        
##' sse.rb.2stage(C=NA, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N=NA, n, rg, se) 
##' sse.rb.2stage(C, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N=NA, n, rg, se) 
##' sse.rb.2stage(C=NA, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N, n, rg, se) 
##' sse.rb.2stage(C, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N, n, rg, se) 
sse.rb.2stage<- function(C=NA, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N=NA, n, rg, se) {
  if (length(se) == 1) se<- rep(se, nrow(n))
  sep<- numeric(nrow(n))
  # calculate sep for all clusters
  if (length(N) == 1)  {
    # cluster sizes not provided so use binomial for all clusters
    for (i in 1:nrow(n)) {
      sep[i]<- sep.rb.bin(pstar.u, rr.u, ppr.u[i,], n[i,], se[i])[[1]]
    } 
  } else {
    # cluster sizes provided so use hypergeometric nless NA for specific clusters
    for (i in 1:nrow(n)) {
      if (is.na(N[i,1])) {
        sep[i]<- sep.rb.bin(pstar.u, rr.u, ppr.u[i,], n[i,], se[i])[[1]]
      } else { 
        sep[i]<- sep.rb.hypergeo(pstar.u, rr.u, N[i,], n[i,], se[i])[[1]]
      }  
    }
  } 
  # calculate system sensitivity
  if (is.na(C)) {  
    # Population size unnown, use binomial
    sse<- sep.rb.bin.varse(pstar.c, rr.c, ppr.c, df=cbind(rg, sep, 1))
  } else {
    sse<- sep.rb.hypergeo.varse(pstar.c, rr.c, C*ppr.c, df=cbind(rg, sep, 1))
  }
  return(list("System sensitivity" = sse[[1]], 
              "Cluster sensitivity" = sep))
}


##' System sensitivity by combining multiple surveillance components
##' @description Calculates overall system sensitivity for 
##'   multiple components, accounting for lack of independence 
##'   (overlap) between components
##' @param C population sizes (number of clusters) for each risk group, 
##' NA or vector of same length as rr
##' @param pstar.c cluster level design prevalence (scalar)
##' @param rr cluster level relative risks (vector, length 
##'   equal to the number of risk strata)
##' @param ppr cluster level population proportions (optional), 
##'   not required if C is specified (NA or vector of same length as rr)
##' @param sep sep values for clusters in each component and 
##'   corresponding risk group. A list with multiple elements, each element 
##'   is a dataframe of sep values from a separate component, 
##'   first column= clusterid, 2nd =cluster-level risk group index, 3rd col = sep
##' @return list of 2 elements, a matrix (or vector if C not specified) 
##'   of population-level (surveillance system) 
##'   sensitivities (binomial and hypergeometric and adjusted vs unadjusted) and 
##'   a matrix of adjusted and unadjusted component sensitivities for each component
##' @keywords methods
##' @export
##' @examples 
##' # example for sse.combined (checked in excel combined components.xlsx)
##' C<- c(300, 1200)
##' pstar<- 0.01
##' rr<- c(3,1)
##' ppr<- c(0.2, 0.8)
##' comp1<- data.frame(id=1:100, rg=c(rep(1,50), rep(2,50)), cse=rep(0.5,100)) 
##' comp2<- data.frame(id=seq(2, 120, by=2), rg=c(rep(1,25), rep(2,35)), cse=runif(60, 0.5, 0.8))
##' comp3<- data.frame(id=seq(5, 120, by=5), rg=c(rep(1,10), rep(2,14)), cse=runif(24, 0.7, 1))
##' sep<- list(comp1, comp2, comp3)
##' sse.combined(C, pstar, rr, sep = sep)
##' sse.combined(C=NA, pstar, rr, ppr, sep = sep)
sse.combined<- function(C = NA, pstar.c, rr, ppr, sep) {
  if (length(C) > 1) ppr<- C/sum(C)
  components<- length(sep)
  epi<- epi.calc(pstar.c, rr, ppr)[[1]]
  # Create master list of clusters sampled
  cluster.list<- sep[[1]]
  i<- 2
  while (i <= components) {
    cluster.list<- merge(cluster.list, sep[[i]], by.x = 1, by.y = 1, all.x=T, all.y=T)                   
    i<- i+1
  }
  # ensure risk group recorded in data
  risk.group<- cluster.list[,2]
  tmp<- which(is.na(risk.group))
  if (length(tmp)>0) {
    for (i in tmp) {
      j<- 2
      while (j<=components && is.na(risk.group[i])) {
        risk.group[i]<- cluster.list[i,(j-1)*2+2]
        j<- j+1
      }
    }
  }
  # Replace NA values with 0
  for (i in 2:ncol(cluster.list)) {
    cluster.list[is.na(cluster.list[,i]), i]<- 0
  }
  # set up arrays for epi  and p.neg (adjusted and unadjusted) for each cluster and each component
  epi.c<- array(0, dim = c(nrow(cluster.list), components))
  epi.c[,1]<- epi[risk.group]
  # dim 3: 1 = adjusted, 2 = unadjusted (independence)
  p.neg<- array(0, dim = c(nrow(cluster.list), components, 2))
  p.neg[,1,1]<- 1-cluster.list[,3]*epi.c[,1]
  p.neg[,1,2]<- p.neg[,1,1]
  for (i in 2:components) {
    for (j in 1:nrow(cluster.list)) {
      epi.c[j,i]<- 1 - pfree.1(cluster.list[j,(i-1)*2+1], 0, 1-epi.c[j,i-1])[,4]
    }
    p.neg[,i,1]<- 1-cluster.list[,(i-1)*2+3]*epi.c[,i]
    p.neg[,i,2]<- 1-cluster.list[,(i-1)*2+3]*epi.c[,1]
  }
  # calculate n, mean sep and mean epi for each risk group and component
  n<- array(0, dim = c(components, length(rr)))
  sep.mean<- array(0, dim = c(components, length(rr)))
  epi.mean<- array(0, dim = c(components, length(rr), 2))
  for (i in 1:components) {
    n[i,]<- table(sep[[i]][2])
    sep.mean[i,]<- sapply(split(sep[[i]][3], sep[[i]][2]), FUN=colMeans)
    epi.mean[i,,1]<- sapply(split(epi.c[cluster.list[,(i-1)*2+2] > 0,i], cluster.list[cluster.list[,(i-1)*2+2] > 0,(i-1)*2+2]), FUN=mean)
    epi.mean[i,,2]<- epi.mean[1,,1]
  }
  # Calculate cse and sse
  cse<- array(0, dim = c(2, components, 2))
  rownames(cse)<- c("Adjusted", "Unadjusted")
  colnames(cse)<- paste("Component", 1:components)
  dimnames(cse)[[3]]<- c("Binomial", "Hypergeometric")
  sse<- array(0, dim = c(2, 2))
  rownames(sse)<- rownames(cse)
  colnames(sse)<- dimnames(cse)[[3]]
  rownames(epi.mean)<- colnames(cse)
  rownames(sep.mean)<- colnames(cse)
  rownames(n)<- colnames(cse)
  colnames(epi.mean)<- paste("RR =",rr)
  colnames(sep.mean)<- paste("RR =",rr)
  colnames(n)<- paste("RR =",rr)
  dimnames(epi.mean)[[3]]<- rownames(cse)
  # rows = adjusted and unadjusted, dim3 = binomial and hypergeometric
  for (i in 1:2) {
    for (j in 1:components) {
      cse[i,j,1]<- 1 - prod(p.neg[,j,i])
      if (length(C) > 1) {
        cse[i,j,2]<- 1 - prod((1 - sep.mean[j,]*n[j,]/C)^(epi.mean[j,,i]*C))
      }
    }
    sse[i,1]<- 1- prod(1 - cse[i,,1])
    sse[i,2]<- 1- prod(1 - cse[i,,2])
  }
  if (length(C) <= 1) {
    sse<- sse[,1]
    cse<- cse[,,1]
  }
  return(list("System sensitivity"= sse, 
                "Component sensitivity" = cse))    
}  



##' Risk-based sample size
##' @description Calculates sample size for risk-based sampling 
##'   for a single risk factor and using binomial method
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector, length equal to the number of risk strata)
##' @param ppr population proportions corresponding to rr values 
##'   (vector of equal length to rr)
##' @param spr planned surveillance proportion for each risk group 
##'   (vector equal length to rr, ppr)
##' @param se unit sensitivity (fixed or vector same length as rr, ppr, n)
##' @param sep required population sensitivity (scalar)
##' @return list of 2 elements, a vector of sample sizes for each risk group
##'   a scalar of total sample size, a vector of EPI values and a vector of
##'   adjusted risks 
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.rb
##' n.rb(0.1, c(5, 3, 1), c(0.1, 0.10, 0.80), c(0.5, 0.3, 0.2), 0.9, 0.95)
##' n.rb(0.01, c(5, 1), c(0.1, 0.9), c(0.8, 0.2), c(0.9, 0.95), 0.95)
n.rb<- function(pstar, rr, ppr, spr, se, sep) {
  epi<- epi.calc(pstar, rr, ppr)
  p.pos<- sum(epi[[1]]*spr*se)
  n.total<- ceiling(log(1-sep)/log(1-p.pos))
  n<- numeric(length(rr))
  for (i in 1:length(rr)) {
    if (i<length(rr)) {
      n[i]<- ceiling(n.total*spr[i])
    } else {
      n[i]<- n.total - sum(n)
    }
  }
  return(list(n=n, total=n.total, epi=epi[[1]], adj.risk=epi[[2]]))
}

##' Risk-based sample size for varying unit sensitivity
##' @description Calculates sample size for risk-based sampling 
##'   for a single risk factor and varying unit sensitivity, 
##'   using binomial method
##' @param pstar design prevalence (scalar)
##' @param rr relative risk values (vector, length equal to the number of risk strata)
##' @param ppr population proportions for each risk group,
##'   vector of same length as rr
##' @param spr planned surveillance proportions for each risk group,
##'   vector of same length as rr
##' @param se unit sensitivities (vector of group values)
##' @param spr.rg proportions of samples for each sensitivity value 
##'   in each risk group (matrix with rows = risk groups, columns = sensitivity values),
##'   row sums must equal 1
##' @param sep required population sensitivity (scalar)
##' @return list of 3 elements, a matrix of sample sizes for each risk 
##'   and sensitivity group, a vector of EPI values and a vector of 
##'   mean sensitivity for each risk group
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.rb.varse
##' m<- rbind(c(0.8, 0.2), c(0.5, 0.5), c(0.7, 0.3))
##' n.rb.varse(0.01, c(5, 3, 1), c(0.1, 0.1, 0.8), c(0.4, 0.4, 0.2), c(0.92, 0.8), m, 0.95)
##' 
##' m<- rbind(c(0.8, 0.2), c(0.6, 0.4))
##' n.rb.varse(0.05, c(3, 1), c(0.2, 0.8), c(0.7, 0.3), c(0.95, 0.8), m, 0.95)
##' 
##' m<- rbind(c(1), c(1))
##' n.rb.varse(0.05, c(3, 1), c(0.2, 0.8), c(0.7, 0.3), c(0.95), m, 0.99)
n.rb.varse<- function(pstar, rr, ppr, spr, se, spr.rg, sep) {
  mean.se<- numeric(length(rr))
  for (r in 1:length(rr)) {
    mean.se[r]<- sum(spr.rg[r,] * se)
  }
  epi<- epi.calc(pstar, rr, ppr)[[1]]
  p.pos<- sum(epi*mean.se*spr)
  n.total<- ceiling(log(1 - sep)/log(1 - p.pos))
  n.rg<- numeric(length(rr))
  n<- array(0, dim = c(nrow(spr.rg), ncol(spr.rg)))
  for (i in 1:length(rr)) {
    if (i<length(rr)) {
      n.rg[i]<- ceiling(n.total*spr[i])
    } else {
      n.rg[i]<- n.total - sum(n.rg)
    }
    for (j in 1:length(se)) {
      if (j<length(se)) {
        n[i, j]<- ceiling(n.rg[i]*spr.rg[i, j])
      } else {
        n[i, j]<- n.rg[i] - sum(n[i,])
      }
    }
  }
  n<- cbind(n, n.rg)
  tmp<- apply(n, FUN=sum, MARGIN=2)
  n<- rbind(n, tmp)
  colnames(n)<- c(paste("Se =", se), "Total")
  rownames(n)<- c(paste("RR =", rr), "Total")
  return(list(n=n, epi=epi, mean.se=mean.se))
}


