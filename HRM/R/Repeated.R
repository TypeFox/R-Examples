#' EEG data of 160 subjects
#' 
#' A dataset containing EEG data (Staffen et al., 2014) of 160 subjects, 4 variables are measured at ten different locations.
#' The columns are as follows:
#' 
#' \itemize{
#'   \item group. Diagnostic group of the subject: Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI), Subject Cognitive Complaints (SCC+, SCC-).
#'   \item value. Measured data of a subject at a specific variable and region.
#'   \item sex. Sex of the subject: Male (M) or Female (W).
#'   \item subject. A unique identification of a subject.
#'   \item variable. The variales measured are activity, complexity, mobility and brain rate coded from 1 to 4.
#'   \item region. Frontal left/right, central left/right, temporal left/right, occipital left/right, parietal left/right coded as 1 to 10. 
#'   \item dimension. Mixing variable and region together, levels range from 1 to 40.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name EEG
#' @usage data(EEG)
#' @format A data frame with 6400 rows and 7 variables.
"EEG"

# functions which define the estimator

#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E1 = function(n,i, M) {
  return ((n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E2 = function(n,i, M) {
  return ((n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E3 = function(M_i, M_j) {
  return (matrix.trace(M_i)*matrix.trace(M_j))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E4 = function(M_i,M_j) {
  return (matrix.trace(M_i%*%M_j))
}

#######################################

#' Function for the output: significant p-values have on or more stars
#' 
#' @param value p-value
#' @export
#' @keywords internal
.hrm.sigcode = function(value) {
  
  char=""
  if(value <= 0.1 & value > 0.05) { char = "."}
  if(value <= 0.05 & value > 0.01) {char = '*'}
  if(value <= 0.01 & value > 0.001) {char = "**"}
  if(value <= 0.001 & value >= 0) {char = "***"}
  return (char)
}

I = function(size){
  return (diag(rep(1,size)))
}

J = function(size){
  return (rep(1,size)%*%t(rep(1,size)))
}

P = function(size){
  return (I(size)-J(size)*1/size)
}


DualEmpirical = function(Data, B){
  n = dim(Data)[1]
  B = B(dim(Data)[2]) # B = e.g. t(J_d)%*%J_d
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}


DualEmpirical2 = function(Data, B){
  n = dim(Data)[1]
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}

#' Test for no main treatment effect (weighted version)
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A.weighted = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0
  if(is.data.frame(X[[1]])){
    X = lapply(X, as.matrix)
  }
  
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  D_a = diag(n)-1/sum(n)*n%*%t(n)
  K_A = kronecker(D_a,J(d))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = J)

  #################################################################################################
  # f = f_1/f_2: numerator degrees of freedom
  f_1 = 0
  f_2 = 0
  
  
  # computation of f_1
  for(i in 1:a){
    f_1 = f_1 + ((1-n[i]/sum(n))^2)*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*((1-n[i]/sum(n)))*((1-n[j]/sum(n)))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  # computation of f_2
  for(i in 1:a){
    f_2 = f_2 + ((1-n[i]/sum(n)))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*(n[i]*n[j])/(d^2*sum(n)^2)*.E4(d/(n[i]-1)*P(n[i])%*%X[[i]],d/(n[j]-1)*J(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%J(d)%*%t(X[[i]])%*%P(n[i]))
      j=j+1
    }
  }
  
  f=f_1/f_2
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0 = f0_1/f0_2: denumerator degrees of freedom
  # f0_1 = f_1 numerator are the same
  f0_2 = 0
  
  # computation of f0_2
  for(i in 1:a){
    f0_2 = f0_2 + ((1-n[i]/sum(n)))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # constructing the  Test
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }

  test = (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="A weighted",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  return (output)
}

# Hypothesis 1 ------------------------------------------------------------





#' Test for no main time effect
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.B = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  D_a = diag(n)-1/sum(n)*n%*%t(n)
  K_A = kronecker(D_a,J(d))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = P)
  
  K_B = kronecker(J(a),P(d))
  
#   V = list(P_d%*%var(X[[1]]),P_d%*%var(X[[2]]))
#   if(a>2){
#     for(i in 3:a){
#       V[[i]] = (P_d%*%var(X[[i]]))
#     }
#   }
  
  #################################################################################################
  
  # f
  f_1 = 0
  f_2 = 0
  
  for(i in 1:a){
    f_1 = f_1 + (1/n[i])^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(1/n[i])*(1/n[j])*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (1/n[i])^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*P(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%P(d)%*%t(X[[i]])%*%P(n[i]))
      j=j+1
    }
  }
  
  f=f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 = f_1
  f0_2 = 0
  
  
  for(i in 1:a){
    f0_2 = f0_2 + (1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # Test
  
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test = (t(X_bar)%*%K_B%*%X_bar)/(t(rep(1,dim(K_B)[1]))%*%(K_B*direct)%*%(rep(1,dim(K_B)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="B",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis 2 ------------------------------------------------------------







#' Test for no interaction between treatment and time
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.AB = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = P)
  
  K_AB = kronecker(P(a),P(d))
  
#   V = list(P_d%*%var(X[[1]]),P_d%*%var(X[[2]]))
#   if(a>2){
#     for(i in 3:a){
#       V[[i]] = (P_d%*%var(X[[i]]))
#     }
#   }
  
  #################################################################################################
  # f
  f_1 = 0
  f_2 = 0
  
  
  # phi = A
  for(i in 1:a){
    f_1 = f_1 + (1/n[i]*(1-1/a))^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(1/n[i]*(1-1/a))*(1/n[j]*(1-1/a))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (1/n[i]*(1-1/a))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2/(n[i]*n[j]*a^2)*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*P(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%P(d)%*%t(X[[i]])%*%P(n[i]))
      j=j+1
    }
  }
  
  f=f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 = f_1
  f0_2 = 0
  
  for(i in 1:a){
    f0_2 = f0_2 + (1/n[i]*(1-1/a))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # Test
  
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test = (t(X_bar)%*%K_AB%*%X_bar)/(t(rep(1,dim(K_AB)[1]))%*%(K_AB*direct)%*%(rep(1,dim(K_AB)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="AB",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis 3 ------------------------------------------------------------






#' Test for no simple treatment effect
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A_B = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = I)
  
  K_A_B = kronecker(P(a),I(d))
  
#   V = list(var(X[[1]]),var(X[[2]]))
#   if(a>2){
#     for(i in 3:a){
#       V[[i]] = (var(X[[i]]))
#     }
#   }
  
  #################################################################################################
  # f
  f_1 = 0
  f_2 = 0
  
  for(i in 1:a){
    f_1 = f_1 + (1/n[i]*(1-1/a))^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(1/n[i]*(1-1/a))*(1/n[j]*(1-1/a))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (1/n[i]*(1-1/a))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2/(n[i]*n[j]*a^2)*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*t(X[[j]])%*%P(n[j])%*%X[[j]]%*%t(X[[i]])%*%P(n[i]))
      j=j+1
    }
  }
  
  f=f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 = f_1
  f0_2 = 0
  
  
  for(i in 1:a){
    f0_2 = f0_2 + (1/n[i]*(1-1/a))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # Test
  
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test = (t(X_bar)%*%K_A_B%*%X_bar)/(t(rep(1,dim(K_A_B)[1]))%*%(K_A_B*direct)%*%(rep(1,dim(K_A_B)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="A|B",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  return (output)
}

# Hypotheses 4 ------------------------------------------------------------




#' Test for no main treatment effect (unweighted version)
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A.unweighted = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  # creating dual empirical covariance matrices, where observations are transformed by a matrix B
  V = lapply(X, DualEmpirical, B = J)
  
  K_A = kronecker(P(a),J(d))
  

  
  #################################################################################################
  # f
  f_1 = 0
  f_2 = 0
  
  
  
  for(i in 1:a){
    f_1 = f_1 + ((1-1/a)/(d*n[i]))^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*((1-1/a)/(d*n[i]))*((1-1/a)/(d*n[j]))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  
  for(i in 1:a){
    f_2 = f_2 + ((1-1/a)/(d*n[i]))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*(1/(a^2*n[i]*n[j]*d^2))*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*J(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%J(d)%*%t(X[[i]])%*%P(n[i]))
      j=j+1                                 
    }
  }
  
  f=f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 = f_1
  f0_2 = 0
  
  for(i in 1:a){
    f0_2 = f0_2 + ((1-1/a)/(d*n[i]))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # Test
  
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test = (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="A unweighted",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(1-pf(test,f,f0)))
  

  
  return (output)
}

# Hypothesis 1/2 ------------------------------------------------------------


#' Test for no main treatment effect, no main time effect, no simple treatment effect and no interaction between treatment and time
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_1.txt
#' @export
hrm.test = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  temp0=hrm.A.weighted(n,a,d,X,alpha)
  temp1 = hrm.A.unweighted(n,a,d,X,alpha)
  temp2=hrm.B(n,a,d,X,alpha)
  temp3=hrm.AB(n,a,d,X,alpha)
  temp4 = hrm.A_B(n,a,d,X,alpha)
  output = rbind(temp0,temp1,temp2,temp3,temp4)
  return (output)
}


#' Test for no main effects and interactino effects of one between-subject factor and two crossed within-subject factors
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.within = function(X, alpha, group , factor1, factor2, subject, data ){
  
  temp4= hrm.1w.2f.A(X, alpha, group , factor1, factor2, subject, data )
  temp5= hrm.1w.2f.B(X, alpha, group , factor1, factor2, subject, data )
  temp6= hrm.1w.2f.C(X, alpha, group , factor1, factor2, subject, data )
  temp0= hrm.1w.2f.AB(X, alpha, group , factor1, factor2, subject, data )
  temp1= hrm.1w.2f.AC(X, alpha, group , factor1, factor2, subject, data )
  temp2= hrm.1w.2f.BC(X, alpha, group , factor1, factor2, subject, data )
  temp3= hrm.1w.2f.ABC(X, alpha, group , factor1, factor2, subject, data )
  
  output = rbind(temp4, temp5, temp6, temp0,temp1,temp2,temp3)
  return (output)
}



#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor column name of the data frame X of within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between = function(X, alpha, group , subgroup, factor, subject, data ){
  
  temp4= hrm.2w.1f.A(X, alpha, group , subgroup, factor, subject, data )
  temp5= hrm.2w.1f.A2(X, alpha, group , subgroup, factor, subject, data )
  temp6= hrm.2w.1f.B(X, alpha, group , subgroup, factor, subject, data  )
  temp0= hrm.2w.1f.AA2(X, alpha, group , subgroup, factor, subject, data  )
  temp1 =  hrm.2w.1f.AB(X, alpha, group , subgroup, factor, subject, data  )
  temp2= hrm.2w.1f.A2B(X, alpha, group , subgroup, factor, subject, data  )
  temp3= hrm.2w.1f.AA2B(X, alpha, group , subgroup, factor, subject, data  )
  
  output = rbind(temp4, temp5, temp6, temp0, temp1, temp2, temp3)
  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between.within = function(X, alpha, group , subgroup, factor1, factor2, subject, data ){
  
  temp0= hrm.2w.2f.A(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp1= hrm.2w.2f.A2(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp2= hrm.2w.2f.B(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp3= hrm.2w.2f.C(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp4= hrm.2w.2f.AA2(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp5= hrm.2w.2f.AB(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp6= hrm.2w.2f.AC(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp7= hrm.2w.2f.A2B(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp8= hrm.2w.2f.A2C(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp9= hrm.2w.2f.BC(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp10= hrm.2w.2f.AA2B(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp11= hrm.2w.2f.AA2C(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp12= hrm.2w.2f.ABC(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp13= hrm.2w.2f.A2BC(X, alpha, group , subgroup, factor1, factor2, subject, data )
  temp14= hrm.2w.2f.AA2BC(X, alpha, group , subgroup, factor1, factor2, subject, data )
  
  
  output = rbind(temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14)
  return (output)
}


#' Test for main effects and interaction effects of one or two between-subject factors and one or two within-subject factors
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor1 column name of the data frame X of within-subject factor
#' @param factor2 column name of the second within-subject factor crossed with factor1
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_2.txt
#' @export
hrm.test.2 = function(X, alpha = 0.05, group , subgroup, factor1, factor2, subject, data ){
  
  if(missing(group) || !is.character(group)){
    print("At least one between-subject factor is needed!")
    stop("group column name not specified ")
  }
  
  if(missing(factor1) || !is.character(factor1)){
    print("At least one between-subject factor is needed!")
    stop("factor1 column name not specified ")
  }
  
  if(missing(X) || !is.data.frame(X)){
    stop("dataframe needed")
  }
  
  if(missing(subject) || !is.character(subject)){
    stop("subject column name not specified")
  }
  
  if(missing(data) || !is.character(data)){
    stop("data column name not specified")
  }
  
  if(sum(is.na(X[,data]))>=1){
    warning("Your data contains missing values!")
  }
  
  if(!is.double(alpha)){
    stop("alpha level needs to be a number between 0 and 1")
  }
  if(is.double(alpha)){
    if(alpha > 1 || alpha < 0){
      stop("alpha level needs to be a number between 0 and 1")
    }
  }
  
  
  if(missing(factor2) & missing(subgroup)){
    stop("The model needs at least two between-subject factors or two within-subject factors. Please use otherwise the function hrm.test.")
  }
  
  if(missing(factor2) & !missing(subgroup)){
    if(!is.character(subgroup)){
      stop("subgroup column name not specified")
    }
    return(hrm.test.2.between(X, alpha, group , subgroup, factor1, subject, data ))
  }
  
  if(!missing(factor2) & missing(subgroup)){
    if(!is.character(factor2)){
      stop("factor2 column name not specified")
    }
    return(hrm.test.2.within(X, alpha, group , factor1, factor2, subject, data ))
  }
  
  if(!missing(factor2) & !missing(subgroup)){
    if(!is.character(factor2)){
      stop("factor2 column name not specified")
    }
    if(!is.character(subgroup)){
      stop("subgroup column name not specified")
    }
    return(hrm.test.2.between.within(X, alpha, group , subgroup, factor1, factor2, subject, data ))
  }
  
}