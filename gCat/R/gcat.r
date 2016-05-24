
gcat.test = function(counts, distmatrix=NULL, C0=NULL, method="C-uMST", Nperm=0){
  if (length(method)==0){
    cat("Specify an appropriate method to use.\n")
    return(0)
  }
  methodid = match(method, c("aMST","C-uMST", "uMST", "C-uNNB", "RC0", "TC0"))
  NAids = which(is.na(methodid))
  if (length(NAids)>0){
    cat(method[NAids], "is/are not in the method list and so removed.\n")
    methodid = methodid[-NAids]
  }
  if (length(methodid)==0){
    cat("Specify an appropriate method to use.\n")
    return(0)
  }
  if ((!is.na(match(1,methodid))) && Nperm==0){
    cat("Specify the number of permutations for calculating p-value for R_aMST\n")
    return(0)
  }
  if (( (!is.na(match(1,methodid))) || (!is.na(match(2,methodid))) || (!is.na(match(3,methodid))) || (!is.na(match(4,methodid)))) && is.null(distmatrix)){
    cat("You need to specify the distance matrix on the categories for calculating R_aMST, R_C-uMST, R_uMST, or R_C-uNNB.\n")
    return(0)
  }
  if (((!is.na(match(5,methodid))) || (!is.na(match(6,methodid)))) && is.null(C0)){
    cat("You need to specify the edge information for calcuating RC0 or TC0.\n")
    return(0)
  }
  K = dim(counts)[1]
  B = Nperm
  permmat = matrix(0,K,2*(B+1))
  permmat[,1:2] = counts
  v1 = counts[,1]
  v2 = counts[,2]
  if (B>0){
    for (j in 1:B){
      temp = perm(v1,v2)
      permmat[,2*j+(1:2)] = t(temp)
    }
  }
  result = list()
  if (!is.na(match(1,methodid))){
    rawresult0 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(0), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(distmatrix), PACKAGE="gCat")
    result0 = rawresult0[[1]]
    r = matrix(c(result0[1], myp_left(result0[2:(B+1)], result0[1])),1)
    colnames(r) = c("R_aMST", "pval.perm")
    result = c(result, list(aMST=r))
  }
  if (!is.na(match(2,methodid)) || (!is.na(match(3,methodid)))){
    rawresult1 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(1), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(distmatrix), PACKAGE="gCat")
    edge = t(matrix(rawresult1[[3]],2)[,1:rawresult1[[2]]]) + 1
    if (!is.na(match(2,methodid))){
      result1 = rawresult1[[1]]
      tmp = getMV(v1, v2, edge)
      if (B<=0){
        r = matrix(c(result1[1], pnorm(result1[1], tmp$Mean, tmp$Sd)),1)
        colnames(r) = c("R_C-uMST", "pval.appr")
      }else{
        r = matrix(c(result1[1], pnorm(result1[1], tmp$Mean, tmp$Sd), myp_left(result1[2:(B+1)], result1[1])),1)
        colnames(r) = c("R_C-uMST", "pval.appr", "pval.perm")
      }
      result = c(result, list(CuMST=r))
    }
    if (!is.na(match(3,methodid))){
      rawresult4 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(4), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(distmatrix), PACKAGE="gCat")
      result4 = rawresult4[[1]]
      tmp = getMV_uMST(v1, v2, edge)
      if (B<=0){
        r = matrix(c(result4[1], pnorm(result4[1], tmp$Mean, tmp$Sd)),1)
        colnames(r) = c("R_uMST", "pval.appr")
      }else{
        r = matrix(c(result4[1], pnorm(result4[1], tmp$Mean, tmp$Sd), myp_left(result4[2:(B+1)], result4[1])),1)
        colnames(r) = c("R_uMST", "pval.appr", "pval.perm")
      }
      result = c(result, list(uMST=r))
    }
  }
  if (!is.na(match(4,methodid))){
    rawresult3 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(3), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(distmatrix), PACKAGE="gCat")
    result3 = rawresult3[[1]]
    edge3 = t(matrix(rawresult3[[3]],2)[,1:rawresult3[[2]]]) + 1
    tmp3 = getMV(v1, v2, edge3)
    if (B<=0){
      r = matrix(c(result3[1], pnorm(result3[1], tmp3$Mean, tmp3$Sd)),1)
      colnames(r) = c("R_C-uNNB", "pval.appr")
    }else{
      r = matrix(c(result3[1], pnorm(result3[1], tmp3$Mean, tmp3$Sd), myp_left(result3[2:(B+1)], result3[1])),1)
      colnames(r) = c("R_C-uNNB", "pval.appr", "pval.perm")
    }
    result = c(result, list(CuNNB=r))
  }
  if (!is.na(match(5,methodid))){
    R = RC(counts, C0)
    tmp = getMV(v1,v2,C0)
    if (B<=0){
      r = matrix(c(R, pnorm(R,tmp$Mean, tmp$Sd)), 1)
      colnames(r) = c("RC0", "pval.appr")
    }else{
      Rv = rep(0,B)
      for (j in 1:B){
        Rv[j] = RC(permmat[,(2*j)+(1:2)], C0)
      }
      r = matrix(c(R, pnorm(R,tmp$Mean, tmp$Sd), myp_left(Rv,R)), 1)
      colnames(r) = c("RC0", "pval.appr", "pval.perm")
    }
    result = c(result, list(RC0=r))
  }
  if (!is.na(match(6,methodid))){
    R = TC(counts, C0)
    tmp = getMV_uMST(v1,v2,C0)
    if (B<=0){
      r = matrix(c(R, pnorm(R,tmp$Mean, tmp$Sd)), 1)
      colnames(r) = c("RC0", "pval.appr")
    }else{
      Rv = rep(0,B)
      for (j in 1:B){
        Rv[j] = TC(permmat[,(2*j)+(1:2)], C0)
      }
      r = matrix(c(R, pnorm(R,tmp$Mean, tmp$Sd), myp_left(Rv,R)), 1)
      colnames(r) = c("TC0", "pval.appr", "pval.perm")
    }
    result = c(result, list(TC0=r))
  }
  return(result)
}


RC = function(counts, edge){
  v1 = counts[,1]
  v2 = counts[,2]
  v = v1+v2
  R = sum(2*v1*v2/v)
  for (i in 1:(dim(edge)[1])){
    k1 = edge[i,1]
    k2 = edge[i,2]
    R = R + (v1[k1]*v2[k2] + v1[k2]*v2[k1])/v[k1]/v[k2]
  }
  R
}

TC = function(counts, edge){
  v1 = counts[,1]
  v2 = counts[,2]
  v = v1+v2
  R = sum(v1*v2)
  for (i in 1:(dim(edge)[1])){
    k1 = edge[i,1]
    k2 = edge[i,2]
    R = R + (v1[k1]*v2[k2] + v1[k2]*v2[k1])
  }
  R
}


