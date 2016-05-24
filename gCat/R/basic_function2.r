eps=1e-10

# calculate the Pearson's statistic given two vectors.
pearson = function(vec1, vec2){
  v = vec1+ vec2
  id = which(v!=0) 
  vec1 = vec1[id]
  vec2 = vec2[id]
  K = length(vec1)
  n1 = sum(vec1)
  n2 = sum(vec2)
  result = 0
  for (i in 1:K){
    m = vec1[i]+vec2[i]
    E1 = n1*m/(n1+n2)
    E2 = n2*m/(n1+n2)
    result = result + (vec1[i]-E1)^2/E1 + (vec2[i]-E2)^2/E2
  }
  result
}

# calculate the likelihood ratio statistic given two vectors.
LR = function(vec1, vec2){
  K = length(vec1)
  n1 = sum(vec1)
  n2 = sum(vec2)
  result = 0
  for (i in 1:K){
    m = vec1[i]+vec2[i]
    E1 = n1*m/(n1+n2)
    E2 = n2*m/(n1+n2)
    if (vec1[i] >0){
      result = result + vec1[i]*log(vec1[i]/E1)
    }
    if (vec2[i]>0){
      result = result + vec2[i]*log(vec2[i]/E2)
    }
  }
  result
}

# calculated the Friedman-Rafsky's statistic for ranked categorical data.
FRranked = function(vec1, vec2){
  v = vec1+ vec2
  id = which(v!=0)
  vec1 = vec1[id]
  vec2 = vec2[id]
  K = length(vec1)
  R = 0
  for (i in 1:K){
    R = R + 2*vec1[i]*vec2[i]/(vec1[i]+vec2[i])
  }
  for (i in 1:(K-1)){
    j = i+1
    R = R + (vec1[i]*vec2[j] + vec2[i]*vec1[j])/((vec1[i]+vec2[i])*(vec1[j]+vec2[j]))
  }
  R
}

# generate 2d distance matrix
getdist_2d = function(vec1){
  n = ceiling(sqrt(length(vec1)))
  distmat = matrix(0, n^2, n^2)
  for (i in 1:n){
    for (j in 1:n){
      for (k in 1:n){
        for (m in 1:n){
          id1 = n*(i-1)+j
          id2 = n*(k-1)+m
          distmat[id1, id2] = abs(i-k)+abs(j-m)
        }
      }
    }
  }
  return(distmat)
}
  
# calculate the nearest-neighbor statistic for ranked categorical data given two vectors.
nb = function(vec1, vec2){
  result = 0
  v = vec1 + vec2
  K = length(vec1)
  for (i in 1:K){
    result = result + vec1[i]*vec2[i]
    if (v[i]==1){
      if (i==1){
        b = min(which(v[2:K]>0))+i
        result = result + vec1[i]*vec2[b]+vec2[i]*vec1[b]
      }else if (i==K){
        a = max(which(v[1:(K-1)]>0))
        result = result + vec1[i]*vec2[a]+vec2[i]*vec1[a]
      }else{
        a = max(which(v[1:(i-1)]>0))
        b = min(which(v[(i+1):K]>0))+i
        if ((i-a)<(b-i)){
          result = result + vec1[i]*vec2[a]+vec2[i]*vec1[a]
        }else if((i-a)>(b-i)){
          result = result + vec1[i]*vec2[b]+vec2[i]*vec1[b]
        }else{
          result = result + vec1[i]*vec2[a]+vec2[i]*vec1[a] + vec1[i]*vec2[b]+vec2[i]*vec1[b]
        }
      }
    }
  }
  result
}

# calculate double factorial
dfac = function(n){
  if (n==0 ||n==1 || n==-1){
    return(1)
  }else{
    return(n*dfac(n-2))
  }
}

# calculate the basic cross-match statistic for one category.
cm0 = function(n1, n2){
  a = min(n1, n2)
  b = max(n1, n2)
  result = 0
  if (a%%2 == 0){
    for (i in seq(0,a,2)){
      result = result + i*factorial(b)/factorial(b-i)*factorial(a)/factorial(a-i)/factorial(i)*dfac(a-i-1)*dfac(b-i-1)
    }
  }else{
    for (i in seq(1,a,2)){
      result = result + i*factorial(b)/factorial(b-i)*factorial(a)/factorial(a-i)/factorial(i)*dfac(a-i-1)*dfac(b-i-1)
    }
  }
  return(result/dfac(a+b-1))
}

# test the cm0 is correct
testcm = function(a,b){
  result1 = result2 = result3 = 0
  if (a%%2 == 0){
    for (i in seq(0,a,2)){
      result1 = result1 + factorial(b)/factorial(b-i)*factorial(a)/factorial(a-i)/factorial(i)*dfac(a-i-1)*dfac(b-i-1)
      result2 = result2 + dfac(a-i-1)*dfac(b-i-1)*factorial(i)
    }
  }else{
    for (i in seq(1,a,2)){
      result1 = result1 + factorial(b)/factorial(b-i)*factorial(a)/factorial(a-i)/factorial(i)*dfac(a-i-1)*dfac(b-i-1)
      result2 = result2 + dfac(a-i-1)*dfac(b-i-1)*factorial(i)
    }
  }
  result3 = dfac(a+b-1)
  return(c(result1, result2, result3))
}

## cm1 = function(vec1, vec2){
##   n = length(vec1)
##   result = matrix(0,2,n)
##   for (i in 1:n){
##     if (vec1[i]==0){
##       result[1,i] = 0
##     }else{
##       result[1,i] = cm0(vec1[i]-1, vec2[i])
##     }
##     if (vec2[i]==0){
##       result[2,i] = 0
##     }else{
##       result[2,i] = cm0(vec1[i], vec2[i]-1)
##     }
##   }
##   return(result)
## }

# cm2: calculate the cross-match for two odd categories. n1A = vec1[1], n2A = vec1[2]
cm2 = function(vec1, vec2){
  v = vec1 + vec2
  result = vec1[1]*vec2[2] + vec2[1]*vec1[2]
  if (vec1[1]>0){
    result = result + vec1[1]*v[2]*cm0(vec1[1]-1,vec2[1])
  }
  if (vec2[1]>0){
    result = result + vec2[1]*v[2]*cm0(vec1[1], vec2[1]-1)
  }
  if (vec1[2]>0){
    result = result + v[1]*vec1[2]*cm0(vec1[2]-1, vec2[2])
  }
  if (vec2[2]>0){
    result = result + v[1]*vec2[2]*cm0(vec1[2], vec2[2]-1)
  }
  return(result/v[1]/v[2])
}


## cmx = function(x){
##   n = length(x)
##   count = 0
##   for (i in 1:(n/2)){
##     j = 2*i-1
##     k = 2*i
##     if (x[j]!=x[k]){
##       count = count + 1
##     }
##   }
##   return(count)
## }


## Ind = function(n){
##   if (n==1){
##     return(matrix(c(1,2), nrow=2))
##   temp = cate(oriv1,oriv2)
##   }else{
##     return(rbind( cbind(rep(1,dim(Ind(n-1))[1]), Ind(n-1)), cbind(rep(2,dim(Ind(n-1))[1]), Ind(n-1))) )
##   }
## }

# cmf: calculate the cross-match for two vectors
cmf = function(vec1, vec2){
  v = vec1 + vec2
  E = which(v%%2==0)
  O = which(v%%2!=0)
  result = 0
  n = length(O)
  if (n>0){
    for (i in 1:(n/2)){
      temp1 = vec1[O[(2*i-1):(2*i)]]
      temp2 = vec2[O[(2*i-1):(2*i)]]
      result = result + cm2(temp1, temp2)
    }
  }
  for (k in E){
    result = result + cm0(vec1[k],vec2[k])
  }
  return(result)
}

# generate a permutation 
perm = function(vec1, vec2){
  m = vec1 + vec2
  N = sum(m)
  X = rep(0, sum(m))
  id = sample(1:N, sum(vec1))
  X[id] = 1
  K = length(vec1)
  newvec1 = rep(0,K)
  newvec2 = rep(0,K)
  for (i in 1:K){
    j1 = sum(m[0:(i-1)])+1
    j2 = sum(m[0:i])
    newvec1[i] = sum(X[j1:j2])
    newvec2[i] = sum(1-X[j1:j2])
  }
  return(rbind(newvec1, newvec2))
}

# calculate right-side p-value
myp_right = function(vec, value){
  return((length(which(vec>=value-1e-10)))/length(vec))
}

# calculate left-side p-value
myp_left = function(vec, value){
  return((length(which(vec<=value+1e-10)))/length(vec))
}

# make continuous inputs to categories. vec1 and vec2 are x values.
cate = function(vec1, vec2){
  v = c(vec1, vec2)
  m = max(2, ceiling(length(v)/5))
  a = min(v)
  b = max(v)
  temp = seq(a,b,(b-a)/m)
  newv1 = rep(0,m)
  newv2 = rep(0,m)
  for (i in 1:(m-1)){
    newv1[i] = length(which(vec1[vec1>=temp[i]]<temp[i+1]))
    newv2[i] = length(which(vec2[vec2>=temp[i]]<temp[i+1]))
  }
  newv1[m] = length(vec1)-sum(newv1[1:(m-1)])
  newv2[m] = length(vec2)-sum(newv2[1:(m-1)])
  return(rbind(newv1,newv2))
}

# make continuous inputs to categories with the lower and upper bounds fixed (a & b) and the number of category (m) fixed.
cate2 = function(vec1, vec2, a, b, m){
  v = c(vec1, vec2)
  temp = seq(a,b,(b-a)/m)
  newv1 = rep(0,m)
  newv2 = rep(0,m)
  for (i in 1:(m-1)){
    newv1[i] = length(which(vec1[vec1>=temp[i]]<temp[i+1]))
    newv2[i] = length(which(vec2[vec2>=temp[i]]<temp[i+1]))
  }
  newv1[m] = length(vec1)-sum(newv1[1:(m-1)])
  newv2[m] = length(vec2)-sum(newv2[1:(m-1)])
  return(rbind(newv1,newv2))
}

# list the index of a vector that is '>=a' and '<b'
extract = function(vec, a, b){
  id1 = which(vec >= a)
  v1  = vec[id1]
  id2 = which(v1 < b)
  id  = id1[id2]
  return(id)
}

# make categories from 2-dimensional data
cate2d = function(x1, y1, x2, y2, catn=5){
  x = c(x1, x2)
  y = c(y1, y2)
  n = length(x)
  m = max(2, ceiling(sqrt(n/catn)))
  xmin = min(x)
  xmax = max(x)
  ymin = min(y)
  ymax = max(y)
  xtemp = seq(xmin, xmax, (xmax-xmin)/m)
  m1 = m2 = matrix(0,m,m)
  for (i in 1:(m-1)){
    v1 = y1[extract(x1, xtemp[i], xtemp[i+1])]
    v2 = y2[extract(x2, xtemp[i], xtemp[i+1])]
    result = cate2(v1, v2, ymin, ymax, m)
    m1[,i] = result[1,]
    m2[,i] = result[2,]
  }
  v1 = y1[extract(x1, xtemp[m], xmax+1)]
  v2 = y2[extract(x2, xtemp[m], xmax+1)]
  result = cate2(v1, v2, ymin, ymax, m)
  m1[,m] = result[1,]
  m2[,m] = result[2,]
  return(list(m1=m1, m2=m2))
}

# calculate the type II error given the p-value matrix and type I error.
roc = function(pmat, p0, N){
  pow = c()
  for (i in 1:dim(pmat)[2]){
    pow = c(pow, length(which(pmat[,i]<=p0))/N)
  }
  return(pow)
}

# generate all possible haplotypes given the length of haplotype
haplotypes = function(haplolength){
  if (haplolength == 1){
    return( c("A", "B"))
  }else{
    temp = haplotypes(haplolength-1)
    m = length(temp)
    temp2 = c()
    for (i in 1:m){
      temp2 = c(temp2, paste(temp[i], "A", sep=""), paste(temp[i], "B", sep=""))
    }
    return(temp2)
  }
}

# randomly generate a haplotype with the probability of having A in each site given.  The parameter 'p' is a vector with the length of the haplotype.
generate = function(p){
  n = length(p)
  result = ""
  for (i in 1:n){
    result = paste(result, generate_sub(p[i]), sep="")
  }
  return(result)
}

# randomly generate 'A' or 'B' given the probability of geting 'A'. Here, 'p' is a number between 0 and 1. 
generate_sub = function(p){
  if (runif(1)<p){
    return("A")
  }else{
    return("B")
  }
}

# randomly generate a vector with each element the count of the different haplotypes.
# 'p' is a vector with each element the probability of having 'A' at that location of the haplotype.
# 'n' is the sum of the vector elements.
generate_vec = function(p, hap, n){
  result = rep(0,length(hap))
  for (i in 1:n){
    temp = generate(p)
    id = which(hap == temp)
    result[id] = result[id] + 1
  }
  return(result)
}

# get the distance matrix given a vector containing haplotypes
getdist = function(hap){
  m = length(hap)
  result = matrix(0,m,m)
  for (i in 1:(m-1)){
    for (j in (i+1):m){
      result[i,j] = result[j,i] = getdist_sub(hap[i], hap[j])
    }
  }
  return(result)
}

# calculate the distance between two haplotypes.
getdist_sub = function(hap1, hap2){
  a = strsplit(hap1, "")[[1]]
  b = strsplit(hap2, "")[[1]]
  m = length(a)
  count = 0
  for (i in 1:m){
    if (a[i] != b[i]){
      count = count +1
    }
  }
  return(count)
}


# simplify distance matrix # this one does not work
simp_dist = function(distmat){
  n = dim(distmat)[1]
  newdist = distmat
  for (i in 1:n){
    idnz = which(distmat[i,] !=0)
    tmp = distmat[i,idnz]
    ids = which(tmp==min(tmp))
    tmp[-ids] = -1
    newdist[i, idnz] = tmp
  }
  return(newdist)
}


# calculate p-values
pfcn = function(control, case, mydist, B=1000, nnb=F){
  v = control + case
  ids = which(v!=0)
  v1 = control[ids]
  v2 = case[ids]
  K = length(v1)
  newdist = mydist[ids,ids]
  oriPear = pearson(v1, v2)
  oriLR = LR(v1, v2)
  submat = matrix(0, B, 2)
  permmat = matrix(0, K, 2*(B+1))
  permmat[,1] = v1
  permmat[,2] = v2
  pmt0 = proc.time()
  for (j in 1:B){
    temp = perm(v1, v2)
    v1new = temp[1,]
    v2new = temp[2,]
    submat[j,1] = pearson(v1new, v2new)
    submat[j,2] = LR(v1new, v2new)
    permmat[,2*j+1] = v1new
    permmat[,2*j+2] = v2new
  }
#  rawresult0 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(0), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(newdist))
  rawresult1 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(1), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(newdist))
  rawresult3 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(3), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(newdist))
  rawresult4 = .C("RG", as.double(rep(0,(B+1))), as.integer(0), as.double(rep(0,K*(K-1))), as.integer(4), as.integer(0), as.integer(K), as.integer(B+1), as.integer(permmat), as.double(newdist))
#  result0 = rawresult0[[1]]
  result1 = rawresult1[[1]]
  result3 = rawresult3[[1]]
  result4 = rawresult4[[1]]
  edge1 = t(matrix(rawresult1[[3]],2)[,1:rawresult1[[2]]]) + 1
  edge3 = t(matrix(rawresult3[[3]],2)[,1:rawresult3[[2]]]) + 1
  tmp1 = getMV(v1, v2, edge1)
  tmp3 = getMV(v1, v2, edge3)
  return(c(myp_left(result1[2:(B+1)], result1[1]), pnorm(result1[1], tmp1$Mean, tmp1$Sd), myp_left(result3[2:(B+1)], result3[1]), pnorm(result3[1], tmp3$Mean, tmp3$Sd), myp_left(result4[2:(B+1)], result4[1]), myp_right(submat[,2], oriLR), myp_right(submat[,1], oriPear)))
}

getMV = function(v1,v2,edge){
  K = length(v1)
  na = sum(v1)
  nb = sum(v2)
  m = v1 + v2
  N = na + nb
  p1 = na*nb/N/(N-1)
  p2 = 4*na*(na-1)*nb*(nb-1)/N/(N-1)/(N-2)/(N-3)
  nodedeg = rep(0,K)
  sumE = dim(edge)[1]
  quan = 0
  for (i in 1:sumE){
    e1 = edge[i,1]
    e2 = edge[i,2]
    nodedeg[e1] = nodedeg[e1]+1
    nodedeg[e2] = nodedeg[e2]+1
    quan = quan + 1/m[e1]/m[e2]
  }
  Mean = (N-K+sumE)*2*p1
  Var = 4*(p1-p2)*(N-K+2*sumE+sum(nodedeg^2/(4*m))-sum(nodedeg/m)) + (6*p2-4*p1)*(K-sum(1/m)) + p2*quan + (N-K+sumE)^2*(p2-4*p1^2)
  return(list(Mean=Mean,Sd=sqrt(Var)))
}

getMV_uMST = function(v1, v2, edge){
  K = length(v1)
  m = v1 + v2
  na = sum(v1)
  nb = sum(v2)
  N = na + nb
  p1 = na*nb/N/(N-1)
  p2 = 4*na*(na-1)*nb*(nb-1)/N/(N-1)/(N-2)/(N-3)
  nE = dim(edge)[1]
  Ebynode = vector("list", K)
  q1 = 0
  for (i in 1:nE){
    e1 = edge[i,1]
    e2 = edge[i,2]
    Ebynode[[e1]] = c(Ebynode[[e1]], e2)
    Ebynode[[e2]] = c(Ebynode[[e2]], e1)
    q1 = q1 + m[e1]*m[e2]
  }
  quan = 0
  for (i in 1:K){
    tmp = sum(m[Ebynode[[i]]])
    quan = quan + m[i]*(m[i]+tmp-1)*(m[i]+tmp-2)
  }
  nGbar = sum(m*(m-1)/2) + q1
  Mean = 2*p1*nGbar
  Var = p1*(2*nGbar+quan) + p2*(nGbar^2-nGbar-quan) - Mean^2
  return(list(Mean=Mean, Sd=sqrt(Var)))
}

## getMV_uMST = function(v1, v2, edge){
##   K = length(v1)
##   m = v1 + v2
##   na = sum(v1)
##   nb = sum(v2)
##   N = na + nb
##   p1 = na*nb/N/(N-1)
##   p2 = 4*na*(na-1)*nb*(nb-1)/N/(N-1)/(N-2)/(N-3)
##   nE = dim(edge)[1]
##   Ebynode = vector("list", K)
##   for (i in 1:nE){
##     e1 = edge[i,1]
##     e2 = edge[i,2]
##     Ebynode[[e1]] = c(Ebynode[[e1]], e2)
##     Ebynode[[e2]] = c(Ebynode[[e2]], e1)
##   }
##   q1 = sum(m*(m-1))
##   q2 = q3 = q4 = q5 = 0
##   q6 = sum(m*(m-1)^2)
##   for (i in 1:nE){
##     e1 = edge[i,1]
##     e2 = edge[i,2]
##     q2 = q2 + m[e1]*m[e2]
##     q3 = q3 + m[e1]*m[e2]*(m[e1]+m[e2])
##   }
##   for (i in 1:K){
##     nu = length(Ebynode[[i]])
##     if (nu > 1){
##       for (j in 1:(nu-1)){
##         for (l in (j+1):nu){
##           e1 = Ebynode[[i]][j]
##           e2 = Ebynode[[i]][l]
##           q4 = q4 + 2*m[i]*m[e1]*m[e2]
##         }
##       }
##     }
##     q5 = q5 + m[i]*(m[i]-1)*sum(m[Ebynode[[i]]])
##   }
##   Mean = p1*(q1+2*q2)
##   Var = (q1+2*q2)^2*(p2-4*p1^2)/4 + (p1-p2)*(q3-q4+2*q5+q6)+p2*(q1/2+q2)
##   return(list(Mean=Mean, Sd=sqrt(Var)))
## }
    
  
