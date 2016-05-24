RobRSVD <-
function(data, irobust = F, huberk = 1.345, 
                   iinituv = F, inits, initu, initv, 
                   niter = 1000, tol = 1e-5, istablize = T, 
                   uspar = 0, vspar = 0, iugcv = F, ivgcv = F,
                   usparmax = 1e4, usparmin = 1e-10, nuspar = 14,iloguspar = T, 
                   vsparmax = 1e4, vsparmin = 1e-10, nvspar = 14, ilogvspar = T){
  
  # Robust Functional Singular Value Decomposition
  
  ugcvmat = as.numeric()
  vgcvmat = as.numeric()
  
  # Check if the paramenters are appropriate. #
  
  size_data = c(dim(data))
  m = size_data[1]
  n = size_data[2]
  if (istablize){
    # The original data is centered by the median of the whole panel/matrix
    myscale = 1.4785*median(abs(c(data) - median(data)))
    localdata = data / myscale
  }else{
    myscale = 1
    localdata = data
  }
  
  if (!iinituv){
    # Classical SVD is applied on the localdata
    # Only the first u vector and v vector are retained 
    svdold = svd(localdata, 1, 1)
    uold = svdold$u
    vold = svdold$v
    sold = svdold$d[1]
  }else{
    uold = initu
    sold = inits
    vold = initv
    if (istablize)
      sold = sold / myscale
  }
  uold = sold * uold
  
  # Use the initial estimate of u and v vectors to reconstruct the data matrix
  Appold = uold %*% t(vold) 
  # Residuals (matrix) after the initial estimate
  Rmat = localdata - Appold
  Rvec = c(Rmat) # stack the Rmat columns
  mysigma = median(abs(Rvec))/0.675
  iter = 1
  localdiff = 9999
  diffvec = as.numeric()
  uspar_current = uspar
  vspar_current = vspar
  
  ugcvscore = as.numeric()
  vgcvscore = as.numeric()
  
  while (localdiff>tol & iter<niter){
    if (irobust){
      Wmat = huberWeightLS(Rmat/mysigma, huberk)
    }else{
      Wmat = matrix(1, m, n)
    }
    # update u first 
    if (!iugcv){
      uterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) + 
        (2*mysigma^2) * 
        (c(t(vold) %*% (diag(n) + vspar*ssmatls(n)$y) %*% vold) * (diag(m) + uspar*ssmatls(m)$y) - diag(sum(vold^2), m))
      uterm2 = (Wmat*localdata) %*% vold
      unew = solve(uterm1) %*% uterm2
    }else{ #choosing the penality parameter by cross-validation according to GCV criteria.
      if (nuspar<0)
        stop("number of smoothing parameter can not be negative")
      else{
        if (iloguspar){
          usparvec = 10^seq(log10(usparmin), log10(usparmax), length.out=nuspar)
        }else{
          usparvec = seq(usparmin, usparmax, length.out=nuspar)
        }
        ugcvvec = as.numeric()
        ugcvmat = as.numeric()
        for (iter_uspar in 1:nuspar){
          u_nsrobust = solve(diag(colSums(diag(c(vold^2)) %*% t(Wmat)))) %*% (Wmat*localdata) %*% vold
          usterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) + 
            (2*mysigma^2) *
            (c(t(vold) %*% (diag(n) + vspar_current*ssmatls(n)$y) %*% vold) * (diag(m)+usparvec[iter_uspar]*ssmatls(m)$y)-diag(sum(vold^2), m))
          usterm2 = (Wmat*localdata) %*% vold
          u_srobust = solve(usterm1) %*% usterm2
          smooth_u = solve(usterm1) %*% diag(colSums(diag(c(vold^2)) %*% t(Wmat)))
          gcv_ut = m * sum((u_nsrobust-u_srobust)^2) / (m - sum(diag(smooth_u)))^2
          ugcvvec = c(ugcvvec, gcv_ut)
          ugcvmat = cbind(ugcvmat, u_srobust/(sqrt(sum(u_srobust^2))))
        }
        uspar_current = usparvec[which.min(ugcvvec)]
        ugcvscore = cbind(usparvec, ugcvvec)
        
        uterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) + 
          (2*mysigma^2) *
          (c(t(vold) %*% (diag(n)+vspar_current*ssmatls(n)$y) %*% vold) * (diag(m)+uspar_current*ssmatls(m)$y) - diag(sum(vold^2), m))
        uterm2 = (Wmat*localdata) %*% vold
        unew = solve(uterm1) %*% uterm2 
      }#nuspar else
    }
    
    # update v secondly
    if (!ivgcv){
      vterm1 =  diag(colSums(diag(c(unew^2))%*%Wmat)) + 
        (2*mysigma^2) * 
        (c(t(unew)%*%(diag(m)+uspar*ssmatls(m)$y)%*%unew) * (diag(n)+vspar*ssmatls(n)$y) - diag(sum(unew^2), n))
      vterm2 = t(Wmat*localdata) %*% unew
      vnew = solve(vterm1) %*% vterm2
    }else{
      if (nvspar<0)
        stop("number of smoothing parameter can not be negative")
      else{
        if (ilogvspar){
          vsparvec = 10^seq(log10(vsparmin), log10(vsparmax), length.out=nvspar)
        }else{
          vsparvec = seq(vsparmin, vsparmax, length.out=nvspar)
        }
        vgcvvec = as.numeric()
        vgcvmat = as.numeric()
        for (iter_vspar in 1:nvspar){
          v_nsrobust = solve(diag(colSums(diag(c(unew^2))%*%Wmat))) %*% t(Wmat*localdata) %*% unew
          vsterm1 = diag(colSums(diag(c(unew^2))%*%Wmat)) + 
            (2*mysigma^2) * 
            (c(t(unew)%*%(diag(m)+uspar_current*ssmatls(m)$y)%*%unew) * (diag(n)+vsparvec[iter_vspar]*ssmatls(n)$y) - diag(sum(unew^2), n))
          vsterm2 = t(Wmat*localdata) %*% unew
          v_srobust = solve(vsterm1) %*% vsterm2
          smooth_v = solve(vsterm1) %*% diag(colSums(diag(c(unew^2))%*%Wmat))
          gcv_vt = n*sum((v_nsrobust-v_srobust)^2) / (n-sum(diag(smooth_v)))^2
          vgcvvec = c(vgcvvec, gcv_vt)
          vgcvmat = cbind(vgcvmat, v_srobust/sqrt(sum(v_srobust^2)))
        }
        vspar_current = vsparvec[which.min(vgcvvec)]
        
        vgcvscore = cbind(vsparvec, vgcvvec)
        
        vterm1 = diag(colSums(diag(c(unew^2))%*%Wmat)) + 
          (2*mysigma^2) * 
          (c(t(unew)%*%(diag(m)+uspar_current*ssmatls(m)$y)%*%unew) * (diag(n)+vspar_current*ssmatls(n)$y) - diag(sum(unew^2), n))
        vterm2 = t(Wmat*localdata) %*% unew
        vnew = solve(vterm1) %*% vterm2
      }# nvspar>=0
    }# ivgcv=T
    
    Appnew = unew %*% t(vnew)
    Rmat = localdata - Appnew
    
    localdiff = max(abs(Appnew-Appold))
    Appold = Appnew
    uold = sqrt(sum(vnew^2)) * unew
    vold = vnew / sqrt(sum(vnew^2)) # v is normalized to length 1
    
    iter = iter + 1
    diffvec = c(diffvec, localdiff)
  }# while
  
  v = vold
  s = myscale * sqrt(sum(uold^2))
  u = uold / sqrt(sum(uold^2))
  
  if (iugcv)
    uspar = uspar_current
  if (ivgcv)
    vspar = vspar_current
  
  diagout = list(ugcvscore=ugcvscore, vgcvscore=vgcvscore, ugcvmat=ugcvmat, vgcvmat=vgcvmat)
  return(list(s=s, v=v, u=u, diagout=diagout))
}
