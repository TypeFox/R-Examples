fitCircle <- function(x, y){ 
  n <- length(x) 
  centroid <- cbind(mean(x),mean(y)) 
  Mxx <- 0; Myy <- 0; Mxy <- 0; Mxz <- 0; Myz <- 0; Mzz <- 0; 
  
  for(i in 1:n){ 
    Xi <- x[[i]] - centroid[1] 
    Yi <- y[[i]] - centroid[2] 
    Zi <- (Xi*Xi) + (Yi*Yi) 
    Mxy  <-  Mxy + Xi*Yi 
    Mxx  <-  Mxx + Xi*Xi 
    Myy  <-  Myy + Yi*Yi 
    Mxz  <-  Mxz + Xi*Zi 
    Myz  <-  Myz + Yi*Zi 
    Mzz  <-  Mzz + Zi*Zi 
  } 
  
  Mxx  <-  Mxx/n 
  Myy  <-  Myy/n 
  Mxy  <-  Mxy/n 
  Mxz  <-  Mxz/n 
  Myz  <-  Myz/n 
  Mzz  <-  Mzz/n 
  
  # computing the coefficients of the characteristic polynomial 
  Mz  <-  Mxx + Myy 
  Cov_xy  <-  Mxx*Myy - Mxy*Mxy
  Mxz2  <-  Mxz*Mxz 
  Myz2  <-  Myz*Myz 
  
  A2  <-  4*Cov_xy - 3*Mz*Mz - Mzz 
  A1  <-  Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz 
  A0  <-  Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy
  A22  <-  A2 + A2 
  epsilon <- 1e-12 
  ynew <- 1e+20 
  IterMax <- 20 
  xnew  <-  0 
  
  # Newton's method starting at x=0 
  epsilon <- 1e-12 
  ynew <- 1e+20 
  IterMax <- 20 
  xnew  <-  0 
  iter <- 1:IterMax 
  
  for (i in 1:IterMax){ 
    yold  <-  ynew 
    ynew  <-  A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew))
    if (abs(ynew) > abs(yold)){ 
      print('Newton-Pratt goes wrong direction: |ynew| > |yold|') 
      xnew  <-  0
      break 
    } 
    Dy  <-  A1 + xnew*(A22 + 16*xnew*xnew)
    xold  <-  xnew; 
    xnew  <-  xold - ynew/Dy; 
    if (abs((xnew-xold)/xnew) < epsilon) {break} 
    if(iter[[i]] >= IterMax){ 
      print('Newton-Pratt will not converge')
      xnew  <-  0 
    } 
    if(xnew < 0.){ 
      print('Newton-Pratt negative root:  x=',xnew)
    } 
  } 
  
  DET  <-  xnew*xnew - xnew*Mz + Cov_xy 
  Center  <-  cbind(Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy)/DET/2
  #    computing the circle parameters 
  DET  <-  xnew*xnew - xnew*Mz + Cov_xy
  Center  <-  cbind(Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy)/DET/2 
  coor <- as.vector(Center+centroid)
  
  result <- list(cx = coor[1], cy = coor[2], radius = sqrt(Center[2]*Center[2]+Mz+2*xnew))
  return(result) 
}