expscore <- function(pers_obj, na_treat=NA){
  # returns a matrix with dims like resp matrix with expected scores and more ...
  # func. by joerg-henrik heine jhheine(at)googlemail.com
  # needs func. \code{pvx} in i.pvx.R and \code{pvx.matrix} in i.pvx.matrix.R
  # Notation and formulas see:  Wright & Masters 1982 p.100
  # in a revised form  (korrigendum) see http://www.rasch.org/rmt/rmt34e.htm
  nitems <- dim(pers_obj$pair$threshold)[1]
  emp_resp <- pers_obj$pair$resp # empirical responses
  N <- dim(pers_obj$pair$resp)[1]
    
  # Expected Mean Eni (of xni)-----------------
  foo1 <- function(theta, pers_obj=pers_obj, nitems=nitems){
    #theta: ein theta wert
    catL <- lapply(1:nitems,function(x){(0:length(na.omit(pers_obj$pair$threshold[x,])))})
    pL <- lapply(1:nitems,function(x){pvx(theta=theta,thres=na.omit(pers_obj$pair$threshold[x,]))})
    exp_scor_theta <- mapply(FUN=function(o1,o2){sum(o1*o2)}, o1=catL, o2=pL)
    return(exp_scor_theta)
  }
  Eni <- t(sapply(pers_obj$pers$WLE, foo1, pers_obj ,nitems))
  dimnames(Eni) <- dimnames(emp_resp)
  ### replace cells missing in emp_resp by NA in Eni new approach!----
  # affects Wni and Cni and all the other _ni's too
  Eni[is.na(emp_resp)] <- NA 
  # Eni[is.na(emp_resp)] <- 0 # nicht nötig passiert weiter unten

  # matplot(Eni[order(pers_obj$pers$WLE),],type="l") # ggf. noch machen
  # return(Eni)
  
  ## Variance Wni (of xni)--------
  foo2 <- function(i, pers_obj=pers_obj, Eni=Eni){
  theta_v <- pers_obj$pers$WLE # empirical theta vector
  thres <- na.omit(pers_obj$pair$threshold[i,]) # thurst. thresholds item i 
  kat <- 0:length(na.omit(pers_obj$pair$threshold[i,])) # categories item i
  xm_v <- pers_obj$pair$resp[,i]+1  # responses item i | +1 ist wichtig es heißt 1 und 2 kategorie nicht 0 und 1
  pnik <- (pvx.matrix(theta_v,thres)); dimnames(pnik)[[2]] <- names(xm_v)
  Lpnik <- as.list(data.frame(t(pnik)))
  Leni <- lapply(kat, function(x){  (x-Eni[,i])^2  })# Eni[,i] # expected scores item i
  wni <- rowSums(mapply(function(p,e){e*p}, p=Lpnik, e=Leni,SIMPLIFY = TRUE))
  return(wni)
  }
  Wni <- sapply(1:nitems, foo2, pers_obj, Eni)
  dimnames(Wni) <- dimnames(emp_resp)
  # return(Wni)
  Wni[is.na(emp_resp)] <- na_treat # new
  
  ## Kurtosis Cni (of xni)--------
  foo3 <- function(i, pers_obj=pers_obj, Eni=Eni){
    theta_v <- pers_obj$pers$WLE # empirical theta vector
    thres <- na.omit(pers_obj$pair$threshold[i,]) # thurst. thresholds item i 
    kat <- 0:length(na.omit(pers_obj$pair$threshold[i,])) # categories item i
    xm_v <- pers_obj$pair$resp[,i]+1  # responses item i | +1 ist wichtig es heißt 1 und 2 kategorie nicht 0 und 1
    pnik <- (pvx.matrix(theta_v,thres)); dimnames(pnik)[[2]] <- names(xm_v)
    Lpnik <- as.list(data.frame(t(pnik)))
    Leni <- lapply(kat, function(x){  (x-Eni[,i])^4  })# Eni[,i] # expected scores item i
    cni <- rowSums(mapply(function(p,e){e*p}, p=Lpnik, e=Leni,SIMPLIFY = TRUE))
    return(cni)
  }
  Cni <- sapply(1:nitems, foo3, pers_obj, Eni)
  dimnames(Cni) <- dimnames(emp_resp)
  # return(Eni)
  
  ###############------------------------------------------------------------
  
  ## Score Residual  Yni--------
  Yni <- emp_resp-Eni
  Yni[is.na(emp_resp)] <- na_treat
  
  ## Standardised Residual  Zni--------
  Zni <- Yni / sqrt(Wni) 
  Zni[is.na(emp_resp)] <- na_treat
  
  ## Score Residual Squared Y2ni--------
  Y2ni <- Wni* (Zni)^2

  ## Standardised Residual Squared Z2ni--------
  Z2ni <- (Zni)^2
  
  return(list(Eni=Eni, Wni=Wni, Cni=Cni, Yni=Yni, Zni=Zni, Y2ni=Y2ni, Z2ni=Z2ni))  
}

