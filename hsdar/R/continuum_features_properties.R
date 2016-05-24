feature_properties <- function(x)
{
  feat_prop <- .area.specfeat(x)
  feat_prop <- cbind(feat_prop, .wlhm.specfeat(x))
  feat_prop <- cbind(feat_prop, .max.specfeat(x))
  attribute(x) <- cbind(attribute(x), feat_prop)
  usagehistory(x) <- "Properties of features calculated"
  return(x)
}

.area.specfeat <- function(
                          x
                         ) 
{
  if (class(x)!="Specfeat") 
    stop("x must be of class 'Specfeat'")
  
  x <- x@features
  
  Area <- matrix(ncol=length(x[[1]]),nrow=length(names(x)))
  
  for (i in 1:length(names(x)))
  {
    for (k in 1:length(x[[1]])) 
    {
      Area[i,k] <- as.numeric(sum(x[[i]][[k]]$y))
    }
  }
  
  rownames(Area) <- paste(names(x)[1:length(names(x))], "_area", sep = "")
  colnames(Area) <- names(x[[i]])
  return(as.data.frame(t(as.matrix(Area))))
}

.max.specfeat <- function(
                          x
                         ) 
{
  if (class(x)!="Specfeat") 
    stop("x must be of class 'Specfeat'")
  
  x <- x@features
  
  maxval <- matrix(ncol=length(x[[1]]),nrow=length(names(x)))
  
  for (i in 1:length(names(x)))
  {
    for (k in 1:length(x[[1]])) 
    {
      maxval[i,k] <- max(x[[i]][[k]]$y)
    }
  }
  
  rownames(maxval) <- paste(names(x)[1:length(names(x))], "_max", sep = "")
  colnames(maxval) <- names(x[[i]])
  return(as.data.frame(t(as.matrix(maxval))))
}  

.maxwl.specfeat <- function(
                          x
                         ) 
{
  if (class(x)!="Specfeat") 
    stop("x must be of class 'Specfeat'")
  wl <- wavelength(x)
  x <- x@features
  
  maxWL <- matrix(ncol=length(x[[1]]),nrow=length(names(x)))
  
  for (i in 1:length(names(x)))
  {
    for (k in 1:length(x[[1]])) 
    {
      tmp <- wl[c(which(wl == x[[i]][[k]]$x1):length(wl))]
      maxWL[i,k] <- tmp[which.max(x[[i]][[k]]$y)]
    }
  }
  
  rownames(maxWL) <- paste(names(x)[1:length(names(x))], "_maxwl", sep = "")
  colnames(maxWL) <- names(x[[i]])
  return(as.data.frame(t(as.matrix(maxWL))))
}  


.wlhm.specfeat <- function(
                          x, return_max = TRUE
                         ) 
{
  rmse <- function(
    obs,                  ## Vector of observed values
    pred,                 ## Vector of predicted values
    percent = FALSE       ## Return normalized root mean square error
  )
  {
    if (percent) sqrt(mean((obs-pred)^2))/(max(obs)-min(obs))*100 else sqrt(mean((obs-pred)^2)) 
  }
  
  if (class(x)!="Specfeat") 
    stop("x must be of class 'Specfeat'")
  m <- t(as.matrix(.maxwl.specfeat(x)))
  wl <- wavelength(x)
  x <- x@features
  
  col_fac <- 5
  
  wlhm <- matrix(ncol=length(x[[1]]),nrow=length(names(x))*col_fac)
  
  for (i in 1:length(names(x)))
  {
    for (k in 1:length(x[[1]])) 
    {
      tmp <- wl[c(which(wl == x[[i]][[k]]$x1):length(wl))]
      tmp <- tmp[1:length(x[[i]][[k]]$y)]
      maxval <- max(x[[i]][[k]]$y)
      tmp_1 <- 1:which(tmp == m[i, k])
      tmp_2 <- which(tmp == m[i, k]):length(tmp)
      tmp_x_1 <- tmp[tmp_1]
      tmp_x_2 <- tmp[tmp_2]
      tmp_y_1 <- x[[i]][[k]]$y[tmp_1]
      tmp_y_2 <- x[[i]][[k]]$y[tmp_2]
      wlhm[(i-1)*col_fac + 1,k] <- tmp_x_1[which.min(abs(tmp_y_1 - maxval/2))]
      wlhm[(i-1)*col_fac + 2,k] <- tmp_x_2[which.min(abs(tmp_y_2 - maxval/2))]
      wlhm[(i-1)*col_fac + 3,k] <- wlhm[(i-1)*col_fac + 2,k] - wlhm[(i-1)*col_fac + 1,k]
      gauss <- dnorm(tmp, mean = mean(c(wlhm[(i-1)*col_fac + 1,k],wlhm[(i-1)*col_fac + 2,k])), 
                     sd = (wlhm[(i-1)*col_fac + 3,k])/2)
      gauss <- (gauss-min(gauss))/(max(gauss)-min(gauss))
      gauss <- gauss * maxval
      wlhm[(i-1)*col_fac + 4,k] <- rmse(gauss[tmp_1], tmp_y_1)
      wlhm[(i-1)*col_fac + 5,k] <- rmse(gauss[tmp_2], tmp_y_2)
    }
  }
 
  rownames(wlhm) <- as.vector(sapply(names(x)[1:length(names(x))], function(i) paste(i, c("_lo", "_up", "_width", "gauss_lo", "gauss_up"), "_wlhm", sep = "")))
  colnames(wlhm) <- names(x[[i]])
  wlhm[!is.finite(wlhm)] <- 0
  if (return_max)
    wlhm <- rbind(m, wlhm)
  return(as.data.frame(t(as.matrix(wlhm))))
} 