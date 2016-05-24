"sample_sel" <-
function (data, sample = c(), sample_time = c(), sample_lambda = c(), 
    sel_time = c(), sel_lambda = c(), sel_time_ab = c(), 
    sel_lambda_ab = c(), sel_special = list() ) 
{
  increasing_x2 <- if(data@x2[1] < data@x2[2]) TRUE else FALSE 
  if (length(sel_time) == 2) {
    if (sel_time[2] > data@nt) 
      sel_time[2] <- data@nt
    data@psi.df <- as.matrix(data@psi.df[sel_time[1]:sel_time[2], ])
    data@x <- data@x[sel_time[1]:sel_time[2]]
    if (data@simdata) 
      data@C2 <- data@C2[sel_time[1]:sel_time[2], ]
  }
  if (length(sel_time_ab) == 2) {
    tmin <- which(data@x >= sel_time_ab[1])[1]
    tmax <- which(data@x <= sel_time_ab[2])[length(
                    which(data@x <= sel_time_ab[2]))]
    data@psi.df <- as.matrix(data@psi.df[tmin:tmax, ])
    data@x <- data@x[tmin:tmax]
    if (data@simdata) 
      data@C2 <- data@C2[tmin:tmax, ]
  }
  if (length(sel_lambda) == 2) {
    lmin <- sel_lambda[1]
    lmax <- sel_lambda[2]
    data@psi.df <- as.matrix(data@psi.df[, lmin:lmax])
    
    data@nl <- ncol(data@psi.df)
    data@x2 <- data@x2[lmin:lmax]
    if (data@simdata) 
      data@E2 <- data@E2[lmin:lmax, ]
  }
  if (length(sel_lambda_ab) == 2) {
    
    if(increasing_x2) {
      lmin <- which(data@x2 >= sel_lambda_ab[1])[1]
      lmax <- which(data@x2 <= sel_lambda_ab[2])[length(
                      which(data@x2 <= sel_lambda_ab[2]))]
    }
    else {
      lmin <- which(data@x2 <= sel_lambda_ab[1])[1]
      lmax <- which(data@x2 >= sel_lambda_ab[2])[length(
                      which(data@x2 >= sel_lambda_ab[2]))]
    }
    data@psi.df <- as.matrix(data@psi.df[,lmin:lmax])
    data@nt <- nrow(data@psi.df)
    
    data@x2 <- data@x2[lmin:lmax]
    if (data@simdata) 
      data@E2 <- data@E2[,lmin:lmax ]
    
  }
  if (length(sel_special) > 0) {
    for(i in 1:length(sel_special)) {
      reg <- sel_special[[i]]
      if(increasing_x2) {
        lmin <- which(data@x2 >= reg[1])[1]
        lmax <- which(data@x2 <= reg[2])[length(which(data@x2 <= reg[2]))]
      }
      else {
        lmin <- which(data@x2 <= reg[1])[1]
        lmax <- which(data@x2 >= reg[2])[length(which(data@x2 >= reg[2]))]
      }
      if(i == 1)
        datanew <- as.matrix(data@psi.df[- (reg[3]:reg[4]),lmin:lmax])
      else
        datanew <- cbind(datanew, as.matrix(data@psi.df[- (reg[3]:reg[4]),
                                                        lmin:lmax]))
    }
    data@psi.df <- datanew
    data@x <- data@x[-(sel_special[[1]][3]:sel_special[[1]][4])]
  }
  data@nl <- length(data@x2)
  data@nt <- length(data@x)
  if (sample != 1) {
    data@psi.df <- as.matrix(data@psi.df[seq(1, data@nt, by = sample), 
                                         seq(1, data@nl, by = sample)])
    data@x <- data@x[seq(1, data@nt, by = sample)]
    data@x2 <- data@x2[seq(1, data@nl, by = sample)]
    data@nl <- length(data@x2)
    data@nt <- length(data@x)
  }
  else {
    if (is.matrix(data@psi.df)) 
      data@psi.df <- as.matrix(data@psi.df[seq(1, data@nt,
                                               by = sample_time), 
                                           seq(1, data@nl,
                                               by = sample_lambda)])
    data@x <- data@x[seq(1, data@nt, by = sample_time)]
    data@x2 <- data@x2[seq(1, data@nl, by = sample_lambda)]
    data@nl <- length(data@x2)
    data@nt <- length(data@x)
  }
  data
}
