vector.alpha <-
function(x, set, type="cor", CI=.95, CItype="xci", minval=-1.0) {
  comb.data <- data.frame(cbind(x, set))
  comp.data <- subset(comb.data, complete.cases(comb.data))
  N <- nrow(comp.data)
    # Alpha for vector of correlations
  if(type=="cor") {
    std.data <- data.frame(scale2(comp.data))
    Zcross.prod <- std.data$x * std.data[2:length(std.data)]
    tZcross.prod <- data.frame(t(Zcross.prod))
    Zcov.mat <- cov(tZcross.prod)
    Zcor.mat <- cor(tZcross.prod)
    Zavg.r <- mean(Zcor.mat[upper.tri(Zcor.mat)])
    Zalpha <- alpha.cov(Zcov.mat)
    Zalpha <- ifelse(Zalpha >= minval, Zalpha, minval)
    if(CItype=="xci") {
      ZCIs <- alpha.xci(Zalpha, k=N, n=ncol(set), CI=CI)
    }
    if(CItype=="aci") {
      ZCIs <- alpha.aci(Zalpha, k=N, n=ncol(set), CI=CI)
    } 
    out <- rbind(N, Zavg.r, Zalpha, ZCIs[1], ZCIs[2])
  }
    # Alpha for vector of covariances
  if(type=="cov") {
    cnt.data <- data.frame(scale2(comp.data, scale=F))
    Ccross.prod <- cnt.data$x * cnt.data[2:length(cnt.data)]
    tCcross.prod <- data.frame(t(Ccross.prod))
    Ccov.mat <- cov(tCcross.prod)
    Ccor.mat <- cor(tCcross.prod)
    Cavg.r <- mean(Ccor.mat[upper.tri(Ccor.mat)])
    Calpha <- alpha.cov(Ccov.mat)
    Calpha <- ifelse(Calpha >= minval, Calpha, minval)
    if(CItype=="xci") {
      CCIs <- alpha.xci(Calpha, k=N, n=ncol(set), CI=CI)
    }
    if(CItype=="aci") {
      CCIs <- alpha.aci(Calpha, k=N, n=ncol(set), CI=CI)
    }
    out <- rbind(N, Cavg.r, Calpha, CCIs[1], CCIs[2])
  }
    # Alpha for X-Y Beta Coefficients (X predicts each Y; same as covariances)
  if(type=="XY") {
    cnt.data <- data.frame(scale2(comp.data, scale=F))
    Ccross.prod <- cnt.data$x * cnt.data[2:length(cnt.data)]
    XYelems <- Ccross.prod / (var(comp.data$x)*(N-1)/N)
    tXYelems <- data.frame(t(XYelems))
    XYcov.mat <- cov(tXYelems)
    XYcor.mat <- cor(tXYelems)
    XYavg.r <- mean(XYcor.mat[upper.tri(XYcor.mat)])
    XYalpha <- alpha.cov(XYcov.mat)
    XYalpha <- ifelse(XYalpha >= minval, XYalpha, minval)
    if(CItype=="xci") {
      XYCIs <- alpha.xci(XYalpha, k=N, n=ncol(set), CI=CI)
    }
    if(CItype=="aci") {
      XYCIs <- alpha.aci(XYalpha, k=N, n=ncol(set), CI=CI)
    }
    out <- rbind(N, XYavg.r, XYalpha, XYCIs[1], XYCIs[2])
  }
    # Alpha for Y-X Beta Coefficients (Each Y predicts X)
  if(type=="YX") {
    cnt.data <- data.frame(scale2(comp.data, scale=F))
    y.vars <- diag(var(cnt.data[2:length(cnt.data)]))*(N-1)/N
    y.mat <- matrix(rep(y.vars, N), nrow=N, ncol=length(y.vars), byrow=T)
    YXelems <- cnt.data$x * cnt.data[2:length(cnt.data)] / y.mat
    tYXelems <- data.frame(t(YXelems))
    YXcov.mat <- cov(tYXelems)
    YXcor.mat <- cor(tYXelems)
    YXavg.r <- mean(YXcor.mat[upper.tri(YXcor.mat)])
    YXalpha <- alpha.cov(YXcov.mat)
    YXalpha <- ifelse(YXalpha >= minval, YXalpha, minval)
    if(CItype=="xci") {
      YXCIs <- alpha.xci(YXalpha, k=N, n=ncol(set), CI=CI)
    }
    if(CItype=="aci") {
      YXCIs <- alpha.aci(YXalpha, k=N, n=ncol(set), CI=CI)
    }    
    out <- rbind(N, YXavg.r, YXalpha, YXCIs[1], YXCIs[2])
  }
  colnames(out) <- c("Results")
  rownames(out) <- c("N", "Average r", "Alpha", "Lower Limit", "Upper Limit")
  return(out)
}
