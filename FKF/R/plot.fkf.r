plot.fkf <- function(x, type = c("state", "resid.qq", "qqchisq", "acf"),
                     CI = 0.95,
                     at.idx = 1:nrow(x$at), att.idx = 1:nrow(x$att), ...)
{
  type <- match.arg(type)

  d <- nrow(x$vt)
  n <- ncol(x$att)

  distance <- sapply(1:n, function(i, Ft, vt)
                     {
                       if(any(is.na(vt[, i])))
                         return(matrix(NA))
                       else
                         return(t(vt[,i]) %*% solve(Ft[,,i]) %*% vt[,i])
                     },
                     Ft = x$Ft, vt = x$vt)

  std.resid <- sapply(1:n, function(i, Ft, vt)
                      {
                        if(any(is.na(vt[, i])))
                          return(matrix(NA, ncol = 1, nrow = nrow(vt)))
                        else
                          return(solve(t(chol(Ft[,,i]))) %*% vt[,i])
                      },
                      Ft = x$Ft, vt = x$vt)

  dim(std.resid) <- c(d, n) # In case d == 1, make std.resid to be a matrix.
  ##< ------------------------------------------------------------ >
  if(type == "state"){
    xlim <- 1:(n + 1)
    xlim.att <- 1:n                     # Nbr. of predictions - 1

    ylim <- range(x$att[att.idx, ], x$at[at.idx, ], na.rm = TRUE)
    plot(xlim, ylim = ylim, type = "n", xlab = "Index",
         ylab = "State variables", ...)

    if(any(is.na(at.idx))){
      at.legend <- character(0)
      col.at <- character(0)
      at.idx <- numeric(0)
    }else{
      col.at <- rainbow(length(at.idx))
      for(i in 1:length(at.idx)){
        lines(xlim, x$at[at.idx[i],],
              col = col.at[i], lty = "dashed")
        if(is.finite(CI)){
          lines(xlim, x$at[at.idx[i],] + qnorm(0.5 - CI / 2) * sqrt(x$Pt[at.idx[i], at.idx[i], ]),
                col = col.at[i], lty = "dotted")
          lines(xlim, x$at[at.idx[i],] + qnorm(0.5 + CI / 2) * sqrt(x$Pt[at.idx[i], at.idx[i], ]),
                col = col.at[i], lty = "dotted")
        }
      }
      at.legend <- paste("at[", at.idx, ",]", sep = "")
    }

    if(any(is.na(att.idx))){
      att.legend <- character(0)
      col.att <- character(0)
      att.idx <- numeric(0)
    }else{
      col.att <- rainbow(length(att.idx))
      for(i in 1:length(att.idx)){
        lines(xlim.att, x$att[att.idx[i],],
              col = col.att[i], lty = "solid")
        if(is.finite(CI)){
          lines(xlim.att, x$att[att.idx[i],] + qnorm(0.5 - CI / 2) * sqrt(x$Ptt[att.idx[i], att.idx[i], ]),
                col = col.att[i], lty = "dotdash")
          lines(xlim.att, x$att[att.idx[i],] + qnorm(0.5 + CI / 2) * sqrt(x$Ptt[att.idx[i], att.idx[i], ]),
                col = col.att[i], lty = "dotdash")
        }

      }
      att.legend <- paste("att[", att.idx, ",]", sep = "")
    }

    legend("topleft", legend = c(at.legend, att.legend),
           lty = c(rep("dashed", length(at.idx)),
             rep("solid", length(att.idx))),
           col = c(col.at, col.att), bg = "white")

    if(is.finite(CI)){
      legend("topright", legend = c(att.legend, at.legend),
             lty = c(rep("dotted", length(at.idx)),
               rep("dotdash", length(att.idx))),
             col = c(col.at, col.att),
             title = paste("Confidence interval:", CI), bg = "white")
    }
    ## < ------------------------------------------------------------ >
  }else if(type == "resid.qq"){
    if(d < 3){
      par(mfrow = c(d, 1))
    }else{
      par(mfrow = c(d %/% 3 + (d %% 3 > 0), 3))
    }
    for(i in 1:d){
      qqnorm(std.resid[i, ],
             main = paste("Normal Q-Q Plot of residuals 'vt[",i,",]'", sep = ""),
             ...)
    }
    ## < ------------------------------------------------------------ >
  }else if(type == "qqchisq"){
    q.distance <- qchisq(ppoints(n), df = d)[order(order(distance))]
    qqplot(q.distance, distance,
           main = expression(paste(chi[2]^n,
               " Q-Q Plot of transformed residuals Dt = t(vt) * Ft^(-1) * vt")),
           xlab = "Theoretical Quantiles",
           ylab = "Sample Quantiles", ...)
    ## < ------------------------------------------------------------ >
  }else if(type == "acf"){
    par(mfrow = c(d, d))
    for(i in 1:d){
      for(j in 1:d){
        if(i == j){
          acf(x$vt[i, ],
              main = paste("ACF: vt[", i,", ] & vt[", i,", ]", sep = ""), ...)
        }else{
          ccf(x$vt[i, ], x$vt[j, ],
              main = paste("CCF: vt[", i,", ] & vt[", j,", ]", sep = ""), ...)
        }
      }
    }

  }
  invisible(list(std.resid = std.resid, distance = distance))
}

