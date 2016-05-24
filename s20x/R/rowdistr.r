prepCrosstabList = function(crosstablist){
  crosstablist  =  list(row.props = sweep(crosstablist, 1, apply(crosstablist, 1, sum), "/"), 
                        col.props = sweep(crosstablist, 2, apply(crosstablist, 2, sum), "/"), 
                        whole.props = crosstablist/sum(apply(crosstablist, 2, sum)), 
                        Totals = rbind(cbind(crosstablist, apply(crosstablist, 1, sum)), 
                                       apply(cbind(crosstablist, apply(crosstablist, 1, sum)), 2, sum)))
  
  if (is.null(dimnames(crosstablist$row.props)))
    dimnames(crosstablist$row.props)  =  list(c(1:nrow(crosstablist$whole.props)),
                                              c(1:ncol(crosstablist$whole.props)))
  
  if (is.null(dimnames(crosstablist$col.props)))
    dimnames(crosstablist$col.props)  =  dimnames(crosstablist$row.props)
  
  if (is.null(dimnames(crosstablist$whole.props)))
    dimnames(crosstablist$whole.props)  =  dimnames(crosstablist$row.props)
  
  if (!is.null(dimnames(crosstablist$Totals)[[1]]) & is.null(dimnames(crosstablist$Totals)[[2]])) {
    dimnames(crosstablist$Totals)  =  list(c(dimnames(crosstablist$Totals)[[1]][1:nrow(crosstablist$whole.props)],
                                             "Total"), c(1:ncol(crosstablist$whole.props),
                                                         "Total"))
    dimnames(crosstablist$row.props)  =  list(c(dimnames(crosstablist$Totals)[[1]][1:nrow(crosstablist$row.props)]),
                                              c(1:ncol(crosstablist$row.props)))
    dimnames(crosstablist$col.props)  =  dimnames(crosstablist$row.props)
  }
  
  if (is.null(dimnames(crosstablist$Totals)))
    dimnames(crosstablist$Totals)  =  list(c(1:nrow(crosstablist$whole.props),
                                             "Total"), c(1:ncol(crosstablist$whole.props),
                                                         "Total"))
  if (is.null(names(dimnames(crosstablist$Totals))))
    names(dimnames(crosstablist$Totals))  =  c("fac1", "fac2")
  
  
  return(crosstablist)
}

drawPlot = function(crosstablist, comp = "basic", conf.level = 0.95){
  
  ad  =  0.5
  rowcol  =  1
  colrect  =  rowcol
  labs  =  "vert"
  
  propmat  =  crosstablist$row.props
  n  =  nrow(propmat)
  opar = NULL
  
  if (comp == "basic" | comp == "within")
    opar = par(mfrow = c(n, 1), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(3, 3, 4, 3))
  
  if (comp == "between")
    opar = par(mfrow = c(1, 1))
  
  chiuse  =  chiadd  =  NULL
  pvaluse  =  pvaladd  =  NULL
  cexval  =  2 / n
  
  for (rowcol in 1:n) {
    maxpl  =  1.17 * max(propmat)
    plusmin  =  0.3
    xuse  =  propmat[rowcol, ]
    len  =  length(xuse)
 
    if (comp == "basic" | comp == "within") {
      plot(c(1 - 2 * plusmin, len + 2 * plusmin), c(0,
                                                    maxpl), type = "n", bty = "n", axes = FALSE, xlab = "",
           ylab = "proportion", cex.lab = 1.15)
      axis(2, las = 2)
      axis(1, at = 1:len, labels = substr(colnames(propmat),
                                          1, 8), cex.axis = 1.15)
      for (i in 1:len) {
        rect(i - plusmin, 0, i + plusmin, xuse[i], col = rowcol)
        box()
      }
    }
    if (comp == "basic") {
      text((len + 1)/2, maxpl - 0.12 * (max(propmat)),
           paste(rownames(propmat)[rowcol], " (n = ", crosstablist$Totals[,
                                                                          ncol(crosstablist$Totals)][rowcol], ")", sep = ""),
           cex = 1.15)
    }
    E  =  crosstablist$Totals[, ncol(crosstablist$Totals)][1:nrow(crosstablist$row.props)][rowcol]/ncol(crosstablist$col.props)
    dfs  =  ncol(crosstablist$col.props) - 1
    chi  =  round(sum((crosstablist$Totals[1:nrow(crosstablist$whole.props),
                                           1:ncol(crosstablist$whole.props)][rowcol, ] - E)^2/E),
                  3)
    chiadd  =  chi
    chiuse  =  cbind(chiuse, chiadd)
    pval  =  round((1 - pchisq(chi, dfs)), 5)
    pvaladd  =  pval
    pvaluse  =  cbind(pvaluse, pvaladd)
   
    if (comp == "within") {
      text((len + 1)/2, maxpl - 0.12 * (max(propmat)),
           paste(rownames(propmat)[rowcol], "(n = ", crosstablist$Totals[,
                                                                         ncol(crosstablist$Totals)][rowcol], ")", ", ",
                 "uniformity p-value=", pval, sep = ""), cex = 1.15)
    }
  }
  
  if (comp == "basic" | comp == "within") {
    mtext(paste(names(dimnames(crosstablist$Totals))[2],
                "distribution for each", sep = " "), outer = TRUE, at = 0.5,
          line = 2, cex = 1.15, adj = ad, font = 2)
    mtext(paste("level of", names(dimnames(crosstablist$Totals))[1],
                "(row proportions)", sep = " "), outer = TRUE, at = 0.5,
          cex = 1.15, adj = ad, font = 2)
  }

  if (comp == "between")
    propslsd.new(crosstablist, conf.level = 1-(1-conf.level)/choose(n,2))
  
  par(opar)
  
}

printOutput = function(crosstablist, comp = "basic", conf.level = 0.95){
  
  cat("Row Proportions\n")

  propmat  =  crosstablist$row.props
  n  =  nrow(propmat)
  
  matprint  =  cbind(round(propmat, 2), apply(propmat, 1, sum),
                     as.numeric(crosstablist$Totals[, ncol(crosstablist$Totals)][1:(nrow(crosstablist$Totals) -
                                                                                      1)]))
  dimnames(matprint)  =  list(dimnames(propmat)[[1]], c(dimnames(propmat)[[2]],
                                                        "Totals", "n"))
  print(matprint)
  
  chiuse  =  chiadd  =  NULL
  pvaluse  =  pvaladd  =  NULL
  cexval  =  2 / n
  
  p  =  propmat
  ns  =  crosstablist$Totals[, ncol(crosstablist$Totals)][1:(nrow(crosstablist$Totals) -
                                                               1)]
  matw  =  matrix(NA, ncol(propmat) - 1, ncol(propmat) - 1)
  mat  =  matrix(NA, n - 1, n - 1)
  colvar  =  names(dimnames(crosstablist$Totals))[2]
  rowvar  =  names(dimnames(crosstablist$Totals))[1]
  name  =  dimnames(p)[[1]]
  dimnames(mat)  =  list(name[-n], name[-1])
  namew  =  dimnames(p)[[2]]
  dimnames(matw)  =  list(namew[-length(namew)], namew[-1])
  
  for (j in 1:ncol(p)) {
    for (i1 in 1:(n - 1)) {
      for (i2 in 2:n) {
        zCrit = abs(qnorm((1 - conf.level)/(2*choose(n,2))))
        seDiff = sqrt(p[i1, j] * (1 - p[i1, j])/ns[i1] + p[i2, j] * (1 - p[i2, j])/ns[i2])
        temp  =  p[i1, j] - p[i2, j] +  zCrit * c(-1, 1) * seDiff
        temp  =  round(temp, 3)
        mat[i1, i2 - 1]  =  ifelse((i1 < i2), paste("(",
                                                    temp[1], ",", temp[2], ")", sep = ""), " ")
        if ((0 <= temp[1] | 0 >= temp[2]) & (i1 < i2)) {
          mat[i1, i2 - 1]  =  paste(mat[i1, i2 - 1], "*",
                                    sep = "")
        }
      }
    }
    
    if (comp == "between") {
      cat("\n")
      cat(paste(paste(100*conf.level,"%",sep=""), "CIs for diffs between proportions with",
                colvar, "=", dimnames(p)[[2]][j]), "\n")
      cat("(rowname-colname)", "\n")
      print(mat, quote = FALSE)
    }
  }
  
  for (k in 1:nrow(p)) {
    for (i1 in 1:(ncol(p) - 1)) {
      for (i2 in 2:ncol(p)) {
        zCrit = abs(qnorm((1 - conf.level)/(2*choose(ncol(p),2))))
        seDiff = sqrt(((p[k, i1] + p[k, i2]) - ((p[k, i1] - p[k, i2])^2))/ns[k])
        tempw  =  p[k, i1] - p[k, i2] + zCrit * c(-1, 1) * seDiff
        tempw  =  round(tempw, 3)
        matw[i1, i2 - 1]  =  ifelse((i1 < i2), paste("(",
                                                     tempw[1], ",", tempw[2], ")", sep = ""), " ")
        if ((0 <= tempw[1] | 0 >= tempw[2]) & (i1 < i2)) {
          matw[i1, i2 - 1]  =  paste(matw[i1, i2 - 1],
                                     "*", sep = "")
        }
      }
    }
    if (comp == "within") {
      print(matprint[k, ])
      cat("\n")
      cat(paste("Chisq test for uniformity:", "chisq = ",
                chiuse[k], ",", "df =", ncol(propmat) - 1, ",",
                "p-value =", pvaluse[k]), "\n")
      cat(paste(paste(100*conf.level,"%",sep=""),"CIs for diffs in propns within the",
                rowvar, "=", dimnames(p)[[1]][k]), "distribution",
          "\n")
      cat("(rowname-colname)", "\n")
      print(matw, quote = FALSE)
      cat("\n")
      cat("----------------------------------------------------------------------")
      cat("\n")
    }
  }

  invisible(propmat)
}


rowdistr = function (crosstablist, comp = "basic", conf.level = 0.95,
                    plot = TRUE, suppressText = FALSE){

      
    if (!is.matrix(crosstablist) && (class(crosstablist) != "ct.20x"))
        stop("check form of crosstablist: input list must be output from crosstabs or a list of similar form")
    
    if (is.matrix(crosstablist) & !is.list(crosstablist)) {
      crosstablist = prepCrosstabList(crosstablist)
    }
    
    propmat = NULL
    
    if(!suppressText){
      propmat = printOutput(crosstablist, comp, conf.level)
    }
    
    if(plot)
      drawPlot(crosstablist, comp, conf.level)
   
    invisible(propmat)
}

