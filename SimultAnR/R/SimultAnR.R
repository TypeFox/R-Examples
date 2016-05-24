################################################################################
#R# #####               SimultAnR (Simultaneous Analysis) Package              #####
#S# #####               SimultAn  (Simultaneous Analysis) Package              #####
################################################################################
### The package includes six functions:
### CorrAn <- function(data, sr = NA, sc = NA, nd = 2, dp = 2)
### SimAn <- function(data, G, acg, weight = 2, nameg = NA, sr = NA, sc = NA,
###    nd = 2, dp = 2)
### summary.CorrAn <- function(object, ...)
### summary.SimAn <- function(object, ...)
### plot.CorrAn <- function(x, s1 = 1, s2 = 2, ...)
### plot.SimAn <- function(x, s1 = 1, s2 = 2, ...)

#R# Cambiar ginverse por solve
#R# Cambiar is.missing por is.na
#R# Cambiar CAres[, g] por CAres[[g]]
#R# Cambiar CAres[, h] por CAres[[h]]
#R# Cambiar CAres[, 1] por CAres[[1]]

################################################################################
#####               CorrAn function: Correspondence Analysis               #####
################################################################################
CorrAn <- function(data, sr = NA, sc = NA, nd = 2, dp = 2, oar = 1, oac = 1, 
   multiple = 0)
{
# -------------------------------------------------------------------
#  The CorrAn function computes the correspondence analysis of the
#  selected data.
#
# Arguments
#   data     Data set
#   sr       Indices of supplementary rows
#   sc       Indices of supplementary columns
#   nd       Number of dimensions in results
#   dp       Number of digits in results
#   oar      Output for active rows
#   oac      Output for active columns
#   multiple     ZZZ 
#
# See also
#   summary.CorrAn
#   plot.CorrAn
#
# -------------------------------------------------------------------
# Reference
#      Zarraga, A. and Goitisolo, B.
#        Simultaneous Analysis is S-PLUS. The SimultAn Package.
#        submitted for publication, 2010.
# -------------------------------------------------------------------
cat("\n \n")
   ### some legends
   if(multiple == 0 ) {
      legend100fi <- "100fi"
      legendF<- "F"
      legend100fj <- "100fj"
      legendG <- "G" 
      legendsr <- "No supplementary rows"
      legendsc <- "No supplementary columns" 
      }
   ### multiple
   if(multiple == 1 ) {
#D      temp <- sr   
#D      sr <- sc
#D      sc <- temp
#D      temp <- oar
#D      oar <- oac 
#D      oac <- temp 
#D      data <- t(data)
      legend100fi <- "100fj"
      legendF<- "G"
      legend100fj <- "100fi"
      legendG <- "F" 
      legendsr <- "No supplementary columns"
      legendsc <- "No supplementary rows" 
      }
   ### input
   datafa <- as.matrix(data)
   ### Isr
   if(is.na(sr[1])) {
      Isr <- 0
      I <- nrow(datafa)   }
   else {
      Isr <- length(sr)
      I <- nrow(datafa) - Isr   }
   ### Jsc
   if(is.na(sc[1])) {
      Jsc <- 0
      J <- ncol(datafa)   }
   else {
      Jsc <- length(sc)
      J <- ncol(datafa) - Jsc   }
   ### Active data
   if(Isr != 0) {
      if(Jsc != 0) {
         datafaA <- datafa[ - sr, - sc]   }
      else {
         datafaA <- datafa[ - sr, ]   }   }
   else {
      if(Jsc != 0) {
         datafaA <- datafa[, - sc]   }
      else {
         datafaA <- datafa   }   }
   namei <- labels(datafaA)[[1]]
   namej <- labels(datafaA)[[2]]
   ### supplementary rows
   nameisr <- labels(datafa)[[1]][sr]
   datafasr <- matrix(0, Isr, J)
   if(Isr != 0) {
      if(Jsc != 0) {
         datafasr[, ] <- as.matrix(datafa[sr, - sc])   }
      else {
         datafasr[, ] <- datafa[sr, ]   }   }
   else {
      datafasr[, ] <- 0   }
   ### supplementary columns
   namejsc <- labels(datafa)[[2]][sc]
   datafasc <- matrix(0, I, Jsc)
   if(Jsc != 0) {
      if(Isr != 0) {
         datafasc[, ] <- datafa[ - sr, sc]   }
      else {
         datafasc[, ] <- datafa[, sc]   }   }
   else {
      datafasc[, ] <- 0   }
   ### active analysis
   ### frequency matrix
   Fd <- datafaA/sum(datafaA)
   fi <- apply(Fd, 1, sum)
   fj <- apply(Fd, 2, sum)
   I <- length(fi)
   J <- length(fj)
   ### Number of axes
   ndim <- min(I, J, nd)
   ########## ########## ########## ########## ########## ########## ###########
   ### Matrix to diagonalize
   temp <- sweep(Fd, 1, fi, FUN = "/")
   temp <- replace(temp, is.na(temp), 0)
   temp <- sweep(temp, 2, fj, FUN = "/")
   temp <- replace(temp, is.na(temp), 0)
   temp <- sweep(temp - 1 , 1, sqrt(fi), FUN = "*")
   X <- sweep(temp , 2, sqrt(fj), FUN = "*")
   if(I >= J) {
      A <- t(X) %*% X   }
   else {
      A <- (X) %*% t(X)   }
   ### eigen values, inertias and percentages
   eigall <- eigen(A)$values
   totalin <- sum(eigall)
   names(totalin) <- "Total inertia"
   porin <- eigall/totalin
   poracuin <- cumsum(porin)
   eig <- eigall[1:ndim]
   fe <- sqrt(solve(diag(abs(eig))))
   resin <- cbind(Re(round(eigall, dp + 5)), Re(round(100 * porin, dp + 5)),
      Re(round(100 * poracuin, dp + 5)))
   dimnames(resin) <- list(paste("s=", 1:length(eigall), sep = ""),
      c("values", "percentage", "cumulated"))
   ### eigen vectors
   temp1 <- eigen(A)$vectors[, 1:ndim]
   temp2 <- abs(diag(t(temp1) %*% temp1))
   if(I >= J) {
      U <- t(t(temp1)/sqrt(temp2))
      V <- X %*% U %*% fe   }
   else {
     V <- t(t(temp1)/sqrt(temp2))
     U <- t(X) %*% V %*% fe  }   
   ########## ########## ########## ########## ########## ########## ###########
   ####################### rows  ###############################################
   ### Projections
   Fs <- sweep(X %*% U, 1, sqrt(fi), FUN = "/")
   Fs <- replace(Fs, is.na(Fs), 0)
   dimnames(Fs) <- list(namei, paste(legendF, 1:ndim, sep = ""))
   ####################### other results for rows (if oar = 1) #################
   if(oar != 0) {
      ### contributions of the points to the dimensions
      ctai <- t(t((Fs * Fs) * fi)/eig)
      ### distances
      temp <- sweep(Fd, 1, fi, FUN = "/")
      temp <- sweep(temp, 2, fj, FUN = "-")
      temp <- temp^2
      temp <- sweep(temp, 2, fj, FUN = "/")
      temp <- replace(temp, is.na(temp), 0)
      d2i <- apply(temp, 1, sum)
      ### squared correlations
      ctri <- (Fs * Fs)/d2i
      ctri <- replace(ctri, is.na(ctri), 0)
      ### results of rows
      resi <- cbind(Re(round(100 * fi, dp)), Re(round(d2i, dp)),
         Re(round(Fs, dp)), Re(round(100 * ctai, dp)),
         Re(round(ctri, dp)))
      dimnames(resi) <- list(namei, 
         c(legend100fi, "d2", paste(legendF, 1:ndim, sep = ""),
         paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
      }
   if(oar == 0) {
      resi = d2i <- "option oar = 0"   }
   	### end rows
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ####################### columns #############################################
   if(oac != 0) {
      ### Projections
      Gs <- sweep(t(X) %*% V, 1, sqrt(fj), FUN= "/")
      Gs <- replace(Gs, is.na(Gs), 0)
      dimnames(Gs) <- list(namej, paste(legendG, 1:ndim, sep = ""))
      ### contributions of the points to the dimensions
      ctaj <- t(t((Gs * Gs) * fj)/eig)
      ### distances
      temp <- sweep(Fd, 2, fj, FUN = "/")
      temp <- sweep(temp, 1, fi, FUN = "-")
      temp <- temp^2
      temp <- sweep(temp, 1, fi, FUN = "/")
      temp <- replace(temp, is.na(temp), 0)
      d2j <- apply(temp, 2, sum)
      ### squared correlations
      ctrj <- (Gs * Gs)/d2j
      ctrj <- replace(ctrj, is.na(ctrj), 0)
      ### results of  columns
      resj <- cbind(Re(round(100 * fj, dp)), Re(round(d2j, dp)),
         Re(round(Gs, dp)), Re(round(100 * ctaj, dp)),
         Re(round(ctrj, dp)))
      dimnames(resj) <- list(namej, 
         c(legend100fj, "d2", paste(legendG, 1:ndim, sep = ""),
         paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
       }
   if (oac == 0) {
	   resj = Gs = d2j <- "Option oac = 0"   }
   ### end columns
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ### supplementary rows
   if(Isr == 0) {
      Fssr = nameisr = fisr = d2isr = Xsr <- 0
      resisr <- legendsr   }
   if(Isr != 0) {
      Fdsr <- as.matrix(datafasr/sum(datafaA))
      fisr <- apply(Fdsr, 1, sum)
      Isr <- length(fisr)
      ### matrix supplementary rows
      Xsr <- matrix(0, Isr, J)
      for(i in 1:Isr) {
         for(j in 1:J) {
            Xsr[i, j] <- (sqrt(fisr[i]) * (Fdsr[i, j] / fisr[i] -
               fj[j])) / sqrt(fj[j])   }   }
      Xsr[!is.finite(Xsr)] <- 0
      dimnames(Xsr) <- list(nameisr, namej)
      ### projections supplementary rows
      Fssr <- sweep(Xsr %*% U, 1, sqrt(fisr), FUN = "/")
      Fssr <- replace(Fssr, is.na(Fssr), 0)
      dimnames(Fssr) <- list(nameisr, paste(legendF, 1:ndim, sep = ""))
      ### distances and squared correlations supplementary rows
      temp <- sweep(Fdsr, 1, fisr, FUN = "/")
      temp <- sweep(temp, 2, fj, FUN = "-")
      temp <- temp^2
      temp <- sweep(temp, 2, fj, FUN = "/")
      temp <- replace(temp, is.na(temp), 0)
      d2isr <- apply(temp, 1, sum)
      ctrisr <- Fssr^2 / d2isr
      ctrisr <- replace(ctrisr, is.na(ctrisr), 0)
      dimnames(ctrisr) <- list(nameisr, NULL)
      ## results supplementary rows
      resisr <- Re(round(cbind(100 * fisr, as.matrix(d2isr), Fssr, ctrisr), dp))
      dimnames(resisr) <- list(nameisr, c(legend100fi, "d2",
         paste(legendF, 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))   }
   ### end supplementary rows
   ########## ########## ########## ########## ########## ########## ###########
   #### supplementary columns
   if(Jsc == 0) {
      Gssc = namejsc = fjsc = d2jsc <- 0
      resjsc <- legendsc   }
   if(Jsc != 0) {
      Fdsc <- datafasc/sum(datafaA)
      fjsc <- apply(Fdsc, 2, sum)
      Jsc <- length(fjsc)
      ### Matrix supplementary columns
      Xsc <- matrix(0, I, Jsc)
      for(i in 1:I) {
         for(j in 1:Jsc) {
            Xsc[i, j] <- (Fdsc[i, j] / fjsc[j]) - fi[i]   }   }
      ### projections supplementary columns
      Gssc <- sweep(t(Xsc) %*% V, 1, sqrt(fjsc), FUN= "/")
      Gssc <- replace(Gssc, is.na(Gssc), 0)
      dimnames(Gssc) <- list(namejsc, paste(legendG, 1:ndim, sep = ""))
      ### distances and squared correlations supplementary columns
      temp <- sweep(Fdsc, 2, fjsc, FUN = "/")
      temp <- sweep(temp, 1, fi , FUN = "-")
      temp <- temp^2
      temp <- sweep(temp, 1, fi, FUN = "/")
      temp <- replace(temp, is.na(temp), 0)
      d2jsc <- apply(temp, 2, sum)
      ctrjsc <- Gssc^2 / d2jsc
      ctrjsc <- replace(ctrjsc, is.na(ctrjsc), 0)
      dimnames(ctrjsc) <- list(namejsc, NULL)
      ## results supplementary columns
      resjsc <- Re(round(cbind(100 * fjsc, as.matrix(d2jsc), Gssc, ctrjsc), dp))
      dimnames(resjsc) <- list(namejsc, c(legend100fj, "d2",
         paste(legendG, 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))   }
   ### end supplementary columns
   ########## ########## ########## ########## ########## ########## ###########
   CorrAn.output <- list(totalin = totalin, eig = eig, resin = resin, resi = resi,
      resj = resj, resisr = resisr, resjsc = resjsc,
      X = X, totalk = sum(datafaA),
      I = I, namei = namei, fi = fi, Fs = Fs, d2i = d2i,
      J = J, namej = namej, fj = fj, Gs = Gs, d2j = d2j,
      Isr = Isr, nameisr = nameisr, fisr = fisr, Fssr = Fssr, d2isr = d2isr,
      Xsr = Xsr,
      Jsc = Jsc, namejsc = namejsc, fjsc = fjsc, Gssc = Gssc, d2jsc = d2jsc)
   class(CorrAn.output) <- "CorrAn"
   return(CorrAn.output)
}
### End CorrAn: Correspondence Analysis
################################################################################




################################################################################
#####               SimAn function: Simultaneous Analysis                  #####
################################################################################
SimAn <- function(data, G, acg, weight = 2, nameg = NA, sr = NA, sc = NA,
   nd = 2, dp = 2, oar = 1, oac = 1, multiple = 0, arg)
{
# -------------------------------------------------------------------
#  Simultaneous analysis is a factorial method developed for the
#  joint treatment of a set of several data tables, especially
#  frequency tables whose row margins are different, for example when
#  the tables are from different samples or different time points,
#  without modifying the internal structure of each table.
#  In the data tables  rows must refer to the same entities, but
#  columns may be different.
#
# Arguments
#   data     Data set
#   G        Number of tables to be jointly analyzed
#   acg      List of number of the active columns for each table
#   weight   Weighting on each table
#   nameg    Prefix for identifying partial rows and tables
#   sr       Indices of supplementary rows
#   sc       Indices of supplementary columns
#   nd       Number of dimensions in results
#   dp       Number of digits in results
#   oar      Output for active rows
#   oac      Output for active columns
#   multiple     ZZZ
#
# See also
#   summary.SimAn
#   plot.SimAn
#
# -------------------------------------------------------------------
# Reference
#      Zarraga, A. and Goitisolo, B.
#         Simultaneous Analysis is S-PLUS. The SimultAn Package.
#         submitted for publication, 2010.
# -------------------------------------------------------------------
cat("\n \n")
   if(multiple != 0 && multiple != 1 ) {
   	cat("\n Only values 0 and 1 for option multiple \n")
   	stop()   }
   ### some legends
   legendpi <- "pi"
   legend100fig <- "100fig"
   legendrows <- "overall rows"
   legendF<- "F"
   legend100fjg <- "100fjg"
   legendG <- "G" 
   legendsr <- "No supplementary rows"
   legendsc <- "No supplementary columns"
   ### multiple
   if(multiple == 1 ) {
      temp <- sr   
      sr <- sc
      sc <- temp
      temp <- oar
      oar <- oac 
      oac <- temp 
      data <- t(data)
      acg <- arg
      legendpi <- "pj"
      legend100fig <- "100fjg"
      legendrows <- "overall columns"
      legendF<- "G"
      legend100fjg <- "100fig"
      legendG <- "F" 
      legendsr <- "No supplementary columns"
      legendsc <- "No supplementary rows" 
      }
   ### input
   datafa <- as.matrix(data)
   ### Isr
   if(is.na(sr[1])) {
      Isr <- 0
      I <- nrow(datafa)
      nameic <- labels(datafa)[[1]]   }
   else {
      Isr <- length(sr)
      I <- nrow(datafa) - Isr
      nameic <- labels(datafa)[[1]][ - sr]
      nameicsr <- labels(datafa)[[1]][sr]   }
   ### names of groups
   if(is.na(nameg[1])) {
      nameg <- paste(legendG, 1:G, sep = "")   }
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ### separate CA
   ### oar = 1 Fs necessary for rCASA
   CAres <- vector("list", G)
   names(CAres) <- paste("Table_", 1:G, sep = "")
   for(g in 1:G) {
      Jg <- length(acg[g][[1]])
      if(!is.na(sc[1])) {
         tablag <- datafa[, c(acg[g][[1]], sc)]
         scg <- c((1 + Jg):(ncol(tablag)))   }
      else {
         tablag <- datafa[, c(acg[g][[1]])]
         scg <- NA   }
      dimnames(tablag)[[1]] <- paste(nameg[g], labels(tablag)[[1]], sep = "")
      CAres[[g]] <- CorrAn(tablag, sr = sr, sc = scg, nd = nd, dp = dp,
         oar = oar, oac = oac, multiple = multiple)
#      names(CAres[[g]]) <- c("totalin", "eig", "resin",
#         "resi", "resj", "resisr", "resjsc", "X", "totalk",
#         "I", "namei", "fi", "Fs", "d2i",
#         "J", "namej", "fj", "Gs", "d2j",
#         "Isr", "nameisr", "fisr", "Fssr", "d2isr", "Xsr",
#         "Jsc", "namejsc", "fjsc", "Gssc", "d2jsc")
   }
   ### end separate CA
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ### start simultaneous analysis
   J <- 0
   for(g in G:1) {
      J <- CAres[[g]]$J + J   }
   ndim <- min(I, J, nd)
   ### number of the columns belonging to each group (for allFs, allGs and ctag)
   cJ <- c()
   for(g in G:1) {
      cJ <- c(CAres[[g]]$J, cJ)   }
   ndimJgm1 <- c(0)
   for(g in G:1) {
      ndimJgm1[g + 1] <- sum(cJ[1:g])   }
   ### max Jg (for allGs)
   maxJg <- 0
   for(g in G:1) {
      maxJg <- max(maxJg, CAres[[g]]$J)   }
   ### names of rows (nameig) and columns (namej)
   nameig <- c()
   for(g in G:1) {
      nameig <- c(CAres[[g]]$namei, nameig)   }
   namej <- c()
   for(g in G:1) {
      namej <- c(CAres[[g]]$namej, namej)   }
   ### weighting of tables: 0=1/1st eig separate; 1=1; 2=1/totalin
   alphag <- c()
   if(weight == 1) {
      alphag <- rep(1, G)   }
   if(weight == 2) {
      for(g in G:1) {
         alphag <- c(1/CAres[[g]]$eig[[1]], alphag)   }   }
   if(weight == 3) {
      for(g in G:1) {
         alphag <- c(1/CAres[[g]]$totalin, alphag)   }   }
   ### weights of overall rows (pic); partial (fiG) and columns (fj)
   rfi <- array(0, c(I, G))
   for(g in G:1) {
      rfi[, g] <- sqrt(CAres[[g]]$fi)   }
   srfi <- apply(rfi, 1, sum)
   pic <- srfi^2
   fiG <- c()
   for(g in G:1) {
      fiG <- c(CAres[[g]]$fi, fiG)   }
   fj <- c()
   for(g in G:1) {
      fj <- c(CAres[[g]]$fj, fj)   }
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ### matrix to diagonalize
   X <- c()
   for(g in G:1) {
      X <- cbind(sqrt(alphag[g]) * CAres[[g]]$X, X)   }
   if(I >= J) {
      A <- t(X) %*% X   }
   else {
      A <- (X) %*% t(X)   }
   ### eigen values, inertias and percentages
   eigall <- eigen(A)$values
   totalin <- sum(eigall)
   porin <- eigall/totalin
   poracuin <- cumsum(porin)
   eig <- eigall[1:ndim]
   fe <- sqrt(solve(diag(abs(eig))))
   resin <- cbind(Re(round(eigall, dp + 5)), Re(round(100 * porin, dp + 5)),
      Re(round(100 * poracuin, dp + 5)))
   dimnames(resin) <- list(paste("s=", 1:length(eigall), sep = ""),
      c("values", "percentage", "cumulated"))
   ### eigen vectors
   temp1 <- eigen(A)$vectors[, 1:ndim]
   temp2 <- abs(diag(t(temp1) %*% temp1))
   if(I >= J) {
      U <- t(t(temp1)/sqrt(temp2))
      V <- X %*% U %*% fe   }
   else {
      V <- t(t(temp1)/sqrt(temp2))
      U <- t(X) %*% V %*% fe   }
   ### end diagonalization
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ####################### rows  ###############################################
   ### projections of overall rows
   Fsic <- sweep(X %*% U, 1, sqrt(pic), FUN = "/")
   Fsic <- replace(Fsic, is.na(Fsic), 0)
   dimnames(Fsic) <- list(nameic, paste(legendF, 1:ndim, sep = ""))
   ### projections of partial rows (Fsig) and allFs for graphs
   allFs <- array(0, c(I, ndim, G + 1))
   dimnames(allFs) <- list(nameic, paste(legendF, 1:ndim, sep = ""),
      c(paste("Table ", nameg, sep = ""), legendrows))
   for(g in 1:G) {
      temp <- X[, (1 + ndimJgm1[g]):ndimJgm1[g + 1]] %*%
         U[(1 + ndimJgm1[g]):ndimJgm1[g + 1], ]
      allFs[, , g] <- sweep(temp, 1, sqrt(CAres[[g]]$fi), FUN= "/")   }
   allFs[, , G + 1] <- Fsic
   Fsig <- c()
   for(g in G:1) {
      Fsig <- rbind(allFs[, , g], Fsig)   }
   dimnames(Fsig) <- list(nameig, paste(legendF, 1:ndim, sep = ""))
   ####################### other results for rows (if oar = 1) #################
   if (oar != 0) {
      ### contributions of the points to the dimensions
      temp <- sweep(Fsic^2, 1, pic, FUN = "*")
      ctaic <- sweep(temp, 2, eig, FUN = "/")
      ### distances
      d2ig <- c()
      for(g in G:1) {
         d2ig <- c((alphag[g]) * CAres[[g]]$d2i, d2ig)   }
      d2ic <- 0
      for(g in G:1) {
         d2ic <- c(((alphag[g]) * (CAres[[g]]$fi/pic) * CAres[[g]]$d2i) +
            d2ic)   }
      d2ic <- replace(d2ic, is.na(d2ic), 0)
      ### squared correlations
      ctric <- (Fsic^2)/d2ic
      ctric <- replace(ctric, is.na(ctric), 0)
      ctrig <- (Fsig^2)/d2ig
      ctrig[!is.finite(ctrig)] <- 0
      ### results
      resic <- cbind(Re(round(pic, dp)), Re(round(d2ic, dp)),
         Re(round(Fsic, dp)), Re(round(100 * ctaic, dp)),
         Re(round(ctric, dp)))
      dimnames(resic) <- list(nameic, c(legendpi, "d2", paste(legendF, 1:ndim, sep = ""),
         paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
      resig <- cbind(Re(round(100 * fiG, dp)), Re(round(d2ig, dp)),
         Re(round(Fsig, dp)), Re(round(ctrig, dp)))
      dimnames(resig) <- list(nameig, c(legend100fig, "d2",
         paste(legendF, 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
      }
   ### end rows
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ####################### columns #############################################
   ### projections of columns
   Gs <- sweep(t(X) %*% V, 1, sqrt(fj), FUN= "/")
   Gs <- replace(Gs, is.na(Gs), 0)
   dimnames(Gs) <- list(namej, paste(legendG, 1:ndim, sep = ""))
   ### contributions of the points to the dimensions
   temp <- sweep(Gs^2, 1, fj, FUN = "*")
   ctaj <- sweep(temp, 2, eig, FUN = "/")
   ### contributions of the groups to the dimensions
   ctag <- c()
   for(g in G:1) {
      ctajg <- ctaj[(1 + ndimJgm1[g]):ndimJgm1[g + 1], ]
      ctag <- rbind(apply(ctajg, 2, sum), ctag)   }
   ctag <- Re(round(ctag, dp))
   dimnames(ctag) <- list(paste("Table", nameg[1:G]),
      paste("Axis ", 1:ndim, sep = ""))
   ####################### other results for columns (if oac = 1) ##############
   if (oac != 0) {
      ### allGs (Projections of columns) for graphs
      allGs <- array(0, c(maxJg, ndim, G))
      for(g in G:1) {
         allGs[(1:CAres[[g]]$J), , g] <- Gs[(1 + ndimJgm1[g]):ndimJgm1[g + 1], ]
         }
         ### without names. Not possible different names for each table in array
      ### distances
      d2j <- c()
      for(g in G:1) {
         d2j <- c((alphag[g]) * matrix(CAres[[g]]$d2j), d2j)   }
      ### squared correlations
      ctrj <- (Gs^2)/d2j
      ctrj <- replace(ctrj, is.na(ctrj), 0)
      ### results 
      resj <- cbind(Re(round(100 * fj, dp)), Re(round(d2j, dp)),
         Re(round(Gs, dp)), Re(round(100 * ctaj, dp)),
         Re(round(ctrj, dp)))
      dimnames(resj) <- list(namej,
         c(legend100fjg, "d2", paste(legendG, 1:ndim, sep = ""),
         paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
      }
   ### end columns
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ####################### groups ##############################################
   ### projections of the groups
   Fsg <- t(eig * t(ctag))
   Fsg <- Re(round(Fsg, dp))
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ####################### relation among factors of CA, SA ####################
   ### relation between the overall rows and the partial rows
   riig <- matrix(0, G, ndim)
   for(is in 1:ndim) {
      temp <- as.matrix(allFs[, is, G + 1])
      temp2 <- sweep(temp, 1, srfi, FUN = "*")
      for(g in 1:G) {
         temp <- as.matrix(allFs[, is, g])
         temp1 <- sweep(temp, 1, rfi[, g], FUN = "*")
         varsg <- t(temp1) %*% temp1
         riig[g,is] <- t((1 / sqrt(varsg[1,1])) * temp1) %*%
            ((1/sqrt(eig[is])) * temp2)   }   }
   riig <- Re(round(riig, dp))
   dimnames(riig) <- list(paste("Table", nameg[1:G]),
      paste("Axis ", 1:ndim, sep = ""))
   ### relation between CA axes and SA axes
   rCASA <- array(0, c(ndim, ndim, G))
   for(is in 1:ndim) {
      temp <- as.matrix(allFs[, is, G + 1])
      temp2 <- sweep(temp, 1, srfi, FUN = "*")
      for(g in 1:G) {
         ndimcs <- min((I-1), ((CAres[[g]]$J)-1), ndim)
         for(sp in 1:ndimcs) {
            if(CAres[[g]]$eig[sp]>0 & eig[is]>0) {
               temp <- as.matrix(CAres[[g]]$Fs[, sp])
               temp1 <- sweep(temp, 1, rfi[, g], FUN = "*")
               rCASA[sp, is, g] <- t((1 / sqrt(CAres[[g]]$eig[sp])) * temp1) %*%
                  ((1/sqrt(eig[is])) * temp2)   }   }   }   }
   dimnames(rCASA) <- list(paste("CA Axis", 1:ndim), paste("SA Axis", 1:ndim),
      paste("Table", nameg[1:G]))
   rCASA <- Re(round(rCASA, dp))
   ### relation between separate CA axes
   rCACA <- array(0, c(ndim, ndim, G * (G - 1) / 2))
   namegh <- c(1:(G * (G - 1) / 2))
   gh <- 1
   for(g in 1:G) {
      ndimis <- min((I-1), ((CAres[[g]]$J)-1), ndim)
      for(h in 1:G) {
         ndimsp <- min((I-1), ((CAres[[h]]$J)-1), ndim)
         if (g < h) {
            for(is in 1:ndimis) {
               temp <- as.matrix(CAres[[g]]$Fs[, is])
               temp2 <- sweep(temp, 1, rfi[, g], FUN = "*")
               for(sp in 1:ndimsp) {
                  if(CAres[[g]]$eig[is]>0 & CAres[[h]]$eig[sp]>0) {
                     temp <- as.matrix(CAres[[h]]$Fs[, sp])
                     temp1 <- sweep(temp, 1, rfi[, h], FUN = "*")
                     rCACA[sp, is, gh] <-
                        t((1/sqrt(CAres[[h]]$eig[sp])) * temp1) %*%
                        ((1/sqrt(CAres[[g]]$eig[is])) * temp2)
                     namegh[gh] <- paste("Table", nameg[g], "Table", nameg[h])
                  }   }   }
            gh <- gh + 1
         }   }   }
   rCACA <- Re(round(rCACA, dp))
   dimnames(rCACA) <- list(paste("CA Axis", 1:ndim), paste("CA Axis", 1:ndim),
      namegh)
   ### end relation among factors of CA, SA
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ### supplementary rows
   if(is.na(sr[1])) {
      resicsr = resigsr = Fsicsr = Fsigsr = allFssr = nameicsr <-
         legendsr
      Isr = 0   }
   if(!is.na(sr[1])) {
      Isr <- length(sr)
      fisr <- c()
      for(g in G:1) {
         fisr <- c(CAres[[g]]$fisr, fisr)   }
      nameigsr <- c()
      ### names of supplementary rows
      for(g in G:1) {
         nameigsr <- c(CAres[[g]]$nameisr, nameigsr)   }
      rfisr <- array(0, c(Isr, G))
      for(g in G:1) {
         rfisr[, g] <- sqrt(CAres[[g]]$fisr)   }
      srfisr <- c(rep(0, Isr))
      for(g in G:1) {
         srfisr <- rfisr[, g] + srfisr   }
      picsr <- srfisr^2
      ### Xsr
      Xsr <- c()
      for(g in G:1) {
         Xsr <- cbind(sqrt(alphag[g]) * CAres[[g]]$Xsr, Xsr)   }
      ### projections
      Fsicsr <- sweep(Xsr %*% U, 1, sqrt(picsr), FUN = "/")
      Fsicsr <- replace(Fsicsr, is.na(Fsicsr), 0)
      dimnames(Fsicsr) <- list(nameicsr, paste(legendF, 1:ndim, sep = ""))
      ### projections of partial rows (Fsigsr) and allFssr for graphs
      allFssr <- array(0, c(Isr, ndim, G + 1))
      dimnames(allFssr) <- list(nameicsr, paste(legendF, 1:ndim, sep = ""),
         c(paste("Table ", nameg, sep = ""), legendrows))
      for(g in 1:G) {
         temp <- Xsr[, (1 + ndimJgm1[g]):ndimJgm1[g + 1]] %*%
            U[(1 + ndimJgm1[g]):ndimJgm1[g + 1], ]
         allFssr[, , g] <- sweep(temp, 1, sqrt(CAres[[g]]$fisr), FUN= "/")   }
      allFssr[, , G + 1] <- Fsicsr
      Fsigsr <- c()
      for(g in G:1) {
         Fsigsr <- rbind(allFssr[, , g], Fsigsr)   }
      dimnames(Fsigsr) <- list(nameigsr, paste(legendF, 1:ndim, sep = ""))
      ### contributions of the points to the dimensions
      temp <- sweep(Fsicsr^2, 1, picsr, FUN = "*")
      ctaicsr <- sweep(temp, 2, eig, FUN = "/")
      ### distances
      d2igsr <- c()
      for(g in G:1) {
         d2igsr <- c((alphag[g]) * CAres[[g]]$d2isr, d2igsr)   }
      d2icsr <- 0
      for(g in G:1) {
         d2icsr <- c(((alphag[g]) * (CAres[[g]]$fisr/picsr) * 
            CAres[[g]]$d2isr) + d2icsr)   }
      d2icsr <- replace(d2icsr, is.na(d2icsr), 0)
      ### squared correlations
      ctricsr <- (Fsicsr^2)/d2icsr
      ctricsr <- replace(ctricsr, is.na(ctricsr), 0)
      ctrigsr <- (Fsigsr^2)/d2igsr
      ctrigsr[!is.finite(ctrigsr)] <- 0
      ### results supplementary rows
      resicsr <- cbind(Re(round(picsr, dp)), Re(round(d2icsr, dp)),
         Re(round(Fsicsr, dp)), Re(round(ctricsr, dp)))
      dimnames(resicsr) <- list(nameicsr, c(legendpi, "d2",
         paste(legendF, 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
      resigsr <- cbind(Re(round(100 * fisr, dp)), Re(round(d2igsr, dp)),
         Re(round(Fsigsr, dp)), Re(round(ctrigsr, dp)))
      dimnames(resigsr) <- list(nameigsr, c(legend100fig, "d2",
         paste(legendF, 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))   }
   ### end supplementary rows
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   #### supplementary columns
   if(is.na(sc[1])) {
      allresjsc = allGssc = Gssc = namejsc <- legendsc
      Jsc = 0   }
   if(!is.na(sc[1])) {
      ### allresjsc, allGssc
      Jsc <- length(sc)
      datafasc <- matrix(0, I, Jsc)
      if(!is.na(sr[1])) {
         datafasc <- as.matrix(datafa[ - sr, sc])   }
      else {
         datafasc <- as.matrix(datafa[, sc])   }
#      namejscorig <- labels(datafasc)[[2]]
      namejscorig <- labels(datafa)[[2]][sc]
      allGssc <- array(0, c(Jsc, ndim, G))
      allresjsc <- 0
      Gssc <- array(0, c(1, ndim))
      for(g in G:1) {
         namejsc <- paste(nameg[g], namejscorig, sep = "")
         Fdsc <- datafasc/CAres[[g]]$totalk
         fjsc <- apply(Fdsc, 2, sum)
         fisc <- CAres[[g]]$fi
         temp <- sweep(Fdsc, 2, fjsc, FUN = "/")
         temp <- replace(temp, is.na(temp), 0)
         temp <- sweep(temp, 1, fisc, FUN = "-")
         temp <- sweep(temp, 1, sqrt(fisc), FUN = "/")
         temp <- replace(temp, is.na(temp), 0)
         Xsc <- sqrt(alphag[g]) * temp
         Gsgsc <- t(Xsc) %*% V
         dimnames(Gsgsc) <- list(namejsc, paste(legendG, 1:ndim, sep = ""))
         temp <- sweep(Fdsc, 2, fjsc, FUN = "/")
         temp <- sweep(temp, 1, fisc , FUN = "-")
         temp <- temp^2
         temp <- sweep(temp, 1, fisc, FUN = "/")
         temp <- replace(temp, is.na(temp), 0)
         d2jsc <- apply(temp, 2, sum)
         ctrjsc <- Gsgsc^2 / d2jsc
         ctrjsc[!is.finite(ctrjsc)] <- 0
         dimnames(ctrjsc) <- list(namejsc, NULL)
         ### results supplementary columns
         resjsc <- Re(round(cbind(100 * fjsc, as.matrix(d2jsc), Gsgsc, ctrjsc),
            dp))
         dimnames(resjsc) <- list(namejsc, c(legend100fjg, "d2",
            paste(legendG, 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
         ### end results supplementary columns
         allresjsc <- rbind(resjsc, allresjsc)
         allGssc[(1:Jsc), , g] <- Gsgsc
         Gssc <- rbind(Gsgsc, Gssc)
         dimnames(allGssc) <- list(CAres[[g]]$namejsc,
            paste(legendG, 1:ndim, sep = ""),
            c(paste("Table ", nameg, sep = "")))   }
      allresjsc <- allresjsc[1:nrow(allresjsc) - 1, ]
      Gssc <- Gssc[1:nrow(Gssc) - 1, ]
   }
   ### end supplementary columns
   ########## ########## ########## ########## ########## ########## ###########
   ########## ########## ########## ########## ########## ########## ###########
   ### oar and oac of elements needed for other things (ctag, riig, ...)
   if(oar == 0) {
      resic = resig = Fsic = Fsig = allFs <- "option oar = 0"   
#      for (g in 1:G) {
#      CAres[[g]]$resi = CAres[[g]]$ctai = 
#         CAres[[g]]$ ctri <- "option oar = 0"   }
      }
   if(oac == 0) {
      resj = Gs = allGs <- "option oac = 0"   }
   ########## ########## ########## ########## ########## ########## ###########
   if(multiple == 1) {
	     for (g in 1:G){
	      CAres[[g]]$X <-t(CAres[[g]]$X)
         names(CAres[[g]]) <- c("totalin", "eig", "resin",
         "resj", "resi", "resjsc", "resisr", "X", "totalk",
         "J", "namej", "fj", "Gs", "d2j",
         "I", "namei", "fi", "Fs", "d2i",
         "Jsc", "namejsc", "fjsc", "Gssc", "d2jsc", "Xsc",
         "Isr", "nameisr", "fisr", "Fssr", "d2isr")   
           }   
   SimAn.output <- list(totalin = totalin, 
      resin = resin, resi = resj, resj = resic, resjg = resig, 
      Fsg = Fsg, ctrg = ctag, rjjg = riig, rCACA = rCACA, rCASA = rCASA,
      Fs = Gs, Gs = Fsic, Gsjg = Fsig, allFs = allGs,  allGs = allFs,      
      J = I, maxIg = maxJg, G = G, namej = nameic, nameg = nameg,
      resisr = allresjsc, resjsc = resicsr, 
      resjgsc = resigsr, 
      Fssr = Gssc, Gssc = Fsicsr,
      Gsjgsc = Fsigsr, 
      allFssr = allGssc, allGssc = allFssr,
      Isr = Jsc, Jsc = Isr, 
      nameisr = namejsc, namejsc = nameicsr, 
      CAres = CAres,
      multiple = multiple)   }
   ########## ########## ########## ########## ########## ########## ###########
   if(multiple == 0) {
      SimAn.output <- list(totalin = totalin, 
         resin = resin, resi = resic, resj = resj, resig = resig, 
         Fsg = Fsg, ctrg = ctag, riig = riig, rCACA = rCACA, rCASA = rCASA,
         Fs = Fsic, Gs = Gs, Fsig = Fsig, allFs = allFs, allGs = allGs,
         I = I, maxJg = maxJg, G = G, namei = nameic, nameg = nameg,
         #namejg = namejg,
         #d2g = d2g, resg = resg, Lgh = Lgh, d2gh = d2gh,
         resisr = resicsr, resjsc = allresjsc, 
         resigsr = resigsr, 
         Fssr = Fsicsr, Gssc = Gssc, 
         Fsigsr = Fsigsr, 
         allFssr = allFssr, allGssc = allGssc, 
         Isr = Isr, Jsc = Jsc,
         nameisr = nameicsr, namejsc = namejsc,
         CAres = CAres,
         multiple = multiple)   }
   class(SimAn.output) <- "SimAn"
   return(SimAn.output)
}
### end SimAn: Simultaneous Analysis
################################################################################






#RyS#
################################################################################
#####               Summary of CorrAn                                      #####
################################################################################
summary.CorrAn <- function(object, oar = 1, oac = 1, ...)
{
# -------------------------------------------------------------------
#  This function summarizes the results of CorrAn.
#
# Arguments
#   object   The output of the correspondence analysis (class "CorrAn")
#   oar      Output for active rows
#   oac      Output for active columns
#
# See also
#   CorrAn
#   plot.CorrAn
#
# -------------------------------------------------------------------
# Reference
#      Zarraga, A. and Goitisolo, B.
#        Simultaneous Analysis is S-PLUS. The SimultAn Package.
#        submitted for publication, 2010.
# -------------------------------------------------------------------
cat("\n \n")
   CAres <- object
   if (oar == 0) { CAres$resi <- "option oar = 0" } 
   if (oac == 0) { CAres$resj <- "option oac = 0" }
   CAR <- list("CA table" = "CA table",
               "Total inertia" = CAres$totalin,
               "Eigenvalues and percentages of inertia" = CAres$resin,
               "Output for rows" = CAres$resi,
               "Output for columns" = CAres$resj,
               "Output for supplementary rows" = CAres$resisr,
               "Output for supplementary columns" = CAres$resjsc)
   return(CAR = CAR)
}
### end summary.CorrAn: Summary of Correspondence Analysis
################################################################################


#RyS
################################################################################
#####               Summary of SimAn                                       #####
################################################################################
summary.SimAn <- function(object, oar = 1, oac = 1, ...)
{
# -------------------------------------------------------------------
#  This function summarizes the results of SimAn.
#
# Arguments
#   object   The output of the simultaneous analysis (class "SimAn")
#   oar      Output for active rows
#   oac      Output for active columns
#
# See also
#   SimAn
#   plot.SimAn
#
# -------------------------------------------------------------------
# Reference
#      Zarraga, A. and Goitisolo, B.
#         Simultaneous Analysis is S-PLUS. The SimultAn Package.
#         submitted for publication, 2010.
# -------------------------------------------------------------------
cat("\n \n")
   ASG <- object
   ### CA
   CAres <- ASG$CAres
   G <- ASG$G
   multiple <- ASG$multiple
   if (oar == 0) {
      for(g in 1:G) {CAres[[g]]$resi <- "option oar = 0"}   }
   if (oac == 0) {
      for(g in 1:G) {CAres[[g]]$resj <- "option oac = 0"}   }
   CAR <- vector("list", G)
   names(CAR) <- paste("CA", ASG$nameg)
   for(g in 1:G) {
         CAR[[g]] <- list(CAres[[g]]$totalin, CAres[[g]]$resin, CAres[[g]]$resi,
         CAres[[g]]$resj, CAres[[g]]$resisr, CAres[[g]]$resjsc)
      names(CAR[[g]]) <- c("Total inertia",
      "Eigenvalues and percentages of inertia", "Output for rows",
      "Output for columns", "Output for supplementary rows",
      "Output for supplementary columns")   }
   ### SA
   if (oar == 0) {
      ASG$resi <- "option oar = 0"
      ASG$resig <- "option oar = 0"   }
   if (oac == 0) {ASG$resj <- "option oac = 0"}
   if (multiple == 1) {
      SAR <- list("Total inertia" = ASG$totalin,
         "Eigenvalues and percentages of inertia" = ASG$resin,
         "Output for rows" = ASG$resi,
         "Output for columns" = ASG$resj,
         "Output for partial columns" = ASG$resjg,
         "Projections of tables" = ASG$Fsg,
         "Contributions of tables to SA" = ASG$ctrg,
         "Relation between overall and partial columns" = ASG$rjjg,
         "Relation between factors of separate CA" = ASG$rCACA,
         "Relation between factors of CA and SA" = ASG$rCASA,
         "Output for supplementary rows" = ASG$resisr,
         "Output for supplementary columns" = ASG$resjsc,
         "Output for supplementary partial columns" = ASG$resjgsc )
      }  
   if (multiple == 0) {
      SAR <- list("Total inertia" = ASG$totalin,
         "Eigenvalues and percentages of inertia" = ASG$resin,
         "Output for rows" = ASG$resi,
         "Output for columns" = ASG$resj,
         "Output for partial rows" = ASG$resig,
         "Projections of tables" = ASG$Fsg,
         "Contributions of tables to SA" = ASG$ctrg,
         "Relation between overall and partial rows" = ASG$riig,
         "Relation between factors of separate CA" = ASG$rCACA,
         "Relation between factors of CA and SA" = ASG$rCASA,
         "Output for supplementary rows" = ASG$resisr,
         "Output for supplementary columns" = ASG$resjsc,
         "Output for supplementary partial rows" = ASG$resigsr )
      }   
   ### output CA y SA
   return(list(CAR = CAR, SAR = SAR))
}
### end summary.SimAn: Summary of Simultaneous Analysis
################################################################################





################################################################################
#####               plot.CorrAn: Graphs for CorrAn                        ######
################################################################################
#R# plot.CorrAn <- function(x, s1 = 1, s2 = 2, screen = TRUE, ...)
#S# plot.CorrAn <- function(x, s1 = 1, s2 = 2, screen = FALSE, ...)
plot.CorrAn <- function(x, s1 = 1, s2 = 2, screen = TRUE, oar = 1, oac = 1, ...)
{
# -------------------------------------------------------------------
#  Graphical display of correspondence analysis results in two
#  dimensions.
#
# Arguments
#   x       The output of the correspondence analysis (class "CorrAn")
#   s1      Dimension to plot on horizontal axis
#   s2      Dimension to plot on vertical axis
#   screen  TRUE for R and FALSE for S-plus
#   oar     Output for active rows
#   oac     Output for active columns
#
# See also
#   CorrAn
#   summary.CorrAn
#
# -------------------------------------------------------------------
# Reference
#      Zarraga, A. and Goitisolo, B.
#         Simultaneous Analysis is S-PLUS. The SimultAn Package.
#         submitted for publication, 2010.
# -------------------------------------------------------------------
cat("\n \n")
   if(oar == 0 && oac == 0) {
      cat ("\n ***************************************************  \n") 
      cat ("\n WARNING \n")
      cat ("\n Options oar=0 and oac=0. No active elements to plot. \n") 
      cat ("\n ***************************************************  \n") 
        }
      
#R# opcion screen=TRUE
#S# quitar opcion screen=TRUE
   SAPplot <- function(x, y)
   {
         plot(x, y, type = "n",
#S#            axes = T,
            axes = T,
            xlab = paste("  "),
            ylab = paste("  "),
            col = 1,
#S#            cex = 0.5)
            cex = 0.5)
   }
   SAPtitle <- function(titulo, titsub, s1, s2, resin)
   {
      title(titulo,
#S#         axes = F,
         col = 1,
#S#          cex = 1.2)
         cex = 1.2)
      title("",
#S#         axes = F,
         xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)   "),
         ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)   "),
         adj = 1,
         col = 1,
#S#         cex = 1)
         cex = 1)
#R#         mtext(titsub, side = 3, line = 0.4, outer = FALSE)
#S#         mtext(titsub, side = 3, line = 1, outer = F)
      mtext(titsub, side = 3, line = 1, outer = F)
   }
   SAPplotpuntos <- function(puntos, gra, s1, s2, etiq,
      colg = 2,
      SAPpch)
      {
         points(c(puntos[gra, s1]), c(puntos[gra, s2]),
            pch = SAPpch,
            col = colg,
#S#            cex = 0.7)
            cex = 0.7)
         text(c(puntos[gra, s1]), c(puntos[gra, s2]),
            c(paste("  ", etiq[gra], sep = "")),
            col = colg,
            adj = 0,
#S#            cex = 0.7)
            cex = 0.7)
      }
   ### Graph options
#R#    par(cex = 1, cex.axis = .6, cex.lab = 1, cex.main = 1.2, cex.sub = 1)
#R#    colg <- c(rep(2:8, 10))
#R#    pchg <- c(rep(c(15, 17, 18, 19), 10))
#R#    pchgs <- c(rep(c(0, 2, 5, 1), 10))
#S#    cex = 1, cex.axis = .6, cex.lab = 1, cex.main = 1.2, cex.sub = 1
#RS    colg.F.CA  y colg.sr.CA    #R 1 #S 16
   colg <- 2                        # color groups
   pchg <- 15                       # pch groups
   pchgs <- 0                       # pch groups supplementary elements
   colg.F.CA <-  1                  # color rows
   colg.sr.CA <- 1                  # color supplementary rows
   pchg.F.CA <- 24                  # pch rows
   pchg.sr.CA <- 4                  # pch supplementary rows
   ### Input
   CAres <- x
   names(CAres) <- c("totalin", "eig", "resin", "resi", "resj",
      "resisr", "resjsc",
      "X", "totalk",
      "I", "namei", "fi", "Fs", "d2i",
      "J", "namej", "fj", "Gs", "d2j",
      "Isr", "nameisr", "fisr", "Fssr", "d2isr",
      "Xsr",
      "Jsc", "namejsc", "fjsc", "Gssc", "d2jsc")
   if (oar != 0 && class(CAres$Fs) == "character" ) {
      cat ("\n ***************************************************  \n") 
      cat ("\n WARNING \n")
      cat ("\n option oar = 1 in plot but oar = 0 in results of CorrAn \n" )   
      cat ("\n ***************************************************  \n")    }                                                                            
   if (oac != 0 && class(CAres$Gs) == "character" ) {
      cat ("\n ***************************************************  \n") 
      cat ("\n WARNING \n")
      cat("\n option oac = 1 in plot but oac = 0 in results of CorrAn \n" )
      cat ("\n ***************************************************  \n")    }
                                                                           
   titulinparc <- "CA"
   ### 
   ### dimensions
   I <- CAres$I
   Isr <- CAres$Isr
   Jsc <- CAres$Jsc
   J <- CAres$J
   ###
   nameic <- CAres$namei
   nameicsr <- CAres$nameisr
   #######################   CA active   #######################################
   if(oar != 0 || oac != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulinparc, sep = "")
      titsub <- "Active elements"
      puntosplot <- c() 
      if (oar != 0 && class(CAres$Fs) != "character" ) {
         namei <- CAres$namei
         puntosplot <- rbind(puntosplot, CAres$Fs)   }
      if (oac != 0 && class(CAres$Gs) != "character" ) {
         puntosplot <- rbind(puntosplot, CAres$Gs)   }
      if (length(puntosplot) != 0) {
         resin <- round(CAres$resin, 2)
         SAPplot(puntosplot[, s1], puntosplot[, s2])
         SAPtitle(titulo, titsub, s1, s2, resin)
         if (oar != 0 && class(CAres$Fs) != "character" ) {
            SAPplotpuntos(puntos = CAres$Fs, gra = c(1:I), s1, s2,
               etiq = CAres$namei, colg = 1, SAPpch = pchg.F.CA)   }
         if (oac != 0 && class(CAres$Gs) != "character" ) {
            SAPplotpuntos(puntos = CAres$Gs,
               gra = c(1:CAres$J), s1, s2, etiq = CAres$namej,
               SAPpch = pchg)   }   }
      }
   #######################   CA active + supplementary   #######################
   if(Isr != 0 || Jsc != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulinparc, sep = "")
      titsub <- "Active and supplementary elements"
      resin <- round(CAres$resin, 2)
      puntosplot <- c()
      if (oar != 0 && class(CAres$Fs) != "character" ) {
         namei <- CAres$namei
         puntosplot <- rbind(puntosplot, CAres$Fs)   }
      if (oac != 0 && class(CAres$Gs) != "character" ) {
         puntosplot <- rbind(puntosplot, CAres$Gs)   }
      if (Isr != 0) {
         puntosplot <- rbind(puntosplot, CAres$Fssr)   }
      if (Jsc != 0) {
         puntosplot <- rbind(puntosplot, CAres$Gssc)   }
      ### plot
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      if (oar != 0 && class(CAres$Fs) != "character" ) {
         SAPplotpuntos(puntos = CAres$Fs, gra = c(1:I), s1, s2,
            etiq = CAres$namei, colg = colg.F.CA,
            SAPpch = pchg.F.CA)   }
      if (oac != 0 && class(CAres$Gs) != "character" ) {
         SAPplotpuntos(puntos = CAres$Gs,
            gra = c(1:CAres$J), s1, s2, etiq = CAres$namej,
            SAPpch = pchg)   }
      if (Isr != 0) {
         puntos <- CAres$Fssr
         SAPplotpuntos(puntos = puntos, gra = c(1:Isr), s1, s2,
            etiq = CAres$nameisr, colg = colg.sr.CA,
            SAPpch = pchg.sr.CA)   }
      if (Jsc != 0) {
         puntos <- CAres$Gssc
         SAPplotpuntos(puntos = puntos,
            gra = c(1:CAres$Jsc), s1, s2,
            etiq = CAres$namejsc, SAPpch = pchgs)   }   }   
   else {   }
}
### end plot.CorrAn
################################################################################



################################################################################
#####               plot.SimAn: Graphs for SimAn                          ######
################################################################################
#R# plot.SimAn <- function(x, s1 = 1, s2 = 2, screen = TRUE, ...)
#S# plot.SimAn <- function(x, s1 = 1, s2 = 2, screen = FALSE, ...)
plot.SimAn <- function(x, s1 = 1, s2 = 2, screen = TRUE, oar = 1, oac = 1, ... )
{
# -------------------------------------------------------------------
#  Graphical representation of simultaneous analysis.
#
# Arguments
#   x       The output of the simultaneous analysis (class "SimAn")
#   s1      Dimension to plot on horizontal axis
#   s2      Dimension to plot on vertical axis
#   screen  TRUE for R and FALSE for S-plus
#   oar     Output for active rows
#   oac     Output for active columns
#
# See also
#   SimAn
#   summary.SimAn
#
# -------------------------------------------------------------------
# Reference
#      Zarraga, A. and Goitisolo, B.
#        Simultaneous Analysis is S-PLUS. The SimultAn Package.
#        submitted for publication, 2010.
# -------------------------------------------------------------------
cat("\n \n")
   if(oar == 0 && oac == 0) {
      cat ("\n ***************************************************  \n") 
      cat ("\n WARNING \n")
      cat ("\n Options oar=0 and oac=0. No active elements to plot. \n") 
      cat ("\n ***************************************************  \n")   }
  SAPplot <- function(x, y)
   {
         plot(x, y, type = "n",
            xlab = paste("  "),
            ylab = paste("  "),
            col = 1,
            cex = 0.5)
   }
   SAPtitle <- function(titulo, titsub, s1, s2, resin)
   {
      title(titulo,
         col = 1,
         cex = 1.2)
      title("",
         xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)   "),
         ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)   "),
         adj = 1,
         col = 1,
         cex = 1)
#R#       mtext(titsub, side = 3, line = 0.4, outer = FALSE)
#S#       mtext(titsub, side = 3, line = 1, outer = F)
      mtext(titsub, side = 3, line = 1, outer = F)
   }
   SAPplotpuntos <- function(g, puntos, gra, s1, s2, etiq,
      colg = c(rep(2:8, 10)),
      SAPpch)
      {
         points(c(puntos[gra, s1]), c(puntos[gra, s2]),
            pch = SAPpch[g],
            col = colg[g],
            cex = 0.7)
         text(c(puntos[gra, s1]), c(puntos[gra, s2]),
            c(paste("  ", etiq[gra], sep = "")),
            col = colg[g],
            adj = 0,
            cex = 0.7)
      }
   ### Input
   ASG <- x
   G <- ASG$G
   multiple <- ASG$multiple
   ########## ########## ########## ########## ########## ########## ###########
   if(multiple == 1) {
	   for (g in 1:G){
	      ASG$CAres[[g]]$X <-t(ASG$CAres[[g]]$X)
         names(ASG$CAres[[g]]) <- c("totalin", "eig", "resin",
            "resi", "resj", "resisr", "resjsc", "X", "totalk",
            "I", "namei", "fi", "Fs", "d2i",
            "J", "namej", "fj", "Gs", "d2j",
            "Isr", "nameisr", "fisr", "Fssr", "d2isr", "Xsr",
            "Jsc", "namejsc", "fjsc", "Gssc", "d2jsc")   }   
         names(ASG) <- c("totalin", "resin", "resj", "resic", "resig", 
            "Fsg", "ctag", "riig", "rCACA", "rCASA",
            "Gs", "Fs", "Fsig", "allGs",  "allFs",      
            "I", "maxJg", "G", "namei", "nameg",
            "allresjsc", "resicsr", 
            "resigsr", 
            "Gssc", "Fssr",
            "Fsigsr", 
            "allGssc", "allFssr",
            "Jsc", "Isr", 
            "namejsc", "nameisr", 
            "CAres",
            "multiple")
      }
   ########## ########## ########## ########## ########## ########## ###########
   if(multiple == 0) {
      if (oar == 1 && class(ASG$Fs) == "character" ) {
         cat("\n ***************************************************  \n") 
         cat("\n WARNING \n")
         cat("\n option oar = 1 in plot but oar = 0 in results of SimAn \n" ) 
         cat("\n ***************************************************  \n")   }
      if (oac == 1 && class(ASG$Gs) == "character" ) {
         cat("\n ***************************************************  \n") 
         cat("\n WARNING \n")
         cat("\n option oac = 1 in plot but oac = 0 in results of SimAn \n" )   
         cat("\n ***************************************************  \n")   }
      }
   if(multiple == 1) {
      if (oar == 1 && class(ASG$Fs) == "character" ) {
         cat("\n ***************************************************  \n") 
         cat("\n WARNING \n")
         cat("\n option oar = 1 in plot but oar = 0 in results of SimAn \n" ) 
         cat("\n ***************************************************  \n")   }
      if (oac == 1 && class(ASG$Gs) == "character" ) {
         cat("\n ***************************************************  \n") 
         cat("\n WARNING \n")
         cat("\n option oac = 1 in plot but oac = 0 in results of SimAn \n" )   
         cat("\n ***************************************************  \n")   }  
      }
   ### Graph options
#R#    par(cex = 1, cex.axis = .6, cex.lab = 1, cex.main = 1.2, cex.sub = 1)
#R#    colg <- c(rep(2:8, 10))
#R#    pchg <- c(rep(c(15, 17, 18, 19), 10))
#R#    pchgs <- c(rep(c(0, 2, 5, 1), 10))
#S#    cex = 1, cex.axis = .6, cex.lab = 1, cex.main = 1.2, cex.sub = 1
#S#    colg <- c(rep(2:5, 10))
#S#    pchg <- c(rep(c(15, 16, 17, 18, 19), 10))
#S#    pchgs <- c(rep(c(0, 1, 2, 5, 6), 10))
#RS    cuidado con cex de cuadrado, circulo y (Fsic+Fsig+Fsicsr+Fsigsr)
#RS   colg.F.CA  y colg.sr.CA    #R 1 #S 16
   colg <- c(rep(2:15, 10))                   # color groups
   pchg <- c(rep(c(15, 16, 17, 18, 19), 10))  # pch groups
   pchgs <- c(rep(c(0, 1, 2, 5, 6), 10))      # pch groups supplem. elements 
   colg.F.CA <-  rep(1, G + 1)                # color columns CA
   colg.sr.CA <- rep (1, G + 1)               # color supplementary columns CA
   pchg.F.CA <- rep (24, G + 1)               # pch Columnas CA
   pchg.sr.CA <- rep (4, G + 1)               # pch supplementary columns 
   ###
   titulin <- "SA"
   titulinparc <- "CA"
   ### CA
   CAres <- ASG$CAres
   J <- 0
   for(g in G:1) {
      J <- CAres[[g]]$J + J}
   ### Dimensions
   I <- ASG$I
   Jg <- ASG$maxJg
   ndim <- ncol(ASG$Fsg)
   Isr <- ASG$Isr
   Jsc <- ASG$Jsc
   ### Names
   ###
   nameic <- ASG$namei
   nameicsr <- ASG$nameisr
   namegroup <- ASG$nameg
   ### Data from separate CA
   if (oar != 0) {
      Fsgroups <- array(0, c(I, ndim, G))
      for(g in G:1) {
         Fsgroups[, , g] <- CAres[[g]]$Fs   }   }
   if (oac != 0) {
      Gsgroups <- array(0, c(Jg, ndim, G))
      for(g in G:1) {
         Gsgroups[1:CAres[[g]]$J, , g] <- CAres[[g]]$Gs   }   }
   if(Isr != 0) {
      Fsgroupssr <- array(0, c(Isr, ndim, G))
      for(g in G:1) {
         Fsgroupssr[, , g] <- CAres[[g]]$Fssr   }   }
   if(Jsc != 0) {
      Gsgroupssc <- array(0, c(Jsc, ndim, G))
      for(g in G:1) {
         Gsgroupssc[1:CAres[[g]]$Jsc, , g] <- CAres[[g]]$Gssc   }   }
   ### Data from SA
   if (oar != 0) {
      Fsic <- ASG$Fs
      Fsig <- ASG$Fsig
      allFs <- ASG$allFs   }
   if (oac != 0) {
      Gs <- ASG$Gs
      allGs <- ASG$allGs   }
   if(Isr != 0) {
      Fsicsr <- ASG$Fssr
      Fsigsr <- ASG$Fsigsr
      allFssr <- ASG$allFssr   }
   if(Jsc != 0) {
      Gssc <- ASG$Gssc
      allGssc <- ASG$allGssc   }
   ### Data from SA for groups (square plot) and correlations CA and SA (circle)
   eig <- ASG$resin[, 1]
   ctag <- ASG$ctrg
   Fsg <- ASG$Fsg
   fact <- ASG$rCASA
   ###################   CA active   ###########################################
   if(oar != 0 || oac != 0) {
      for(g in 1:G) {
         if(screen == TRUE) {dev.new()}
         titulo <- paste(titulinparc, ": Table ", namegroup[g], sep = "")
         titsub <- "Active elements"
         puntosplot <- c() 
         if (oar != 0 && class(ASG$Fs) != "character" ) {
            namei <- CAres[[g]]$namei
            puntosplot <- rbind(puntosplot, Fsgroups[, , g])   }
         if (oac != 0 && class(ASG$Gs) != "character" ) {
            puntosplot <- rbind(puntosplot, Gsgroups[1:CAres[[g]]$J, , g])   }
         if (length(puntosplot) != 0) {
            resin <- round(CAres[[g]]$resin, 2)
            SAPplot(puntosplot[, s1], puntosplot[, s2])
            SAPtitle(titulo, titsub, s1, s2, resin)
            if (oar != 0 && class(ASG$Fs) != "character" ) {
            SAPplotpuntos(g, puntos = Fsgroups[, , g], gra = c(1:I), s1, s2,
               etiq = CAres[[g]]$namei, colg = colg.F.CA, SAPpch = pchg.F.CA)   }
            if(oac != 0 && class(ASG$Gs) != "character" ) {
               SAPplotpuntos(g, puntos = Gsgroups[1:CAres[[g]]$J, , g],
                  gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej,
                  SAPpch = pchg)   }   }
      }   }
   ###################   CA active + supplementary   ###########################
   if(Isr != 0 || Jsc != 0) {
      for(g in 1:G) {
         if(screen == TRUE) {dev.new()}
         titulo <- paste(titulinparc, ": Table ", namegroup[g], sep = "")
         titsub <- "Active and supplementary elements"
         resin <- round(CAres[[g]]$resin, 2)
         puntosplot <- c() 
         if (oar != 0 && class(ASG$Fs) != "character" ) {
            namei <- CAres[[g]]$namei
            puntosplot <- rbind(puntosplot, Fsgroups[, , g])   }
         if (oac != 0 && class(ASG$Gs) != "character" ) {
            puntosplot <- rbind(puntosplot, Gsgroups[1:CAres[[g]]$J, , g])   }   
         if(Isr != 0) {
            puntosplot <- rbind(puntosplot, Fsgroupssr[ , , g])   }
         if(Jsc != 0) {
            puntosplot <- rbind(puntosplot,
               Gsgroupssc[1:CAres[[g]]$Jsc, , g])   }
         ### Plot
         SAPplot(puntosplot[, s1], puntosplot[, s2])
         SAPtitle(titulo, titsub, s1, s2, resin)
         if (oar != 0 && class(ASG$Fs) != "character" ) {
            SAPplotpuntos(g, puntos = Fsgroups[, , g], gra = c(1:I), s1, s2,
               etiq = CAres[[g]]$namei, colg = colg.F.CA, SAPpch = pchg.F.CA)   }
         if(oac != 0 && class(ASG$Gs) != "character" ) {
            SAPplotpuntos(g, puntos = Gsgroups[1:CAres[[g]]$J, , g],
               gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej,
               SAPpch = pchg)   }   
         if(Isr != 0) {
            puntos <- matrix(0, Isr, ndim)
            puntos[1:Isr, ] <- Fsgroupssr[1:Isr, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:Isr), s1, s2,
               etiq = CAres[[g]]$nameisr, colg = colg.sr.CA,
               SAPpch = pchg.sr.CA)   }
         if(Jsc != 0) {
            puntos <- matrix(0, CAres[[g]]$Jsc, ndim)
            puntos[1:CAres[[g]]$Jsc, ] <- Gsgroupssc[1:CAres[[g]]$Jsc, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:CAres[[g]]$Jsc),
               s1, s2, etiq = CAres[[g]]$namejsc, SAPpch = pchgs)   }   }   }
   else {   }
   ###################   SA   ##################################################
   resin <- round(ASG$resin, 2)
   ###################   SA Fsic   #############################################
   if (oar != 0 && class(ASG$Fs) != "character" ) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall rows", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Overall columns", sep = "")   }        
      titsub <- "Active elements"
      puntosplot <- allFs[, , (G + 1)]
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      SAPplotpuntos(g = G + 1, puntos = allFs[, , (G + 1)], gra = c(1: I),
         s1, s2, etiq = nameic, SAPpch = pchg)   }
   ###################   SA Gsj   ##############################################
   if (oac != 0 && class(ASG$Gs) != "character" ) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Columns", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Rows", sep = "")   }     
      titsub <- "Active elements"
      puntosplot <- Gs
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      for(g in 1:G) {
         SAPplotpuntos(g, puntos = allGs[1:CAres[[g]]$J, , g],
            gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej,
            SAPpch = pchg)   }   }
   ###################   SA Fsic + Gsj   #######################################
   if (oar != 0 && class(ASG$Fs) != "character" && 
      oac != 0 && class(ASG$Gs) != "character") {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall rows and columns", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Rows and overall columns", sep = "")   }
      titsub <- "Active elements"
      puntosplot <- rbind(allFs[, , (G + 1)], Gs)
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      SAPplotpuntos(g = G + 1, puntos = allFs[, , (G + 1)], gra = c(1:I),
         s1, s2, etiq = nameic, SAPpch = pchg)
      if(oac == 1) {
         for(g in 1:G) {
            SAPplotpuntos(g, puntos = allGs[1:CAres[[g]]$J, , g],
               gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej,
               SAPpch = pchg)   }   }   }   
   ###################   SA Fsic + Fsig   ######################################
   if (oar != 0 && class(ASG$Fs) != "character") {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall and partial rows", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Overall and partial columns", sep = "")   }
      titsub <- "Active elements"
      puntosplot <- allFs
      SAPplot(puntosplot[, s1, ], puntosplot[, s2, ])
      SAPtitle(titulo, titsub, s1, s2, resin)
      SAPplotpuntos(G + 1, puntos = allFs[, , (G + 1)], gra = c(1:I), s1, s2,
         etiq = nameic, SAPpch = pchg)
      for(g in 1:G) {
         SAPplotpuntos(g, puntos = allFs[, , g], gra = c(1:I), s1, s2,
            etiq = c(paste(namegroup[g], nameic, sep = "")), SAPpch = pchg)   }
      }
   ###################   end SA active   #######################################
   ###################   SA Fsic + Fsicsr   ####################################
   if(Isr != 0) {
      if(screen == TRUE) {dev.new()}
      ### 
      titulo <- paste(titulin, ": Overall rows", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Overall columns", sep = "")   }
      titsub <- "Active and supplementary elements"
      g <- G + 1
      puntosplot <- c()
      if (oar != 0 && class(ASG$Fs) != "character") {
         puntosplot <- rbind(puntosplot, allFs[, , (G + 1)])   }
      puntosplot <- rbind(puntosplot, allFssr[, , (G + 1)])
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      ### Fsic
      if (oar != 0 && class(ASG$Fs) != "character") {
         puntos <- allFs[, , (G + 1)]
         etiq <- c(nameic)
         SAPplotpuntos(G + 1, puntos = allFs[, , (G + 1)], gra = c(1:I),
            s1, s2, etiq = nameic, SAPpch = pchg)   }
      ### Fsicsr
      puntos <- matrix(0, Isr, ndim)
      puntos[1:Isr, ] <- allFssr[1:Isr, , (G + 1)]
      SAPplotpuntos(G + 1, puntos = puntos, gra = c(1:Isr),
         s1, s2, etiq = nameicsr, SAPpch = pchgs)   }
   else {   }
   ###################   SA Gsj + Gsjsc   ######################################
   if(Jsc != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Columns", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Rows", sep = "")   }
      titsub <- "Active and supplementary elements"
      puntosplot <- c()
      if (oac != 0 && class(ASG$Gs) != "character") {
         puntosplot <- rbind(puntosplot, Gs)   }
      puntosplot <- rbind(puntosplot, Gssc)
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      if (oac != 0 && class(ASG$Gs) != "character") {
         for(g in 1:G) {
            SAPplotpuntos(g, puntos = allGs[1:CAres[[g]]$J, , g],
               gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej,
               SAPpch = pchg)   }   }
      for(g in 1:G) {
         puntos <- matrix(0, Jsc, ndim)
         puntos[1:Jsc, ] <- allGssc[1:Jsc, , g]
         SAPplotpuntos(g, puntos = puntos, gra = c(1:Jsc), s1, s2,
            etiq = paste(namegroup[g], CAres[[g]]$namejsc, sep = ""),
            SAPpch = pchgs)   }   }
   else {   }
   ###################   SA Fsic + Fsicsr + Gsj + Gsjsc   ######################
   if(Isr != 0 || Jsc != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall rows and columns", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Rows and overall columns", sep = "")   }
      titsub <- "Active and supplementary elements"
      puntosplot <- c()
      if (oar != 0 && class(ASG$Fs) != "character") {
         puntosplot <- rbind(puntosplot, allFs[, , (G + 1)])   }
      if (oac != 0 && class(ASG$Gs) != "character") {
         puntosplot <- rbind(puntosplot, Gs)   }
      if(Isr != 0) {
         puntosplot <- rbind(puntosplot, allFssr[, , (G + 1)])   }
      if(Jsc != 0) {
         puntosplot <- rbind(puntosplot, Gssc)   }
      ### Plot
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      if (oar != 0 && class(ASG$Fs) != "character") {
         SAPplotpuntos(g = G + 1, puntos = allFs[, , (G + 1)], gra = c(1:I),
            s1, s2, etiq = nameic, SAPpch = pchg)   }
      if (oac != 0 && class(ASG$Gs) != "character") {
         for(g in 1:G) {
            SAPplotpuntos(g, puntos = allGs[1:CAres[[g]]$J, , g],
               gra = c(1:CAres[[g]]$J), s1, s2,
               etiq = CAres[[g]]$namej, SAPpch = pchg)   }   }
      if(Isr != 0) {
         puntos <- matrix(0, Isr, ndim)
         puntos[1:Isr, ] <- allFssr[1:Isr, , (G + 1)]
         SAPplotpuntos(g = G + 1, puntos = puntos,
            gra = c(1:Isr), s1, s2, etiq = nameicsr, SAPpch = pchgs)   }
      if(Jsc != 0) {
         for(g in 1:G) {
            puntos <- matrix(0, Jsc, ndim)
            puntos[1:Jsc, ] <- allGssc[1:Jsc, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:Jsc), s1, s2,
               etiq = paste(namegroup[g], CAres[[g]]$namejsc, sep = ""),
               SAPpch = pchgs)   }   }   }
   else {   }
   ###################   SA Fsic + Fsig + Fsicsr + Fsigsr   ####################
   if(Isr != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall and partial rows", sep = "")
      if(multiple == 1) {
         titulo <- paste(titulin, ": Overall and partial columns", sep = "")   }
      titsub <- "Active and supplementary elements"
      ### puntosplot rbind(allFs, allFssr) don't work
      puntosplot <- rbind(Fsicsr, Fsigsr)
      if (oar != 0 && class(ASG$Fs) != "character") {
         puntosplot <- rbind(Fsic, Fsig, puntosplot)   }
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      ### Active
      if (oar != 0 && class(ASG$Fs) != "character") {
         puntos <- allFs
         etiq <- c(nameic)
         SAPplotpuntos(G + 1, puntos = allFs[, , (G + 1)], gra = c(1:I),
            s1, s2, etiq = nameic, SAPpch = pchg)
         for(g in 1:G) {
            SAPplotpuntos(g, puntos = allFs[, , g], gra = c(1:I), s1, s2,
               etiq = c(paste(namegroup[g], nameic, sep = "")),
               SAPpch = pchg)   }   }
      ### Supplementary
      puntos <- allFssr
      etiq <- c(nameicsr)
      ### Fcsr
      points(c(puntos[, s1, (G + 1)]), c(puntos[, s2, (G + 1)]),
         pch =pchgs[G + 1], col = colg[G + 1], cex = 0.7)
      text(c(puntos[, s1, (G + 1)]), c(puntos[, s2, (G + 1)]),
         c(paste("  ", etiq, sep = "")), col = colg[G + 1], adj = 0, cex = 0.7)
      ### Fgsr
      for(g in 1:G) {
         points(c(puntos[, s1, g]), c(puntos[, s2, g]), pch =pchgs[g],
            col = colg[g], cex = 0.7)
         text(c(puntos[, s1, g]), c(puntos[, s2, g]),
            c(paste("  ", namegroup[g], etiq, sep = "")),
            col = colg[g], adj = 0, cex = 0.7 )   }   }
   else {   }
   ###################   end SA active + suplemmentary   #######################
   ###################   Fsg y SA + CA (sqauare & circle)   ####################
   ###################   SA Fsg (square)   #####################################
   if(screen == TRUE) {dev.new()}
   titulo <- paste(titulin, ": Tables", sep = "")
   puntosplot <- Fsg
   puntos <- Fsg
   etiq <- namegroup
   oldpar <- par(pin = c(5, 5), las = 1)
   plot(puntosplot[, s1], puntosplot[, s2], type = "n",
      xlab = paste("  "), ylab = paste("  "),
      col = 1, xlim = c(0, 1), ylim = c(0, 1), cex = 0.5)
   title(titulo, col = 1, cex = 1.2)
   title("",
      xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)    "),
      ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)    "),
      adj = 1, col = 1)
   for(g in 1:G) {
      points(c(puntos[g, s1]), c(puntos[g, s2]), pch =pchg[g], col = colg[g],
         cex = 0.7)
      text(c(puntos[g, s1]), c(puntos[g, s2]),
         c(paste("  ", etiq[g], sep = "")), col = colg[g], adj = 0,
         cex = 0.7)   }
   ### Volver a pin originales (no grafico cuadrado)
   par(oldpar)
   ###################   SA + CA (circle)   ####################################
   if(screen == TRUE) {dev.new()}
   titulo <- "Relation between factors of CA and SA"
   puntosplot <- fact
   puntos <- fact
   etiq <- namegroup
   ### Circle
   z <- seq(0, 2 * pi, length = 1000)
   x <- sin(z)
   y <- cos(z)
   oldpar <- par(pin = c(5, 5), las = 1)
#RS FALSE F   plot(x, y, type = "l", axes = FALSE, xlab = "", ylab = "")
   plot(x, y, type = "l", axes = FALSE, xlab = "", ylab = "")
#RS FALSE F   axis(1, pos = 0, labels = FALSE)
   axis(1, pos = 0, labels = FALSE)
#RS FALSE F   axis(2, pos = 0, labels = FALSE)
   axis(2, pos = 0, labels = FALSE)
   title(titulo, col = 1, cex = 1.2)
   title("",
      xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)    "),
      ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)    "),
      adj = 1, col = 1)
   for(sac in 1:2) {
      for(g in 1:G) {
         points(c(puntos[sac, s1, g]), c(puntos[sac, s2, g]), pch =pchg[g],
            col = colg[g], cex = 0.7)
#S#          arrows(x1 = 0, y1 = 0, x2 = c(puntos[sac, s1, g]), y2 = c(puntos[sac, s2, g]), col = colg[g], size = 0.2)
#R#          arrows(x0 = 0, y0 = 0, x1 = c(puntos[sac, s1, g]), y1 = c(puntos[sac, s2, g]), col = colg[g], lwd=1.5)
#SR          arrows(0, 0, c(puntos[sac, s1, g]), c(puntos[sac, s2, g]), col = colg[g], lwd = 1.5)
         arrows(0, 0, c(puntos[sac, s1, g]), c(puntos[sac, s2, g]),
            col = colg[g],
            lwd = 1.5)
         text(c(puntos[sac, s1, g]), c(puntos[sac, s2, g]),
            c(paste("  F", sac, "(", etiq[g], ")", sep = "")), col = colg[g],
            adj = 0, cex = 0.7)   }   }
   ### Volver a pin originales (no grafico redondo)
   par(oldpar)
   cat ("\n")
   ###################   fin Fsg y SA + CA (square & circle)   #################
}
### End plot.SimAn
################################################################################



################################################################################
#####               End of SimultAnR (Simultaneous Analysis) Package       #####
################################################################################








