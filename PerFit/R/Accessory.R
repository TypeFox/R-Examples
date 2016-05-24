# final.PFS(): Compute final PFS vector ----
final.PFS <- function(res.red, all.0s, all.1s, N)
{
  if (length(c(all.0s, all.1s)) > 0)
  {
    res <- rep(NA, N)
    res[(1:N)[-c(all.0s, all.1s)]] <- res.red
  } else
  {
    res <- res.red
  }
  return(res)
}

# res.process(): Process the PFS scores from vector to data frame ----
res.process <- function(matrix, N, res)
{
  res <- data.frame(PFscores = round(res, 4))
  if (is.null(row.names(matrix)))
  {
    row.names(res) <- rep(paste0("Resp.", 1:N))
  } else
  {
    row.names(res) <- row.names(matrix)
  }
  return(res)
}

# export.res.NP(): Export results (nonparametric) ----
export.res.NP <- function(matrix, N, res, PFStatistic, part.res, Ncat, NAs, 
                          IRT.PModel, IP, Ability.PModel, Ability, IP.NA, Ability.NA, NAs.imp)
{
  res <- list(PFscores  = res.process(matrix, N, res), PFStatistic = PFStatistic, 
              PerfVects = part.res$perfect.vectors, ID.all0s = part.res$all.0s, ID.all1s = part.res$all.1s, 
              Matrix = matrix, Ncat=Ncat, 
              IRT.PModel = if(IP.NA & (NAs == "PModel")) {IRT.PModel} else {NULL}, 
              IP = IP, 
              Ability.PModel = if(Ability.NA & (NAs == "PModel")) {Ability.PModel} else {NULL}, 
              Ability = Ability, 
              NAs.method = if(NAs.imp) {NAs} else {NULL})
  class(res) <- "PerFit"
  return(res)
}

# export.res.P(): Export results (parametric) ----
export.res.P <- function(matrix, N, res, PFStatistic, part.res, Ncat, NAs, 
                         IRT.PModel, IP, Ability.PModel, Ability, IP.NA, Ability.NA, NAs.imp)
{
  res <- list(PFscores  = res.process(matrix, N, res), PFStatistic = PFStatistic, 
              PerfVects = part.res$perfect.vectors, ID.all0s = part.res$all0s, ID.all1s = part.res$all.1s, 
              Matrix = matrix, Ncat=Ncat, 
              IRT.PModel = if(IP.NA) {IRT.PModel} else {NULL}, 
              IP = IP, 
              Ability.PModel = if(Ability.NA) {Ability.PModel} else {NULL}, 
              Ability = Ability, 
              NAs.method = if(NAs.imp) {NAs} else {NULL})
  class(res) <- "PerFit"
  return(res)
}

# estIP(): Estimate item parameters if not provided (using 'irtoys') ----
estIP <- function(matrix, ip, model)
{
  I <- dim(matrix)[2]
  if (is.null(ip)) 
  {
    # Sanity check - IP model:
    Sanity.IPm(model)
    ip <- est(matrix, model, engine="ltm", rasch=TRUE, nqp=20)$est
  } else 
  {
    ip <- as.matrix(ip)
    # Sanity check - IP matrix adequacy:
    Sanity.IPa(ip, I)
  }
  ip
}

# estIP.poly(): Estimate item parameters if not provided (polytomous) ----
estIP.poly <- function(matrix, Ncat, ip, model)
{
  I       <- dim(matrix)[2]
  matrix2 <- data.frame(apply(matrix + 1, 2, as.factor)) # eliminates item levels with no answers
  # Sanity check - IP model (polytomous):
  Sanity.IPm.poly(model)
  if (is.null(ip)) 
  {
    ip <- switch(model,
                 PCM  = gpcm(matrix2, constraint = "rasch", IRT.param = TRUE),
                 GPCM = gpcm(matrix2, constraint = "gpcm" , IRT.param = TRUE),
                 GRM  = grm (matrix2, constrained = FALSE , IRT.param = TRUE))
    ip.coef <- coef(ip)
  } else 
  {
#     ip <- as.matrix(ip)
#     # Sanity check - IP matrix adequacy (polytomous):
#     Sanity.IPa.poly(ip, I, Ncat)
#     ip.coef <- ip
    
    # Sanity check - IP matrix adequacy (polytomous):
    Sanity.IPa.poly(ip, I, Ncat)
    ip     <- as.matrix(ip)
    ip.ltm <- cbind(ip[, -ncol(ip)] * ip[, ncol(ip)], ip[, ncol(ip)])
    ip     <- switch(model,
                 PCM  = gpcm(matrix2, constraint = "rasch", IRT.param = TRUE, 
                             start.val = unlist(apply(ip.ltm, 1, list), recursive = FALSE), 
                             control = list(iter.qN = 0, optimizer = "optim")),
                 GPCM = gpcm(matrix2, constraint = "gpcm" , IRT.param = TRUE, 
                             start.val = unlist(apply(ip.ltm, 1, list), recursive = FALSE), 
                             control = list(iter.qN = 0, optimizer = "optim")),
                 GRM  = grm (matrix2, constrained = FALSE , IRT.param = TRUE, 
                             start.val = unlist(apply(ip.ltm, 1, list), recursive = FALSE), 
                             control = list(iter.qN = 0)))
    ip.coef <- coef(ip)
  }
  # In case NOT all answer categories of all items were used:
  if (is.list(ip.coef)) 
  {
    abs.freqs <- apply(matrix, 2, table)
    abs.freqs <- lapply(abs.freqs, function(vect) as.numeric(names(vect)))
    tmp       <- matrix(NA, nrow = I, ncol = Ncat)
    for (i in 1:I) 
    {
      tmp[i,abs.freqs[[i]][-length(abs.freqs[[i]])]+1] <- ip.coef[[i]][-length(abs.freqs[[i]])]
      tmp[i,Ncat] <- ip.coef[[i]][length(ip.coef[[i]])]
    }
    ip.coef <- tmp
  }
  # 
  list(ip.coef, ip)
}

# estAb(): Estimate ability parameters if not provided (using 'ltm') ----
estAb <- function(matrix, ip, ability, method, mu, sigma)
{
  N <- dim(matrix)[1]
  if (is.null(ability))
  {
    # Sanity check - Ability method:
    Sanity.Abm(method)
    ability <- switch(method,
                      ML = mlebme(matrix, ip, mu, sigma, method = "ML")[, 1],
                      BM = mlebme(matrix, ip, mu, sigma, method = "BM")[, 1],
                      WL = wle(matrix, ip)[, 1])
  } else
  {
    ability <- as.vector(ability)
    # Sanity check - Ability matrix adequacy:
    Sanity.ABa(ability, N)
  }
  ability
}

# estAb.poly(): Estimate ability parameters if not provided (using 'ltm') (polytomous) ----
estAb.poly <- function(matrix, ip.ltm, ability, method)
{
  N       <- dim(matrix)[1]
  matrix2 <- data.frame(apply(matrix + 1, 2, as.factor)) # eliminates item levels with no answers
  if (is.null(ability)) 
  {
    # Sanity check - Ability method:
    Sanity.Abm.poly(method)
    ability <- ltm::factor.scores(ip.ltm, resp.patterns = matrix2, method = method)
    ability <- ability$score.dat[, ncol(ability$score.dat) - 1]
  } else
  {
    ability <- as.vector(ability)
    # Sanity check - Ability matrix adequacy:
    Sanity.ABa(ability, N)
  }
  ability
}

# estP.CRF(): Compute P.CRF (polytomous) ----
estP.CRF <- function(I, Ncat, model, ip.coef, ability)
{
  N <- length(ability)
  M <- Ncat - 1
  #  Based on GRM:
  if (model == "GRM") 
  {
    # P.ISRF is N x I*M:
    P.ISRF <- t(
      sapply(ability, function(x)
      {
        as.vector(t(1 / (1 + exp(-ip.coef[, ncol(ip.coef)]*(x - ip.coef[, -ncol(ip.coef)])))))
      })
    )
    # Fix for datasets with non-chosen answer options (the NA entries):
    #   1s for NAs in the first item steps
    #   0s for NAs in the last item steps
    #   entry (x+1) for NAs in entry x
    if (sum(is.na(ip.coef)) > 0) 
    {
      first.cols  <- (which(is.na(ip.coef[, 1]))-1) * M + 1; P.ISRF[, first.cols][is.na(P.ISRF[, first.cols])] <- 1
      last.cols   <- which(is.na(ip.coef[, M]))*M; P.ISRF[, last.cols][is.na(P.ISRF[, last.cols])] <- 0
      middle.cols <- sort(which(is.na(t(cbind(rep(0, I), ip.coef[, -c(1, M, Ncat)], rep(0, I))))), decreasing = TRUE)
      for (i in 1:length(middle.cols)){P.ISRF[, middle.cols] <- P.ISRF[, middle.cols + 1]}
    }
    P.CRF <- matrix(, nrow = N, ncol = I * Ncat)
    for (i in 1:I) 
    {
      P.ISRF.item <- -cbind(rep(1, N), P.ISRF[, ((i-1)*M+1):(i*M)], rep(0, N))
      P.CRF[, ((i-1)*Ncat+1):(i*Ncat)] <- P.ISRF.item[, (2:(M+2))] - P.ISRF.item[, (1:(M+1))]
    }
  }
  
  # Based on PCM or GPCM:
  if (model == "PCM" | model == "GPCM") 
  {
    lin.it <- t(sapply(ability, function(x){as.vector(t(ip.coef[, ncol(ip.coef)] * (x - ip.coef[, -ncol(ip.coef)])))}))
    lin.it[, is.na(as.vector(t(ip.coef[, -ncol(ip.coef)])))] <- 0 # NAs -> 0 to eliminate these terms from the sums
    P.CRF  <- matrix(, nrow = N, ncol = I * Ncat)
    tri    <- matrix(1, nrow = M + 1, ncol = M + 1)
    tri    <- upper.tri(tri, diag = TRUE)
    for (i in 1:I) 
    {
      num <- exp(cbind(rep(0, N), lin.it[, ((i-1)*M+1):(i*M)]) %*% tri)
      P.CRF[, ((i-1) * Ncat+1):(i*Ncat)] <- num / rowSums(num)
    }  
  }
  return(P.CRF)
}
