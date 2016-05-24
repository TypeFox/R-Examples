MalikPvalue <-
function (hfobj, N = 500,pnote=TRUE) 
{
  tallset <- hfobj$tall
  y <- tallset$y
  rows <- tallset$rows
  cols <- tallset$cols
  a <- max(as.numeric(as.character(rows)))
  b <- max(as.numeric(as.character(cols)))
  mod <- lm(y ~ rows + cols)
  r <- resid(mod)
  names(r) <- NULL
  rmat <- matrix(r, nrow = length(r), ncol = 1)
  kmean <- kmeans(x = rmat, centers = 3, nstart = 100)
  assn <- kmean$cluster
  modclus <- lm(y ~ rows + cols + as.factor(assn))
  amodclus <- anova(modclus)
  Tc <- (amodclus[3, 2]/amodclus[3, 1])/(amodclus[4, 2]/amodclus[4,1])
  Tcsim <- c()
  for (i in 1:N) {
    ysim <- rnorm((a * b), 0, 1)
    modsim <- lm(ysim ~ rows + cols)
    rsim <- resid(modsim)
    rmatsim <- matrix(rsim, nrow = length(rsim), ncol = 1)
    kmeansim <- kmeans(x = rmatsim, centers = 3, nstart = 100)
    assnsim <- kmeansim$cluster
    modclussim <- lm(ysim ~ rows + cols + as.factor(assnsim))
    amodclussim <- anova(modclussim)
    Tcsim[i] <- (amodclussim[3, 2]/amodclussim[3, 1])/(amodclussim[4,2]/amodclussim[4, 1])
  }
  malik.p <- mean(Tcsim > Tc)
  # list(pvalue = malik.p)
if(pnote){
  cat(paste("(Pvalue from Malik's test estimated with N=", 
            N, " Monte Carlo datasets) \n",sep=""))}
  list(pvalue = malik.p, Tc=Tc, N=N)
}
