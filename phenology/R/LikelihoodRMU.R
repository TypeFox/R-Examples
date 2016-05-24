.LikelihoodRMU <- function(x, fixed, model.trend, RMU.data, colname.year=NULL, RMU.names=NULL, index=NULL) {
  
  x<-c(x, fixed)
  
  if (is.null(index)) {
    nm <- colnames(RMU.data)
    index.year <- which(nm==colname.year)
    index.mean <- match(RMU.names$mean, nm)
    index.se <- match(RMU.names$se, nm)
    d <- RMU.data[,index.mean]
    nabeach <- colnames(d)
    nbeach <- length(nabeach)
    nyear <- dim(RMU.data)[1]
    nayear <- RMU.data[,index.year]
    maxL <- 1E9

    index <- list(year=index.year, mean=index.mean, se=index.se, colnames=nabeach, 
                  nyear=nyear, nbeach=nbeach)
  } else {
    nabeach <- index$colnames
    nbeach <- length(nabeach)
    nyear <- dim(RMU.data)[1]
    maxL <- index$maxL
  }
  
  SD <- x[substr(names(x), 1, 3)=="SD_"]
  if (length(SD)==0) {
    SD <- rep(0, nbeach)
    names(SD) <- nabeach
  } else {
    if (length(SD)==1) {
      SD <- rep(abs(SD), nbeach)
      names(SD) <- nabeach
    } else {
      SD <- abs(SD)
    }
  }
  
  dtaL_obs <- RMU.data[,index$mean]
  dtaL_SD <- RMU.data[,index$se]+matrix(rep(SD, nyear), nrow=nyear, byrow=TRUE)

  #___________________________________________________________
  # Je crée le vecteur avec les proportions de chaque site
  #___________________________________________________________
  
  La0 <- x[paste0("a0_", nabeach[paste0("a0_", nabeach) %in% names(x)])]
  La1 <- x[paste0("a1_", nabeach[paste0("a1_", nabeach) %in% names(x)])]
  La2 <- x[paste0("a2_", nabeach[paste0("a2_", nabeach) %in% names(x)])]
  
  if (any(is.na(La2))) {
    La2 <- La0
    La2[] <- 0
    names(La2) <- gsub("^a0_", "a2_", names(La0))
  }
  if (any(is.na(La1))) {
    La1 <- La0
    La1[] <- 0
    names(La1) <- gsub("^a0_", "a1_", names(La0))
  }
  
  map <- matrix(rep(NA, nyear*(nbeach-1)), ncol=nbeach-1)
  for (i in 1:(nbeach-1)) {
    map[,i] <- 1/(1+exp(La2[i]*(1:nyear)^2+La1[i]*(1:nyear)+La0[i]))
  }
  mapp <- matrix(rep(NA, nyear*nbeach), ncol=nbeach, 
                 dimnames=list(NULL, beach=nabeach))
  for (j in 1:nyear) {
    mapp[j,] <- getFromNamespace(".p.multinomial", ns="phenology")(map[j,])
  }
  
  if (model.trend=="year-specific") {
    Tot <- abs(x[paste0("T_", RMU.data[,index$year])])
    dtaL_theo <- matrix(rep(Tot, nbeach) , ncol=nbeach, byrow = FALSE, 
                        dimnames=list(NULL, beach=nabeach))
  }
  if (model.trend=="constant") {
    Tot <- abs(x["T_"])
    dtaL_theo <- matrix(rep(Tot,nbeach*nyear) , ncol=nbeach, byrow = FALSE, 
                        dimnames=list(NULL, beach=nabeach))
  }
  if (model.trend=="exponential") {
    Tot <- abs(x["T_"])*exp(x["r"]*(1:nyear))
    dtaL_theo <- matrix(rep(Tot,nbeach) , ncol=nbeach, byrow = FALSE, 
                        dimnames=list(NULL, beach=nabeach))
  }

  # dtaL_theo <- dtaL_theo[, names(pp)]
  dtaL_obs <- dtaL_obs[, nabeach]
  dtaL_theo <- dtaL_theo*mapp
  valide <- !is.na(dtaL_obs)
  
#  if ((any(is.na(dtaL_theo))) | (min(colMeans(dtaL_theo, na.rm=TRUE), na.rm=TRUE)<1)) {
  if ((any(is.na(dtaL_theo)))) {
      L <- maxL
  } else {
  
  L <- sum(-dnorm(x=dtaL_obs[valide], 
                mean=dtaL_theo[valide], 
                sd=dtaL_SD[valide], log=TRUE))
}
  #	if (!is.finite(L)) {print(x)}
  return(L)
  
}

# je donne les N probabilités et j'ai les N-1 probabilités conditionelles 
.inv.p.multinomial <- function(x) {
  p <- x[1]
  for (i in 2:(length(x)-1)) 
    p <- c(p, x[i]/(prod(1-p[1:(i-1)])))
  return(p)
}

# je donne les N-1 probabilités conditionelles et j'ai les N probabilités
.p.multinomial <- function(p) c(p, 1)*c(1, sapply(seq_along(p), function(i) prod(1-p[1:i])))

