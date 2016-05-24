"elo" <- function(x, status=NULL, init=2200, gamma=0, kfac=27, history=FALSE, sort=TRUE, ...)
{
  if(length(init) != 1) stop("the length of 'init' must be one")
  if(ncol(x) != 4) stop("'x' must have four variables")
  if(nrow(x) == 0) {
    if(is.null(status)) stop("'x' is empty and 'status' is NULL")
	lout <- list(ratings = status, history = NULL, gamma = gamma, kfac=kfac, type = "Elo")
    class(lout) <- "rating"
    return(lout)
  }
  gammas <- rep(gamma, length.out = nrow(x))         
  names(x) <- c("Month","White","Black","Score")
  if(!is.numeric(x$Month)) 
    stop("Time period must be numeric")
  if(!is.numeric(x$White) && !is.character(x$White))
    stop("Player identifiers must be numeric or character")
  if(!is.numeric(x$Black) && !is.character(x$Black))
    stop("Player identifiers must be numeric or character")	
  if(!is.numeric(x$Score) || any(x$Score > 1) || any(x$Score < 0))
    stop("Game scores must be in the interval [0,1]")
  
  play <- sort(unique(c(x$White,x$Black)))
  np <- length(play)
  x$White <- match(x$White, play)
  x$Black <- match(x$Black, play)

  if(!is.null(status)) {
    npadd <- play[!(play %in% status$Player)]
	zv <- rep(0, length(npadd))
    npstatus <- data.frame(Player = npadd, Rating = rep(init,length(npadd)), Games = zv, 
      Win = zv, Draw = zv, Loss = zv, Lag = zv)
	if(!("Games" %in% names(status))) status <- cbind(status, Games = 0)
	if(!("Win" %in% names(status))) status <- cbind(status, Win = 0)
	if(!("Draw" %in% names(status))) status <- cbind(status, Draw = 0)
	if(!("Loss" %in% names(status))) status <- cbind(status, Loss = 0)
    if(!("Lag" %in% names(status))) status <- cbind(status, Lag = 0)
    status <- rbind(status[,c("Player","Rating","Games","Win","Draw","Loss","Lag")], npstatus)
    rinit <- status[,2]
    ngames <- status[,3]
    nwin <- status[,4]
    ndraw <- status[,5]
    nloss <- status[,6]
    nlag <- status[,7]
    names(rinit) <- names(ngames) <- status$Player
  }
  else {
    rinit <- rep(init, length.out=np)
    ngames <- nwin <- ndraw <- nloss <- nlag <- rep(0, length.out=np)
    names(rinit) <- names(ngames) <- names(nlag) <- play
  }
 
  if(!all(names(rinit) == names(ngames)))
    stop("names of ratings and ngames are different")
  if(!all(play %in% names(rinit))) 
    stop("Payers in data are not within current status")

  nm <- length(unique(x$Month))
  curplay <- match(play, names(rinit))
  orats <- rinit[-curplay] 
  ongames <- ngames[-curplay]
  onwin <- nwin[-curplay]
  ondraw <- ndraw[-curplay]
  onloss <- nloss[-curplay]
  olag <- nlag[-curplay]
  olag[ongames != 0] <- olag[ongames != 0] + nm
  crats <- rinit[curplay] 
  ngames <- ngames[curplay] 
  nwin <- nwin[curplay]
  ndraw <- ndraw[curplay]
  nloss <- nloss[curplay]
  nlag <- nlag[curplay]

  gammas <- split(gammas, x$Month)
  x <- split(x, x$Month)
  if(history) {
    histry <- array(NA, dim=c(np,nm,3), dimnames=list(play,1:nm,c("Rating","Games","Lag")))
  }

  for(i in 1:nm) {
    traini <- x[[i]]
	gammai <- gammas[[i]] 
    nr <- nrow(traini)
    dscore <- .C("elo_c",
      as.integer(np), as.integer(nr), as.integer(traini$White-1), as.integer(traini$Black-1),
      as.double(traini$Score), as.double(crats), as.double(gammai), dscore = double(np))$dscore
    if(!is.function(kfac)) {
      crats <- crats + kfac * dscore
    }
    else {
      crats <- crats + kfac(crats, ngames, ...) * dscore
    }
    trainipl <- c(traini$White,traini$Black)
    trainiplw <- c(traini$White[traini$Score==1],traini$Black[traini$Score==0])
    trainipld <- c(traini$White[traini$Score==0.5],traini$Black[traini$Score==0.5])
    trainipll <- c(traini$White[traini$Score==0],traini$Black[traini$Score==1])
    ngames <- ngames + tabulate(trainipl, np)
    nwin <- nwin + tabulate(trainiplw, np)
    ndraw <- ndraw + tabulate(trainipld, np)
    nloss <- nloss + tabulate(trainipll, np)
    playi <- unique(trainipl)
    nlag[ngames!=0] <- nlag[ngames!=0] + 1
    nlag[playi] <- 0

     if(history) {
      histry[,i,1] <- crats
      histry[,i,2] <- ngames
      histry[,i,3] <- nlag
    }
  }
  if(!history) histry <- NULL
  player <- suppressWarnings(as.numeric(names(c(crats,orats))))
  if (any(is.na(player))) player <- names(c(crats,orats))
  dfout <- data.frame(Player=player, Rating=c(crats,orats), Games=c(ngames,ongames), 
    Win=c(nwin,onwin), Draw=c(ndraw,ondraw), Loss=c(nloss,onloss), Lag=c(nlag,olag),
	stringsAsFactors = FALSE)
  if(sort) dfout <- dfout[order(dfout$Rating,decreasing=TRUE),] else dfout <- dfout[order(dfout$Player),]
  row.names(dfout) <- 1:nrow(dfout)  

  lout <- list(ratings = dfout, history = histry, gamma = gamma, kfac=kfac, type = "Elo")
  class(lout) <- "rating"
  lout
}

"fide" <- function(x, status=NULL, init=2200, gamma=0, kfac=kfide,
  history=FALSE, sort=TRUE, ...)
{  
  if(length(init) != 1) stop("the length of 'init' must be one")
  if(ncol(x) != 4) stop("'x' must have four variables")
  if(nrow(x) == 0) {
    if(is.null(status)) stop("'x' is empty and 'status' is NULL")
	lout <- list(ratings = status, history = NULL, gamma = gamma, kfac=kfac, type = "Elo")
    class(lout) <- "rating"
    return(lout)
  }
  gammas <- rep(gamma, length.out = nrow(x)) 
  names(x) <- c("Month","White","Black","Score")
  if(!is.numeric(x$Month)) 
    stop("Time period must be numeric")
  if(!is.numeric(x$White) && !is.character(x$White))
    stop("Player identifiers must be numeric or character")
  if(!is.numeric(x$Black) && !is.character(x$Black))
    stop("Player identifiers must be numeric or character")	
  if(!is.numeric(x$Score) || any(x$Score > 1) || any(x$Score < 0))
    stop("Game scores must be in the interval [0,1]")
	
  play <- sort(unique(c(x$White,x$Black)))
  np <- length(play)
  x$White <- match(x$White, play)
  x$Black <- match(x$Black, play)

  if(!is.null(status)) {
    npadd <- play[!(play %in% status$Player)]
	zv <- rep(0, length(npadd))
	ev <- rep(as.numeric(init > 2400), length.out=length(npadd))
    npstatus <- data.frame(Player = npadd, Rating = rep(init,length(npadd)), Games = zv, 
      Win = zv, Draw = zv, Loss = zv, Lag = zv, Elite = ev, Opponent = zv)
	if(!("Games" %in% names(status))) status <- cbind(status, Games = 0)
	if(!("Win" %in% names(status))) status <- cbind(status, Win = 0)
	if(!("Draw" %in% names(status))) status <- cbind(status, Draw = 0)
	if(!("Loss" %in% names(status))) status <- cbind(status, Loss = 0)
    if(!("Lag" %in% names(status))) status <- cbind(status, Lag = 0)
    if(!("Elite" %in% names(status))) status <- cbind(status, Elite = as.numeric(status$Player >= 2400))
    if(!("Opponent" %in% names(status))) status <- cbind(status, Opponent = status$Rating)
    status <- rbind(status[,c("Player","Rating","Games","Win","Draw","Loss","Lag","Elite","Opponent")], npstatus)
    rinit <- status[,2]
    ngames <- status[,3]
    nwin <- status[,4]
    ndraw <- status[,5]
    nloss <- status[,6]
    nlag <- status[,7]
    elite <- status[,8]
    opponent <- status[,9]
    names(rinit) <- names(ngames) <- status$Player
  }
  else {
    rinit <- rep(init, length.out=np)
    ngames <- nwin <- ndraw <- nloss <- nlag <- elite <- opponent <- rep(0, length.out=np)
    names(rinit) <- names(ngames) <- names(nlag) <- play
  }
 
  if(!all(names(rinit) == names(ngames)))
    stop("names of ratings and ngames are different")
  if(!all(play %in% names(rinit))) 
    stop("Payers in data are not within current status")

  nm <- length(unique(x$Month))
  curplay <- match(play, names(rinit))
  orats <- rinit[-curplay] 
  ongames <- ngames[-curplay]
  onwin <- nwin[-curplay]
  ondraw <- ndraw[-curplay]
  onloss <- nloss[-curplay]
  olag <- nlag[-curplay]
  olag[ongames != 0] <- olag[ongames != 0] + nm
  oelite <- elite[-curplay]
  oopponent <- opponent[-curplay]
  crats <- rinit[curplay] 
  ngames <- ngames[curplay] 
  nwin <- nwin[curplay]
  ndraw <- ndraw[curplay]
  nloss <- nloss[curplay]
  nlag <- nlag[curplay]
  elite <- elite[curplay]
  opponent <- opponent[curplay]

  gammas <- split(gammas, x$Month)
  x <- split(x, x$Month)
  if(history) {
    histry <- array(NA, dim=c(np,nm,3), dimnames=list(play,1:nm,c("Rating","Games","Lag")))
  }

  for(i in 1:nm) {
    traini <- x[[i]]
	gammai <- gammas[[i]]
    nr <- nrow(traini)
    trainiW <- traini$White; trainiB <- traini$Black; trainiS <- traini$Score

    dscore <- .C("elo_c",
      as.integer(np), as.integer(nr), as.integer(trainiW-1), as.integer(trainiB-1),
      as.double(trainiS), as.double(crats), as.double(gammai), dscore = double(np))$dscore
    if(!is.function(kfac)) {
      crats <- crats + kfac * dscore
    }
    else {
      crats <- crats + kfac(crats, ngames, elite, ...) * dscore
    }

    trainipl <- c(trainiW,trainiB)
    trainiplw <- c(trainiW[traini$Score==1],trainiB[traini$Score==0])
    trainipld <- c(trainiW[traini$Score==0.5],trainiB[traini$Score==0.5])
    trainipll <- c(trainiW[traini$Score==0],trainiB[traini$Score==1])
    ngamesi <- tabulate(trainipl, np)
    ngames <- ngames + ngamesi
    nwin <- nwin + tabulate(trainiplw, np)
    ndraw <- ndraw + tabulate(trainipld, np)
    nloss <- nloss + tabulate(trainipll, np)
    playi <- unique(trainipl)
    nlag[ngames!=0] <- nlag[ngames!=0] + 1
    nlag[playi] <- 0
    elite[crats >= 2400] <- 1
    
    opponentiw <- sapply(split(crats[trainiB], trainiW), sum) 
    opponentib <- sapply(split(crats[trainiW], trainiB), sum) 
    opponenti <- numeric(np)
    opponenti[as.numeric(names(opponentiw))] <- opponentiw
    opponenti[as.numeric(names(opponentib))] <- opponenti[as.numeric(names(opponentib))] + opponentib
    opponent[ngames!=0] <- (((ngames-ngamesi)/ngames)*opponent + opponenti/ngames)[ngames!=0]

    if(history) {
      histry[,i,1] <- crats
      histry[,i,2] <- ngames
      histry[,i,3] <- nlag
    }
  }
  if(!history) histry <- NULL
  player <- suppressWarnings(as.numeric(names(c(crats,orats))))
  if (any(is.na(player))) player <- names(c(crats,orats))
  dfout <- data.frame(Player=player, Rating=c(crats,orats), Games=c(ngames,ongames), 
    Win=c(nwin,onwin), Draw=c(ndraw,ondraw), Loss=c(nloss,onloss), Lag=c(nlag,olag),
    Elite=c(elite,oelite), Opponent=c(opponent,oopponent),
	stringsAsFactors = FALSE)
  if(sort) dfout <- dfout[order(dfout$Rating,decreasing=TRUE),] else dfout <- dfout[order(dfout$Player),]
  row.names(dfout) <- 1:nrow(dfout)  

  lout <- list(ratings = dfout, history = histry, gamma = gamma, kfac=kfac, type = "Elo")
  class(lout) <- "rating"
  lout
}

"glicko" <- function(x, status=NULL, init=c(2200,300), gamma=0, cval=15, history=FALSE, sort=TRUE, ...)
{ 
  if(length(init) != 2) stop("the length of 'init' must be two")
  if(ncol(x) != 4) stop("'x' must have four variables")
  if(nrow(x) == 0) {
    if(is.null(status)) stop("'x' is empty and 'status' is NULL")
	lout <- list(ratings=status, history=NULL, gamma=gamma, cval=cval, type = "Glicko")
    class(lout) <- "rating"
    return(lout)
  }
  gammas <- rep(gamma, length.out = nrow(x)) 
  names(x) <- c("Month","White","Black","Score")
  if(!is.numeric(x$Month)) 
    stop("Time period must be numeric")
  if(!is.numeric(x$White) && !is.character(x$White))
    stop("Player identifiers must be numeric or character")
  if(!is.numeric(x$Black) && !is.character(x$Black))
    stop("Player identifiers must be numeric or character")	
  if(!is.numeric(x$Score) || any(x$Score > 1) || any(x$Score < 0))
    stop("Game scores must be in the interval [0,1]")
	
  play <- sort(unique(c(x$White,x$Black)))
  np <- length(play)
  x$White <- match(x$White, play)
  x$Black <- match(x$Black, play)

  if(!is.null(status)) {
    npadd <- play[!(play %in% status$Player)]
	zv <- rep(0, length(npadd))
    npstatus <- data.frame(Player = npadd, Rating = rep(init[1],length(npadd)), 
      Deviation = rep(init[2],length(npadd)), Games = zv, Win = zv, Draw = zv, 
	  Loss = zv, Lag = zv)
	if(!("Games" %in% names(status))) status <- cbind(status, Games = 0)
	if(!("Win" %in% names(status))) status <- cbind(status, Win = 0)
	if(!("Draw" %in% names(status))) status <- cbind(status, Draw = 0)
	if(!("Loss" %in% names(status))) status <- cbind(status, Loss = 0)
    if(!("Lag" %in% names(status))) status <- cbind(status, Lag = 0)
    status <- rbind(status[,c("Player","Rating","Deviation","Games","Win","Draw","Loss","Lag")], npstatus)
    rinit <- status[,2]
    dinit <- status[,3]
    ngames <- status[,4]
    nwin <- status[,5]
    ndraw <- status[,6]
    nloss <- status[,7]
    nlag <- status[,8]
    names(rinit) <- names(dinit) <- names(ngames) <- status$Player
  }
  else {
    rinit <- rep(init[1], length.out=np)
    dinit <- rep(init[2], length.out=np)
    ngames <- nwin <- ndraw <- nloss <- rep(0, length.out=np)
    nlag <- rep(0,np)
    names(rinit) <- names(dinit) <- names(ngames) <- names(nlag) <- play
  }
 
  if(!all(names(rinit) == names(ngames)))
    stop("names of ratings and deviations are different")
  if(!all(names(rinit) == names(ngames)))
    stop("names of ratings and ngames are different")
  if(!all(play %in% names(rinit))) 
    stop("Payers in data are not within current status")

  nm <- length(unique(x$Month))
  curplay <- match(play, names(rinit))
  orats <- rinit[-curplay] 
  odevs <- dinit[-curplay]^2
  ongames <- ngames[-curplay]
  onwin <- nwin[-curplay]
  ondraw <- ndraw[-curplay]
  onloss <- nloss[-curplay]
  olag <- nlag[-curplay]
  olag[ongames != 0] <- olag[ongames != 0] + nm
  crats <- rinit[curplay] 
  cdevs <- dinit[curplay]^2
  ngames <- ngames[curplay] 
  nwin <- nwin[curplay]
  ndraw <- ndraw[curplay]
  nloss <- nloss[curplay]
  nlag <- nlag[curplay]

  qv <- log(10)/400; qv2 <- qv^2; qip3 <- 3*(qv/pi)^2 
  gammas <- split(gammas, x$Month)
  x <- split(x, x$Month)
  players <- lapply(x, function(y) unique(c(y$White, y$Black)))
  if(history) {
    histry <- array(NA, dim=c(np,nm,4), dimnames=list(play,1:nm,c("Rating","Deviation","Games","Lag")))
  }

  for(i in 1:nm) {
    traini <- x[[i]]
	gammai <- gammas[[i]]
    nr <- nrow(traini)
    playi <- players[[i]]

    cdevs[playi] <- pmin(cdevs[playi] + (nlag[playi]+1)*(cval^2), 122500)
    gdevs <- 1/sqrt(1 + qip3*cdevs) 
    ngamesi <- tabulate(c(traini$White,traini$Black), np)
    dscore <- .C("glicko_c",
      as.integer(np), as.integer(nr), as.integer(traini$White-1), as.integer(traini$Black-1),
      as.double(traini$Score), as.double(crats), as.double(gdevs), as.double(gammai),
      dscore = double(2*np))$dscore
    dval <- dscore[(np+1):(2*np)]; dscore <- dscore[1:np]
    cdevs <- 1/(1/cdevs + dval)
    crats <- crats + cdevs * qv * dscore

    trainiplw <- c(traini$White[traini$Score==1],traini$Black[traini$Score==0])
    trainipld <- c(traini$White[traini$Score==0.5],traini$Black[traini$Score==0.5])
    trainipll <- c(traini$White[traini$Score==0],traini$Black[traini$Score==1])
    ngames <- ngames + ngamesi
    nwin <- nwin + tabulate(trainiplw, np)
    ndraw <- ndraw + tabulate(trainipld, np)
    nloss <- nloss + tabulate(trainipll, np)
    nlag[ngames!=0] <- nlag[ngames!=0] + 1
    nlag[playi] <- 0

    if(history) {
      histry[,i,1] <- crats
      histry[,i,2] <- sqrt(cdevs)
      histry[,i,3] <- ngames
      histry[,i,4] <- nlag
    }
  }
  if(!history) histry <- NULL
  player <- suppressWarnings(as.numeric(names(c(crats,orats))))
  if (any(is.na(player))) player <- names(c(crats,orats))
  dfout <- data.frame(Player=player, Rating=c(crats,orats), Deviation=sqrt(c(cdevs,odevs)), 
    Games=c(ngames,ongames), Win=c(nwin,onwin), Draw=c(ndraw,ondraw), Loss=c(nloss,onloss), 
    Lag=c(nlag,olag),
	stringsAsFactors = FALSE)
  if(sort) dfout <- dfout[order(dfout$Rating,decreasing=TRUE),] else dfout <- dfout[order(dfout$Player),]
  row.names(dfout) <- 1:nrow(dfout)

  lout <- list(ratings=dfout, history=histry, gamma=gamma, cval=cval, type = "Glicko")
  class(lout) <- "rating"
  lout
}

"steph" <- function(x, status=NULL, init=c(2200,300), gamma=0, cval=10, hval=10, 
  bval=0, lambda = 2, history=FALSE, sort=TRUE, ...)
{
  if(length(init) != 2) stop("the length of 'init' must be two")
  if(ncol(x) != 4) stop("'x' must have four variables")
  if(nrow(x) == 0) {
    if(is.null(status)) stop("'x' is empty and 'status' is NULL")
    lout <- list(ratings=status, history=NULL, gamma=gamma, cval=cval, hval=hval, 
    bval=bval, lambda=lambda, type = "Stephenson")
    class(lout) <- "rating"
    return(lout)
  }
  gammas <- rep(gamma, length.out = nrow(x)) 
  names(x) <- c("Month","White","Black","Score")
  if(!is.numeric(x$Month)) 
    stop("Time period must be numeric")
  if(!is.numeric(x$White) && !is.character(x$White))
    stop("Player identifiers must be numeric or character")
  if(!is.numeric(x$Black) && !is.character(x$Black))
    stop("Player identifiers must be numeric or character")	
  if(!is.numeric(x$Score) || any(x$Score > 1) || any(x$Score < 0))
    stop("Game scores must be in the interval [0,1]")
	
  play <- sort(unique(c(x$White,x$Black)))
  np <- length(play)
  x$White <- match(x$White, play)
  x$Black <- match(x$Black, play)

  if(!is.null(status)) {
    npadd <- play[!(play %in% status$Player)]
	zv <- rep(0, length(npadd))
    npstatus <- data.frame(Player = npadd, Rating = rep(init[1],length(npadd)), 
      Deviation = rep(init[2],length(npadd)), Games = zv, Win = zv, Draw = zv, 
	  Loss = zv, Lag = zv)
	if(!("Games" %in% names(status))) status <- cbind(status, Games = 0)
	if(!("Win" %in% names(status))) status <- cbind(status, Win = 0)
	if(!("Draw" %in% names(status))) status <- cbind(status, Draw = 0)
	if(!("Loss" %in% names(status))) status <- cbind(status, Loss = 0)
    if(!("Lag" %in% names(status))) status <- cbind(status, Lag = 0)
    status <- rbind(status[,c("Player","Rating","Deviation","Games","Win","Draw","Loss","Lag")], npstatus)
    rinit <- status[,2]
    dinit <- status[,3]
    ngames <- status[,4]
    nwin <- status[,5]
    ndraw <- status[,6]
    nloss <- status[,7]
    nlag <- status[,8]
    names(rinit) <- names(dinit) <- names(ngames) <- status$Player
  }
  else {
    rinit <- rep(init[1], length.out=np)
    dinit <- rep(init[2], length.out=np)
    ngames <- nwin <- ndraw <- nloss <- rep(0, length.out=np)
    nlag <- rep(0,np)
    names(rinit) <- names(dinit) <- names(ngames) <- names(nlag) <- play
  }

  if(!all(names(rinit) == names(ngames)))
    stop("names of ratings and deviations are different")
  if(!all(names(rinit) == names(ngames)))
    stop("names of ratings and ngames are different")
  if(!all(play %in% names(rinit))) 
    stop("Payers in data are not within current status")

  nm <- length(unique(x$Month))
  curplay <- match(play, names(rinit))
  orats <- rinit[-curplay] 
  odevs <- dinit[-curplay]^2
  ongames <- ngames[-curplay]
  onwin <- nwin[-curplay]
  ondraw <- ndraw[-curplay]
  onloss <- nloss[-curplay]
  olag <- nlag[-curplay]
  olag[ongames != 0] <- olag[ongames != 0] + nm
  crats <- rinit[curplay] 
  cdevs <- dinit[curplay]^2
  ngames <- ngames[curplay] 
  nwin <- nwin[curplay]
  ndraw <- ndraw[curplay]
  nloss <- nloss[curplay]
  nlag <- nlag[curplay]

  qv <- log(10)/400; qv2 <- qv^2; qip3 <- 3*(qv/pi)^2
  gammas <- split(gammas, x$Month)  
  x <- split(x, x$Month)
  players <- lapply(x, function(y) unique(c(y$White, y$Black)))
  if(history) {
    histry <- array(NA, dim=c(np,nm,4), dimnames=list(play,1:nm,c("Rating","Deviation","Games","Lag")))
  }
    
  for(i in 1:nm) {
    traini <- x[[i]]
	gammai <- gammas[[i]]
    nr <- nrow(traini)
    playi <- players[[i]]

    cdevs[playi] <- pmin(cdevs[playi] + (nlag[playi]+1)*(cval^2), 122500)
    gdevs <- 1/sqrt(1 + qip3*cdevs) 
    ngamesi <- tabulate(c(traini$White,traini$Black), np)
    
    dscore <- .C("stephenson_c",
      as.integer(np), as.integer(nr), as.integer(traini$White-1), as.integer(traini$Black-1),
      as.double(traini$Score), as.double(crats), as.double(gdevs), as.double(gammai),
      as.double(bval/100), dscore = double(3*np))$dscore
  
    l1t <- dscore[(2*np+1):(3*np)]
    dval <- dscore[(np+1):(2*np)]
    dscore <- dscore[1:np]
    
    cdevs <- 1/(1/(cdevs + ngamesi*(hval^2)) + dval)
    crats <- crats + cdevs * qv * dscore 
    crats[playi] <- crats[playi] + (lambda/100)*l1t[playi]/(ngamesi[playi])

    trainiplw <- c(traini$White[traini$Score==1],traini$Black[traini$Score==0])
    trainipld <- c(traini$White[traini$Score==0.5],traini$Black[traini$Score==0.5])
    trainipll <- c(traini$White[traini$Score==0],traini$Black[traini$Score==1])
    ngames <- ngames + ngamesi
    nwin <- nwin + tabulate(trainiplw, np)
    ndraw <- ndraw + tabulate(trainipld, np)
    nloss <- nloss + tabulate(trainipll, np)
    nlag[ngames!=0] <- nlag[ngames!=0] + 1
    nlag[playi] <- 0

    if(history) {
      histry[,i,1] <- crats
      histry[,i,2] <- sqrt(cdevs)
      histry[,i,3] <- ngames
      histry[,i,4] <- nlag
    }
  }
  if(!history) histry <- NULL
  player <- suppressWarnings(as.numeric(names(c(crats,orats))))
  if (any(is.na(player))) player <- names(c(crats,orats))
  dfout <- data.frame(Player=player, Rating=c(crats,orats), Deviation=sqrt(c(cdevs,odevs)), 
    Games=c(ngames,ongames), Win=c(nwin,onwin), Draw=c(ndraw,ondraw), Loss=c(nloss,onloss), 
    Lag=c(nlag,olag),
	stringsAsFactors = FALSE)
  if(sort) dfout <- dfout[order(dfout$Rating,decreasing=TRUE),] else dfout <- dfout[order(dfout$Player),]
  row.names(dfout) <- 1:nrow(dfout)

  lout <- list(ratings=dfout, history=histry, gamma=gamma, cval=cval, hval=hval, 
    bval=bval, lambda=lambda, type = "Stephenson")
  class(lout) <- "rating"
  lout
}

"metrics" <- function(act, pred, cap = c(0.01,0.99), which = 1:3, na.rm = TRUE,
  sort = TRUE, digits = 3, scale = TRUE)
{
  if(!is.numeric(pred)) stop("'pred' must be numeric")
  if(!is.numeric(act)) stop("'act' must be numeric")
  pred <- as.matrix(pred)
  np <- ncol(pred); nr <- nrow(pred)
  mets <- matrix(NA, ncol=3, nrow=np, 
    dimnames=list(colnames(pred),c("bdev","mse","mae")))
  for(i in 1:np) {
    predc <- pmax.int(pmin.int(pred[,i], cap[2]), cap[1])
    mets[i,1] <- -mean(act*log(predc) + (1-act)*log(1-predc), na.rm = na.rm)
    if(scale) mets[i,1] <- mets[i,1]/(-mean(act*log(0.5) + (1-act)*log(0.5), na.rm = na.rm))
    mets[i,2] <- sqrt(mean((pred[,i]-act)^2, na.rm = na.rm))
    if(scale) mets[i,2] <- mets[i,2]/sqrt(mean((0.5-act)^2, na.rm = na.rm))
    mets[i,3] <- mean(abs(pred[,i]-act), na.rm = na.rm)
    if(scale) mets[i,3] <- mets[i,3]/mean(abs(0.5-act), na.rm = na.rm)
  }
  mets <- 100*mets[,which]
  if(sort && is.matrix(mets)) mets <- mets[order(mets[,1]),]
  round(drop(mets), digits)
}

"kfide" <- function(rating, games, elite = NULL, kv = c(10,15,30)) 
{
  if(any(is.na(rating))) stop("missing values in 'ratings' vector")
  if(any(is.na(games))) stop("missing values in 'games' vector")
  if(length(rating) != length(games)) 
    stop("lengths of 'ratings' and 'games' must be the same")
  kfac <- rep(NA, length(rating))
  if(is.null(elite)) elite <- (rating >= 2400) else elite <- as.logical(elite)
  kfac[!elite & games < 30] <- kv[3]
  kfac[!elite & games >= 30] <- kv[2]
  kfac[elite] <- kv[1]
  if(any(is.na(kfac))) stop("missing values in K factor")
  kfac
}

"krating" <- function(rating, games, elite = NULL, rv = 2300, kv = c(32,26)) 
{
  if(any(is.na(rating))) stop("missing values in 'ratings' vector")
  if(any(is.na(games))) stop("missing values in 'games' vector")
  if(length(rating) != length(games)) 
    stop("lengths of 'ratings' and 'games' must be the same")
  if(length(rv) != (length(kv)-1)) 
    stop("length of 'kv' must be one more than 'gv'")

  rv <- c(-Inf, rv, Inf)
  rind <- as.numeric(cut(rating, rv))
  kfac <- kv[rind]
  if(any(is.na(kfac))) stop("missing values in K factor")
  kfac
}

"kgames" <- function(rating, games, elite = NULL, gv = 30, kv = c(32,26)) 
{
  if(any(is.na(rating))) stop("missing values in 'ratings' vector")
  if(any(is.na(games))) stop("missing values in 'games' vector")
  if(length(rating) != length(games)) 
    stop("lengths of 'ratings' and 'games' must be the same")
  if(length(gv) != (length(kv)-1)) 
    stop("length of 'kv' must be one more than 'gv'")

  gv <- c(-Inf, gv, Inf)
  gind <- as.numeric(cut(games, gv))
  kfac <- kv[gind]
  if(any(is.na(kfac))) stop("missing values in K factor")
  kfac
}

"print.rating" <-  function(x, digits = 0, ...) 
{
    rdf <- x$ratings
    rdf$Rating <- round(rdf$Rating, digits)
    if(x$type == "Glicko" || x$type == "Stephenson") 
	  rdf$Deviation <- round(rdf$Deviation, digits+2)
    np <- nrow(rdf); ng <- round(sum(rdf$Games)/2)
    cat(paste("\n",x$type," Ratings For ",np," Players Playing ",ng," Games\n\n", sep=""))
    print(rdf[1:(min(1000,np)),])
    if(np > 1000) cat("\nOutput Tructated To First 1000 Players \n")
    cat("\n")
    invisible(0)
}

"summary.rating" <-  function(object, ...) 
{
  obj <- object$ratings
  obj$Games <- factor(obj$Games)
  obj$Lag <- factor(obj$Lag)
  summary(obj)
}

"predict.rating" <- function(object, newdata, tng=15, trat=NULL, gamma=30, 
  thresh, ...)
{
  if(nrow(newdata) == 0) stop("'newdata' must have non-zero rows")
  obj <- object$ratings
  wmat <- match(newdata[,2], obj$Player)
  bmat <- match(newdata[,3], obj$Player)

  if(!is.null(trat)) obj$Rating[obj$Games < tng] <- trat[1]
  else is.na(obj$Rating[obj$Games < tng]) <- TRUE
  if(!is.null(trat)) obj$Deviation[obj$Games < tng] <- trat[2]
  else is.na(obj$Deviation[obj$Games < tng]) <- TRUE
  wrat <- obj$Rating[wmat]; brat <- obj$Rating[bmat]
  if(!is.null(trat)) wrat[is.na(wrat)] <- trat[1]
  if(!is.null(trat)) brat[is.na(brat)] <- trat[1] 

  qv <- log(10)/400; qip3 <- 3*(qv/pi)^2
  if(object$type == "Glicko" || object$type == "Stephenson") {
    if(!is.null(trat) && length(trat) != 2) 
	  stop("'trat' must be vector of length two")
    wdev <- obj$Deviation[wmat]; bdev <- obj$Deviation[bmat]
    if(!is.null(trat)) wdev[is.na(wdev)] <- trat[2]
    if(!is.null(trat)) bdev[is.na(bdev)] <- trat[2]
  }
  if(object$type == "Elo")
    preds <- 1/(1+10^((brat-wrat-gamma)/400))
  if(object$type == "Glicko" || object$type == "Stephenson") {
    vec <- 1/sqrt(1 + qip3*(wdev^2 + bdev^2))
    preds <- 1/(1+10^(vec * (brat-wrat-gamma)/400))
  }
  if(!missing(thresh)) preds <- as.numeric(preds >= thresh)
  preds
}

"hist.rating" <- function(x, which = "Rating", tng=15, history = FALSE, log = FALSE, 
  xlab = which, main = paste(x$type," Ratings System"), density=FALSE, add=FALSE, ...)
{
   if(!history) {
     obj <- x$ratings
     obj <- obj[obj$Games >= tng,]
     obj <- obj[[which]]
     if(log) obj <- log(obj+1) 
     if(density) {
       if(add) lines(density(obj), xlab=xlab, main=main, ...)
       else plot(density(obj), xlab=xlab, main=main, ...)
     } else hist(obj, xlab=xlab, main=main, ...)
   } else {
     if(is.null(x$history)) stop("Need Full History For Plotting")
     obj <- x$history[,,which]
     ngm <- x$history[,,"Games"]
     nt <- ncol(obj)
     old <- par(ask = TRUE)
     for(i in 1:nt) {
       ngi <- (ngm[,i] >= tng)
	   if(all(!ngi)) next
       if(density) {
         if(add) lines(density(obj[ngi,i]), xlab=xlab, main=main, ...)
         else plot(density(obj[ngi,i]), xlab=xlab, main=main, ...)
       } else hist(obj[ngi,i], xlab=xlab, main=main, ...)
     }
     par(old)
   }
   invisible(obj)
}

"plot.rating" <- function(x, which = "Rating", players = NULL, t0 = 1, tv = NULL, 
  npl = 10, random = FALSE, xlab = "Time Period", ylab = paste(x$type," Ratings"), 
  main = paste(x$type," Ratings System"), inflation = FALSE, add=FALSE, ...) 
{
  if(is.null(x$history)) stop("Need Full History For Plotting")
  dmh <- dim(x$history)
  np <- dmh[1]
  if(length(t0) == 2) {
    nt <- t0[2]
    t0 <- t0[1] 
  } else nt <- dmh[2]
  if(nt > dmh[2] || nt==1) stop("Not enough history available")
  obj <- x$history[,t0:nt,which]
  ngm <- x$history[,t0:nt,"Games"]
  if(is.null(tv)) tv <- t0:nt 
  if(inflation == TRUE) {
    obj <- x$history[,t0:nt,"Rating"]
    objG <- x$history[,t0:nt,"Games"]
    objL <- x$history[,t0:nt,"Lag"]
    is.na(obj) <- (objL > 11 | objG < 25)
    obj <- apply(obj, 2, function(x) mean(sort(x, decreasing = TRUE)[1:npl]))
    if(!add) plot(tv, obj, type="l", xlab=xlab, ylab=ylab, main=main, ...)
    else lines(tv, obj, ...)
    return(invisible(0))
  }
  if(!is.null(players)) {
    obj <- t(as.matrix(obj[as.character(players),]))
    matplot(tv, obj, type="l", xlab=xlab, ylab=ylab, main=main, add=add, ...)
  } else {
    if(!random) players <- order(ngm[,1],decreasing=TRUE)[1:npl]
    else players <- sample(1:npl, 10)
    obj <- t(as.matrix(obj[players,]))
    matplot(tv, obj, type="l", xlab=xlab, ylab=ylab, main=main, add=add, ...)
  }
  invisible(0)
}


