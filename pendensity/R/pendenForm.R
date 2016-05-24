pendenForm <- function(penden.env) {
  pf <- get("frame",penden.env)
  char.vec <- as.character(get("form",penden.env))
  response.name <- char.vec[2]
  response.val <- eval(parse(text=response.name),envir=pf)
  val <- paste(string.help(char.vec[3]),collapse= "")
  val <- string.help(val, "+")
  no.cov <- list()
  parcov <- list()
  nonparcov <- list()
  no.cov$names <- NULL
  no.cov$x <- NULL
  parcov$names <- NULL
  parcov$x <- NULL
  parcov$contrasts <- NULL
  parcov$length.how <- 0
  parcov$combi <- c()
  parcov$combi.how <- 1
  nonparcov$names <- NULL
  nonparcov$x <- NULL
  nonparcov$knots <- NULL
  nonparcov$num.knots <- NULL
  
  if (val[1] == "1") {
    no.cov <- 1
    parcov <- NULL
    nonparcov <- NULL
    parcov$length.how <- 1
    parcov$combi.how <- 1
  }
  if (val[1] != "1") {
    col.nam <- c()
    for (i in 1:length(val)) {
      term <- val[i]
      parcov$name <- c(parcov$name, term)
      if (is.factor(eval(parse(text = term),envir=pf))) {#penden.env
        parcov$x <- cbind(parcov$x, eval(parse(text = term),envir=pf))#penden.env
        parcov$contrasts[[i]] <- contrasts(eval(parse(text = term),envir=pf))#penden.env
        parcov$levels[[i]] <- levels(eval(parse(text = term),envir=pf))#penden.env
        parcov$length.how <- parcov$length.how+1
        parcov$len.cov[i] <- length(parcov$levels[[i]])
        parcov$combi.how <- parcov$combi.how * parcov$len.cov[i]
      }
      else {
        stop("Covariates have to be factors!",call.=FALSE)
      }
      col.nam <- c(col.nam,paste(parcov$name[i],"=",parcov$levels[[i]],sep=""))
    }
    help <- c()
    help2 <- c()
    for(j in 1:parcov$length.how) {
    help <- seq(1,parcov$len.cov[j],by=1)
    help2 <- c(help2,help)
  }
    ind <- which(help2==1)
    col.nam <- col.nam[-ind]
    
    parcov$x.mat <- matrix(0,length(parcov$x[,1]),length(col.nam))
    
    colnames(parcov$x.mat) <- col.nam
    N <- matrix(1:length(parcov$x[,1]))
    k <- 1
    for(j in 2:parcov$len.cov[1]) {
    parcov$x.mat[,k] <- apply(N,1,function(i,levels,level,x) if(levels[x[i]]==as.numeric(level)) return(1) else return(0),parcov$levels[[1]],parcov$levels[[1]][j],parcov$x[,1])
    k <- k+1
  }   
    if(parcov$length.how>1) {
    for(i in 2:parcov$length.how) {
      for(j in 2:parcov$len.cov[i]) {
        parcov$x.mat[,k] <- apply(N,1,function(i,levels,level,x) if(levels[x[i]]==as.numeric(level)) return(1) else return(0),parcov$levels[[i]],parcov$levels[[i]][j],parcov$x[,i])
        k <- k+1
      }
    }
  }
    cons <- (parcov$contrasts[[1]])
    if(parcov$length.how > 1) {
      for(j in 2:parcov$length.how) {
        help <- length(cons[,1])
        help.j <- parcov$len.cov[j]
        help.mat1 <- kronecker(matrix(1,help.j,1),cons)
        help.mat2 <- kronecker(diag(1,help.j),matrix(1,help,1))[,-1]
        cons <- cbind(help.mat1,help.mat2)
      }
    }
    parcov$cons <- cons
  }
  info <- list(formula = get("form",penden.env), y = response.val, no.cov=no.cov, parcov=parcov)
  return(info)
}
