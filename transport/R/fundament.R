wasserstein <- function(a, b, p=1, tplan=NULL, prob=TRUE, ...) {
  stopifnot(compatible(a,b))
  if (a$dimension < 2) stop("pixel grids of dimension >= 2 required")
  
  if (is.null(tplan)) {
  	tplan <- transport(a,b,p=p,...)
  }
  
  K <- dim(tplan)[1]

  # computes (sum(m * x^pp))^(1/pp)
  wpsum <- function(x, m=rep(1,length(x)), pp) {
  	if (pp == 1) {
      mmax <- max(m)
      xmax <- max(abs(x))
      return(ifelse(mmax == 0, 0, mmax * xmax * sum((m/mmax)*(x/xmax))))
  	} else if (length(unique(m)) == 1 && unique(m) == 1) {
  	  xmax <- max(abs(x))
  	  return(ifelse(xmax == 0, 0, xmax * sum((x/xmax)^pp)^(1/pp)))
  	} else {
      mmax <- max(m)
      xmax <- max(abs(x))
      return(ifelse(mmax == 0 || xmax == 0, 0, (mmax)^(1/pp) * xmax * sum((m/mmax)*(x/xmax)^pp)^(1/pp)))
    }
  }
  
  if (class(a) == "pgrid" && class(b) == "pgrid") {
    gg <- expand.grid(a$generator)
    orig <- gg[tplan$from,]
    dest <- gg[tplan$to,]
  } else if (class(a) == "pp" && class(b) == "pp") {
  	orig <- a$coordinates[tplan$from,]
    dest <- b$coordinates[tplan$to,]
  } else {
  	stop("a and b must be either both of class 'pgrid' or both of class 'pp'")
  }
  
  dd <- apply(orig-dest,1,wpsum,pp=2)
  res <- wpsum(dd, tplan$mass, p)

  if (prob) {
  	if (class(a) == "pgrid") {
      res <- res/(a$totmass^(1/p))
      # note: a$totmass might be larger than sum(tplan$mass)    	
      # due to static mass
    } else {
      res <- res/(a$N^(1/p))	
    }
  }

  return(res)
}


transport <- function(a, b, ...) {
  stopifnot(class(a) == class(b))	
  UseMethod("transport")
}

trcontrol <- function(method = c("shortsimplex", "revsimplex", "primaldual", "aha", "auction", "auctionbf"),
  para=list(), start = c("auto", "modrowmin", "nwcorner", "russell"), nscales = 1, scmult = 2, returncoarse = FALSE,
  a=NULL, b=NULL, M=NULL, N=NULL) {
# a,b,M,N serve for computing parameters or start solutions automatically
# by default M=a$N, N=b$N, which are overridden if M and/or N are specified
  method <- match.arg(method)
  start <- match.arg(start)
  if (is.null(M) && !is.null(a)) {
  	M <- ifelse(class(a) %in% c("pgrid","pp"), a$N, length(a))
  }
  if (is.null(N) && !is.null(b)) {
  	N <- ifelse(class(b) %in% c("pgrid","pp"), b$N, length(b))
  }
  if (is.null(M) && !is.null(N)) {M <- N}
  if (is.null(N) && !is.null(M)) {N <- M}
  if (is.null(N) && method %in% c("shortsimplex","auction")) {
  	stop("At least M or N must be specified for method ", method)
  }
    
  if (method == "aha") {

    if (is.numeric(para) && length(para) == 2) {
      para = list(factr=para[1], maxit=para[2])
      message("Parameter vector interpreted by trcontrol as: \n factr = ",para[1],"; maxit = ",para[2],".")
    }
  	
  	propnames <- c("factr","maxit")
    done <- pmatch(names(para), propnames)
    if (!is.list(para) || (length(para) > 0 && length(done) == 0) || any(is.na(done))) {
      stop("expected for 'para' with method='aha' either an empty list or a named list with 1-2 components out of 'factr', 'maxit' (after partial matching)")
    }
  
    newpara <- list(factr=0, maxit=0)
    
    if (1 %in% done) {
      newpara$factr <- para[[which(done == 1)]]
      if (newpara$factr < 1) {
      	stop("parameter 'factr' for aha algorithm has to be >= 1 (and is typically >= 1e3)")
      }
    } else {
      newpara$factr <- 1e5 
      
          # /400 since after reduction about half of the N producers disappear
      #newpara$slength <- min(b$N, newpara$slength)
          # Can't choose the shortlist longer then the number of consumers
          # This is caught now in transport functions
    }
    
    if (2 %in% done) {
      newpara$maxit <- para[[which(done == 2)]]
      if (newpara$maxit < 1000) {
      	warning("parameter 'maxit' for aha algorithm < 1000 is not recommended")
      }
    } else {
      newpara$maxit <- 3000
    }
  
    para <- newpara
  }
  # fi (method == "shortsimplex")


  if (method == "shortsimplex") {

    if (is.numeric(para) && length(para) == 3) {
      para = list(slength=para[1], kfound=para[2], psearched=para[3])
      message("Parameter vector interpreted by trcontrol as: \n slength = ",para[1],"; kfound = ",para[2],"; psearched = ",para[3],".")
    }

  	propnames <- c("slength","kfound","psearched")
    done <- pmatch(names(para), propnames)
    if (!is.list(para) || (length(para) > 0 && length(done) == 0) || any(is.na(done))) {
      stop("expected for 'para' with method='shortlist' either an empty list or a named list with 1-3 components out of 'slength', 'cfound', or 'psearched' (after partial matching)")
    }
 
    newpara <- list(slength=0, kfound=0, psearched=0)

    if (1 %in% done) {
      newpara$slength <- para[[which(done == 1)]]
#      if (newpara$slength > b$N) {
#      	warning("parameter 'slength' for shortlist algorithm was larger then no. of target points... Fixed.")
#      	newpara$slength <- b$N
      if (newpara$slength < 1) {
      	stop("parameter 'slength' for shortlist algorithm has to be >= 1")
      }
    } else {
      newpara$slength <- min(N,15 + max(0, floor(15 * log(N/400)/log(2)))) 
      
          # /400 since after reduction about half of the N producers disappear
      #newpara$slength <- min(b$N, newpara$slength)
          # Can't choose the shortlist longer then the number of consumers
          # This is caught now in transport functions
    }
    
    if (2 %in% done) {
      newpara$kfound <- para[[which(done == 2)]]
      if (newpara$kfound < 1) {
      	stop("parameter 'kfound' for shortlist algorithm has to be >= 1")
      }
    } else {
      newpara$kfound <- newpara$slength
    }

    if (3 %in% done) {
      newpara$psearched <- para[[which(done == 3)]]
      if (newpara$psearched > 1 || newpara$psearched <= 0) {
      	stop("parameter 'psearched' for shortlist algorithm has to be > 0 and <= 1")
      }
    } else {
      newpara$psearched <- 0.05
    }
    # Note that at least one row is searched anyway (even if psearch was 0 or negative)
    
    para <- newpara
  }
  # fi (method == "shortsimplex")


  if (method == "auction" || method == "auctionbf") {

  # lasteps Wert muss < 1/n sein fuer exaktes Ergebnis (an gewissen Stellen in altem Code steht = 1/n, aber ich sehe
  # keinen guten Grund dafuer)
  # lasteps und epsfac sind nur fuer auction und auctionbf relevant
  # epsfac = NA bedeutet kein eps-scaling verwenden, habe experimentiert mit epsfac = 10 und = 80 
  # Bertsekas empfiehlt 4-10, ist aber nach meiner Erfahrung fuer unsere Probleme zu klein
  # ich faende einen eps-power natuerlicher als einen epsfac, aber das mit dem epsfac habe ich von Bertsekas

    if (is.numeric(para) && length(para) == 2) {
      para = list(lasteps=para[1], epsfac=para[2])
      message("Parameter vector interpreted by trcontrol as: \n lasteps = ",para[1],"; epsfac = ",para[2],".")
    }
  	
  	propnames <- c("lasteps","epsfac")
    done <- pmatch(names(para), propnames)
    if (!is.list(para) || (length(para) > 0 && length(done) == 0) || any(is.na(done))) {
      stop("expected for 'para' with method='auction'/'auctionbf' either an empty list or a named list with 1-2 components out of 'lasteps', 'epsfac' (after partial matching)")      
    }
 
    newpara <- list(lasteps=0, epsfac=0)

    if (1 %in% done) {
      newpara$lasteps <- para[[which(done == 1)]]
#      if (newpara$lasteps > 1) {
#      	stop("parameter 'slength' for shortlist algorithm has to be >= 1")
#      }
    } else {
      stopifnot(M==N)
      if (is.null(N)) { stop("Not enough information available to determine lasteps.") }
      newpara$lasteps <- 1/(N+1)
    }
    
    if (2 %in% done) {
      newpara$epsfac <- para[[which(done == 2)]]
#      if (newpara$kfound <= 1) {
#      	stop("parameter 'kfound' for shortlist algorithm has to be >= 1")
#      }
    } else {
      newpara$epsfac <- 10
    }
  	
    para <- newpara
  }
  # fi (method == "auction" || method == "auctionbf")

  res <- list(method=method, para=para, start=start, nscales=nscales, scmult=scmult, returncoarse=returncoarse)
  class(res) <- "trc"

  return(res)
  
}



transport.pgrid <- function(a, b, p = NULL, method = c("revsimplex", "shortsimplex", "primaldual", "aha"), control = list(), ...) {
  # returncoarse: gibt im Falle von nscales >= 2 an, ob groebere Probleme und deren Loesungen auch ausgegeben werden sollen;
  # Anderenfalls wird nur die feinste Loesung (ohne Problem) zurueckgegeben

  # Check inputs
  # ======================================================================                             	
  stopifnot(class(a) == "pgrid" && class(b) == "pgrid")
  stopifnot(compatible(a,b))
  if (a$dimension < 2) stop("pixel grids of dimension >= 2 required")
#  if (a$dimension > 2) warning("transport.pgrid for pixel grids of dimension > 2 is still somewhat experimental")
  if (!(a$structure %in% c("square", "rectangular")))
    stop("transport.pgrid is currently only implemented for rectangular pixel grids")
  ngrid <- a$n
  Ngrid <- a$N
  
  method <- match.arg(method)
    
  if (is.null(p))
  	p <- ifelse(method == "aha",2,1)
  if (method == "aha" && p != 2)
    stop("method 'aha' works only for p = 2")
  if (method == "aha" && a$dimension != 2)
    stop("method 'aha' works only in two dimensions")    
  if (p < 1) {
  	stop("p has to be >= 1")
  }  

  if (class(control) != "trc") {
  	control$method <- method
  	control$a <- a
  	control$b <- b
  	control = do.call(trcontrol, control)
  }
  
  start <- control$start
  nscales <- control$nscales
  scmult <- control$scmult
  returncoarse <- control$returncoarse
  is.natural <-
    function(x, tol = .Machine$double.eps^0.5)  all((abs(x - round(x)) < tol) & x > 0.5)
  # aus der Hilfe zu is.integer, checks for a vector whether all entries are approximately natural numbers
  if (nscales != 1 && (length(unique(ngrid)) != 1 || length(ngrid) != 2)) {
  	stop("multiscale approach is currently only implemented for quadratic grids of dimension 2")
  }    
  if (nscales != 1 && !(p == 1 && method=="revsimplex") && !(p == 2 && method=="aha")) {
  	stop('the multiscale approach is currently only implemented for p = 1 and method = "revsimplex" or p = 2 and method = "aha".
  	Note that for p != 1 and method = "revsimplex" the default starting solution is computed by a 1-step multiscale approach')
  }      
  if (nscales != 1 && !(is.natural(scmult) && is.natural(nscales) && is.natural(ngrid/scmult^(nscales-1)))) {
  	stop("grid of size ", ngrid, " cannot be scaled down ", nscales, " times by a factor of ", scmult)
  }
  
  gg <- expand.grid(a$generator)
  dd <- as.matrix(dist(gg))^p
  
# we exclude aha, because currently only grid approach is implemented
  if (method != "aha") {
  	# if costs are based on metric, 
  	# moving mass that is needed at a site never pays off
  	if (p == 1) {
      minab <- pmin(a$mass,b$mass)
      ared <- a$mass - minab  
      bred <- b$mass - minab
      wha <- ared > 0
      whb <- bred > 0
      apos <- ared[wha]
      bpos <- bred[whb]
      dd <- dd[wha,whb]
    } else {
      ared <- a$mass
      bred <- b$mass
      wha <- rep(TRUE,Ngrid)
      whb <- rep(TRUE,Ngrid)
      apos <- ared
      bpos <- bred	
    }
    m <- length(apos)
    n <- length(bpos)
    asum <- sum(apos)
    bsum <- sum(bpos)
    if (!isTRUE(all.equal(asum,bsum))) {
      warning("total mass in a and b differs. Normalizing a and b to probability measures.")
  	  apos <- apos/asum
  	  bpos <- bpos/bsum
  	  asum <- bsum <- 1
    }
    fudgeN <- fudgesum <- 1 
    if (!is.natural(apos) || !is.natural(bpos)) {
  	  fudgeN <- 1e9
  	  fudgesum <- asum
  	  apos <- round(apos/asum * fudgeN)
  	  bpos <- round(bpos/bsum * fudgeN)
  	  apos <- fudge(apos,fudgeN)
  	  bpos <- fudge(bpos,fudgeN)
    }
  } 
    
  # Interesting part starts here
  # ======================================================================
  if (method == "primaldual") {
    precision=9
    dd <- dd/max(dd)
    dd <- dd*(10^precision)
    #cat(as.integer(ddpos), as.integer(apos), as.integer(bpos), as.integer(m), as.integer(n),
             # flowmatrix = integer(m*n), sep="\n")
    #stop("primaldual soon")
  	ti <- proc.time()
    res1 <- .C("primaldual", as.integer(dd), as.integer(apos), as.integer(bpos), as.integer(m), as.integer(n),
              flowmatrix = integer(m*n), DUP=TRUE, PACKAGE="transport")
    print(proc.time()-ti)
    temp <- list(assignment=res1$flowmatrix, basis=as.numeric(res1$flowmatrix > 0))
     # make pretty output
    nbasis <- sum(temp$basis)
    res <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
    ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
    res$from <- which(wha)[ind[,1]]
    res$to <- which(whb)[ind[,2]]
    res$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]
    res$mass <- res$mass * fudgesum / fudgeN
    return(res)
  }

  if (method == "aha") {
    #cat(ddpos, apos, bpos, sep="\n")
    #stop("aha soon")
    # braucht noch Verfeinerung (wir moechten idealerweise auch nur apos, bpos verwenden, was direkt nicht geht;
    # die Verwendung von ared, bred verlangsamt extrem (bei 32x32 ca. 5 sek versus 12 sek), was ok ist, denke ich)
    ti <- proc.time()
  	res <- aha(a$mass,b$mass,nscales=1,scmult=2,maxit=control$para$maxit,factr=control$para$factr,wasser=FALSE,wasser.spt=NA)
  	print(proc.time()-ti)
  	return(res)
  	# Bring in right form
  }

  if (method == "revsimplex") {
  	
  	# Short-cut if nscale = 1
  	# (saves extremely little time)
  	# ----------------------------
    if (nscales == 1) {
      #
      # Compute starting solution
      if (start == "auto" && (p == 1 || min(ngrid) < 8 || max(ngrid) < 16 || scmult == 1 || a$dimension > 2)) {
      	start <- "modrowmin"
      } else if (start == "auto" && !is.natural(ngrid/scmult)) {
      	warning('grid dimensions are not divisible by ", scmult, "for multiscale starting solution. Using uniscale starting solution instead.
      	Note that the computation *might* be much more efficient when setting control = list(scmult=k), where k is a small common divisor
      	of the numbers of pixels in each direction')
      	start <- "modrowmin"
      }
  	  if (start == "russell") {
        temp <- russell(apos,bpos,dd)
        initassig <- temp$assignment
        initbasis <- temp$basis
        startgiven <- 1
      } else if (start == "nwcorner") {
        temp <- northwestcorner(apos,bpos)
        initassig <- temp$assignment
        initbasis <- temp$basis	
        startgiven <- 1
      } else if (start == "modrowmin"){   # modrowmin in C-Code
      	initassig <- rep(0L,m*n)
      	initbasis <- rep(0L,m*n)
      	startgiven <- 0
      } else {     # scalestart (is usually the best choice if p!=1 except for problems with very local optimal solutions)
      	           # for p==1 it is suboptimal only, because the current implementation does not take reduction of problem
      	           # due to static mass into account; use nscales = 2, scmult >= 2 for (essentially) same behaviour with p==1
      	           # Note 1: scmult also works for scalestart
      	           # Note 2: nscales *only* works for p=1 because it always removes static mass.
        temp <- scalestart(ared,bred,a$generator,ngrid[1],ngrid[2],p=p,scmult=scmult)
        initassig <- temp$assignment
        initbasis <- temp$basis	
        startgiven <- 1
      }
    ti <- proc.time()
    res <- .C("revsimplex", as.integer(m), as.integer(n), as.integer(apos), as.integer(bpos),
	          as.double(dd), assignment = as.integer(initassig), basis = as.integer(initbasis), startgiven = as.integer(startgiven),
	          DUP=TRUE, PACKAGE="transport")
	print(proc.time()-ti)           
	temp <- list(assignment=res$assignment, basis=res$basis)
    } else { 

      # if (a$dimension != 2) stop("Multi-scale approach only implemented for pixel grids of dimension 2")
      x <- a$generator[[1]]
      y <- a$generator[[2]]
  	  # Compute coarser problems
  	  # (note: coarsening includes subtracting pixelwise minimum on coarser grid!!)
  	  # ----------------------------
  	  ngridvec <- ngrid[1]/scmult^(0:(nscales-1))
  	  Ngridvec <- ngridvec^2
      problem <- vector("list", nscales)
      problem[[1]] <- list(ared=ared, bred=bred, x=x, y=y)
      # x, y are the grid sequences in x and y direction (currently always x=y)
      grmat <- matrix(unlist(lapply(as.list(1:ngridvec[2]), function(x) {rep(((x-1)*ngridvec[2]+1):(x*ngridvec[2]), 
      	              times = scmult, each = scmult)})), ngridvec[1], ngridvec[1])
      for (k in 2:nscales) {
      	problem[[k]] <- list(ared=numeric(0), bred=numeric(0)) 
      	problem[[k]]$ared <- 
      	  matrix( aggregate(as.vector(problem[[k-1]]$ared), by=list(as.vector(grmat[1:ngridvec[k-1],1:ngridvec[k-1]])),
      	    FUN="sum")[,2], ngridvec[k], ngridvec[k])
      	problem[[k]]$bred <- 
      	  matrix( aggregate(as.vector(problem[[k-1]]$bred), by=list(as.vector(grmat[1:ngridvec[k-1],1:ngridvec[k-1]])),
      	    FUN="sum")[,2], ngridvec[k], ngridvec[k])
      	  if (p == 1) {
       	    minabnew <- pmin(problem[[k]]$ared, problem[[k]]$bred)
            problem[[k]]$ared <- problem[[k]]$ared - minabnew  
            problem[[k]]$bred <- problem[[k]]$bred - minabnew
          }
          # deal with grid sequences
          problem[[k]]$x <- aggregate(problem[[k-1]]$x, by=list(as.vector(grmat[1,1:ngridvec[k-1]])), FUN="mean")[,2]
          problem[[k]]$y <- aggregate(problem[[k-1]]$y, by=list(as.vector(grmat[1:ngridvec[k-1],1])), FUN="mean")[,2]
      }
     
  	  # Solve them starting with the coarsest
  	  # ----------------------------
  	  if (returncoarse) { 	 
      	sol <- vector("list",nscales)
      }
  	  for (k in nscales:1) {
  	  	# inputs
  	  	# print(problem[[k]])
  	  	wha <- problem[[k]]$ared>0
  	    whb <- problem[[k]]$bred>0
        apos <- problem[[k]]$ared[wha]
        bpos <- problem[[k]]$bred[whb]
        gg <- expand.grid(problem[[k]]$x, problem[[k]]$y)
        dd <- as.matrix(dist(gg))   # of course everything here is also already p=1
        dd <- dd[wha,whb]
        mcoarse <- length(apos)
        ncoarse <- length(bpos)
        if (k == nscales) {
          # Compute starting solution
  	      if (start == "russell") {
            temp <- russell(apos,bpos,dd)
            initassig <- temp$assignment
            initbasis <- temp$basis
            startgiven <- 1
          } else if (start == "nwcorner") {
            temp <- northwestcorner(apos,bpos)
            initassig <- temp$assignment
            initbasis <- temp$basis	
            startgiven <- 1
          } else {   # modrowmin in C-Code
      	    initassig <- rep(0L,mcoarse*ncoarse)
      	    initbasis <- rep(0L,mcoarse*ncoarse)
      	    startgiven <- 0
          }
        } else {
            
          if (p == 1) {  # <================================ 15/10/12: for p != 1 the strict separation in producers
                                                # and consumers in refinesol is wrong and the removal of static mass above also
            newtemp <- refinesol(problem[[k+1]]$ared, problem[[k+1]]$bred, problem[[k]]$ared, problem[[k]]$bred,
                                 temp$assignment, temp$basis, mult=scmult)  
          } else {
          	stop("the multiscale approach for p != 1 is out of order.")
          }
          initassig <- newtemp$assig2
          initbasis <- newtemp$basis2 
          startgiven <- 1      
        }  
        
        ti <- proc.time()
        res <- .C("revsimplex", as.integer(mcoarse), as.integer(ncoarse), as.integer(apos), as.integer(bpos),
	              as.double(dd), assignment = as.integer(initassig), basis = as.integer(initbasis), startgiven = as.integer(startgiven),
	              DUP=TRUE, PACKAGE="transport")
	    print(proc.time()-ti)           
	    temp <- list(assignment=res$assignment, basis=res$basis)
        if (returncoarse) {
      	  nbasis <- sum(temp$basis)
          sol[[k]] <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
  	  	  ind <- which(matrix(as.logical(temp$basis),mcoarse,ncoarse), arr.ind=TRUE) 
          sol[[k]]$from <- which(wha)[ind[,1]]
          sol[[k]]$to <- which(whb)[ind[,2]]
          sol[[k]]$mass <- temp$assignment[(ind[,2]-1)*mcoarse + ind[,1]]
        }
      }
    }
    
    if (nscales > 1 && returncoarse) {
      return(list(sol=sol, prob=problem))
    } else {
      nbasis <- sum(temp$basis)
      res <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
      ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
      res$from <- which(wha)[ind[,1]]
      res$to <- which(whb)[ind[,2]]
      res$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]
      res$mass <- res$mass * fudgesum / fudgeN
      return(res)
    }
  }
  # fi (method == "revsimplex") 
  
  if (method == "shortsimplex") {
  	
  	# violating the first condition consistently tosses segfaults 
  	# didn't check non-squares, but the error clearly does not come from the number of sources/targets alone
  	#if (a$N <= 25 || any(a$n <= 5)) {
  	#  stop("Execution halted. There is currently a bug in method 'shortsimplex' that causes segfaults on small
  	#     grids (<= 5 points in any one dimension). Current workaround: choose method='revsimplex'")
  	#}
  	
  	# Works currently only if nscale = 1
  	# (other values of nscale are ignored)
  	# we would have to adapt C-Program so that shortlist
  	# is only used for phase 3 not for phase 2 --> do it.
  	# ----------------------------
    if (nscales || TRUE) {
      initassig <- rep(0,m*n)
      initbasis <- rep(0,m*n)
      if (control$para$slength > n) {
      	control$para$slength <- n
      	control$para$kfound <- n
      	warning("Shortlist parameter 'slength' too large...  decreased to maximal value.")
      }
      #
      #Starting solution not needed
      #   cat(as.integer(control$para$slength), as.integer(control$para$kfound), as.double(control$para$psearched),
      #   as.integer(m), as.integer(n), 
      #   as.integer(apos), as.integer(bpos), as.double(ddpos), assignment = length(initassig), 
      #   basis = length(initbasis), sep="\n")
         # stop("test")
      ti <- proc.time()
      res <- .C("shortsimplex",
          as.integer(control$para$slength), as.integer(control$para$kfound), as.double(control$para$psearched),
          as.integer(m), as.integer(n), 
          as.integer(apos), as.integer(bpos), as.double(dd), assignment = as.integer(initassig), 
          basis = as.integer(initbasis), DUP = TRUE, PACKAGE="transport")
	  print(proc.time()-ti)           
	  temp <- list(assignment=res$assignment, basis=res$basis)
    }    
    nbasis <- sum(temp$basis)
    res <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
    ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
    res$from <- which(wha)[ind[,1]]
    res$to <- which(whb)[ind[,2]]
    res$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]
    res$mass <- res$mass * fudgesum / fudgeN
    return(res)    
  }
  # fi (method == "shortsimplex") 
}



# (der Vollstaendigkeit halber solllte spaeter auch der Fall m != n wieder dazu kommen)
# vorlaeufig ohne Beobachtungsfenster
transport.pp <- function(a, b, p = 1, method = c("auction", "auctionbf", "shortsimplex", "revsimplex", "primaldual"),
                           control = list(), ...) {
  # Check inputs
  # ====================================================================== 
  if (class(a) != "pp")  a <- pp(a)
  if (class(b) != "pp")  b <- pp(b)                           	
  stopifnot(compatible(a,b))
  if (a$dimension < 2) stop("dimension must be >=2")
  N <- a$N
  
  method <- match.arg(method)

  if (class(control) != "trc") {
  	control$method <- method
  	control$a <- a
  	control$b <- b
  	control = do.call(trcontrol, control)
  }

  if (control$start != "auto") {
  	warning("control$start = ", sQuote(control), " is ignored for function transport.pp")
  }
  # nwcorner gibt Identitaet (mit erster Nebendiag fuer Basis), Russell wohl was aehnlich Degeneriertes
  # der Einfachheit halber Testen wir nur mal mit nwcorner 

  x <- a$coordinates
  y <- b$coordinates
  
  dd <- as.matrix(dist(rbind(x,y)))[1:N,(N+1):(2*N)]
  dd <- dd^p
    
  if (method != "shortsimplex" && method != "revsimplex") {
    precision=9
    dd <- dd/max(dd)
    dd <- round(dd*(10^precision))
  }
    # wir sollten mit unseren Berechnungen .Machine$integer.max nicht ueberschreiten, gemaess R-Hilfe ist dies 
    # *auf jedem System* 2147483647 (4 Bytes)
    #
    # Beachte: wenn wir Distanz zurueckgeben wollen, muessen wir natuerlich mit urspruenglichem dd^p rechnen

  if (method == "auction" || method == "auctionbf") {  
    dupper <- max(dd)/10
    lasteps <- control$para$lasteps
    epsvec <- lasteps
    # Bertsekas von dupper/2 bis 1/(n+1) durch fortgesetzt konstante Zahl teilen
    while (lasteps < dupper) {
      lasteps <- lasteps*control$para$epsfac
      epsvec <- c(epsvec,lasteps)
    }
    epsvec <- rev(epsvec)[-1]
    neps <- length(epsvec)
    stopifnot(neps >= 1)
  } 

  if (method == "auction") {
  	desirem <- max(dd)-dd
  	ti <- proc.time()
    temp <- .C("auction", as.integer(desirem), as.integer(N), pers_to_obj = as.integer(rep(-1,N)),
               price = as.double(rep(0,N)), as.integer(neps), as.double(epsvec), DUP=TRUE, PACKAGE="transport")
    print(proc.time()-ti)           
    # make pretty output
    res <- data.frame(from = 1:N, to = temp$pers_to_obj+1, mass = rep(1,N))
    return(res)
  }

  if (method == "auctionbf") {
  	desirem <- max(dd)-dd
  	ti <- proc.time()
    temp <- .C("auctionbf", as.integer(desirem), as.integer(N), pers_to_obj = as.integer(rep(-1,N)),
               price = as.double(rep(0,N)), profit = as.double(rep(0,N)), as.integer(neps), as.double(epsvec),
               DUP=TRUE, PACKAGE="transport")
    print(proc.time()-ti)           
    # make pretty output
    res <- data.frame(from = 1:N, to = temp$pers_to_obj+1, mass = rep(1,N))
    return(res)
  }

  if (method == "primaldual") {
  	ti <- proc.time()
  	temp <- .C("primaldual", as.integer(dd), as.integer(rep.int(1,N)), as.integer(rep.int(1,N)),
  	           as.integer(N), as.integer(N), flowmatrix = as.integer(integer(N^2)), 
  	           DUP=TRUE, PACKAGE="transport")
  	print(proc.time()-ti) 
  	# flowmatrix is the old term for assignment   
  	nassig <- sum(temp$flowmatrix)  # nassig sollte natuerlich gleich N sein, das ist nur zur Kontrolle
    res <- data.frame(from = rep(0,nassig), to = rep(0,nassig), mass = rep(1,nassig))
    ind <- which(matrix(as.logical(temp$flowmatrix),N,N), arr.ind=TRUE) 
    res$from <- ind[,1]
    res$to <- ind[,2]
    return(res)
  }

  if (method == "shortsimplex") {
    initassig <- initbasis <- rep(0,N*N)
    if (control$para$slength > N) {
      control$para$slength <- N
      control$para$kfound <- N
      warning("Shortlist parameter 'slength' too large...  decreased to maximal value.")
    }
    ti <- proc.time()
    temp <- .C("shortsimplex",
                as.integer(control$para$slength), as.integer(control$para$kfound), as.double(control$para$psearched),
                as.integer(N), as.integer(N), as.integer(rep.int(1,N)), as.integer(rep.int(1,N)),
	            as.double(dd), assignment = as.integer(initassig), basis = as.integer(initbasis),
	            DUP=TRUE, PACKAGE="transport")
	print(proc.time()-ti)           
    nassig <- sum(temp$assignment)  # nassig sollte natuerlich gleich N sein, das ist nur zur Kontrolle
    res <- data.frame(from = rep(0,nassig), to = rep(0,nassig), mass = rep(1,nassig))
    ind <- which(matrix(as.logical(temp$assignment),N,N), arr.ind=TRUE) 
    res$from <- ind[,1]
    res$to <- ind[,2]
    return(res)
  }

  if (method == "revsimplex") {
  	# this is nwcorner, at least modrowmin should be feasible and faster
    initassig <- initbasis <- diag(1,N,N)
    initbasis[cbind(2:N,1:(N-1))] <- 1 
    ti <- proc.time()
    temp <- .C("revsimplex", as.integer(N), as.integer(N), as.integer(rep.int(1,N)), as.integer(rep.int(1,N)),
	            as.double(dd), assignment = as.integer(initassig), basis = as.integer(initbasis),
	            DUP=TRUE, PACKAGE="transport")
	print(proc.time()-ti)           
    nassig <- sum(temp$assignment)  # nassig sollte natuerlich gleich N sein, das ist nur zur Kontrolle
    res <- data.frame(from = rep(0,nassig), to = rep(0,nassig), mass = rep(1,nassig))
    ind <- which(matrix(as.logical(temp$assignment),N,N), arr.ind=TRUE) 
    res$from <- ind[,1]
    res$to <- ind[,2]
    return(res)
  }
  
}


# Achtung: laengerfristig unbedingt so aendern, dass auch eine Startloesung uebergeben werden kann!! 
#
# # Input ist ein m - Vektor von Produktionsmengen a, ein n - Vektor von Konsumationsmengen b,
# # und eine m x n - Matrix von Transportkosten
# # primal-dual gibt in der Regel keine Basisloesung zurueck (>= m+n-1 assignments, evtl. weniger bei Degeneriertheit)
# # rev-simplex ergibt Basisloesung, immer (= m+n-1 assignments), bei Degeneriertheit ist es theoretisch in
# # extrem seltenen Spezialfaellen moeglich, dass Endlosschleife entsteht (siehe Luenberger), sonst auch m+n-1 assignments 
transport.default <- function(a, b, costm, method=c("shortsimplex", "revsimplex", "primaldual"),
                             control = list(), ...) {
  # maxmass=1e6, precision=9, 	
  # wir sollten mit unseren Berechnungen .Machine$integer.max nicht ueberschreiten, gemaess R-Hilfe ist dies 
  # *auf jedem System* 2147483647 (4 Bytes)
  # evtl. kann man bis maxmass=1e9 gehen
  method <- match.arg(method)
  M <- length(a)
  N <- length(b)
  costm <- as.matrix(costm)
  if (!all.equal(sum(a),sum(b))) {
  	warning("Sums of a and b differ substantially. sum(a)-sum(b) = ", sum(a)-sum(b), ". Scaling to probability vectors.")
  	a <- a/sum(a)
  	b <- b/sum(b)
  	# note: in general sum(a) == sum(b) will still be all FALSE (but that's ok)
  }
  stopifnot(all(dim(costm) == c(M,N)))
  # init_given = ifelse(!(is.null(initassig) || is.null(initbasis)), TRUE, FALSE)
  # if (init_given && method == "revsimplex") {
  	# stopifnot(all(dim(initassig) == c(m,n)))
    # stopifnot(all(dim(initbasis) == c(m,n)))
    # if (!all(apply(initassig, 1, sum) == a) || !all(apply(initassig, 2, sum) == b)) {
      # stop("Initial transference plan doesn't have the correct marginals.")
    # }
    # if (sum(initbasis) != m + n - 1) {
      # stop("Not the correct number of basis vectors in initbasis")	
    # }
    # if (any((initassig > 0) & (initbasis == 0))) {
      # stop("Positiv mass transfer via non-basis entry in initassig")
    # }
  # }

  if (class(control) != "trc") {
  	control$method <- method
  	control$a <- a
  	control$b <- b
  	control = do.call(trcontrol, control)
  }

  start <- control$start

  wha <- a > 0
  whb <- b > 0
  apos <- a[wha]
  bpos <- b[whb]
  dd <- costm[wha,whb]
  m <- length(apos)
  n <- length(bpos)
  asum <- sum(apos)
  bsum <- sum(bpos)
  if (!isTRUE(all.equal(asum,bsum))) {
    warning("total mass in a and b differs. Normalizing a and b to probability measures.")
    apos <- apos/asum
    bpos <- bpos/bsum
  	asum <- bsum <- 1
  }
  fudgeN <- fudgesum <- 1 
  is.natural <- function(x, tol = .Machine$double.eps^0.5)  all((abs(x - round(x)) < tol) & x > 0.5)
  if (!is.natural(apos) || !is.natural(bpos)) {
    fudgeN <- 1e9
  	fudgesum <- asum
  	apos <- round(apos/asum * fudgeN)
  	bpos <- round(bpos/bsum * fudgeN)
  	apos <- fudge(apos,fudgeN)
  	bpos <- fudge(bpos,fudgeN)
  } 
        
  if (method == "primaldual") {
    precision=9
    dd <- dd/max(dd)
    dd <- dd*(10^precision)
    #cat(as.integer(ddpos), as.integer(apos), as.integer(bpos), as.integer(m), as.integer(n),
             # flowmatrix = integer(m*n), sep="\n")
    #stop("primaldual soon")
  	ti <- proc.time()
    res1 <- .C("primaldual", as.integer(dd), as.integer(apos), as.integer(bpos), as.integer(m), as.integer(n),
              flowmatrix = integer(m*n), DUP=TRUE, PACKAGE="transport")
    print(proc.time()-ti)
    temp <- list(assignment=res1$flowmatrix, basis=as.numeric(res1$flowmatrix > 0))
     # make pretty output
    nbasis <- sum(temp$basis)
    res <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
    ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
    res$from <- which(wha)[ind[,1]]
    res$to <- which(whb)[ind[,2]]
    res$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]
    res$mass <- res$mass * fudgesum / fudgeN
    return(res)
  }

  if (method == "revsimplex") {
      #
      # Compute starting solution
  	  if (start == "russell") {
        temp <- russell(apos,bpos,dd)
        initassig <- temp$assignment
        initbasis <- temp$basis
        startgiven <- 1
      } else if (start == "nwcorner") {
        temp <- northwestcorner(apos,bpos)
        initassig <- temp$assignment
        initbasis <- temp$basis	
        startgiven <- 1
      } else {   # modrowmin in C-Code
      	initassig <- rep(0L,m*n)
      	initbasis <- rep(0L,m*n)
      	startgiven <- 0
      }
#        print(startgiven)
#        cat(as.integer(m), as.integer(n), as.integer(apos), as.integer(bpos),
#	    as.double(dd), "hello", assignment = as.integer(initassig),
#       basis = as.integer(initbasis), as.integer(startgiven), sep="\n")
#       stop("revsimplex soon")
    ti <- proc.time()
    res <- .C("revsimplex", as.integer(m), as.integer(n), as.integer(apos), as.integer(bpos),
	          as.double(dd), assignment = as.integer(initassig), basis = as.integer(initbasis), startgiven = as.integer(startgiven),
	          DUP=TRUE, PACKAGE="transport")
	print(proc.time()-ti)           
	temp <- list(assignment=res$assignment, basis=res$basis)
	
    nbasis <- sum(temp$basis)
    res <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
    ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
    res$from <- which(wha)[ind[,1]]
    res$to <- which(whb)[ind[,2]]
    res$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]
    res$mass <- res$mass * fudgesum / fudgeN
    return(res)
  }
  
  if (method == "shortsimplex") {
  	#
    initassig <- rep(0,m*n)
    initbasis <- rep(0,m*n)
    if (control$para$slength > n) {
      control$para$slength <- n
      control$para$kfound <- n
      warning("Shortlist parameter 'slength' too large...  decreased to maximal value.")
    }
    #
    # Starting solution not needed
    #cat(as.integer(control$para$slength), as.integer(control$para$kfound), as.double(control$para$psearched),
    #as.integer(m), as.integer(n), 
    #as.integer(apos), as.integer(bpos), as.double(ddpos), assignment = as.integer(initassig), basis = as.integer(initbasis), sep="\n")
    #stop("shortsimplex soon")
    ti <- proc.time()
    res <- .C("shortsimplex",
        as.integer(control$para$slength), as.integer(control$para$kfound), as.double(control$para$psearched),
        as.integer(m), as.integer(n), 
        as.integer(apos), as.integer(bpos), as.double(dd), assignment = as.integer(initassig), 
        basis = as.integer(initbasis), DUP = TRUE, PACKAGE="transport")
	print(proc.time()-ti)           
	temp <- list(assignment=res$assignment, basis=res$basis)
     
    nbasis <- sum(temp$basis)
    res <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
    ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
    res$from <- which(wha)[ind[,1]]
    res$to <- which(whb)[ind[,2]]
    res$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]
    res$mass <- res$mass * fudgesum / fudgeN
    return(res)    
  }
  # fi (method == "shortsimplex") 
}



# sum(a) == sum(b) has to be checked in external function 
# this function should deal correctly with any zeros in a and b
# always producing exactly m+n-1 basis vectors
northwestcorner <- function(a,b) {
  m <- length(a)
  n <- length(b)
  assignment <- matrix(0,m,n)
  basis <- matrix(0,m,n)
#  rowrules <- TRUE
  icur <- jcur <- 1
  aleft <- a
  bleft <- b
  massleft <- (sum(aleft)+sum(bleft) > 0)
  while (massleft) {
    mm <- min(aleft[icur],bleft[jcur])
    assignment[icur,jcur] <- mm
    aleft[icur] <- aleft[icur] - mm
    bleft[jcur] <- bleft[jcur] - mm
    if (aleft[icur] == 0) {
      basis[icur,jcur] <- 1  
      # this is the best place for assigning the basis
      # avoids out-of-bounds index an deals correctly with the
      # degenerate case
      icur <- icur+1
    } 
    if (bleft[jcur] == 0 && icur < m+1) {
      # icur < m+1 removes the only possibility to still get an out-of-bounds index
      # this happens only (provided sum(a)==sum(b)) if icur=m+1, jcur=n
      basis[icur,jcur] <- 1
      jcur <- jcur+1
    }
    massleft <- (sum(aleft)+sum(bleft) > 0)
  }
  # the following lines deal with trailing zeros in a and b
  if (icur == m+1) {   # no trailing zeros in a
  	while (jcur < n+1) {  # deal with trailing zeros in b if any
  	  basis[icur-1,jcur] <- 1
  	  jcur <- jcur+1
  	}
  } else {  # icur < m+1 meaning trailing zeros in a
  	while (icur < m) {  # deal with them if more than one (first one was alredy dealt with in final bleft-step above)
  	  icur <- icur+1
  	  basis[icur,jcur-1] <- 1
  	}
  	while (jcur < n+1) {
  	  basis[icur,jcur] <- 1
  	  jcur <- jcur+1
  	}
  }
  rownames(assignment) <- paste(a,"*",sep="")
  colnames(assignment) <- paste("*",b,sep="")
  if (!all(as.logical(basis)==(assignment > 0))) {
    warning("Starting solution is degenerate. Nothing to worry about!")
  }
  if (sum(basis) != m+n-1) {
    stop("Basis contains only ", sum(basis), " != m+n-1 = ", m+n-1, " vectors")
  }
  return(list(assignment=assignment,basis=basis))
}




# Russell 1969
russell <- function(a,b,costm) {
  stopifnot(all(a>0,b>0))
  mm <- m <- length(a)
  nn <- n <- length(b)
  assignment <- matrix(0,m,n)
  basis <- matrix(0,m,n)
  rownames(assignment) <- paste(a,"*",sep="")
  colnames(assignment) <- paste("*",b,sep="")
  
  actrow <- rep(TRUE,m)
  actcol <- rep(TRUE,n)
  
  totalmass <- sum(a)
  count <- 0
  
  while (m>0) {
  	count <- count+1
    costsub <- costm[actrow,actcol,drop=FALSE] 
  
    w <- apply(costsub,1,max)
    y <- apply(costsub,2,max)
    ind <- which.min(costsub - matrix(w,m,n) - matrix(y,m,n, byrow=TRUE))
    i0 <- matrix(1:m,m,n)[ind]
    i0 <- which(actrow)[i0]
    j0 <- matrix(1:n,m,n,byrow=TRUE)[ind]
    j0 <- which(actcol)[j0]
    
    mass <- min(a[i0],b[j0])
    if (mass == 0) {
      warning("mass 0 has been assigned!??")
    }
    assignment[i0,j0] <- mass
    basis[i0,j0] <- 1
    a[i0] <- a[i0] - mass
    b[j0] <- b[j0] - mass
    if (a[i0] == 0) {
      actrow[i0] <- FALSE
      m <- m-1
    }
    if (b[j0] == 0) {
      actcol[j0] <- FALSE
      n <- n-1
    }
  }

  if (!all(as.logical(basis)==(assignment > 0))) {
    warning("Starting solution is degenerate. Nothing to worry about!")
  }  

  if (sum(basis) < mm+nn-1) {
  	cat("Russell produced too few basis vectors. Is ", sum(basis), ", should be ", mm+nn-1,
  	        ".\n Trying to fix this...  ", sep="")
  	basis <- dedegenerate(basis)
  	cat("Done.\n")
  } else if (sum(basis) > mm+nn-1) {
  	stop("Russell produced too many basis vectors. Is ", sum(basis), ", should be ", mm+nn-1)
  	# this should not be possible
  }
      
  return(list(assignment=assignment,basis=basis))
}


# THIS FUNCTION NEEDS MUCH IMPROVEMENT!!!!
# a1 ist "a_old", a2 is "a_new"
refinesol <- function(a1,b1,a2,b2, assig1, basis1, mult=2) {
  a1 <- as.vector(a1)   # Laenge 16
  b1 <- as.vector(b1)   # Laenge 16
  a2 <- as.vector(a2)   # Laenge 64
  b2 <- as.vector(b2)
  Ngrid1 <- length(a1)
  ngrid1 <- sqrt(Ngrid1)
  Ngrid2 <- length(a2)
  ngrid2 <- sqrt(Ngrid2)
  stopifnot(ngrid2 == mult*ngrid1)
  m1 <- sum(a1 > 0)
  n1 <- sum(b1 > 0)
  m2 <- sum(a2 > 0)
  n2 <- sum(b2 > 0)
  stopifnot(length(assig1) == m1*n1)
  stopifnot(length(basis1) == m1*n1)
  assig1 <- matrix(assig1,m1,n1)
  basis1 <- matrix(basis1,m1,n1)
  assig2 <- matrix(0,m2,n2)
  basis2 <- matrix(0,m2,n2)

  # print(paste("REFINESOL IN:", "bvecs", sum(basis1), " ---- ", "expected", m1+n1-1))

  # Hilfsgroessen  
  multiple <- function(v,mult2) {
    n <- sqrt(length(v))
    stopifnot(isTRUE(all.equal(round(n), n)))
    grvec <- unlist(lapply(as.list(1:n), function(x) {rep(((x-1)*n+1):(x*n), times = mult2, each = mult2)}))
    return(v[grvec]) 
  }
  whbox <- unlist(lapply(as.list(1:ngrid1), function(x) {rep(((x-1)*ngrid1+1):(x*ngrid1), times = mult, each = mult)}))
  # tells us for a number in 1,..,Ngrid2, which box the corresponding index in terms of a2/b2 lies in 
  bybox <- order(whbox)
  # gives us a list of indices from 1,..,Ngrid2 that evaluate vectors like a2/b2 boxwise 

  isprod1 <- (a1 > 0)    # "is producer"
  isprod2old <- multiple(a1 > 0, mult)
  isprod2 <- (a2 > 0)

  pnum1 <- isprod1 
  pnum1[isprod1] <- 1:m1
  cnum1 <- !isprod1 
  cnum1[!isprod1] <- 1:n1
  pcnum1 <- pnum1 + cnum1
  # For numbers from 1 to Ngrid1: which producer or consumer number is the corresponding index in the a1/b1 matrix
  # isprod1 tells us what it is (producer or consumer)
  # for p != 1 we usually don't have the separation in producers and consumers
  # isprod1 is all TRUE, and pcnum1 = 1:m1 (could be a problem)

  pnum2 <- isprod2 
  pnum2[isprod2] <- 1:m2
  cnum2 <- !isprod2 
  cnum2[!isprod2] <- 1:n2
  pcnum2 <- pnum2 + cnum2
  # For numbers from 1 to Ngrid2: which producer or consumer number is the corresponding index in the a2/b2 matrix
  # isprod2 tells us what it is (producer or consumer)
  # for p != 1 we usually don't have the separation in producers and consumers
  # isprod1 is all TRUE, and pcnum1 = 1:m2 (could be a problem)

  # --------------------------------
  # ab hier: fuelle assig2 / basis2
  # -------------------------------- 

  # loese das +- Problem in einzelnen Feldern der groben Matrix
  a2left <- a2
  b2left <- b2

  Mult <- mult^2
  for (i in 1:Ngrid1) {
  for (j in 1:Mult) {
  k <- bybox[(i-1)*Mult + j]
  # current index in the finer grid
  # zuerst lokalen Bedarf stillen
  if (isprod2old[k] && !isprod2[k]) {
  	# k needs stuff
  	for (jj in 1:Mult) {
  	  l <- bybox[(i-1)*Mult + jj]
  	  needs <- b2left[k]  # belongs here
  	  gives <- a2left[l]
  	  amount <- min(gives, needs)
  	  if (amount > 0) {
  	  	a2left[l] <- a2left[l]-amount
  	  	b2left[k] <- b2left[k]-amount
  	    basis2[ pcnum2[l] , pcnum2[k] ] <- 1
  	    assig2[ pcnum2[l] , pcnum2[k] ] <- amount
  	  }	
  	}
  # dann lokalen UEberschuss abbauen	
  } else if (!isprod2old[k] && isprod2[k]) {
  	# k gives stuff
  	for (jj in 1:Mult) {
  	  l <- bybox[(i-1)*Mult + jj]
  	  gives <- a2left[k]  # belongs here
  	  needs <- b2left[l]
  	  amount <- min(gives, needs)
  	  if (amount > 0) {
   	  	a2left[k] <- a2left[k]-amount
  	  	b2left[l] <- b2left[l]-amount
  	    basis2[ pcnum2[k] , pcnum2[l] ] <- 1
  	    assig2[ pcnum2[k] , pcnum2[l] ] <- amount
  	  }	
  	}
  }  
  }
  }

  # Emuliere grobe Loesung auf feiner Matrix
  for (i in 1:m1) {
    transfrom1 <- which(isprod1)[i]
    transto1 <- which(!isprod1)[basis1[i,] == 1]
    # in the coarser solution there is mass transfer from transfrom1 to transto1 (in a/b terms)
    # (careful: transto1 is a vector)
    # and the following are the amounts
    transmass <- assig1[i,][basis1[i,] == 1]
    #cat("i =",i,"\n")
    #cat(transfrom1, "---", transto1, "---", transmass, "\n")
    for (r in 1:length(transto1)) {
      if (transmass[r] == 0) {
      	print(paste("input to refinesol degenerate:",i,r))
      	# degenerierter Fall
      	# noch nicht getestet, auch nicht 100%ig klar, ob die Mult+j unten genau das richtige sind
      	for (j in 1:Mult) {
      	  k <- bybox[(transfrom1-1)*Mult+j]
  	      l <- bybox[(transto1[r]-1)*Mult+j]
      	  basis2[ pcnum2[k] , pcnum2[l] ] <- 1
      	  stopifnot(assig2[ pcnum2[k] , pcnum2[l] ] == 0)
      	  # sanity check, remove asap!!
      	}
      } else {
        for (j in 1:Mult) {
        for (jj in 1:Mult) {
          # das ist im Prinzip nix anderes als lokale Northwest-Corner-Rule (mit Loechern
          # wegen internen Verschiebungen innerhalb der Pixel), geht vermutlich eleganter
          # beides die NW-Corner-Rule und ohne diese
  	      k <- bybox[(transfrom1-1)*Mult+j]
  	      l <- bybox[(transto1[r]-1)*Mult+jj]
  	      #cat(k, l, a2left[k],b2left[l],"\n")
  	      if (a2left[k] > 0 && b2left[l] > 0 && transmass[r] > 0) {  	      	
  	        amount <- min(a2left[k],b2left[l],transmass[r])
  	        a2left[k] <- a2left[k]-amount
  	        b2left[l] <- b2left[l]-amount
  	        transmass[r] <- transmass[r]-amount
  	        basis2[ pcnum2[k] , pcnum2[l] ] <- 1
  	        assig2[ pcnum2[k] , pcnum2[l] ] <- amount
  	      }    
        }  
        }
      }
    }
  }
  
  # Das folgende ist ein rechtes Gebastel und es sollte besser im Hauptteil von refinesol dagegen vorgebeugt werden
  # dedegenerate ist wahrsch. nicht noetig, Basisvektoren zufaellig hinzufuegen, sollte fast so gut sein...?
  if (sum(basis2) < m2+n2-1) {
  	cat("Refinesol produced too few basis vectors. Is ", sum(basis2), ", should be ", m2+n2-1,
  	        ".\n Trying to fix this...  ", sep="")
  	basis2 <- dedegenerate(basis2)
  	cat("Done.\n")
  } else if (sum(basis2) > m2+n2-1) {
  	stop("Refinesol produced too many basis vectors. Is ", sum(basis2), ", should be ", m2+n2-1)
  }
  return(list(basis2=basis2,assig2=assig2))
}


# refinesol for p != 1 
# currently not implemented
refinesol2 <- function(a1,b1,a2,b2, assig1, basis1, mult=2) {
  assig2 <- NULL
  basis2 <- NULL  
  return(list(basis2=basis2,assig2=assig2))
}


# vorlaeufig nur fuer quadratische Matrizen
# das folgende ist eh alles viel zu muehsam, koennte aber mal noch nuetzlich sein
triangulate <- function(basis) {
  n <- dim(basis)[1]
  stopifnot(dim(basis)[2] == n)
  stopifnot(sum(basis) <= 2*n -1)
  stopifnot(min(apply(basis,1,sum)) > 0 && min(apply(basis,2,sum)) > 0)
  # spaetere Zeilen/Spaltensummen duerfen 0 sein
  roworder <- rep(0,n)
  colorder <- rep(0,n)
  rowleft <- 1:n
  colleft <- 1:n
  for (k in 1:(n-1)) {
  	# print(k)
    rsum <- apply(basis[rowleft,colleft],1,sum)
    # wir legen hier Wert auf "moeglichst viele" Einsen auf der Diagonalen
    # rein, weil es uebersichtlicher aussieht
    wh <- which(rsum == 1)[1]
    target <- 1
    if (is.na(wh)) {
      wh <- which(rsum == 0)[1]
      target <- 0
      if (is.na(wh)) stop("Cannot triangulate basis")
    }
    roworder[k] <- rowleft[wh]
    colorder[k] <- colleft[which(basis[roworder[k],colleft] == target)[1]]
    rowleft <- rowleft[rowleft != roworder[k]]
    colleft <- colleft[colleft != colorder[k]]    	
  }
  roworder[n] <- rowleft
  colorder[n] <- colleft
  #print(basis[roworder,colorder])
  return(list(roworder=roworder,colorder=colorder,tbasis=basis[roworder,colorder]))
}

# sollte auch ohne Triangulierung funktionieren
# Vor allem koennen wir hier quadratisch gar nicht brauchen
findblocks <- function(tbasis) {
  blocks <- vector("list",0)
  bb <- 0
  m <- dim(tbasis)[1]
  n <- dim(tbasis)[2]
  rows2go <- 1:m
  cols2go <- 1:n
  while (length(rows2go) > 0 || length(cols2go) > 0) {
  	searchforrows <- TRUE
  	rowlist <- numeric(0)
    collist <- min(cols2go)
    cols2go <- cols2go[cols2go != collist]
  	bb <- bb+1
  	blocks[[bb]] <- list(row=numeric(0),col=collist)
    while (length(rowlist) > 0 || length(collist) > 0) {
  	  if (searchforrows) {
  	  	for (k in collist) {
  	      rowlist <- union(rowlist,intersect(which(tbasis[,k] == 1),rows2go))
  	      blocks[[bb]]$row <- union(blocks[[bb]]$row, rowlist)
  	      rows2go <- setdiff(rows2go,rowlist)
  	      searchforrows <- FALSE  	  
  	    }
  	    collist <- numeric(0)
  	  } else {
  	  	for (k in rowlist) {
  	  	  collist <- union(collist,intersect(which(tbasis[k,] == 1),cols2go))
  	  	  blocks[[bb]]$col <- union(blocks[[bb]]$col, collist)
  	  	  cols2go <- setdiff(cols2go,collist)
  	  	  searchforrows <- TRUE
  	  	}
  	  	rowlist <- numeric(0)	
  	  }
    }
  }
  return(blocks)
}


dedegenerate <- function(basis) {
  res <- findblocks(basis)
  ll <- length(res)
  m <- dim(basis)[1]
  n <- dim(basis)[2]
  stopifnot(sum(basis) + ll == m+n)
  basis2 <- basis
  for (i in 1:(ll-1)) {
  	k <- res[[i+1]]$row[1]
  	l <- res[[i]]$col[length(res[[i]]$col)]
  	# print(paste(k,l))
  	basis2[k,l] <- 1
  }
  return(basis2)
}


# the following function *should* be able to deal with all kinds of degeneracies;
# even zeros in aredf,bredf, but this is not (well) tested 
scalestart <- function(aredf,bredf,genf,n1f,n2f,p=1,scmult=2) { 
#  stopifnot(all(aredf>0, bredf>0))   # no longer necessary I leave it for the moment
  stopifnot(all(dim(aredf) == c(n1f,n2f) && dim(bredf) == c(n1f,n2f)))
  is.natural <- function(x, tol = .Machine$double.eps^0.5)  all((abs(x - round(x)) < tol) & x > 0.5)
  stopifnot(is.natural(scmult) && is.natural(n1f) && is.natural(n2f) && is.natural(n1f/scmult) && is.natural(n2f/scmult))
  n1c <- round(n1f/scmult)
  n2c <- round(n2f/scmult)   # f for fine, c for coarse
  Nf <- n1f * n2f
  Nc <- n1c * n2c
  grmat <- matrix(unlist(lapply(as.list(1:n2c), function(x) {rep(((x-1)*n1c+1):(x*n1c), 
      	          times = scmult, each = scmult)})), n1f, n2f)
  aredc <- matrix( aggregate(as.vector(aredf), by=list(as.vector(grmat[1:n1f,1:n2f])),
      	    FUN="sum")[,2], n1c, n2c)
  bredc <- matrix( aggregate(as.vector(bredf), by=list(as.vector(grmat[1:n1f,1:n2f])),
      	    FUN="sum")[,2], n1c, n2c)
  genc <- list(aggregate(genf[[1]], by=list(as.vector(grmat[1:n1f,1])), FUN="mean")[,2],
               aggregate(genf[[2]], by=list(as.vector(grmat[1,1:n2f])), FUN="mean")[,2])
  # preliminary version where we do not distinguish between p=1 and p!=1
  gg <- expand.grid(genc[[1]], genc[[2]])
  costmc <- as.matrix(dist(gg))^p
  assigc <- rep(0L,Nc*Nc)
  basisc <- rep(0L,Nc*Nc)
  startgiven <- 0
  stopifnot(sum(aredc)==sum(bredc))
  resc <- .C("revsimplex", as.integer(Nc), as.integer(Nc), as.integer(aredc), as.integer(bredc),
	        as.double(costmc), assignment = as.integer(assigc), basis = as.integer(basisc), startgiven = as.integer(startgiven),
	        DUP=TRUE, PACKAGE="transport")
  assigc <- matrix(resc$assignment,Nc,Nc)
  basisc <- matrix(resc$basis,Nc,Nc)
  nbasis <- sum(resc$basis)
  res2 <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
  ind <- which(matrix(as.logical(resc$basis),Nc,Nc), arr.ind=TRUE) 
  ofrom <- order(ind[,1])
  res2$from <- ind[,1][ofrom]
  res2$to <- ind[,2][ofrom]
  res2$mass <- (resc$assignment[(ind[,2]-1)*Nc + ind[,1]])[ofrom]
  stopifnot(length(res2$from) == 2*Nc-1, length(res2$to) == 2*Nc-1, length(res2$mass) == 2*Nc-1)

  #######################
  ##  refine solution  ##
  #######################
  assigf <- matrix(0,Nf,Nf)
  basisf <- matrix(0,Nf,Nf)
  qmult <- scmult^2
  
  # new strategy row-wise in the fine grid; this way it seems to be easier to treat degeneracies 
  whblock <- as.vector(grmat)
  # tells us for a number in 1,..,Nf (colmajor pos in matrix), which of the Nc blocks the corresponding index lies in 
  byblock <- order(whblock)
  # gives us a list of colmajor positions (from 1, ..., Nf) in matrix (nf x nf) that evaluates it (the matrix) by Nc-blocks  
  block <- rep(1:Nc,each=qmult)
  # same as whbox[byblock] (vector of the block numbers we are in with byblock)
  # firstpos <- unlist(lapply(seq(1,1+(n2c-1)*qmult*n1c,qmult*n1c), function(k) {seq(k,k+(n1c-1)*scmult,scmult)}))
  # "coarse to fine, first index"; vector of first indices in blocks
  acleft <- aredc
  bcleft <- bredc
  afleft <- aredf
  bfleft <- bredf
  aremfpos <- rep(1,Nc)  # if-position to start with in blocks of a
  bremfpos <- rep(1,Nc)  # jf-position to start with in blocks of b
  # in what follows we basically use nwcorner rule for each row of blocks (later improvement: use rowmostneg)
  # i is the block row number, icur is the row number (in the fine grid)
  for (i in 1:Nc) {   
#  	cat("i: ", i, "\n", sep="")
  	allifs <- byblock[qmult*(i-1) + 1:qmult]  # the qmult if-values that make up block i
  	icurpos <- aremfpos[i]
  	icur <- allifs[icurpos]
    alljcs <- which(basisc[i,] == 1)   # numbers of the to-blocks in row i
    lenjcs <- length(alljcs)
    for (jpos in 1:lenjcs) {    # ipos = i = 1 automatically, but we need icurs because enumeration along i not linear 
#      cat("jpos: ", jpos, "\n", sep="")   	
      j <- alljcs[jpos]
      alljfs <- which(whblock == j)  # all jf-indices to which we can transport
      jcurpos <- bremfpos[j]   
      jcur <- alljfs[jcurpos]
#    bcrowleft <- rep(0,Nc)   # how much mass is left for the various blocks in the i-th row
#    bcrowleft[alljcs] <- assigc[1,alljcs] 
      massleft <- assigc[i,j]
      repeat {   # break condition massleft == 0
        mass <- min(afleft[icur],bfleft[jcur],massleft)
        # whenever two of the three above terms coincide as minima, we degeneration occurs and we have to add
        # a basic vector with zero transport (if all three coincide, we have to add two basis vectors)
        assigf[icur,jcur] <- mass
        basisf[icur,jcur] <- 1
        afleft[icur] <- afleft[icur] - mass
        bfleft[jcur] <- bfleft[jcur] - mass
        massleft <- massleft - mass
        if (massleft == 0) {
          break
        } else {
          if (afleft[icur] == 0 && bfleft[jcur] == 0) {
      	    icurpos <- icurpos+1
            icur <- allifs[icurpos]
            basisf[icur,jcur] <- 1
            jcurpos <- jcurpos+1
            jcur <- alljfs[jcurpos]
          } else if (afleft[icur] == 0) {
            icurpos <- icurpos+1
            icur <- allifs[icurpos]
          } else if (bfleft[jcur] == 0) {
            jcurpos <- jcurpos+1
            jcur <- alljfs[jcurpos]          	
          }
        }
      }  # repeat loop  
      acleft[i] <- acleft[i] - assigc[i,j]
      bcleft[j] <- bcleft[j] - assigc[i,j]
      # do finishing work for block (i,j); note: we just assigned last bit of mass
      if (acleft[i] == 0 && bcleft[j] == 0) {
        # we do not have to check if there are any more blocks in this row or in this col
        # because any such blocks can only have total transport zero --> code below works out  
        # First we fill L-shape basis vectors for trailing zeros in a and b         	                                            
        while (icurpos < qmult) {
          icurpos <- icurpos+1
          icur <- allifs[icurpos]
          basisf[icur,jcur] <- 1
        }
        while (jcurpos < qmult) {
          jcurpos <- jcurpos+1
          jcur <- alljfs[jcurpos]
          basisf[icur,jcur] <- 1
        }
        # Then we make sure that the potential zero-transport blocks to come fill in basis-1s
        # from the right position
        aremfpos[i] <- qmult   # = icurpos
        bremfpos[j] <- qmult   # = jcurpos
      } else if (acleft[i] == 0) {  # hence bcleft[j] > 0, hence we know that
        	                        # there must be another block in the same col
        if (bfleft[jcur] == 0) {
          jcurpos <- jcurpos+1
          jcur <- alljfs[jcurpos]
          basisf[icur,jcur] <- 1
        }
        while (icurpos < qmult) {
          icurpos <- icurpos+1
          icur <- allifs[icurpos]
          basisf[icur,jcur] <- 1
        }
        aremfpos[i] <- qmult    # = icurpos
        bremfpos[j] <- jcurpos
      } else if (bcleft[j] == 0) {  # hence acleft[i] > 0, hence we know that
        	                        # there must be another block in the same row
        if (afleft[icur] == 0) {
          icurpos <- icurpos+1
          icur <- allifs[icurpos]
          basisf[icur,jcur] <- 1 
        }
        while (jcurpos < qmult) {
          jcurpos <- jcurpos+1
          jcur <- alljfs[jcurpos]
          basisf[icur,jcur] <- 1
        }
        aremfpos[i] <- icurpos
        bremfpos[j] <- qmult    # = jcurpos      	
      } else {   # the case acleft[i] > 0 && bcleft[j] > 0
        if (afleft[icur] == 0 && bfleft[jcur] == 0) {
          icurpos <- icurpos+1
          icur <- allifs[icurpos]
          basisf[icur,jcur] <- 1
          jcurpos <- jcurpos+1
          jcur <- alljfs[jcurpos]
          basisf[icur,jcur] <- 1
        } else if (afleft[icur] == 0) {
          icurpos <- icurpos+1
          icur <- allifs[icurpos]
          basisf[icur,jcur] <- 1
        } else if (bfleft[jcur] == 0) {
          jcurpos <- jcurpos+1
          jcur <- alljfs[jcurpos]
          basisf[icur,jcur] <- 1          	
        }
        aremfpos[i] <- icurpos
        bremfpos[j] <- jcurpos
      }	
    }   # for jpos loop
  }     # for i loop 
  rownames(assigf) <- paste(aredf,"*",sep="")
  colnames(assigf) <- paste("*",bredf,sep="")
  if (!all(as.logical(basisf)==(assigf > 0))) {
    warning("Starting solution is degenerate. Nothing to worry about!")  
  }  
  if (sum(basisf) != 2*Nf-1) {
    stop("Basis contains only ", sum(basisf), " != 2*Nf-1 = ", 2*Nf-1, " vectors")
  }
  nbasis <- sum(basisf)
  res3 <- data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis))
  ind <- which(matrix(as.logical(basisf),Nf,Nf), arr.ind=TRUE) 
  ofrom <- order(ind[,1])
  res3$from <- ind[,1][ofrom]
  res3$to <- ind[,2][ofrom]
  res3$mass <- (assigf[(ind[,2]-1)*Nf + ind[,1]])[ofrom]
  return(list(assignment=assigf,basis=basisf,res=res3))
}
