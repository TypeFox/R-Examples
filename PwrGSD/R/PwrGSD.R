# Nodes and weights for Gauss-Legendre quadrature on [-1,1] for n=24

glegx24 <- c(-0.995187219997, -0.974728555971, -0.938274552003, -0.886415527004,
             -0.820001985974, -0.740124191579, -0.648093651937, -0.545421471389,
             -0.433793507626, -0.315042679696, -0.191118867474, -0.0640568928626,
              0.0640568928626, 0.191118867474,  0.315042679696,  0.433793507626,
              0.545421471389,  0.648093651937,  0.740124191579,  0.820001985974,
              0.886415527004,  0.938274552003,  0.974728555971,  0.995187219997)

glegw24 <- c( 0.0123412297185, 0.0285313860439, 0.0442774398676, 0.0592985850337,
              0.0733464816501, 0.0861901617117, 0.0976186522496, 0.107444270284,
              0.115505668234,  0.121670473118,  0.125837456543,  0.127938195546,
              0.127938195546,  0.125837456543,  0.121670473118,  0.115505668234,
              0.107444270284,  0.0976186522496, 0.0861901617117, 0.0733464816501,
              0.0592985850337, 0.0442774398676, 0.0285313860439, 0.0123412297185)

"%,%" <- function(x,y) paste(x,y,sep="")

DX <- function(x)c(x[1],diff(x))

"PwrGSD" <- 
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         NonBindingFutility=TRUE,
         sided =c("2>", "2<", "1>", "1<"),method=c("S","A"),accru,accrat,tlook,
         tcut0 = NULL,h0 = NULL,s0 = NULL,tcut1 = NULL,rhaz = NULL,
         h1 = NULL,s1 = NULL,tcutc0 = NULL,hc0 = NULL,sc0 = NULL,tcutc1 = NULL,hc1 = NULL,
         sc1 = NULL,tcutd0A = NULL,hd0A = NULL,sd0A = NULL,tcutd0B = NULL,hd0B = NULL,sd0B = NULL,
         tcutd1A = NULL,hd1A = NULL,sd1A = NULL,tcutd1B = NULL,hd1B = NULL,sd1B = NULL,
         tcutx0A = NULL,hx0A = NULL,sx0A = NULL,tcutx0B = NULL,hx0B = NULL,sx0B = NULL,
         tcutx1A = NULL,hx1A = NULL,sx1A = NULL,tcutx1B = NULL,hx1B = NULL,sx1B = NULL,
         noncompliance = c("none","crossover","mixed","user"),gradual = FALSE,
         WtFun = c("FH","SFH","Ramp"),ppar = cbind(c(0,0)),
         Spend.Info=c("Variance","Events","Hybrid(k)","Calendar"),
         RR.Futility = NULL,qProp.one.or.Q = c("one","Q"),Nsim=NULL,
         detail = FALSE,StatType = c("WLR","ISD"),doProj=FALSE)
{
  .call. <- match.call()
  prfx <- c("Sim","Asy")
  if(missing(method)) method <- "A"
  .call.[[1]] <- as.name(prfx[1+(method=="A")] %,% as.character(.call.[[1]]))
  .call.$method <- NULL
  if(method=="A") .call.$Nsim <- .call.$detail <- .call.$StatType <- NULL
  if(method=="S"){
    if(missing(Nsim)) stop("Argument 'Nsim' is required for method==\"S\"")
    if(missing(StatType)) StatType <- "WLR"
  }
  ans <- eval(.call.)
  ans$call[[1]] <- as.name("PwrGSD")
  ans$call$method <- method
  ans
}

"AsyPwrGSD" <- 
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         NonBindingFutility=TRUE,
         sided=c("2>","2<", "1>","1<"),accru,accrat,tlook,tcut0=NULL,h0=NULL,s0=NULL,tcut1=NULL,
         rhaz=NULL, h1=NULL, s1=NULL, tcutc0=NULL, hc0 = NULL, sc0=NULL, 
         tcutc1=NULL, hc1 = NULL, sc1 = NULL, tcutd0A=NULL, hd0A=NULL, 
         sd0A=NULL, tcutd0B=NULL, hd0B=NULL, sd0B=NULL, tcutd1A=NULL, hd1A=NULL, 
         sd1A=NULL, tcutd1B=NULL, hd1B=NULL, sd1B=NULL, tcutx0A=NULL, hx0A=NULL, 
         sx0A=NULL, tcutx0B=NULL, hx0B=NULL, sx0B=NULL, tcutx1A=NULL, hx1A=NULL, 
         sx1A=NULL, tcutx1B=NULL, hx1B=NULL, sx1B=NULL, 
         noncompliance=c("none", "crossover", "mixed", "user"), 
         gradual = FALSE, WtFun=c("FH","SFH","Ramp"), ppar=c(0, 0),
         Spend.Info=c("Variance", "Events", "Hybrid(k)", "Calendar"),
         RR.Futility=NULL, V.end = NULL,qProp.one.or.Q = c("one","Q"))
{
    .call. <- match.call()
    normut <- 8.20953615160139
    psimin <- 1.0e-15

    if(missing(WtFun)) {
      WtFun <- "FH"
      wttyp <- 0
      ppar <- c(0,0)
      nstat <- 1
    }
    else{
      if(missing(ppar) && any(WtFun!="FH"))
        stop("You must specify parameters for the chosen weight function(s) in the 'ppar' argument.")
      nstat <- length(WtFun)
      wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1
    }
    nlook <- length(tlook)
    tend <- tlook[nlook]
    if (missing(sided)) sided <- "2>"
    if(!(sided%in%c("2>","2<","1>","1<")))
      stop("Argument 'sided' must be " %,% 
           "equal to \"2>\", \"2<\", \"1>\" or \"1<\"")
    sided <- c(2,-2, 1,-1)[grep(sided, c("2>","2<", "1>","1<"))]

    call.PBS <- .call.
    call.PBS[[1]] <- as.name("ParseBoundarySelection")

    call.PBS$sided <- call.PBS$accru <- call.PBS$accrat <- call.PBS$tlook <- 
    call.PBS$tcut0 <- call.PBS$h0 <- call.PBS$s0 <- call.PBS$tcut1 <- 
    call.PBS$rhaz <- call.PBS$h1 <- call.PBS$s1 <- call.PBS$tcutc0 <-
    call.PBS$hc0  <- call.PBS$sc0 <- call.PBS$tcutc1 <- call.PBS$hc1 <-
    call.PBS$sc1 <- call.PBS$tcutd0A <- call.PBS$hd0A <- call.PBS$sd0A <-
    call.PBS$tcutd0B <- call.PBS$hd0B <- call.PBS$sd0B <- call.PBS$tcutd1A <-
    call.PBS$hd1A <- call.PBS$sd1A <- call.PBS$tcutd1B <- call.PBS$hd1B <- 
    call.PBS$sd1B<- call.PBS$tcutx0A<- call.PBS$hx0A<- call.PBS$sx0A <-
    call.PBS$tcutx0B <- call.PBS$hx0B<- call.PBS$sx0B<- call.PBS$tcutx1A <-
    call.PBS$hx1A <- call.PBS$sx1A <- call.PBS$tcutx1B <- call.PBS$hx1B <-
    call.PBS$sx1B <- call.PBS$noncompliance <- call.PBS$gradual <-
    call.PBS$WtFun <- call.PBS$ppar <- call.PBS$Spend.Info <-
    call.PBS$RR.Futility <- call.PBS$V.end <- call.PBS$qProp.one.or.Q <-
    call.PBS$Nsim <- call.PBS$doProj <- NULL
    call.PBS$n.looks <- nlook
    call.PBS$check.drift <- FALSE
    
    eval.PBS <- eval(call.PBS)

    do.efficacy <- TRUE
    do.futility <- !is.null(eval.PBS$frontend$FutilityBoundary)

    Alpha.Efficacy <- eval.PBS$backend$Alpha.Efficacy
    Alpha.Futility <- eval.PBS$backend$Alpha.Futility
    if(is.null(Alpha.Futility)) Alpha.Futility <- 0
    nbnd.e <- eval.PBS$backend$nbnd.e
    nbnd.f <- eval.PBS$backend$nbnd.f
    nsf.e <- eval.PBS$backend$nsf.e
    nsf.f <- eval.PBS$backend$nsf.f
    rho.Efficacy <- eval.PBS$backend$rho.Efficacy
    rho.Futility <- eval.PBS$backend$rho.Futility
    b.Haybittle.e <- eval.PBS$backend$b.Haybittle.e
    b.Haybittle.f <- eval.PBS$backend$b.Haybittle.f
    drift.end <- eval.PBS$backend$drift.end
    be.end <- eval.PBS$backend$be.end
    prob.e <- eval.PBS$backend$prob.e
    prob.f <- eval.PBS$backend$prob.f
    my.Efficacy <- eval.PBS$backend$my.Efficacy
    my.Futility <- eval.PBS$backend$my.Futility
    is.myE <- !all(my.Efficacy==0)
    is.myF <- !all(my.Futility==0)                                          

    b.e <- rep(0, nlook)
    if (is.myE) b.e <- my.Efficacy
    
    b.f <- rep(0, nlook)
    if (is.myF) b.f <- my.Futility

    Alpha.Efficacy <- Alpha.Efficacy/2^(abs(sided) == 2)
    
    no.SpndInfo <- missing(Spend.Info)
    if(no.SpndInfo) Spend.Info <- "Variance"
    if(!no.SpndInfo){
      Spend.Info <- as.character(.call.$Spend.Info)
      if(!(Spend.Info[1] %in% c("Variance", "Events", "Hybrid", "Calendar")))
        stop("Argument 'Spend.Info' is an expression of the form 'Variance', 'Events', 'Hybrid(k)', or 'Calendar'")
    }
    spend.info <- grep(Spend.Info[1], c("Variance", "Events", "Hybrid", "Calendar")) - 1
    spend.info.k <- 0
    if(spend.info==2 && length(Spend.Info)>1) spend.info.k <- as.integer(as.numeric(Spend.Info[2])) - 1

    user.V.end <- !missing(V.end)
    if(!user.V.end) V.end <- 0

    if(missing(qProp.one.or.Q))
      qProp.one.or.Q <- 0
    else{
      if(!(qProp.one.or.Q %in% c("one","Q")))
        stop("Argument 'qProp.one.or.Q' must be either \"one\" or \"Q\"")
      qProp.one.or.Q <- grep(qProp.one.or.Q, c("one", "Q")) - 1
    }

    tcut.u <- sort(unique(c(tcut0, tcut1, tcutc0, tcutc1, tcutd0A, tcutd0B, tcutd1A, tcutd1B, tcutx0A,
                                 tcutx0B, tcutx1A, tcutx1B)))
    nunq <- length(tcut.u)

    glegx <- glegx24
    glegw <- glegw24

    NGaussQ <- length(glegx)
    stoh <- function(tcut, s)
    {
    	ncut <- length(tcut)
    	Sold <- c(1, s[ - (ncut - 1)])
    	dt <- diff(tcut)
        log.0 <-
          function(x)
          {
            y <- x*0
            y[x>0] <- log(x[x>0])
            y
          }
    	h <- log.0(s/Sold)/dt
    	h
    }

    no.t0 <- missing(tcut0)
    no.s0 <- missing(s0)
    no.h0 <- missing(h0)
    if((no.s0 && no.h0) || no.t0)
    	stop("Must specify 'tcut0' and ('h0' or 's0').")
    if(no.h0)
    	h0 <- stoh(tcut0, s0)
    ncut0 <- length(tcut0)

    no.t1 <- missing(tcut1)
    no.s1 <- missing(s1)
    no.h1 <- missing(h1)
    no.rhaz <- missing(rhaz)
    if((no.s1 && no.rhaz && no.h1) || no.t1)
    	stop("Must specify 'tcut1' and ('rhaz' or 'h1' or 's1').")
    if(!no.s1)
    	h1 <- stoh(tcut1, s1)
    if(!no.rhaz){
    	tcut.01 <- support(c(tcut0, tcut1))
    	h0.01 <- approx(tcut0, h0, tcut.01, method="constant", f=0, yleft=0, yright=h0[ncut0])$y
    	rhaz.01 <- approx(tcut1, rhaz, tcut.01, method="constant", f=0, yleft=0, yright=rhaz[length(rhaz)])$y
    	tcut1 <- tcut.01
    	h1 <- rhaz.01 * h0.01
    }
    ncut1 <- length(tcut1)

    tlook.not <- tlook[sapply(tlook, FUN=function(x,Y){all(abs(x - Y)>1e-8)}, Y=tcut0)]
    tcut0.new <- sort(unique(c(tcut0, tlook.not)))
    ncut0.new <- length(tcut0.new)
    h0.new <- approx(tcut0, h0, tcut0.new, method="const", f=0, yleft=0, yright=h0[ncut0])$y
    tcut0 <- tcut0.new
    h0 <- h0.new
    ncut0 <- ncut0.new

    tlook.not <- tlook[sapply(tlook, FUN=function(x,Y){all(abs(x - Y)>1e-8)}, Y=tcut1)]    
    tcut1.new <- sort(unique(c(tcut1, tlook.not)))
    ncut1.new <- length(tcut1.new)
    h1.new <- approx(tcut1, h1, tcut1.new, method="const", f=0, yleft=0, yright=h1[ncut1])$y
    tcut1 <- tcut1.new
    h1 <- h1.new
    ncut1 <- ncut1.new
    
    
    use.rhaz.fu <- 0
    if(missing(RR.Futility)){
      RR.Futility <- rep(1,ncut0)
      if(do.futility==1) use.rhaz.fu <- 1
    }
    if(length(RR.Futility)<ncut0) RR.Futility <- rep(RR.Futility[1], ncut0)
    
    no.tc0 <- missing(tcutc0)
    no.sc0 <- missing(sc0)
    no.hc0 <- missing(hc0)
    if((no.sc0 && no.hc0) || no.tc0)
    	stop("Must specify 'tcutc0' and ('hc0' or 'sc0').")
    if(no.hc0)
    	hc0 <- stoh(tcutc0, sc0)
    ncutc0 <- length(tcutc0)

    
    no.tc1 <- missing(tcutc1)
    no.sc1 <- missing(sc1)
    no.hc1 <- missing(hc1)
    if((no.sc1 && no.hc1) || no.tc1)
    	stop("Must specify 'tcutc1' and ('hc1' or 'sc1').")
    if(no.hc1)
    	hc1 <- stoh(tcutc1, sc1)
    ncutc1 <- length(tcutc1)
    noncompliance <- as.character(.call.$noncompliance)
    no.noncomp <- (length(noncompliance) == 0)
    if(no.noncomp)
    	noncompliance <- "none"
    switch(noncompliance,
      none =
      {
          tcutd0A <- c(0,tend)
          hd0A <- c(1e-7,0)
          ncutd0A <- 2

          tcutd0B <- c(0,tend)
    	  hd0B <- c(1e-7,0)
    	  ncutd0B <- 2
          
          tcutd1A <- c(0,tend)
          hd1A <- c(1e-7,0)
          ncutd1A <- 2

          tcutd1B <- c(0,tend)
          hd1B <- c(1e-7,0)
          ncutd1B <- 2

          tcutx0A <- tcut0
          hx0A <- h0
          ncutx0A <- ncut0

          tcutx0B <- tcut0
          hx0B <- h0
          ncutx0B <- ncut0

          tcutx1A <- tcut1
          hx1A <- h1
          ncutx1A <- ncut1

          tcutx1B <- tcut1
          hx1B <- h1
          ncutx1B <- ncut1
      },
      crossover =
      {
          no.td0B <- missing(tcutd0B)
          no.sd0B <- missing(sd0B)
          no.hd0B <- missing(hd0B)

          no.td1B <- missing(tcutd1B)
          no.sd1B <- missing(sd1B)
          no.hd1B <- missing(hd1B)

          no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) || 
                      no.td0B || no.td1B
          if(no.dB)
              stop("crossover option requires specification of \n'tcutd0B', 'tcutd1B', ('hd0B' or 'sd0B') and " %,%
                   "('hd1B' or 'sd1B').\n")

          if(no.hd0B) hd0B <- stoh(tcutd0B, sd0B)
          ncutd0B <- length(tcutd0B)

          if(no.hd1B) hd1B <- stoh(tcutd1B, sd1B)
          ncutd1B <- length(tcutd1B)

          tcutd0A <- c(0,tend)
          hd0A <- c(1e-7,0)
          ncutd0A <- 2

          tcutd1A <- c(0,tend)
          hd1A <- c(1e-7,0)
          ncutd1A <- 2

          tcutx0A <- tcut0
          hx0A <- h0
          ncutx0A <- ncut0

          tcutx1A <- tcut1
          hx1A <- h1
          ncutx1A <- ncut1

          tcutx0B <- tcut1
          hx0B <- h1
          ncutx0B <- ncut1

          tcutx1B <- tcut0
          hx1B <- h0
          ncutx1B <- ncut0
      },
      mixed =
      {
          no.td0A <- missing(tcutd0A)
          no.sd0A <- missing(sd0A)
          no.hd0A <- missing(hd0A)

          no.td0B <- missing(tcutd0B)
          no.sd0B <- missing(sd0B)
          no.hd0B <- missing(hd0B)

          no.td1A <- missing(tcutd1A)
          no.sd1A <- missing(sd1A)
          no.hd1A <- missing(hd1A)

          no.td1B <- missing(tcutd1B)
          no.sd1B <- missing(sd1B)
          no.hd1B <- missing(hd1B)

          no.tx0A <- missing(tcutx0A)
          no.sx0A <- missing(sx0A)
          no.hx0A <- missing(hx0A)

          no.tx1A <- missing(tcutx1A)
          no.sx1A <- missing(sx1A)
          no.hx1A <- missing(hx1A)

          no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) || 
                               no.td0A || no.td1A
          no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) || 
                               no.td0B || no.td1B
          no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) || 
                               no.tx0A || no.tx1A
          if(no.dA || no.dB || no.xA) 
              stop("mixed option requires specification of \n'tcutd0A', 'tcutd1A', ('hd0A' or 'sd0A')," %,%
                   " ('hd1A' or 'sd1A'),\n'tcutd0B', 'tcutd1B', ('hd0B' or 'sd0B')," %,%
                   " ('hd1B' or 'sd1B'),\n'tcutx0A', 'tcutx1A', ('hx0A' or 'sx0A')," %,%
                   " and ('hx1A' or 'sx1A').\n")

          if(no.hd0A) hd0A <- stoh(tcutd0A, sd0A)
          ncutd0A <- length(tcutd0A)

          if(no.hd1A) hd1A <- stoh(tcutd1A, sd1A)
          ncutd1A <- length(tcutd1A)

          if(no.hd0B) hd0B <- stoh(tcutd0B, sd0B)
          ncutd0B <- length(tcutd0B)

          if(no.hd1B) hd1B <- stoh(tcutd1B, sd1B)
          ncutd1B <- length(tcutd1B)

          if(no.hx0A) hx0A <- stoh(tcutx0A, sx0A)
          ncutx0A <- length(tcutx0A)

          if(no.hx1A) hx1A <- stoh(tcutx1A, sx1A)
          ncutx1A <- length(tcutx1A)

          tcutx0B <- tcut1
          hx0B <- h1
          ncutx0B <- ncut1

          tcutx1B <- tcut0
          hx1B <- h0
          ncutx1B <- ncut0
      },
      user =
      {
          no.td0A <- missing(tcutd0A)
          no.sd0A <- missing(sd0A)
          no.hd0A <- missing(hd0A)

          no.td0B <- missing(tcutd0B)
          no.sd0B <- missing(sd0B)
          no.hd0B <- missing(hd0B)

          no.td1A <- missing(tcutd1A)
          no.sd1A <- missing(sd1A)
          no.hd1A <- missing(hd1A)

          no.td1B <- missing(tcutd1B)
          no.sd1B <- missing(sd1B)
          no.hd1B <- missing(hd1B)

          no.tx0A <- missing(tcutx0A)
          no.sx0A <- missing(sx0A)
          no.hx0A <- missing(hx0A)

          no.tx1A <- missing(tcutx1A)
          no.sx1A <- missing(sx1A)
          no.hx1A <- missing(hx1A)

          no.tx0B <- missing(tcutx0B)
          no.sx0B <- missing(sx0B)
          no.hx0B <- missing(hx0B)

          no.tx1B <- missing(tcutx1B)
          no.sx1B <- missing(sx1B)
          no.hx1B <- missing(hx1B)

          no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) || 
                               no.td0A || no.td1A
          no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) || 
                               no.td0B || no.td1B
          no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) || 
                               no.tx0A || no.tx1A
          no.xB <- (no.sx0B && no.hx0B) || (no.sx1B && no.hx1B) || 
                               no.tx0B || no.tx1B

          if((no.dA || no.xA) && (no.dB || no.xB)) 
              stop("user option requires specification of \n('tcutd0A', 'tcutd1A', ('hd0A' or 'sd0A')," %,%
                   " ('hd1A' or 'sd1A') and\n'tcutx0A', 'tcutx1A', ('hx0A' or 'sx0A')," %,%
                   " ('hx1A' or 'sx1A')) or \n('tcutd0B', 'tcutd1B', ('hd0B' or 'sd0B')," %,%
                   " and ('hd1B' or 'sd1B') and\n 'tcutx0B', 'tcutx1B', ('hx0B' or 'sx0B')," %,%
                   " and ('hx1B' or 'sx1B')).\n")

          if(!no.dA) {
              if(no.hd0A) hd0A <- stoh(tcutd0A, sd0A)
                 ncutd0A <- length(tcutd0A)

              if(no.hd1A) hd1A <- stoh(tcutd1A, sd1A)
                  ncutd1A <- length(tcutd1A)

              if(no.hx0A) hx0A <- stoh(tcutx0A, sx0A)
                 ncutx0A <- length(tcutx0A)

              if(no.hx1A) hx1A <- stoh(tcutx1A, sx1A)
                 ncutx1A <- length(tcutx1A)

              if(no.dB){
                  tcutd0B <- c(0,tend)
                  hd0B <- c(1e-7,0)
                  ncutd0B <- 2

                  tcutd1B <- c(0,tend)
                  hd1B <- c(1e-7,0)
                  ncutd1B <- 2

                  tcutx0B <- tcut0
                  hx0B <- h0
                  ncutx0B <- ncut0

                  tcutx1B <- tcut1
                  hx1B <- h1
                  ncutx1B <- ncut1
              }
          }
          if(!no.dB){
              if(no.hd0B) hd0B <- stoh(tcutd0B, sd0B)
              ncutd0B <- length(tcutd0B)

              if(no.hd1B) hd1B <- stoh(tcutd1B, sd1B)
              ncutd1B <- length(tcutd1B)

              if(no.hx0B) hx0B <- stoh(tcutx0B, sx0B)
              ncutx0B <- length(tcutx0B)

              if(no.hx1B) hx1B <- stoh(tcutx1B, sx1B)
              ncutx1B <- length(tcutx1B)

              if(no.dA){
                  tcutd0A <- c(0,tend)
                  hd0A <- c(1e-7,0)
                  ncutd0A <- 2

                  tcutd1A <- c(0,tend)
                  hd1A <- c(1e-7,0)
                  ncutd1A <- 2

                  tcutx0A <- tcut0
                  hx0A <- h0
                  ncutx0A <- ncut0

                  tcutx1A <- tcut1
                  hx1A <- h1
                  ncutx1A <- ncut1
              }
          }
      })
    nmax <- max(ncut0,ncut1,nlook)
    cumsum.nppar <- 0
    stat.nms <- NULL
    for(j in 1:nstat){
      prfx <- c("FH-", "SFH-", "R-")[1+wttyp[j]]
      nppar <- c(2,3,1)[1+wttyp[j]]
      ppar.j <- ppar[(cumsum.nppar+1):(cumsum.nppar+nppar)]
      par.string <- glue.2.string(ppar.j, sep="|")
      stat.nms <- c(stat.nms, prfx %,% par.string)
      cumsum.nppar <- cumsum.nppar + nppar
    }
    nbetyp <- length(nbnd.e)
    nbftyp <- length(nbnd.f)
    nsum <- ncut0 + ncut1 + nlook - 2
    wttyp.nms <- c("FH","SFH","Ramp")[1+wttyp]
    ints <- c(nlook,nstat,NGaussQ,ncut0,ncut1,ncutc0,ncutc1,ncutd0A,ncutd0B,ncutd1A,
              ncutd1B,ncutx0A,ncutx0B,ncutx1A,ncutx1B,gradual,nbnd.e,nbnd.f,
              nsf.e,nsf.f,do.futility,use.rhaz.fu,spend.info,user.V.end,is.myE,is.myF,spend.info.k,
              qProp.one.or.Q, sided, NonBindingFutility,wttyp)
    ints.nms <- c("nlook","nstat","NGaussQ","ncut0","ncut1","ncutc0",
                  "ncutc1","ncutd0A","ncutd0B","ncutd1A","ncutd1B","ncutx0A","ncutx0B",
                  "ncutx1A","ncutx1B","gradual","nbnd.e." %,% (1:nlook),"nbnd.f." %,% (1:nlook),
                  "nsf.e." %,% (1:nlook),"nsf.f." %,% (1:nlook),"do.futility","use.rhaz.fu",
                  "spend.info","user.V.end","is.myE","is.myF","spend.info.k","qis1orQ","sided","nbf",
                  wttyp.nms)

    dbls <- c(accru, accrat, rho.Efficacy, rho.Futility, prob.e, prob.f)
    ntrial <- floor(accru*accrat)
    dbls.nms <- c("accru", "accrat", "rho.Efficacy." %,%(1:nlook), "rho.Futility." %,% (1:nlook),
                  "prob.e." %,%(1:nlook), "prob.f." %,%(1:nlook))
    
    names(ints) <- ints.nms
    ans <- .C("AsyPwrGSD",
    	ints = as.integer(ints),
        dbls = as.double(dbls),
        pttlook = as.double(tlook),
    	palphatot = as.double(c(Alpha.Efficacy,Alpha.Futility)),
    	lrrf = as.double(log(RR.Futility)),
        bHay = as.double(c(b.Haybittle.e, b.Haybittle.f)),
    	ppar = as.double(ppar),
    	pgqxw = as.double(c(glegx,glegw)),
    	tcut0 = as.double(tcut0),
    	h0 = as.double(h0),
    	tcut1 = as.double(tcut1),
    	h1 = as.double(h1),
    	tcutc0 = as.double(tcutc0),
    	hc0 = as.double(hc0),
    	tcutc1 = as.double(tcutc1),
    	hc1 = as.double(hc1),
    	tcutd0A = as.double(tcutd0A),
    	hd0A = as.double(hd0A),
    	tcutd0B = as.double(tcutd0B),
    	hd0B = as.double(hd0B),
    	tcutd1A = as.double(tcutd1A),
    	hd1A = as.double(hd1A),
    	tcutd1B = as.double(tcutd1B),
    	hd1B = as.double(hd1B),
    	tcutx0A = as.double(tcutx0A),
    	hx0A = as.double(hx0A),
    	tcutx0B = as.double(tcutx0B),
    	hx0B = as.double(hx0B),
    	tcutx1A = as.double(tcutx1A),
    	hx1A = as.double(hx1A),
    	tcutx1B = as.double(tcutx1B),
    	hx1B = as.double(hx1B),
        V.end = as.double(V.end),
        pinffrac = double(nstat*nlook),
        pinffrac.ii= double(nstat*nlook),
    	pbounds = as.double(c(rep(b.e, nstat), rep(b.f, nstat))),  
    	mufu = double(nstat*nlook),
    	mu = double(nstat*nlook),
        betastar = double(nstat),
    	palpha0vec = double(2*nstat*nlook),
    	palpha1vec = double(2*nstat*nlook),
    	RR = double(5*nsum),
        pnnjmp = integer(1),
        E.NT = double(nstat*nsum),
        Var = double(nstat*nsum),
        Eta = double(nstat*nsum),
        tidx = integer(nlook),
        betabdry = double(2*nstat*nlook),
        bstar = double(2*nstat),
    	PACKAGE = "PwrGSD")
    detail <- ans
    detail$ntrial <- ntrial
    Pwr <- t(matrix(detail$palpha1vec,nlook,2*nstat))
    dErrII <- Pwr[nstat + (1:nstat),,drop=FALSE]
    Pwr <- Pwr[1:nstat,,drop=FALSE]
    dimnames(dErrII) <- dimnames(Pwr) <- list(stat.nms,tlook)
    njmp <- detail$pnnjmp
    detail$pinffrac <- matrix(detail$pinffrac, nlook, nstat)
    detail$pinffrac.ii <- matrix(detail$pinffrac.ii, nlook, nstat)
    detail$pbounds <- matrix(detail$pbounds,nlook,2*nstat)
    detail$betabdry <- matrix(detail$betabdry, nlook, 2*nstat)
    detail$mu <- matrix(detail$mu, nlook, nstat)
    detail$mufu <- matrix(detail$mufu, nlook, nstat)
    detail$Var <- matrix(detail$Var[1:(njmp*nstat)], njmp, nstat)
    detail$Eta <- matrix(detail$Eta[1:(njmp*nstat)], njmp, nstat)
    detail$palpha0vec <- matrix(detail$palpha0vec, nlook, 2*nstat)
    detail$palpha1vec <- matrix(detail$palpha1vec, nlook, 2*nstat)
    detail$RR <- matrix(detail$RR[1:(5*njmp)],njmp,5)
    detail$E.NT <- 4*matrix(detail$E.NT[1:(njmp*nstat)], njmp, nstat)
    detail$tidx <- detail$tidx + 1
    detail$bstar <- matrix(detail$bstar[1:(2*nstat)], nstat, 2)
    dimnames(detail$pinffrac) <-
    dimnames(detail$pinffrac.ii) <-
    dimnames(detail$mu) <- list(tlook, stat.nms)
    dimnames(detail$mufu) <- list(tlook, stat.nms)
    dimnames(detail$Var) <-
    dimnames(detail$Eta) <-
    dimnames(detail$E.NT) <- list(round(detail$RR[,1],7), stat.nms)
    dimnames(detail$pbounds) <- 
    dimnames(detail$betabdry) <- 
    dimnames(detail$palpha0vec) <- 
    dimnames(detail$palpha1vec) <- list(tlook, c(outer(stat.nms,c("-e","-f"),
    				    FUN="%,%")))
    dimnames(detail$RR) <- list(rep("",njmp), c("t","htlde0","RRtlde","h0","RR"))
    dimnames(detail$bstar) <- list(stat.nms, c("bstartlde","bstar"))
    names(detail$ints) <- ints.nms
    names(detail$dbls) <- dbls.nms

    detail$Var.nlook <- with(detail, Var[tidx,])
    detail$Eta.nlook <- with(detail, Eta[tidx,])
    detail$RR.nlook <- with(detail, RR[tidx,])
    detail$E.NT.nlook <- with(detail, E.NT[tidx,])
    
    out <- list(dPower = Pwr, dErrorII = dErrII, detail = detail,call=.call.)
    class(out) <- "PwrGSD"
    out
}

"print.PwrGSD" <- function (x, ...)
{
    print(x$call)
    nlook <- x$detail$ints["nlook"]
    nstat <- x$detail$ints["nstat"]
    pwr <- c(x$dPower %*% rep(1, nlook))
    dErrII <- 0*x$dPower
    do.futility <- (x$detail$ints["do.futility"]==1)
    if(do.futility) dErrII <- x$dErrorII
    ed <- c(apply(x$dPower + dErrII, 1, 
	  FUN = function(a, b, do.futility) 
		{
        	     n <- length(a)
		     y <- a
        	     if(!do.futility) y <- c(a[-n], 1 - sum(a[-n]))
        	     sum(y * b)
    		}, b = eval(x$detail$pttlook), do.futility=do.futility))
    eII <- c(dErrII %*% rep(1, nlook))
    stat.nms <- dimnames(x$dPower)[[1]]
    ans <- cbind(pwr, eII, ed)
    dimnames(ans) <- list(stat.nms, c("Power", "Type.II.Err", "ExpDur"))
    print(ans)
    invisible(x)
}

"summary.PwrGSD" <- 
function (object, ...) 
{
    nlook <- object$detail$ints["nlook"]
    nstat <- object$detail$ints["nstat"]
    pwr <- c(object$dPower %*% rep(1, nlook))
    dErrII <- 0 * object$dPower
    do.futility <- (object$detail$ints["do.futility"] == 1)
    if (do.futility) 
        dErrII <- object$dErrorII
    ed <- c(apply(object$dPower + dErrII, 1, FUN = function(a, b, 
        do.futility) {
        n <- length(a)
        y <- a
        if (!do.futility) y <- c(a[-n], 1 - sum(a[-n]))
        sum(y * b)
    }, b = eval(object$detail$pttlook), do.futility = do.futility))
    eII <- c(dErrII %*% rep(1, nlook))
    stat.nms <- dimnames(object$dPower)[[1]]
    ans <- cbind(pwr, eII, ed)
    dimnames(ans) <- list(stat.nms, c("Power", "Type.II.Err", 
        "ExpDur"))
    out <- object
    out$Tbl <- ans
    out
}

"print.cpd.PwrGSD" <- 
function(x, ...)
{
        print(x$call)
        sapply(x$Elements, FUN = print.PwrGSD)
}

"summary.cpd.PwrGSD" <- 
function(object, ...)
{
        sapply(object$Elements, FUN = summary.PwrGSD)
}

"cpd.PwrGSD" <-
function(descr)
{
  res <- list()
  d.descr <- dim(descr)
  n <- d.descr[1]
  p <- d.descr[2]
  length(res) <- n
  if(!("index" %in% names(descr)))
    descr$index <- 1:n
  
  ans <- list(date=format(Sys.time(), "%Y.%m.%d.%H:%M:%S"), Elements=res, descr=descr)
  class(ans) <- "cpd.PwrGSD"
  ans
}

"Elements" <-
function(object, subset, na.action=na.pass)
{
  m <- .call. <- match.call()
  descr <- object$descr
  n <- dim(descr)[1]
  index <- rep(TRUE, n)
  if(!missing(subset))
    index <- with(descr, eval(.call.$subset))
  
  if(length(index>1)){
    ans <- list(date=object$date, Elements=object$Elements[index], descr=descr[index, ])
    class(ans) <- "cpd.PwrGSD"
  }
  else{
    ans <- object$Elements[index]
    class(ans) <- "PwrGSD"
  }
  ans
}

"as.data.frame.boundaries" <-
function(x, row.names, optional, ...)
{
  nr <- nrow(x$table)
  as.data.frame(cbind(frac=x$frac[1:nr], frac.ii=x$frac.ii[1:nr], x$table[,-1], drift=x$drift[1:nr]))
}

"as.data.frame.cpd.PwrGSD" <- 
function(x, row.names, optional, ...) 
{
    N <- length(x$Elements)
    idx <- x$descr$index
    det.1 <- x$Elements[[1]]$detail
    det.a <- c(det.1$ntrial, det.1$ints["nlook"], det.1$ints["nstat"])
    names(det.a) <- c("ntrial", "nlook", "nstat")

    nstat <- det.a["nstat"]
    nlook <- det.a["nlook"]
    nms <- outer(1:nstat, 1:nlook, FUN = function(x, y) {
        "nstat[" %,% x %,% "].nlook[" %,% y %,% "]"
    })
    ans <- NULL
    for (k in 1:nstat) {
      tmp.b <- sapply(x$Elements, FUN = function(x, stat) as.boundaries(x, 
                      stat)$table, stat = k)
      d.tmp.b <- dim(tmp.b)
      col.gsbt <- d.tmp.b[1]/nlook
      tmp.b <- matrix(aperm(array(tmp.b, c(nlook, col.gsbt, N))[,c(1,2,5),],c(3,1,2)), N, 3*nlook)
      nstatis <- "stat[" %,% k %,% "]"
      col.nms <- c(t(outer(c("frac.","bFut.","bEff."), nstatis %,% ".nlook[" %,% 1:nlook %,% "]",
                           FUN = "%,%")))
      dimnames(tmp.b) <- list(idx, col.nms)
      tmp.p <- t(sapply(x$Elements, FUN=function(x, stat)c(x$dE[stat,],x$dP[stat,]), stat=k))
      col.nms <- c(t(outer(c("dErrII.","dP."),nstatis %,% ".nlook[" %,% 1:nlook %,% "]",FUN="%,%")))
      dimnames(tmp.p) <- list(idx, col.nms)
      ans <- cbind(ans, tmp.b, tmp.p)
    }

    ans <- as.data.frame(ans)
    ans$index <- x$descr$index
    ans <- merge(x$descr, ans, by = "index")
    ans <- mystack(ans, fu.vars=c("frac","bEff","bFut","dP","dE"))
    ans$stat <- floor(ans$fu/nlook)+1
    ans$nlook <- (ans$fu %% nlook) + 1
    ans$fu <- NULL
    p <- dim(x$descr)[2]
    pp <- dim(ans)[2]
    ans <- ans[,c(1:p, pp-1, pp, (p+1):(pp-2))]
    attr(ans, "detail") <- det.a
    ans
}

"Power" <- 
function (object, subset, nlook.ind=NULL) 
{
    object.df <- as.data.frame(object)
    
    m <- .call. <- match.call()
    m1 <- .call.
    m1[[1]] <- as.name("model.frame")
    m1$object <- m1$nlook.ind <- NULL
    m1$data <- as.name("object.df")
    m1$formula <- ~.
    m1 <- eval(m1)
    nc.m1 <- ncol(m1)
    nr.m1 <- nrow(m1)
    descr <- with(m1, m1[nlook==1,1:(nc.m1-6)])    
    n.look <- length(unique(m1$nlook))
    n.stat <- length(unique(m1$stat))
    n.scen <- nr.m1/(n.look*n.stat)
    Pow <- NULL
    if(missing(nlook.ind)) nlook.ind <- 1:n.look
    n.ind <- length(nlook.ind)
    for (k in 1:n.scen)
      for (j in 1:n.stat)
      {
        Pow <- rbind(Pow, rep(1, n.ind) %*% as.matrix(m1[n.look*n.stat*
          (k - 1) + n.look*(j-1) + nlook.ind, c("dP", "dE")]))
      }
    dimnames(Pow)[[2]] <- c("Power", "Type II error")
    ans <- list(table = cbind(descr, Pow), call = .call.)
    ans
}

"PwrGSDcall" <- 
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         sided =c("2>", "2<", "1>", "1<"),method=c("S","A"),accru,accrat,tlook,
         tcut0 = NULL,h0 = NULL,s0 = NULL,tcut1 = NULL,rhaz = NULL,
         h1 = NULL,s1 = NULL,tcutc0 = NULL,hc0 = NULL,sc0 = NULL,tcutc1 = NULL,hc1 = NULL,
         sc1 = NULL,tcutd0A = NULL,hd0A = NULL,sd0A = NULL,tcutd0B = NULL,hd0B = NULL,sd0B = NULL,
         tcutd1A = NULL,hd1A = NULL,sd1A = NULL,tcutd1B = NULL,hd1B = NULL,sd1B = NULL,
         tcutx0A = NULL,hx0A = NULL,sx0A = NULL,tcutx0B = NULL,hx0B = NULL,sx0B = NULL,
         tcutx1A = NULL,hx1A = NULL,sx1A = NULL,tcutx1B = NULL,hx1B = NULL,sx1B = NULL,
         noncompliance = c("none","crossover","mixed","user"),gradual = FALSE,
         WtFun = c("FH","SFH","Ramp"),ppar = cbind(c(0,0)),
         Spend.Info=c("Variance","Events","Hybrid(k)","Calendar"),
         RR.Futility = NULL,qProp.one.or.Q = c("one","Q"),Nsim=NULL,
         detail = FALSE,StatType = c("WLR","ISD"))
{
  match.call()
}

"GrpSeqBnds"<-
function (EfficacyBoundary = LanDemets(alpha=0.05, spending=ObrienFleming),
          FutilityBoundary = LanDemets(alpha=0.10, spending=ObrienFleming),
          NonBindingFutility=TRUE, frac, frac.ii = NULL, drift = NULL)
{
    .call. <- match.call()
    
    is.frac.ii <- TRUE
    if (missing(frac.ii)) {
        frac.ii <- frac
        is.frac.ii <- FALSE
    }
    if (is.frac.ii && length(frac.ii) != length(frac)) 
        stop("Lengths of 'frac' and 'frac.ii' must agree")
    
    normut <- 8.20953615160139
    psimin <- 1e-15
    nlooks <- length(frac)

    call.PBS <- .call.
    call.PBS[[1]] <- as.name("ParseBoundarySelection")
    call.PBS$frac <- call.PBS$frac.ii <- NULL
    call.PBS$check.drift <- TRUE
    call.PBS$n.looks <- nlooks
    eval.PBS <- eval(call.PBS)
    do.efficacy <- TRUE
    do.futility <- !is.null(eval.PBS$frontend$FutilityBoundary)
    is.drift <- !missing(drift)
    nbf <- do.futility && NonBindingFutility
    if(!do.futility)
      drift <- rep(0, nlooks)

    Alpha.Efficacy <- eval.PBS$backend$Alpha.Efficacy
    Alpha.Futility <- eval.PBS$backend$Alpha.Futility
    if(is.null(Alpha.Futility)) Alpha.Futility <- 0
    nbnd.e <- eval.PBS$backend$nbnd.e
    nbnd.f <- eval.PBS$backend$nbnd.f
    nsf.e <- eval.PBS$backend$nsf.e
    nsf.f <- eval.PBS$backend$nsf.f
    rho.Efficacy <- eval.PBS$backend$rho.Efficacy
    rho.Futility <- eval.PBS$backend$rho.Futility
    b.Haybittle.e <- eval.PBS$backend$b.Haybittle.e
    b.Haybittle.f <- eval.PBS$backend$b.Haybittle.f
    drift.end <- eval.PBS$backend$drift.end
    be.end <- eval.PBS$backend$be.end
    prob.e <- eval.PBS$backend$prob.e
    prob.f <- eval.PBS$backend$prob.f
    my.Efficacy <- eval.PBS$backend$my.Efficacy
    my.Futility <- eval.PBS$backend$my.Futility
    is.myE <- !all(my.Efficacy==0)
    is.myF <- !all(my.Futility==0)

    nunq <- nlooks
    glegx <- glegx24
    glegw <- glegw24
    ngqnodes <- length(glegx)

    b.e <- alpha.e <- rep(0, nlooks)
    b.f <- alpha.f <- rep(0, nlooks)
    l <- 1
    l.act.e <- 1
    l.act.f <- 1
    fracold <- rep(0, 2)
    fracold.ii <- rep(0, 2)
    fracnew <- frac[1]
    fracnew.ii <- frac.ii[1]

    sc.drift.factor <- 1
    mu.end <- drift.end
  
    if (nbnd.e[1] != 2)
        bold.e <- normut
    if (nbnd.e[1] == 2)
        bold.e <- b.Haybittle.e[1]
    if (nbnd.f[1] != 2)
        bold.f <- -normut
    if (nbnd.f[1] == 2) 
        bold.f <- b.Haybittle.f[1]
    y.e <- tmp.e <- rep(0, ngqnodes)
    y.f <- tmp.f <- rep(0, ngqnodes)
    x.e <- intgrndx.e <- rep(0, ngqnodes)
    x.f <- intgrndx.f <- rep(0, ngqnodes)
    if (nbnd.e[1] == 3)
        my.Efficacy <- rep(0, nlooks)
    if (nbnd.f[1] == 3) 
        my.Futility <- rep(0, nlooks)
    while (l <= nlooks) {
      mu <- drift[l]
        if (nbnd.e[l] == 3) {
            my.Efficacy[l] <-
              .C("StCu2Bnds",
                 pmu = double(2),
                 pfrac = as.double(fracnew.ii),
                 pzcrit = as.double(be.end),
                 prho = as.double(prob.e[l]),
                 pef = as.integer(0),
                 b = double(1),
                 PACKAGE = "PwrGSD")$b
            is.myE <- TRUE
            nsf.e[l] <- 1
        }
        if (nbnd.f[l] == 3) {
            my.Futility[l] <-
              .C("StCu2Bnds",
                 pmu = as.double(c(mu,mu.end)*sc.drift.factor),
                 pfrac = as.double(fracnew.ii),
                 pzcrit = as.double(be.end),
                 prho = as.double(prob.f[l]), 
                 pef = as.integer(1),
                 b = double(1),
                 PACKAGE = "PwrGSD")$b
            is.myF <- TRUE
            nsf.f[l] <- 1
        }
        bnew <- c(normut, -normut)
        if (nbnd.e[l] == 3 || nbnd.e[l] == 4) 
            bnew[1] <- my.Efficacy[l]
        if ((nbnd.f[l] == 3 || nbnd.f[l] == 4) && (1 - fracnew.ii >= 
            1e-06)) 
            bnew[2] <- my.Futility[l]

        ans <- .C("grpseqbnds",
                  do.futility = as.integer(do.futility),
                  nbf = as.integer(nbf),
                  nbnd = as.integer(c(nbnd.e[l], nbnd.f[l])),
                  nsf = as.integer(c(nsf.e[l],nsf.f[l])),
                  rho = as.double(c(rho.Efficacy[l], rho.Futility[l])), 
                  pnlook = as.integer(c(l.act.e, l.act.f)),
                  palphtot = as.double(c(Alpha.Efficacy, Alpha.Futility)),
                  palpha = double(2),
                  psimin = as.double(psimin), 
                  dlact = integer(2),
                  pfracold = as.double(fracold), 
                  pfracnew = as.double(fracnew),
                  pfracold.ii = as.double(fracold.ii), 
                  pfracnew.ii = as.double(fracnew.ii),
                  x = as.double(c(x.e, x.f)),
                  y = as.double(c(y.e, y.f)),
                  tmp = as.double(c(tmp.e, tmp.f)),
                  intgrndx = as.double(c(intgrndx.e, intgrndx.f)), 
                  gqxw = as.double(c(glegx, glegw)),
                  pngqnodes = as.integer(ngqnodes), 
                  mu = as.double(mu),
                  bold = as.double(c(bold.e, bold.f)), 
                  bnew = as.double(bnew),
                  mybounds = as.integer(c(is.myE, is.myF)),
                  PACKAGE = "PwrGSD")
        dlact <- ans$dlact
        if (dlact[1] == 1)
        {
            if (nbnd.e[l] == 1 || nbnd.e[l] == 3 || nbnd.e[l] == 4) 
                b.e[l] <- bold.e <- ifelse(is.myE, bnew[1], ans$bnew[1])
            if (nbnd.e[l] == 2)
            {
                if (1 - fracnew.ii >= 1e-06) {
                  b.e[l] <- bold.e <- b.Haybittle.e[l]
                  Alpha.Efficacy <- Alpha.Efficacy - ans$palpha[1]
                }
              else b.e[l] <- ans$bnew[1]
            }
            alpha.e[l] <- ans$palpha[1]
            x.e <- ans$x[1:ngqnodes]
            intgrndx.e <- ans$intgrndx[1:ngqnodes]
            fracold[1] <- ans$pfracnew
            fracold.ii[1] <- ans$pfracnew.ii
            l.act.e <- l.act.e + 1
        }
        else
        {
            b.e[l] <- bold.e <- ifelse(is.myE, bnew[1], normut)
            alpha.e[l] <- ifelse(is.myE, ans$palpha[1], psimin)
            fracold.ii[1] <- ans$pfracnew.ii
        }
        if (do.futility && dlact[2] == 1) {
            if (nbnd.f[l] == 1 || nbnd.f[l] == 3 || nbnd.f[l] == 4) 
                b.f[l] <- bold.f <- ifelse((nbnd.f[l] == 3|| nbnd.f[l] == 4) && (1 - fracnew.ii >= 1e-06),
                                           bnew[2], ans$bnew[2])
            if (nbnd.f[l] == 2) {
                if (l < nlooks) {
                  b.f[l] <- bold.f <- b.Haybittle.f[l]
                  Alpha.Futility <- Alpha.Futility - ans$palpha[2]
                }
                else b.f[l] <- ans$bnew[2]
            }
            alpha.f[l] <- ans$palpha[2]
            x.f <- ans$x[ngqnodes + (1:ngqnodes)]
            intgrndx.f <- ans$intgrndx[ngqnodes + (1:ngqnodes)]
            fracold[2] <- ans$pfracnew
            fracold.ii[2] <- ans$pfracnew.ii
            l.act.f <- l.act.f + 1
        }
        if (do.futility && dlact[2] == 0) {
            b.f[l] <- bold.f <- ifelse(is.myF || (l<nlooks), bnew[2], 
                -normut)
            alpha.f[l] <- ifelse(is.myF, ans$palpha[2], psimin)
            fracold.ii[2] <- ans$pfracnew.ii
        }
        l <- l + 1
        fracnew <- frac[l]
        fracnew.ii <- frac.ii[l]
        mu <- drift[l]
    }
    out <- frac.ii
    nms <- "frac"
    if (do.futility) {
        out <- cbind(out, b.f, alpha.f, cumsum(alpha.f))
        nms <- c(nms, "b.f", "alpha.f", "cum-alpha.f")
    }
    out <- cbind(out, b.e, alpha.e, cumsum(alpha.e))
    nms <- c(nms, "b.e", "alpha.e", "cum-alpha.e")
    dimnames(out) <- list(1:nlooks, nms)
    ans <- list(table = out, frac=frac, frac.ii=frac.ii, drift=drift, call = .call.)
    class(ans) <- "boundaries"
    ans
}

"ParseBoundarySelection" <- 
function(EfficacyBoundary,  n.looks, check.drift=FALSE, NonBindingFutility=NULL,
         FutilityBoundary=NULL, drift=NULL)
{
  m <- m.sv <- match.call()
  if(missing(NonBindingFutility)) NonBindingFutility <- TRUE
  do.efficacy <- !missing(EfficacyBoundary)
  if(!do.efficacy)
    stop("The argument 'EfficacyBoundary' which specifies the Efficacy Boundary construction " %,%
         "method, must be specified")
  do.futility <- !missing(FutilityBoundary)
  is.drift <- !missing(drift)
  n.drift <- 0
  if(is.drift) n.drift <- length(eval(drift))
  if(do.futility && check.drift)
  {
    if(!is.drift)
      stop("In order to construct a futility boundary, you must specify the 'drift' argument. " %,%
           "This is the expected value under the design alternative of the statistic, normalized " %,%
           "to have variance equal to information fraction")
    if(is.drift && n.drift != n.looks)
      stop("The argument 'drift' must be equal in length to that of the argument 'frac'")        
  }

  EB <- eval(m$EfficacyBoundary)
  n.types <- 1
  mode.EB <- mode(EB)
  if(mode.EB=="numeric")
  {
    m <- eval(parse(text="as.call(expression(c," %,% paste(EB, collapse=", ") %,% "))"))
    .eb. <- list(type="numeric", value=EB, len=length(EB), call=m)
    class(.eb.) <- "boundary.construction.method"
  }
  if(mode.EB!="numeric")
  {
    if(length(EB[[1]])>1)
    {
      n.types <- length(EB)
      .eb. <- list()
      for(k in 1:n.types)
      {
        mode.k <- mode(EB[[k]])
        if(mode.k=="numeric")
        {
          m <- eval(parse(text="as.call(expression(c," %,% paste(EB[[k]], collapse=", ") %,% "))"))
          .eb.[[k]] <- list(type="numeric", value=EB[[k]], len=length(EB[[k]]), call=m)
          class(.eb.[[k]]) <- "boundary.construction.method"
        }
        if(mode.k!="numeric")
        {
          msg <- "Boundary type must be one of the following: 'LanDemets()'," %,%
                  "'Haybittle()' or 'SC()'"
          if(mode(EB[[k]])!="call") stop(msg)
          type.EB.k <- as.character(EB[[k]][[1]])
          if(!(type.EB.k %in% c("LanDemets", "Haybittle", "SC"))) stop(msg)
          .eb.[[k]] <- EB[[k]]
          .eb.[[k]]$type <- type.EB.k
        }
      }
      if((is.null(.eb.[[k]]$from) || is.null(.eb.[[k]]$to)) && is.null(.eb.[[k]]$len))
          stop("Split boundary definition require the specifification of either" %,% 
               "the pair of arguments 'from' and 'to' or the 'len' argument")
    }
    else
    {
      msg <- "Boundary type must be one of the following: 'LanDemets()'," %,%
             "'Haybittle()' or 'SC()'"
      if(mode(m$EfficacyBoundary)!="call") stop(msg)
      type.EB <- as.character(m$EfficacyBoundary[[1]])
      if(!(type.EB %in% c("LanDemets", "Haybittle", "SC"))) stop(msg)
      .eb. <- EB
      .eb.$type <- type.EB
    }
  }
  if(n.types==1) .eb. <- list(.eb.)
  EB <- .eb.

  m <- m.sv
  
  if(do.futility)
  {
    FB <- eval(m$FutilityBoundary)
    n.types <- 1
    mode.FB <- mode(FB)
    if(mode.FB=="numeric")
    {
      m <- eval(parse(text="as.call(expression(c," %,% paste(FB, collapse=", ") %,% "))"))
      .fb. <- list(type="numeric", value=FB, len=length(FB), call=m)
      class(.fb.) <- "boundary.construction.method"
    }
    if(mode.FB!="numeric")
    {
      if(length(FB[[1]])>1)
      {
        n.types <- length(FB)
        .fb. <- list()
        for(k in 1:n.types)
        {
          mode.k <- mode(FB[[k]])
          if(mode.k=="numeric")
          {
            m <- eval(parse(text="as.call(expression(c," %,% paste(FB[[k]], collapse=", ") %,% "))"))
            .fb.[[k]] <- list(type="numeric", value=FB[[k]], len=length(FB[[k]]), call=m)
            class(.fb.[[k]]) <- "boundary.construction.method"
          }
          if(mode.k!="numeric")
          {
            type.FB.k <- as.character(FB[[k]][[1]])
            if(!(type.FB.k %in% c("LanDemets", "Haybittle", "SC")))
            stop("Boundary type must be one of the following: 'LanDemets()'," %,%
                 "'Haybittle()' or 'SC()'")
            .fb.[[k]] <- FB[[k]]
            .fb.[[k]]$type <- type.FB.k
          }
          if((is.null(.fb.[[k]]$from) || is.null(.fb.[[k]]$to)) && is.null(.fb.[[k]]$len))
            stop("Split boundary definition require the specifification of either" %,% 
                 "the pair of arguments 'from' and 'to' or the 'len' argument")
        }
      }
      else
      {
        type.FB <- as.character(m$FutilityBoundary[[1]])
        if(!(type.FB %in% c("LanDemets", "Haybittle", "SC")))
          stop("Boundary type must be one of the following: 'LanDemets()'," %,%
               "'Haybittle()' or 'SC()'")
        .fb. <- FB
        .fb.$type <- type.FB
      }
    }
    if(n.types==1) .fb. <- list(.fb.)
    FB <- .fb.
  }
  ans <- list(frontend=list(EfficacyBoundary=EB), backend=list())
  if(do.futility) ans$frontend$FutilityBoundary <- FB

  ans$backend$Alpha.Efficacy <- 0
  ae <- ans$frontend$EfficacyBoundary[[1]]$alpha
  if(!is.null(ae)) ans$backend$Alpha.Efficacy <- ae
  n.EBtypes <- length(ans$frontend$EfficacyBoundary)
  if(n.EBtypes==1)
  {
    ans$frontend$EfficacyBoundary[[1]]$from <- 1
    ans$frontend$EfficacyBoundary[[1]]$to <- n.looks
  }
  nbnd.e <- nsf.e <- rho.Efficacy <- prob.e <- my.Efficacy <- b.Haybittle.e <- rep(0, n.looks)
  be.end <- 0

  coveredEB <- NULL
  to.k <- 0
  for (k in 1:n.EBtypes)
  {
    EB.k <- ans$frontend$EfficacyBoundary[[k]]
    len.k <- EB.k$len
    if(EB.k$type=="numeric")
    {
      from.k <- to.k + 1
      to.k <- to.k + len.k
      nsf.e[from.k:to.k] <- 1
      my.Efficacy[from.k:to.k] <- EB.k$value
    }
    if(EB.k$type!="numeric")
    {
      from.k <- EB.k$from
      to.k <- EB.k$to
    }
    disjoint <- length(intersect(coveredEB, from.k:to.k))==0
    if(!disjoint) stop("There is overlap in your specification of split boundaries")
    coveredEB <- c(coveredEB, from.k:to.k)
    nbnd.e[from.k:to.k] <- grep(EB.k$type, c("LanDemets", "Haybittle", "SC", "numeric"))
    if(EB.k$type=="LanDemets")
    {
      nsf.e[from.k:to.k] <- grep(EB.k$spending$type, c("ObrienFleming", "Pocock", "Pow"))
      if(EB.k$spending$type=="Pow")
        rho.Efficacy[from.k:to.k] <- EB.k$spending$rho
    }
    if(EB.k$type=="Haybittle")
    {
      b.Haybittle.e[from.k:to.k] <- EB.k$b.Haybittle
    }
    if(EB.k$type=="SC")
    {
       nsf.e[from.k:to.k] <- 1
       prob.e[from.k:to.k] <- EB.k$prob
       be.end <- EB.k$be.end
    }
  }
  coveredEB <- sort(coveredEB)
  if(sum(abs(coveredEB - (1:n.looks)))>0)
    stop("Your split Efficacy Boundary definition does not specifiy all " %,% n.looks %,% " analyses")
  ans$backend$nbnd.e <- nbnd.e
  ans$backend$nsf.e <- nsf.e
  ans$backend$rho.Efficacy <- rho.Efficacy
  ans$backend$b.Haybittle.e <- b.Haybittle.e
  ans$backend$be.end <- be.end
  ans$backend$prob.e <- prob.e
  ans$backend$my.Efficacy <- my.Efficacy

  nbnd.f <- nsf.f <- rho.Futility <- prob.f <- my.Futility <- b.Haybittle.f <- rep(0, n.looks)
  be.end <- drift.end <- 0
  
  ans$backend$Alpha.Futility <- 0
  if(do.futility)
  {
    ans$backend$Alpha.Futility <- ans$frontend$FutilityBoundary[[1]]$alpha
    n.FBtypes <- length(ans$frontend$FutilityBoundary)
    if(n.FBtypes==1)
    {
      ans$frontend$FutilityBoundary[[1]]$from <- 1
      ans$frontend$FutilityBoundary[[1]]$to <- n.looks
    }

    coveredFB <- NULL
    to.k <- 0
    for (k in 1:n.FBtypes)
    {
      FB.k <- ans$frontend$FutilityBoundary[[k]]
      len.k <- FB.k$len
      if(FB.k$type=="numeric")
      {
        from.k <- to.k + 1
        to.k <- to.k + len.k
        nsf.f[from.k:to.k] <- 1
        my.Futility[from.k:to.k] <- FB.k$value
      }
      if(FB.k$type!="numeric")
      {
        from.k <- FB.k$from
        to.k <- FB.k$to
      }

      disjoint <- length(intersect(coveredFB, from.k:to.k))==0
      if(!disjoint) stop("There is overlap in your specification of split boundaries")
      coveredFB <- c(coveredFB, from.k:to.k)
      nbnd.f[from.k:to.k] <- grep(FB.k$type, c("LanDemets", "Haybittle", "SC", "numeric"))
      if(FB.k$type=="LanDemets")
      {
        nsf.f[from.k:to.k] <- grep(FB.k$spending$type, c("ObrienFleming", "Pocock", "Pow"))
        if(FB.k$spending$type=="Pow")
          rho.Futility[from.k:to.k] <- FB.k$spending$rho
      }
      if(FB.k$type=="Haybittle")
      {
        b.Haybittle.f[from.k:to.k] <- FB.k$b.Haybittle
      }
      if(FB.k$type=="SC")
      {
         nsf.f[from.k:to.k] <- 1
         prob.f[from.k:to.k] <- FB.k$prob
         if(be.end==0) be.end <- FB.k$be.end
         if(drift.end==0) drift.end <- FB.k$drift.end
      }      
    }
    coveredFB <- sort(coveredFB)
    if(sum(abs(coveredFB - (1:n.looks)))>0)
      stop("Your split Futility Boundary definition does not specifiy all " %,% n.looks %,% " analyses")
  }
  ans$backend$nbnd.f <- nbnd.f
  ans$backend$nsf.f <- nsf.f
  ans$backend$rho.Futility <- rho.Futility
  ans$backend$b.Haybittle.f <- b.Haybittle.f
  ans$backend$drift.end <- drift.end
  if(ans$backend$be.end==0) ans$backend$be.end <- be.end
  ans$backend$prob.f <- prob.f
  ans$backend$my.Futility <- my.Futility
  #-------------------------------#
  # check for valid combinations: #
  #-------------------------------#
  #  0=none                       #
  #  1=LD                         #
  #  2=H                          #
  #  3=SC                         #
  #  4=U                          #
  #-------------------------------#
  nbnd <- 25*NonBindingFutility + 5*nbnd.e + nbnd.f
  ans$backend$nbnd.pair <- nbnd
  
  # No Efficacy specified--already checked above, but just for completenes
  if(any(nbnd %in% c(0:4, 25:29)))
  {
    indx <- which(nbnd %in% c(0:4, 25:29))
    indx.str <- paste(indx, collapse=",")
    stop("The argument 'EfficacyBoundary' which specifies the Efficacy Boundary construction " %,%
         "method, must be specified (" %,% indx.str %,% ")")
  }
  # We don't do Haybittle Futility
  if(any((nbnd-2) %% 5 ==0))
  {
    indx <- which((nbnd-2) %% 5 ==0)
    indx.str <- paste(indx, collapse=",")
    stop("Haybittle futility boundary not currently supported--and its doubtful whether or not it " %,%
         "makes sense (" %,% indx.str %,% ")")
  }
  ans
}

"LanDemets" <-
function(alpha, spending, from=NULL, to=NULL)
{
  if(missing(alpha)||missing(spending))
    stop("'alpha' and 'spending' are required arguments")
  m <- match.call()
  sf <- as.call(expression(fff))
  n.sf <- length(m$spending)
  if(n.sf==1) sf[[1]] <- m$spending
  if(n.sf>1)
    for(k in 1:n.sf)
      sf[[k]] <- m$spending[[k]]
  sf.name <- as.character(sf[[1]])
  if(!(sf.name %in% c("ObrienFleming", "Pow", "Pocock")))
    stop("Argument 'spending' must be 'ObrienFleming', 'Pow' or 'Pocock'")
  sf <- eval(sf)
  ans <- list(type="LanDemets", alpha=alpha, spending=sf, from=from, to=to, call=m)
  class(ans) <- "boundary.construction.method"
  ans
}
  
"Haybittle" <-
function(alpha, b.Haybittle, from=NULL, to=NULL)
{
  m <- match.call()
  if(missing(alpha)||missing(b.Haybittle))
    stop("'alpha' and 'b.Haybittle' are required arguments")
  ans <- list(type="Haybittle",alpha=alpha, b.Haybittle=b.Haybittle, from=from, to=to, call=m)
  class(ans) <- "boundary.construction.method"
  ans
}

"SC" <-
function(be.end, prob, drift.end=NULL, from=NULL, to=NULL)
{
  m <- match.call()
  if(missing(be.end)||missing(prob))
    stop("'be.end' and 'prob' are required arguments")
  ans <- list(type="SC", be.end=be.end, prob=prob, from=from, to=to, call=m)
  if(!missing(drift.end)) ans$drift.end <- drift.end
  class(ans) <- "boundary.construction.method"
  ans
}

"SCtoBdry" <- 
function(prob, frac, be.end, drift=NULL, drift.end=NULL)
{
  ef <- !missing(drift)
  if(ef && missing(drift.end))
    stop("You must specify both arguments 'drift' and 'drift.end' for a futility boundary " %,%
         "using stochastic curtailment\n")
  
  nlook <- length(frac)
  if(ef)
  {
    if(length(drift)>0 && length(drift) < nlook)
      stop("Argument 'drift' must be either unspecified (efficacy) or " %,%
           "of the same length as 'frac': length(drift)=" %,% length(drift) %,% ", " %,%
           "length(frac)=" %,% length(frac) %,% ".\n")
  }
  if(!ef)
  {
    drift <- rep(0, nlook)
    drift.end <- 0
  }
  
  if(length(prob)!=1 && length(prob) < nlook)
    stop("Argument 'prob' must be either of length 1 or of the same length as " %,%
         "'frac': length(prob)=" %,% length(prob) %,% ", length(frac)=" %,% length(frac) %,% ".\n")

  
  if(length(prob)==1) prob <- rep(prob, nlook)

  b <- NULL
  frac.old <- 0
  for(l in 1:nlook)
  {
    frac.new <- frac[l]
    b.l <- .C("StCu2Bnds",
              pmu = as.double(c(drift[l], drift.end)),
              pfrac = as.double(frac.new),
              pzcrit = as.double(be.end),
              prho = as.double(prob[l]), 
              pef = as.integer(ef),
              b = double(1),
              PACKAGE = "PwrGSD")$b
    b <- c(b, b.l)
    frac.old <- frac.new
  }
  b
}


"ObrienFleming" <-
function()
{
  m <- match.call()
  ans <- list(type="ObrienFleming", call=m)
  class(ans) <- "spending.function"
  ans
}

"Pow" <-
function(rho)
{
  m <- match.call()
  if(missing(rho))
    stop("'rho' is a required argument")
  ans <- list(type="Pow", rho=rho, call=m)
  class(ans) <- "spending.function"
  ans
}

"Pocock" <-
function()
{
  m <- match.call()
  ans <- list(type="Pocock", call=m)
  class(ans) <- "spending.function"
  ans
}

"print.boundary.construction.method" <-
function(x, ...)
{
  print(x$call)
}

"print.spending.function" <-
function(x, ...)
{
  print(x$call)
}

"print.boundaries" <- function(x, ...)
{
  dots <- c(...)
  nms.dots <- names(dots)
  nocall <- FALSE
  if("nocall" %in% nms.dots)
    nocall <- dots["nocall"]  
  ans <- x$table
  if(!nocall){
    cat("call:\n")
    print(x$call)
  }
  print(ans)
  invisible(ans)
}

"plot.boundaries" <- 
function (x, yrng=NULL, ...) 
{
    .call. <- match.call(expand.dots=TRUE)
    m <- length(.call.) - 2
    nms <- names(.call.)[-(1:2)]
    xtra <- list()
    if(m>0) for(k in 1:m) xtra[[names(.call.)[k]]] <- .call.[[2+k]]
    names(xtra) <- nms
    tbl <- x$table
    d.tbl <- dim(tbl)
    is.fut <- (d.tbl[2] == 7)
    x <- tbl[, "frac"]
    eff <- tbl[, "b.e"]
    fut <- NULL
    if (is.fut) 
        fut <- tbl[, "b.f"]
    is.yrng <- TRUE
    if(missing(yrng))
    {
      yrng <- c(min(c(eff, fut)), max(c(eff, fut)))
      is.yrng <- FALSE
    }
    plot.cmd <- as.call(expression(plot))
    plot.cmd$x <- c(0, 1)
    plot.cmd$y <- yrng
    plot.cmd$type <- "n"
    plot.cmd$xlab <- "Information Fraction"
    plot.cmd$ylab <- "Normalized Statistic"
    lines.cmd.e <- lines.cmd.f <- as.call(expression(lines))
    lines.cmd.e$x <- lines.cmd.f$x <- as.name("x")
    lines.cmd.e$y <- as.name("eff")
    lines.cmd.f$y <- as.name("fut")
    if(m>0){
      l.pl <- l.li <- 1
      for(k in 1:m)
      {
        if(nms[k] %in% c("type", "main", "sub", "xlab", "ylab", "asp"))
        {
          plot.cmd[[6+l.pl]] <- xtra[[k]]
          names(plot.cmd)[6+l.pl] <- nms[k]
          l.pl <- l.pl + 1
        }
        if(nms[k] %in% c("type", "lty", "lwd", "pch", "lend", "ljoin", "lmitre"))
        {
          lines.cmd.e[[3+l.li]] <- lines.cmd.f[[3+l.li]] <- xtra[[k]]
          names(lines.cmd.e)[3+l.li] <- names(lines.cmd.f)[3+l.li] <- nms[k]
          l.li <- l.li + 1
        }
        if(nms[k] == "col")
        {
          if(!is.character(xtra[[k]]))
            stop("Specify color as a character string, please")
          col.e <- xtra[[k]]
          col.f <- c(c(16^4,16^2,1)%*%(cbind(c(255,255,255)) -col2rgb(col.e)))
          class(col.f) <- "hexmode"
          col.f <- "#" %,% col.f
          lines.cmd.e[[3+l.li]] <- col.e
          names(lines.cmd.e)[3+l.li] <- "col"
          lines.cmd.f[[3+l.li]] <- col.f
          names(lines.cmd.f)[3+l.li] <- "col"
        }
      }
    }
    eval(plot.cmd)
    eval(lines.cmd.e)
    points(x, eff)
    if (is.fut) {
        eval(lines.cmd.f)
        points(x, fut)
    }
    invisible(x)
}

"SimGSB" <-
function(object, nsim=1e+05, ...)
{
  UseMethod("SimGSB")
}

"SimGSB.boundaries" <-
function (object, nsim = 1e+05, ...)
{
    .call. <- match.call()
    tab <- object$table
    nms <- dimnames(tab)[[2]]
    .fr. <- eval(object$call$frac)
    .fr.ii <- eval(object$call$frac.ii)
    if(is.null(.fr.ii)) .fr.ii <- .fr.
    n <- length(.fr.)
    do.fu <- length(grep("b.f", nms)) > 0
    a <- rep(-Inf, n)
    al.f <- rep(0, n)
    if (do.fu) {
        .dr. <- eval(object$call$drift)
        a <- .fr.^0.5 * tab[, "b.f"]
        al.f <- tab[, "cum-alpha.f"]
    }
    b <- .fr.^0.5 * tab[, "b.e"]
    al.e <- tab[, "cum-alpha.e"]
    sided <- eval(object$call$sided)
    if (is.null(sided))
        sided <- 1
    dW <- matrix(rnorm(nsim * n), n, nsim)
    W0 <- NULL
    W.i <- 0 * dW[1, ]
    f.old <- 0
    for (i in 1:n) {
        f.new <- .fr.[i]
        df <- f.new - f.old
        W.i <- W.i + df^0.5 * dW[i, ]
        W0 <- rbind(W0, W.i)
        f.old <- f.new
    }
    if (do.fu)
        W1 <- W0 + .dr.
    acc.H0.yet <- rej.H0.yet <- rep(FALSE, nsim)    
    if (sided == 1) {
        cont0.old <- cont1.old <- rep(TRUE, nsim)
        rej.H0 <- acc.H0 <- NULL
        for (i in 1:n) {
            cont0.new <- cont0.old & (W0[i, ] > a[i] & W0[i,
                ] < b[i])
            rej.H0.i <- cont0.old & (!acc.H0.yet) & (W0[i, ] >= b[i])
            rej.H0 <- rbind(rej.H0, rej.H0.i)
            cont0.old <- cont0.new
            if (do.fu) {
                cont1.new <- cont1.old & (W1[i, ] > a[i] & W1[i, ] < b[i])
                acc.H0.i <- cont1.old & (!rej.H0.yet) & (W1[i, ] <= a[i])
                acc.H0 <- rbind(acc.H0, acc.H0.i)
                acc.H0.yet <- acc.H0.yet | acc.H0.i
                cont1.old <- cont1.new
            }
            rej.H0.yet <- rej.H0.yet | rej.H0.i
        }
    }
    if (sided == 2) {
        cont0.old <- cont1.old <- rep(TRUE, nsim)
        rej.H0 <- acc.H0 <- NULL
        for (i in 1:n) {
            cont0.new <- cont0.old & (((W0[i, ] > a[i]) & (W0[i, ] < b[i])) | ((W0[i, ] > -b[i]) & (W0[i, ] < -a[i])))
            rej.H0.i <- cont0.old & (!acc.H0.yet) & ((W0[i, ] >= b[i]) | (W0[i, ] <= -b[i]))
            rej.H0 <- rbind(rej.H0, rej.H0.i)
            cont0.old <- cont0.new
            if (do.fu && a[i] > 0) {
                cont1.new <- cont1.old & (((W1[i, ] > a[i]) & (W1[i, ] < b[i])) | ((W1[i, ] > -b[i]) & (W1[i, ] < -a[i])))
                acc.H0.i <- cont1.old & (!rej.H0.yet) & (W1[i, ] >= -a[i]) & (W1[i, ] <= a[i])
                acc.H0 <- rbind(acc.H0, acc.H0.i)
                acc.H0.yet <- acc.H0.yet | acc.H0.i
                cont1.old <- cont1.new
            }
            rej.H0.yet <- rej.H0.yet | rej.H0.i
        }
    }
    ErrI <- cumsum(c(rej.H0 %*% rep(1/nsim, nsim)))
    if (do.fu)
        ErrII <- cumsum(c(acc.H0 %*% rep(1/nsim, nsim)))
    ans <- cbind(ErrI, al.e)
    nms <- c("eI.est", "eI.act")
    if (do.fu) {
        n.diff <- n - length(ErrII)
        zeros <- rep(0, n.diff)
        ans <- cbind(ans, c(zeros, ErrII), al.f)
        nms <- c(nms, "eII.est", "eII.act")
    }
    dimnames(ans) <- list(signif(.fr.ii, 4), nms)
    ans
}

"SimGSB.PwrGSD" <-
function(object, nsim=1e+05, ...)
{
  object <- as.boundaries(object)
  SimGSB(object)
}

"as.boundaries" <-
function(object, ...){
  UseMethod("as.boundaries")
}

"as.boundaries.boundaries" <-
function(object, ...){
  object
}

"as.boundaries.PwrGSD" <-
function(object, ...)
{
  dots <- c(...)
  nms.dots <- names(dots)
  stat <- 1
  if("stat" %in% nms.dots)
    stat <- dots["stat"]
  nlook <- object$detail$ints["nlook"]
  dofu <- object$detail$ints["do.futility"]
  nstat <- object$detail$ints["nstat"]
  if(dofu==0){
    cl <- as.call(expression(GrpSeqBnds)) 
    cl$frac <- frac <- c(object$detail$pinffrac[,stat])
    if(sum(abs(c(object$detail$pinffrac[,stat])-c(object$detail$pinffrac.ii[,stat]))) > 1e-8)
      cl$frac.ii <- frac.ii <- c(object$detail$pinffrac.ii[,stat])
    else frac.ii <- frac
    cl$EfficacyBoundary <- object$call$EfficacyBoundary
    out <- cbind(object$detail$pinffrac.ii[,stat], object$detail$pbounds[,stat], object$detail$palpha0vec[,stat],
                  cumsum(object$detail$palpha0vec[,stat]))
    dimnames(out) <- list(1:nlook, c("frac","b.e", "alpha.e", "cum-alpha.e"))
    ans <- list(table=out, frac=frac, frac.ii=frac.ii, call=cl)
  }
  if(dofu==1){
    cl <- as.call(expression(GrpSeqBnds))
    cl$frac <- frac <- c(object$detail$pinffrac[,stat])
    if(sum(abs(c(object$detail$pinffrac[,stat])-c(object$detail$pinffrac.ii[,stat]))) > 0)
      cl$frac.ii <- frac.ii <- c(object$detail$pinffrac.ii[,stat])
    else frac.ii <- frac
    cl$EfficacyBoundary <- object$call$EfficacyBoundary
    cl$FutilityBoundary <- object$call$FutilityBoundary
    cl$drift <- drift <- c(object$detail$mufu[,stat])
    if(as.character(cl$FutilityBoundary[[1]])=="SC")
      cl$FutilityBoundary$drift.end <- object$detail$mufu[nlook, stat]
    out <- cbind(object$detail$pinffrac.ii[,stat], object$detail$pbounds[,nstat + stat], object$detail$palpha0vec[,nstat + stat],
                  cumsum(object$detail$palpha0vec[,nstat + stat]), object$detail$pbounds[,stat],
                  object$detail$palpha0vec[,stat], cumsum(object$detail$palpha0vec[,stat]))
    dimnames(out) <- list(1:nlook, c("frac","b.f","alpha.f","cum-alpha.f","b.e", "alpha.e", "cum-alpha.e"))
    ans <- list(table=out, frac=frac, frac.ii=frac.ii, drift=drift, call=cl)
  }
  ans <- update(ans)
  ans
}

"SimPwrGSD" <-
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         NonBindingFutility=TRUE,
         sided =c("2>","2<","1>","1<"),accru,accrat,tlook,
         tcut0 = NULL,h0 = NULL,s0 = NULL,tcut1 = NULL,rhaz = NULL,
         h1 = NULL,s1 = NULL,tcutc0 = NULL,hc0 = NULL,sc0 = NULL,tcutc1 = NULL,hc1 = NULL,
         sc1 = NULL,tcutd0A = NULL,hd0A = NULL,sd0A = NULL,tcutd0B = NULL,hd0B = NULL,sd0B = NULL,
         tcutd1A = NULL,hd1A = NULL,sd1A = NULL,tcutd1B = NULL,hd1B = NULL,sd1B = NULL,
         tcutx0A = NULL,hx0A = NULL,sx0A = NULL,tcutx0B = NULL,hx0B = NULL,sx0B = NULL,
         tcutx1A = NULL,hx1A = NULL,sx1A = NULL,tcutx1B = NULL,hx1B = NULL,sx1B = NULL,
         noncompliance = c("none","crossover","mixed","user"),gradual = FALSE,
         WtFun = c("FH","SFH","Ramp"),ppar = cbind(c(0,0)),
         Spend.Info=c("Variance","Events","Hybrid(k)","Calendar"),
         RR.Futility = NULL,qProp.one.or.Q = c("one","Q"),Nsim=NULL,
         detail = FALSE,StatType = c("WLR","ISD"),doProj=FALSE)
{

    if(!missing(StatType) && any(!(StatType %in% c("WLR","ISD"))))
      stop("Elements of the vector argument 'StatType' be either \"WLR\" or \"ISD\"")
    if(!missing(StatType) && !missing(WtFun) && length(StatType)!=length(WtFun))
      stop("Vector arguments 'StatType' and 'WtFun' must be of the same length")

    if(missing(WtFun)) {
      WtFun <- "FH"
      wttyp <- 0
      ppar <- c(0,0)
      nstat <- 1
    }
    else{
      if(missing(ppar) && any(WtFun!="FH"))
        stop("You must specify parameters for the chosen weight function(s) in the 'ppar' argument.")
      nstat <- length(WtFun)
      wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1
    }
    if(missing(StatType)) StatType <- rep("WLR", nstat)

    stattype <- 1*(StatType == "ISD")
    is.detail <- !missing(detail)
    
    .call. <- match.call()
    ppar <- t(ppar)
    ndbl <- floor(accru * accrat)
    ndbl <- ndbl + (ndbl%%2)
    n <- ndbl/2
    nlook <- length(tlook)
    if (missing(sided))
        sided <- "2>"
    if (!all(sided %in% c("2>","2<","1>", "1<")))
        stop("Elements of argument \"sided\" must be " %,%
             "equal to \"2>\", \"2<\", \"1>\" or \"1<\"")
    sided <- c(2, -2, 1, -1)[grep(sided, c("2>", "2<", "1>", "1<"))]
    
    call.PBS <- .call.
    call.PBS[[1]] <- as.name("ParseBoundarySelection")

    call.PBS$sided <- call.PBS$accru <- call.PBS$accrat <- call.PBS$tlook <- 
    call.PBS$tcut0 <- call.PBS$h0 <- call.PBS$s0 <- call.PBS$tcut1 <- 
    call.PBS$rhaz <- call.PBS$h1 <- call.PBS$s1 <- call.PBS$tcutc0 <-
    call.PBS$hc0  <- call.PBS$sc0 <- call.PBS$tcutc1 <- call.PBS$hc1 <-
    call.PBS$sc1 <- call.PBS$tcutd0A <- call.PBS$hd0A <- call.PBS$sd0A <-
    call.PBS$tcutd0B <- call.PBS$hd0B <- call.PBS$sd0B <- call.PBS$tcutd1A <-
    call.PBS$hd1A <- call.PBS$sd1A <- call.PBS$tcutd1B <- call.PBS$hd1B <- 
    call.PBS$sd1B<- call.PBS$tcutx0A<- call.PBS$hx0A<- call.PBS$sx0A <-
    call.PBS$tcutx0B <- call.PBS$hx0B<- call.PBS$sx0B<- call.PBS$tcutx1A <-
    call.PBS$hx1A <- call.PBS$sx1A <- call.PBS$tcutx1B <- call.PBS$hx1B <-
    call.PBS$sx1B <- call.PBS$noncompliance <- call.PBS$gradual <-
    call.PBS$WtFun <- call.PBS$ppar <- call.PBS$Spend.Info <-
    call.PBS$RR.Futility <- call.PBS$V.end <- call.PBS$qProp.one.or.Q <-
    call.PBS$detail <- call.PBS$Nsim <- call.PBS$doProj <- NULL
    call.PBS$n.looks <- nlook
    call.PBS$check.drift <- FALSE
        
    eval.PBS <- eval(call.PBS)

    do.efficacy <- TRUE
    do.futility <- !is.null(eval.PBS$frontend$FutilityBoundary)
       
    Alpha.Efficacy <- eval.PBS$backend$Alpha.Efficacy
    Alpha.Futility <- eval.PBS$backend$Alpha.Futility
    if(is.null(Alpha.Futility)) Alpha.Futility <- 0
    nbnd.e <- eval.PBS$backend$nbnd.e
    nbnd.f <- eval.PBS$backend$nbnd.f
    nsf.e <- eval.PBS$backend$nsf.e
    nsf.f <- eval.PBS$backend$nsf.f
    rho.Efficacy <- eval.PBS$backend$rho.Efficacy
    rho.Futility <- eval.PBS$backend$rho.Futility
    b.Haybittle.e <- eval.PBS$backend$b.Haybittle.e
    b.Haybittle.f <- eval.PBS$backend$b.Haybittle.f
    drift.end <- eval.PBS$backend$drift.end
    prob.e <- eval.PBS$backend$prob.e
    prob.f <- eval.PBS$backend$prob.f
    my.Efficacy <- eval.PBS$backend$my.Efficacy
    my.Futility <- eval.PBS$backend$my.Futility
    is.myE <- !all(my.Efficacy==0)
    is.myF <- !all(my.Futility==0)
   
    b.e <- rep(0, nlook)
    if (is.myE) b.e <- my.Efficacy

    b.f <- rep(0, nlook)
    if (is.myF) b.f <- my.Futility

    Alpha.Efficacy <- Alpha.Efficacy/2^(abs(sided) == 2)

    no.SpndInfo <- missing(Spend.Info)
    if(no.SpndInfo) Spend.Info <- "Variance"
    if(!no.SpndInfo){
      Spend.Info <- as.character(.call.$Spend.Info)
      if(!(Spend.Info[1] %in% c("Variance", "Events", "Hybrid", "Calendar")))
        stop("Argument 'Spend.Info' is an expression of the form 'Variance', 'Events', 'Hybrid(k)', or 'Calendar'")
    }
    spend.info <- grep(Spend.Info[1], c("Variance", "Events", "Hybrid", "Calendar")) - 1
    spend.info.k <- 0
    if(spend.info==2 && length(Spend.Info)>1) spend.info.k <- as.integer(as.numeric(Spend.Info[2])) - 1

    if(missing(qProp.one.or.Q))
      qProp.one.or.Q <- 0
    else{
      if(!(qProp.one.or.Q %in% c("one","Q")))
        stop("Argument 'qProp.one.or.Q' must be either \"one\" or \"Q\"")
      qProp.one.or.Q <- grep(qProp.one.or.Q, c("one", "Q")) - 1
    }

    stoh <- function(tcut, s) {
        ncut <- length(tcut)
        Sold <- c(1, s[-(ncut - 1)])
        dt <- diff(tcut)
        log.0 <-
        function(x)
        {
          y <- 0*x
          y[x>0] <- log(x[x>0])
          y
        }
        h <- log.0(s/Sold)/dt
        h
    }
    no.t0 <- missing(tcut0)
    no.s0 <- missing(s0)
    no.h0 <- missing(h0)
    if ((no.s0 && no.h0) || no.t0)
        stop("Must specify 'tcut0' and ('h0' or 's0').")
    if (no.h0)
        h0 <- stoh(tcut0, s0)
    ncut0 <- length(tcut0)
    
    no.t1 <- missing(tcut1)
    no.s1 <- missing(s1)
    no.h1 <- missing(h1)
    no.rhaz <- missing(rhaz)
    if ((no.s1 && no.rhaz && no.h1) || no.t1)
        stop("Must specify 'tcut1' and ('rhaz' or 'h1' or 's1').")
    if (!no.s1)
        h1 <- stoh(tcut1, s1)
    if (!no.rhaz) {
        tcut.01 <- support(c(tcut0, tcut1))
        h0.01 <- lookup(tcut0, h0, tcut.01)$y
        rhaz.01 <- lookup(tcut1, rhaz, tcut.01)$y
        tcut1 <- tcut.01
        h1 <- rhaz.01 * h0.01
    }
    ncut1 <- length(tcut1)

    use.rhaz.fu <- 0
    if(missing(RR.Futility)) {
      RR.Futility <- rep(1,ncut0)
      if(do.futility==1) use.rhaz.fu <- 1
    }
    if(length(RR.Futility)<ncut0) RR.Futility <- rep(RR.Futility[1], ncut0)
    
    no.tc0 <- missing(tcutc0)
    no.sc0 <- missing(sc0)
    no.hc0 <- missing(hc0)
    if ((no.sc0 && no.hc0) || no.tc0)
        stop("Must specify 'tcutc0' and ('hc0' or 'sc0').")
    if (no.hc0)
        hc0 <- stoh(tcutc0, sc0)
    ncutc0 <- length(tcutc0)
    no.tc1 <- missing(tcutc1)
    no.sc1 <- missing(sc1)
    no.hc1 <- missing(hc1)
    if ((no.sc1 && no.hc1) || no.tc1)
        stop("Must specify 'tcutc1' and ('hc1' or 'sc1').")
    if (no.hc1)
        hc1 <- stoh(tcutc1, sc1)
    ncutc1 <- length(tcutc1)
    noncompliance <- as.character(.call.$noncompliance)
    no.noncomp <- (length(noncompliance) == 0)
    if (no.noncomp)
        noncompliance <- "none"
    switch(noncompliance, none = {
        tcutd0A <- 0
        hd0A <- 0
        ncutd0A <- 1
        tcutd0B <- 0
        hd0B <- 0
        ncutd0B <- 1
        tcutd1A <- 0
        hd1A <- 0
        ncutd1A <- 1
        tcutd1B <- 0
        hd1B <- 0
        ncutd1B <- 1
        tcutx0A <- tcut0
        hx0A <- h0
        ncutx0A <- ncut0
        tcutx0B <- tcut0
        hx0B <- h0
        ncutx0B <- ncut0
        tcutx1A <- tcut1
        hx1A <- h1
        ncutx1A <- ncut1
        tcutx1B <- tcut1
        hx1B <- h1
        ncutx1B <- ncut1
    }, crossover = {
        no.td0B <- missing(tcutd0B)
        no.sd0B <- missing(sd0B)
        no.hd0B <- missing(hd0B)
        no.td1B <- missing(tcutd1B)
        no.sd1B <- missing(sd1B)
        no.hd1B <- missing(hd1B)
        no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) ||
            no.td0B || no.td1B
        if (no.dB)
            stop("crossover option requires specification of\n 'tcutd0B', 'tcutd1B', ('hd0B' or " %,%
                "'sd0B') and ('hd1B' or 'sd1B').\n")
        if (no.hd0B)
            hd0B <- stoh(tcutd0B, sd0B)
        ncutd0B <- length(tcutd0B)
        if (no.hd1B)
            hd1B <- stoh(tcutd1B, sd1B)
        ncutd1B <- length(tcutd1B)
        tcutd0A <- 0
        hd0A <- 0
        ncutd0A <- 1
        tcutd1A <- 0
        hd1A <- 0
        ncutd1A <- 1
        tcutx0A <- tcut0
        hx0A <- h0
        ncutx0A <- ncut0
        tcutx1A <- tcut1
        hx1A <- h1
        ncutx1A <- ncut1
        tcutx0B <- tcut1
        hx0B <- h1
        ncutx0B <- ncut1
        tcutx1B <- tcut0
        hx1B <- h0
        ncutx1B <- ncut0
    }, mixed = {
        no.td0A <- missing(tcutd0A)
        no.sd0A <- missing(sd0A)
        no.hd0A <- missing(hd0A)
        no.td0B <- missing(tcutd0B)
        no.sd0B <- missing(sd0B)
        no.hd0B <- missing(hd0B)
        no.td1A <- missing(tcutd1A)
        no.sd1A <- missing(sd1A)
        no.hd1A <- missing(hd1A)
        no.td1B <- missing(tcutd1B)
        no.sd1B <- missing(sd1B)
        no.hd1B <- missing(hd1B)
        no.tx0A <- missing(tcutx0A)
        no.sx0A <- missing(sx0A)
        no.hx0A <- missing(hx0A)
        no.tx1A <- missing(tcutx1A)
        no.sx1A <- missing(sx1A)
        no.hx1A <- missing(hx1A)
        no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) ||
            no.td0A || no.td1A
        no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) ||
            no.td0B || no.td1B
        no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) ||
            no.tx0A || no.tx1A
        if (no.dA || no.dB || no.xA)
            stop("mixed option requires specification of \n'tcutd0A', 'tcutd1A', ('hd0A' or " %,%
                 "'sd0A'), ('hd1A' or 'sd1A'),\n'tcutd0B', 'tcutd1B', ('hd0B' or " %,%
                 "'sd0B'), ('hd1B' or 'sd1B'),\n'tcutx0A', 'tcutx1A', ('hx0A' or " %,%
                 "'sx0A'), and ('hx1A' or 'sx1A').\n")
        if (no.hd0A)
            hd0A <- stoh(tcutd0A, sd0A)
        ncutd0A <- length(tcutd0A)
        if (no.hd1A)
            hd1A <- stoh(tcutd1A, sd1A)
        ncutd1A <- length(tcutd1A)
        if (no.hd0B)
            hd0B <- stoh(tcutd0B, sd0B)
        ncutd0B <- length(tcutd0B)
        if (no.hd1B)
            hd1B <- stoh(tcutd1B, sd1B)
        ncutd1B <- length(tcutd1B)
        if (no.hx0A)
            hx0A <- stoh(tcutx0A, sx0A)
        ncutx0A <- length(tcutx0A)
        if (no.hx1A)
            hx1A <- stoh(tcutx1A, sx1A)
        ncutx1A <- length(tcutx1A)
        tcutx0B <- tcut1
        hx0B <- h1
        ncutx0B <- ncut1
        tcutx1B <- tcut0
        hx1B <- h0
        ncutx1B <- ncut0
    }, user = {
        no.td0A <- missing(tcutd0A)
        no.sd0A <- missing(sd0A)
        no.hd0A <- missing(hd0A)
        no.td0B <- missing(tcutd0B)
        no.sd0B <- missing(sd0B)
        no.hd0B <- missing(hd0B)
        no.td1A <- missing(tcutd1A)
        no.sd1A <- missing(sd1A)
        no.hd1A <- missing(hd1A)
        no.td1B <- missing(tcutd1B)
        no.sd1B <- missing(sd1B)
        no.hd1B <- missing(hd1B)
        no.tx0A <- missing(tcutx0A)
        no.sx0A <- missing(sx0A)
        no.hx0A <- missing(hx0A)
        no.tx1A <- missing(tcutx1A)
        no.sx1A <- missing(sx1A)
        no.hx1A <- missing(hx1A)
        no.tx0B <- missing(tcutx0B)
        no.sx0B <- missing(sx0B)
        no.hx0B <- missing(hx0B)
        no.tx1B <- missing(tcutx1B)
        no.sx1B <- missing(sx1B)
        no.hx1B <- missing(hx1B)
        no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) ||
            no.td0A || no.td1A
        no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) ||
            no.td0B || no.td1B
        no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) ||
            no.tx0A || no.tx1A
        no.xB <- (no.sx0B && no.hx0B) || (no.sx1B && no.hx1B) ||
            no.tx0B || no.tx1B
        if ((no.dA || no.xA) && (no.dB || no.xB))
            stop("user option requires specification of \n('tcutd0A', 'tcutd1A', ('hd0A' or " %,%
                 "'sd0A'), ('hd1A' or 'sd1A') and\n'tcutx0A', 'tcutx1A', ('hx0A' or " %,%
                 "'sx0A'), ('hx1A' or 'sx1A')) or \n('tcutd0B', 'tcutd1B', ('hd0B' or " %,%
                 "'sd0B'), and ('hd1B' or 'sd1B') and\n 'tcutx0B', 'tcutx1B', ('hx0B' or " %,%
                 "'sx0B'), and ('hx1B' or 'sx1B')).\n")
        if (!no.dA) {
            if (no.hd0A)
                hd0A <- stoh(tcutd0A, sd0A)
            ncutd0A <- length(tcutd0A)
            if (no.hd1A)
                hd1A <- stoh(tcutd1A, sd1A)
            ncutd1A <- length(tcutd1A)
            if (no.hx0A)
                hx0A <- stoh(tcutx0A, sx0A)
            ncutx0A <- length(tcutx0A)
            if (no.hx1A)
                hx1A <- stoh(tcutx1A, sx1A)
            ncutx1A <- length(tcutx1A)
            if (no.dB) {
                tcutd0B <- 0
                hd0B <- 0
                ncutd0B <- 1
                tcutd1B <- 0
                hd1B <- 0
                ncutd1B <- 1
                tcutx0B <- tcut0
                hx0B <- h0
                ncutx0B <- ncut0
                tcutx1B <- tcut1
                hx1B <- h1
                ncutx1B <- ncut1
            }
        }
        if (!no.dB) {
            if (no.hd0B)
                hd0B <- stoh(tcutd0B, sd0B)
            ncutd0B <- length(tcutd0B)
            if (no.hd1B)
                hd1B <- stoh(tcutd1B, sd1B)
            ncutd1B <- length(tcutd1B)
            if (no.hx0B)
                hx0B <- stoh(tcutx0B, sx0B)
            ncutx0B <- length(tcutx0B)
            if (no.hx1B)
                hx1B <- stoh(tcutx1B, sx1B)
            ncutx1B <- length(tcutx1B)
            if (no.dA) {
                tcutd0A <- 0
                hd0A <- 0
                ncutd0A <- 1
                tcutd1A <- 0
                hd1A <- 0
                ncutd1A <- 1
                tcutx0A <- tcut0
                hx0A <- h0
                ncutx0A <- ncut0
                tcutx1A <- tcut1
                hx1A <- h1
                ncutx1A <- ncut1
            }
        }
    })
    nunq <- length(unique(sort(c(tcut0, tcut1, tcutc0, tcutc1, tcutd0A, tcutd0B, tcutd1A, tcutd1B, tcutx0A,
                                 tcutx0B, tcutx1A, tcutx1B))))
    glegx <- glegx24
    glegw <- glegw24
    
    NGaussQ <- length(glegx)
    
    cumsum.nppar <- 0
    stat.nms <- NULL
    for(j in 1:nstat){
      stat.<- c("WLR[","ISD[")[1 + stattype[j]]
      wt. <- c("FH(", "SFH(", "R(")[1+wttyp[j]]
      nppar <- c(2,3,1)[1+wttyp[j]]
      ppar.j <- ppar[(cumsum.nppar+1):(cumsum.nppar+nppar)]
      par.string <- glue.2.string(ppar.j, sep=",")
      stat.nms <- c(stat.nms, stat. %,% wt. %,% par.string %,% ")]")
      cumsum.nppar <- cumsum.nppar + nppar
    }
    
    nbetyp <- length(nbnd.e)
    nbftyp <- length(nbnd.f)

    stattype.nms <- c("WLR","ISD")[1+stattype]
    wttyp.nms <- c("FH","SFH","Ramp")[1+stattype]
    ints <- c(nlook,nstat,NGaussQ,ncut0,ncut1,ncutc0,ncutc1,ncutd0A,ncutd0B,ncutd1A,
              ncutd1B,ncutx0A,ncutx0B,ncutx1A,ncutx1B,gradual,nbnd.e,nbnd.f,
              nsf.e,nsf.f,do.futility,use.rhaz.fu,spend.info,Nsim,is.myE,is.myF,spend.info.k,
              qProp.one.or.Q,sided,NonBindingFutility,stattype,wttyp,doProj)

    ints.nms <- c("nlook","nstat","NGaussQ","ncut0","ncut1","ncutc0","ncutc1",
                  "ncutd0A","ncutd0B","ncutd1A","ncutd1B","ncutx0A","ncutx0B","ncutx1A",
                  "ncutx1B","gradual","nbnd.e."%,%(1:nlook),"nbnd.f."%,%(1:nlook),
                  "nsf.e."%,%(1:nlook),"nsf.f."%,%(1:nlook),"do.futility",
                  "use.rhaz.fu","spend.info","Nsim","is.myE","is.myF","spend.info.k",
                  "qis1orQ", "sided","NonBindingFutility",stattype.nms,wttyp.nms,"doProj")
    
    dbls <- c(accru,accrat, rho.Efficacy, rho.Futility, prob.e, prob.f)

    dbls.nms <- c("accru","accrat", "rho.Efficacy."%,%(1:nlook), "rho.Futility."%,%(1:nlook),
                  "prob.e."%,%(1:nlook), "prob.f."%,%(1:nlook))
    
    logRR.F <- log(RR.Futility)
    
    details <- .C("SimPwrGSD", 
	ints = as.integer(ints),
        dbls = as.double(dbls),
        pttlook = as.double(tlook), 
	palphatot = as.double(c(Alpha.Efficacy,Alpha.Futility)), 
	lrrf = as.double(logRR.F),
	bHay = as.double(c(b.Haybittle.e, b.Haybittle.f)),
	ppar = as.double(ppar), 
	pgqxw = as.double(c(glegx,glegw)),
	tcut0 = as.double(tcut0), 
	h0 = as.double(h0),
        tcut1 = as.double(tcut1), 
	h1 = as.double(h1), 
	tcutc0 = as.double(tcutc0),
        hc0 = as.double(hc0), 
	tcutc1 = as.double(tcutc1), 
	hc1 = as.double(hc1),
        tcutd0A = as.double(tcutd0A), 
	hd0A = as.double(hd0A),
        tcutd0B = as.double(tcutd0B), 
	hd0B = as.double(hd0B),
        tcutd1A = as.double(tcutd1A), 
	hd1A = as.double(hd1A),
        tcutd1B = as.double(tcutd1B), 
	hd1B = as.double(hd1B),
        tcutx0A = as.double(tcutx0A), 
	hx0A = as.double(hx0A),
        tcutx0B = as.double(tcutx0B), 
	hx0B = as.double(hx0B),
        tcutx1A = as.double(tcutx1A), 
	hx1A = as.double(hx1A),
        tcutx1B = as.double(tcutx1B), 
	hx1B = as.double(hx1B),
        t0 = double(n), 
	t1 = double(n),
        tc0 = double(n), 
	tc1 = double(n),
        td0A = double(n), 
	td0B = double(n), 
	td1A = double(n), 
	td1B = double(n), 
	code = integer(2 + ndbl), 
	u = double(2 + ndbl), 
	TT = double(2 + ndbl), 
	delta = integer(2 + ndbl), 
	z = integer(2 + ndbl),
	time = double(2 + ndbl), 
	nrisk = integer(2 + 2*ndbl), 
	nevent = integer(2 + 2*ndbl), 
	pndths = integer(1), 
	avg.inffrac = double(nstat * nlook),
        avg.inffrac.ii = double(nstat * nlook),
	avg.bounds = as.double(c(rep(b.e, nstat), rep(b.f, nstat))),
        mufu = double(nstat*nlook),
	palphavec = double(2 * nstat * nlook), 
	pRejAcc = integer(2 * nstat * Nsim),
        kstop = integer(nstat * Nsim),
        duration = double(nstat * Nsim),
	pStatTend = double(nstat * Nsim), 
        pVarTend = double(nstat * Nsim),
        pmTend1 = double(nstat * Nsim),
	pstatK = double(nstat * Nsim),
        pvarK = double(nstat * Nsim),
        v.Tend.proj = double(nstat * Nsim),
        m.Tend.proj = double(nstat * Nsim),
	PACKAGE = "PwrGSD")

    ndths <- details$pndths
    names(details$ints) <- ints.nms
    names(details$dbls) <- dbls.nms
    details$palpha0vec <- matrix(details$palphavec, nlook, 2*nstat)
    details$palphavec <- NULL
    details$mufu <- matrix(details$mufu, nlook, nstat)
    details$RejAcc <- matrix(details$pRejAcc, Nsim, 2*nstat)
    details$time <- details$time[1:ndths]
    details$nrisk <- details$nrisk[1:ndths]
    details$nevent <- details$nevent[1:ndths]
    details$time1 <- details$time1[1:ndths]
    details$nrisk1 <- details$nrisk1[1:ndths]
    details$nevent1 <- details$nevent1[1:ndths]
    details$StatK<- matrix(details$pstatK, Nsim, nstat)
    details$pstatK <- NULL
    details$VarK <- matrix(details$pvarK, Nsim, nstat)
    details$pvarK <- NULL
    details$pinffrac <- matrix(details$avg.inffrac, nlook, nstat)
    details$avg.inffrac <- NULL
    details$pinffrac.ii <- matrix(details$avg.inffrac.ii, nlook, nstat)
    details$avg.inffrac.ii <- NULL
    details$pbounds <- matrix(details$avg.bounds, nlook, (1 + do.futility)*nstat)
    details$avg.bounds <- NULL   

    RejNull <- details$RejAcc[,1:nstat,drop=FALSE]
    AccNull <- details$RejAcc[,nstat + (1:nstat)]
    duration <- matrix(details$duration, Nsim, nstat)
    StatTend <- matrix(details$pStatTend, Nsim, nstat)
    VarTend <- matrix(details$pVarTend, Nsim, nstat)
    mTend1 <- matrix(details$pmTend1, Nsim, nstat)
    details$RejAcc <- details$duration <- details$StatTend <- details$VarTend <- details$mTend1 <- NULL
      
    kstop <- matrix(details$kstop, Nsim, nstat)
    sum.kstop <- matrix(0, nlook, nstat)
    for(k in 1:nlook) sum.kstop[k,] <- c(rep(1,Nsim) %*% (kstop==k))
    Nstop <- rbind(apply(sum.kstop[nlook:1,,drop=FALSE], 2, FUN=cumsum))[nlook:1,]

    details$pinffrac <- details$pinffrac/Nstop
    details$pinffrac.ii <- details$pinffrac.ii/Nstop
    for(ef in 1:2) details$pbounds[, nstat*(ef-1) + (1:nstat)] <- details$pbounds[, nstat*(ef-1) + (1:nstat)]/Nstop

    dimnames(RejNull) <- list(1:Nsim, stat.nms)
    dimnames(duration) <- list(1:Nsim, stat.nms)
    dimnames(details$StatK) <- list(1:Nsim, stat.nms)
    dimnames(details$VarK) <- list(1:Nsim, stat.nms)
    dimnames(details$pinffrac) <- list(1:nlook, stat.nms)
    dimnames(details$pinffrac.ii) <- list(1:nlook, stat.nms)
    
    bstat.nms <- "Eff." %,% stat.nms
    if(do.futility==1) bstat.nms <- c(bstat.nms, "Fut." %,% stat.nms)
    dimnames(details$pbounds) <- list(1:nlook, bstat.nms)
    dimnames(sum.kstop) <- list(1:nlook, stat.nms)

    two.more <- 0
    v.Tend.proj <- m.Tend.proj <- vTp.nm <- mTp.nm <- NULL
    if(doProj)
    {
      v.Tend.proj <- matrix(details$v.Tend.proj, Nsim, nstat)
      m.Tend.proj <- matrix(details$m.Tend.proj, Nsim, nstat)
      vTp.nm <- "VarTendProj"
      mTp.nm <- "mTend1Proj"
      two.more <- 2
    }

    rslts <- cbind(details$StatK, details$VarK^0.5, StatTend, duration, v.Tend.proj, VarTend,
                   m.Tend.proj, mTend1, RejNull)[,outer(nstat*(0:(6+two.more)), 1:nstat, FUN="+")]
    dimnames(rslts) <- list(1:Nsim, c(outer(c("Stat", "SE", "StatTend","Duration", vTp.nm, "VarTend", mTp.nm, "mTend1",
                                              "RejNull"), stat.nms, FUN=paste, sep="-")))
    details$StatK <- details$VarK <- details$pStatTend <- details$pVarTend <- details$pmTend1 <- NULL
    dPower <- sapply(tlook, FUN=function(x,tt,rr)apply(rr*(tt==x),2,FUN=mean),tt=duration,rr=RejNull)
    dErrorII <- sapply(tlook, FUN=function(x,tt,rr)apply(rr*(tt==x),2,FUN=mean),tt=duration,rr=AccNull)
    se.dPower <- ((dPower * (1 - dPower))/Nsim)^0.5
    mu.d <- apply(duration, 2, FUN = mean)
    se.d <- (apply(duration, 2, FUN = var)/Nsim)^0.5
    
    out <- list(dPower=dPower, se.dPower=se.dPower, dErrorII=dErrorII, Exp.Dur=mu.d, se.Exp.Dur=se.d,
                Nsim = Nsim, StatLast = details$StatLast, ndeaths = details$pndeaths, Results = rslts,
                fail = details$fail,stat.nms=stat.nms, RejNull=RejNull, duration=duration,
                InfFrac=details$pinffrac, InfFrac.ii=details$pinffrac.ii, Bounds=details$pbounds,
                Nstop=sum.kstop, detail=list(ints=details$ints, dbls=details$dbls), call = .call.)
    if (is.detail)
        out$detail <- details
    class(out) <- c("SimPwrGSD", "PwrGSD")
    out
}

"print.SimPwrGSD" <-
function (x, ...)
{
    print(x$call)
    stat.nms <- x$stat.nms
    nstat <- x$detail$ints["nstat"]
    nlook <- x$detail$ints["nlook"]
    ans <- cbind(x$dPower %*% rep(1,nlook), (x$se.dPower^2) %*% rep(1, nlook), x$Exp.Dur, x$se.Exp.Dur)
    dimnames(ans) <- list(stat.nms, c("Power","se(Power)", "Exp.Dur", "se(Exp.Dur)"))
    print(ans)
    invisible(x)
}

"summary.SimPwrGSD" <- 
function(object, ...)
{
	out <- object
	stat.nms <- object$stat.nms
	nstat <- object$detail$ints["nstat"]
        nlook <- object$detail$ints["nlook"]
        ans <- cbind(object$dPower %*% rep(1,nlook), (object$se.dPower^2) %*% rep(1, nlook), object$Exp.Dur, object$se.Exp.Dur)
        dimnames(ans) <- list(stat.nms, c("Power","se(Power)", "Exp.Dur", "se(Exp.Dur)"))
	out$Tbl <- ans
        out
}

"SimPwrGSDcall" <-
function (Nsim, accru, accrat, tlook,
    Alpha.Efficacy=0, Alpha.Futility=0, sided = c("2>", "2<",
    "1>","1<"), tcut0 = NULL, h0 = NULL, s0 = NULL, tcut1 = NULL, 
    rhaz = NULL, h1 = NULL, s1 = NULL, tcutc0 = NULL, hc0 = NULL, 
    sc0 = NULL, tcutc1 = NULL, hc1 = NULL, sc1 = NULL, tcutd0A = NULL, 
    hd0A = NULL, sd0A = NULL, tcutd0B = NULL, hd0B = NULL, 
    sd0B = NULL, tcutd1A = NULL, hd1A = NULL, sd1A = NULL, 
    tcutd1B = NULL, hd1B = NULL, sd1B = NULL, tcutx0A = NULL, 
    hx0A = NULL, sx0A = NULL, tcutx0B = NULL, hx0B = NULL, 
    sx0B = NULL, tcutx1A = NULL, hx1A = NULL, sx1A = NULL,
    tcutx1B = NULL, hx1B = NULL, sx1B = NULL, 
    noncompliance = c("none", "crossover", "mixed", "user"), 
    gradual = FALSE, detail = FALSE, WtFun =c("FH","SFH","Ramp"), ppar = cbind(c(0, 0)), 
    Boundary.Efficacy = c("Lan-Demets", "Haybittle"), 
    Boundary.Futility = c("Lan-Demets", "Haybittle"),
    RR.Futility = NULL,                          
    Spending.Efficacy = c("Obrien-Fleming", "Pocock", "Power"),
    Spending.Futility = c("Obrien-Fleming", "Pocock", "Power"), 
    rho.Efficacy=0, rho.Futility=0, b.Haybittle.e = 3, b.Haybittle.f)  
{
	match.call()
}

"gsd.dens" <-
function(x, frac=NULL, scale="Standard")
{
  if(class(x)=="boundaries")
  {
    frac <- x$frac
    x.eff <- x$table[,"b.e"]*frac^0.5
  }
  if(class(x)=="numeric")
  {
    if(missing(frac))
      stop("You must either specify an object of class 'boundaries' or specify both the " %,%
           "efficacy boundary, 'b.e' and information fraction, 'frac'")
    x.eff <- x
    if(scale=="Standard") x.eff <- x*frac^0.5
    if(length(x.eff)!=length(frac))
      stop("Informtion fraction and efficacy boundary must be numeric vectors of the same length.")
  }
  nlook <- length(frac)
  
  gqx <- glegx24
  gqw <- glegw24
  ngq <- length(gqx)

  ans <-
  .C("gsd_dens",
     frac=as.double(frac),
     xeff=as.double(x.eff),
     gqxw=as.double(c(gqx, gqw)),
     pngq=as.integer(ngq),
     pnlooks=as.integer(nlook),
     x=double(nlook*ngq),
     dF=double(nlook*ngq),
     x1c=double(ngq),
     dF1c=double(ngq),
     PACKAGE="PwrGSD")

  x <- matrix(ans$x, ngq, nlook)
  dF <- matrix(ans$dF, ngq, nlook)
  x1c <- ans$x1c
  dF1c <- ans$dF1c
  dimnames(x) <- dimnames(dF) <- list(1:ngq, 1:nlook)
  list(x=x, dF=dF, x1c=x1c, dF1c=dF1c)
}

"EX1gXK" <-
function(xk, b.eff, frac)
{
    nlook <- length(b.eff)

    x.eff <- frac^0.5 * b.eff
    fr.1 <- frac[1]
    fr.k <- frac[nlook]
    dP.Xk <- gsd.dens(x.eff, frac, scale="Brownian")
    ngq <- length(dP.Xk$x[,1])
    xk.vals <- with(dP.Xk, x[,nlook])
    dP.Xk.vals <- with(dP.Xk, dF[,nlook])
    dP.Xk.at.xk <- approx(xk.vals, dP.Xk.vals, xk)$y

    x1.vals <- dP.Xk$x1c
    dP.X1.vals <- dP.Xk$dF1c

    x.eff.dr1 <- x.eff[-1]
    frac.dr1 <- frac[-1] - fr.1
    dP.Xk.gX1.at.xkmx1 <- rep(0, ngq)
    for(i in 1:ngq)
    {
      dP.Xk.gX1 <- gsd.dens(x.eff.dr1 - x1.vals[i], frac.dr1, scale="Brownian")
      xk.vals <- with(dP.Xk.gX1, x[,nlook-1])
      dP.Xk.gX1.vals <-  with(dP.Xk.gX1, dF[,nlook-1])
      dP.Xk.gX1.at.xkmx1[i] <- approx(xk.vals, dP.Xk.gX1.vals, xk - x1.vals[i], yleft=0, yright=0)$y
    }
    ans <- sum(x1.vals * dP.Xk.gX1.at.xkmx1*dP.X1.vals)/sum(dP.Xk.gX1.at.xkmx1*dP.X1.vals)
    ans
}

"lookup" <-
function (xgrid, ygrid, x, y0 = 0)
{
        nx <- length(x)
        ngrid <- length(xgrid)
        ans <- .C("lookup",
        xgrid = as.double(xgrid),
        ygrid = as.double(ygrid),
        pngrid = as.integer(ngrid),
        x = as.double(x),
        pnx = as.integer(nx),
        py0 = as.double(y0),
        yatx = as.double(rep(0,nx)),
        index = as.integer(rep(0,nx)),
        PACKAGE = "PwrGSD")
        data.frame(x=ans$x,y=ans$yatx,index=ans$index)
}

"support" <- 
function(x)
{
	sort(unique(x))
}

"glue.2.string" <-
function(x, sep="")
{
	str <- x[1]
	n <- length(x)
	if (n>1)
	for(i in 2.:length(x))
		str <- str %,% sep %,% x[i]
	str
}

"html.table" <- 
function(x, file = "", append = FALSE, main = "", center = FALSE, html.tag = TRUE,
        table.attributes = "BORDER", column.attributes = c(
        row.header.column.attributes, rep(data.column.attributes, length = ncol(
        x))), cell.attributes = "", row.header.attributes = "ALIGN=RIGHT",
        row.header.column.attributes = "", data.column.attributes = "", ...)
{
        # formats data.frames, matrices, and vectors as html tables
        # also accepts a list of data.frames, matrices, and/or vectors
        # does not handle more sophisticated structures
        out <- c()
        # add <HTML> tag
        if(html.tag) out <- c(out, "<HTML>")
        # add <CENTER> tag
        if(center) out <- c(out, "<CENTER>")
        # if x is a list, iterate over elements
        # print list component names for each element
        # pad between elements with empty paragraph
        if(is.list(x) && !inherits(x, "data.frame")) {
                if(!is.null(main) && main != "")
                        out <- c(out, paste("<H2>", main, "</H2>"))
                n <- length(x)
                for(i in seq(length = n)) {
                        out <- c(out, html.table(x[[i]], main = names(x)[i]))
                        if(i < n)
                                out <- c(out, "<P> </P>")
                }
        }
        else {
                # if it's not a list, format it as a character matrix
                if(inherits(x, "data.frame")) {
                        dimnames.x <- dimnames(x)
                        dim.x <- dim(x)
                        x <- sapply(x, format, ...)
                        if(dim.x[1] <= 1) {
                                x <- matrix(x, nrow = dim.x[1], ncol = dim.x[
                                        2], byrow = TRUE)
                                dimnames(x) <- dimnames.x
                        }
                }
                else if(is.matrix(x)) {
                        dimnames.x <- dimnames(x)
                        dim.x <- dim(x)
                        x <- sapply(split(x, col(x)), format, ...)
                        if(dim.x[1] <= 1) {
                                x <- matrix(x, nrow = dim.x[1], ncol = dim.x[
                                        2], byrow = TRUE)
                                dimnames(x) <- dimnames.x
                        }
                }
                else {
                        dimnames.x <- list(names(x), NULL)
                        x <- as.matrix(format(x, ...))
                }
                out <- c(out, paste(collapse = " ", "<TABLE", paste(collapse =
                        " ", table.attributes), ">"))
                if(!is.null(main) && main != "")
                        out <- c(out, paste("<CAPTION> <H3>", main,
                                "</H3> </CAPTION>", collapse = " "))
                out <- c(out, paste("<TR>", paste(paste("<TH",
                        column.attributes, ">"), c(" ", dimnames.x[[2]]),
                        "</TH>", collapse = " "), "</TR>", collapse = " "))
                # data rows
                row.header.attributes <- rep(row.header.attributes, len = nrow(
                        x))
                if(is.matrix(cell.attributes)) {
                        if(!identical(dim(cell.attributes), dim(x)))
                                stop("If cell.attributes is a matrix it must " %,% 
                                     "have same dimensions as x")
                        for(i in seq(length = nrow(x))) {
                                out <- c(out, paste(paste("\n   <TR",
                                        row.header.attributes[i], ">"), "<TH>",
                                        dimnames.x[[1]][i], "</TH>", paste(
                                        sep = "", paste("\n      <TD",
                                        cell.attributes[i,  ], ">", sep = " "),
                                        x[i,  ], "</TD>", collapse = " "),
                                        "</TR>", collapse = " "))
                        }
                }
                else {
                        for(i in seq(length = nrow(x))) {
                                out <- c(out, paste(paste("\n   <TR",
                                        row.header.attributes[i], ">"), "<TH>",
                                        dimnames.x[[1]][i], "</TH>", paste(
                                        sep = "", paste("\n      <TD",
                                        cell.attributes, ">", sep = " "), x[
                                        i,  ], "</TD>", collapse = " "),
                                        "</TR>", collapse = " "))
                        }
                }
                out <- c(out, "</TABLE>")
        }
        if(center)
                out <- c(out, "</CENTER>")
        if(html.tag)
                out <- c(out, "</HTML>")
        # return vector of strings or write to file and return file name
        if(file == "") return(out) else {
                write(out, file = file, append = append)
                invisible(file)
        }
}

"encode" <-
function (x, basis)
{
        n <- length(basis)
        if(length(x)!=length(basis)) stop("lengths of 'x' and 'basis' must agree")
        sum(cumprod(c(1,basis[-n]))*x)
}

"decode" <- 
function (x, basis)
{
        n <- length(basis)
        if(length(x)!=1) stop("'x' must be of length 1")
        if(any(basis<=1)|length(basis)<=1) 
                stop("'basis' must be of length greater than 1, " %,%
        "and consist of elements greater than '1'")
        cpb <- cumprod(c(1,basis[-n]))
        ans <- NULL
        res <- x
        for (i in 1:n){
                new <- floor(res/cpb[n-i+1])
                res <- res - new*cpb[n-i+1]
                ans <- c(ans,new)
        }
        ans[n:1]
}

"Min" <- 
function (x, y)
{
    if (length(y) == 1)
        y <- rep(y, length(x))
    ans <- x
    ans[y < x] <- y[y < x]
    ans
}

"Max" <- 
function (x, y)
{
    if (length(y) == 1)
        y <- rep(y, length(x)) 
    ans <- x
    ans[y > x] <- y[y > x]
    ans
}

"RR2RCM" <- 
function (tlook, tcut.i, tcut.ii, h, rr, hOth, accru)
{
        h.i <- h
        h.ii <- h*rr
        mArmz(tlook, tcut.ii, h.ii, hOth, accru)/
        mArmz(tlook,  tcut.i,  h.i, hOth, accru)
}

"RCM2RR" <- 
function (tlook,tcut.i,h.i,hOth,accru,rcm)
{
        nlook <- length(tlook)
        fit <- list(objective=0)
        h.ii <- attr(fit$objective,"h") <- numeric(0)
        obj <-  function(theta,tlk,tcut.i,h.i,tcut.ii,h.ii,hOth,accru){
                nlk <- length(tlk)
                ans <- ((mArmz(tlk,tcut.ii,c(h.ii,exp(theta)),hOth,accru)/
                        mArmz(tlk,tcut.i,h.i,hOth,accru))[nlk] - rcm[nlk])^2
                attr(ans,"h") <- c(h.ii,exp(theta))
                ans
        }
        ans <- rep(0,nlook)
        for(i in 1:nlook){
                theta <- log(1e-12)
                tlk <- tlook[1:i]
                tcut.ii <- tcut.i[1]
                if(i>1) tcut.ii <- c(tcut.ii,tlook[1:(i-1)])
                h.ii <- attr(fit$objective,"h")
                fit <- optimize(f=obj,interval=c(log(1e-12),log(1e10)),tlk=tlk,
                        tcut.i=tcut.i,h.i=h.i,tcut.ii=tcut.ii,h.ii=h.ii,hOth=hOth,accru=accru)
                ans[i] <- exp(fit$min)
        }
        ans
}

"CRRtoRR" <- 
function(CRR, DT, h = NULL)
{
        RR <- CRR[1]
        m <- length(CRR)
        if(length(DT) != m)
                stop("lengths of 'CRR' and 'DT' must agree")
        if(missing(h))
                h <- rep(1, m)
        for(i in 2:m)
                RR <- c(RR, CRR[i] + sum(h[1:(i - 1)] * DT[1:(i - 1)] * (CRR[
                        i] - RR[1:(i - 1)]))/(h[i] * DT[i]))
        RR
}

"CY2TOShaz" <-
function(tcut, t.eor, m, verbose=FALSE)
{

  ti <- c(0, tcut)
  "Psi.1x" <-
  function(x, ti, h){
    L <- length(ti)
    H <- cumsum(DX(ti)*c(0,h))
    k.x <- max(which(x >= ti))
    ans <- 0
    if(x>ti[k.x])
      ans <- exp(-H[k.x])/h[k.x] * (1 - exp(-h[k.x] * (x - ti[k.x])))
    if(k.x >1)
      for(j in 1:(k.x-1)){
        ans <- ans + exp(-H[j])/h[j] * (1 - exp(-h[j] * (ti[j+1] -
        ti[j])))
    }
    ans
  }
  "Psi" <-
    function(x, ti, h){
      if(length(x)==1)
        ans <- Psi.1x(x, ti, h)
      if(length(x)>1)
        ans <- sapply(x, FUN=Psi.1x, ti=ti, h=h)
    ans
  }
  obj <-
  function(theta, ti, t.eor, m){
    x <- ti[-1]
    h <- exp(theta)
    CDF <- NULL
    if(any(x<=t.eor)){
      x.le.t <- x[x<=t.eor]
      CDF <- c(CDF, (x.le.t - Psi(x.le.t, ti, h))/t.eor)
    }
    if(any(x>t.eor)){
      x.gt.t <- x[x>t.eor]
      CDF <- c(CDF, 1 - (Psi(x.gt.t, ti, h) - Psi(x.gt.t - t.eor, ti,
               h))/t.eor)
    }
    n <- length(CDF)
    m.est <- 1 - (1-CDF)/(1-c(0,CDF[-n]))
    ans <- sum((m - m.est)^2)^0.5
    if(verbose) cat("\n objective: ",ans,"h: ",h,"\n")
    attr(ans,"tbl") <- cbind(x=x, CDF=CDF, m=m, m.est=m.est)
    ans
  }
  th0 <- log(-log(1-m)/DX(ti[-1]))
#  return(Psi.1x(tcut[1], ti, exp(th0)))
  fit <- optim(par=th0, fn=obj, method = "Nelder-Mead", ti=ti, t.eor=t.eor,
               m=m, control=list(maxit=10000))
  h <- exp(fit$par)
  obj. <- obj(fit$par, ti, t.eor, m)
  list(hazard=h, table=attr(obj., "tbl"))
}

"CDFOR2LRR" <-
function(tcut, tmax, h0, CDFOR)
{
  ti <- c(tcut, tmax)
  ind <- which(ti>0)
  ti <- ti[ind]
  m <- length(tcut)
  if(length(CDFOR)!=m)
    stop("lengths of 'tcut' and 'CDFR' must agree")
  DX <- function(x)c(x[1],diff(x))
  H0 <- cumsum(DX(ti) * h0)
  CDF.0 <- 1-exp(-H0)
  Odds.0 <- CDF.0/(1-CDF.0)
  Odds.1 <- CDFOR*Odds.0
  CDF.1 <- Odds.1/(1+Odds.1)
  beta <- log((DX(CDF.1)/(1-CDF.1))/(DX(CDF.0)/(1-CDF.0)))
  ans <- cbind(ti, beta)
  dimnames(ans) <- list(rep("", length(ti)), rep("", 2))
  ans
}

"mArmz" <- 
function (tlook, tcut, h, hOth, accru)
{
        fff <- function(tlk, tcut, h, hOth, accru){
                Min <- function (x, y)
                {
                    if (length(y) == 1)
                        y <- rep(y, length(x))
                    ans <- x
                    ans[y < x] <- y[y < x]
                    ans
                }
                Max <- function (x, y)
                {
                    if (length(y) == 1)
                        y <- rep(y, length(x))
                    ans <- x
                    ans[y > x] <- y[y > x]
                    ans
                }
                tER <- accru
                if(tlk - tER %in% tcut) tER <- tER-1e-10
                mu <- hOth
                ncut <- length(tcut)
                if(tcut[ncut]>tlk-tER) JtlkmtER <- max(min(which(tcut>tlk-tER))-1,0)
                else JtlkmtER <- ncut
                H <- cumsum(c(0,h[-ncut]*diff(tcut)))
                H.tlkmtER <- 0
                if(JtlkmtER>0) {
                        DT <- max(tlk-tER,0)-tcut[JtlkmtER]
                        H.tlkmtER <- H[JtlkmtER] + h[JtlkmtER]*DT
                }
                H. <- Max(H,H.tlkmtER)
                tcut. <- c(tcut,Inf)[-1]
                ans1 <- 0
                if(JtlkmtER>=1){
                        .b0. <- Min(tcut.,tlk-tER)
                        .a0. <- tcut
                        ans1 <- sum((h*exp(-(H+mu*.a0.))/(h+mu)*(1-exp(-(h+mu)*(.b0.-.a0.))))[1:JtlkmtER])
                }
                if(tcut[ncut]>tlk) Jtlk <- max(min(which(tcut>tlk))-1,0)
                else Jtlk <- ncut
                .c. <- Min(tcut.,tlk)
                .b. <- Max(tcut,tlk-tER)
                .a. <- tcut
                factor1 <- h*exp(-(H.+mu*.a.))/tER
                t1 <- (tlk - .b.)*exp(-(h+mu)*(.b.-.a.))/(h+mu)
                t2 <- (tlk - .c.)*exp(-(h+mu)*(.c.-.a.))/(h+mu)
                t3 <- exp(-(h+mu)*(.c.-.a.))/(h+mu)^2
                t4 <- exp(-(h+mu)*(.b.-.a.))/(h+mu)^2
                ans2 <- sum((factor1*(t1-t2+t3-t4))[JtlkmtER:Jtlk])
                ans1+ans2
        }
        sapply(tlook, FUN=fff, tcut=tcut, h=h, hOth=hOth, accru=accru)
}

"EX" <- function(x, h){
	nh <- length(h)
	dx <- diff(x)
	H <- cumsum(c(0,h[-nh]*dx))
	sum(exp(-H) * (1.0 - exp(-h*c(dx,1))*c(rep(1,nh-1),0))/h)
}

"thtilde.const" <- function(x,h,hA,hB,lA,lB)
{
        nx <- length(x)
	ans <- .C("htildeConst",
	x = as.double(x),
	nx = as.integer(nx),
	h = as.double(h),
	hA = as.double(hA),
	hB = as.double(hB),
	lA = as.double(lA),
	lB = as.double(lB),
	Stilde = as.double(rep(0,nx)),
	htlde = as.double(rep(0,nx)),
	PACKAGE = "PwrGSD")
}

"thtilde" <- 
function(x,xh,h,xhA,hA,xhB,hB,xlA,lA,xlB,lB,gradual=FALSE)
{
	.call. <- match.call()
	nx <- length(x)
	nh <- length(xh)
	nhA <- length(xhA)
	nhB <- length(xhB)
	nlA <- length(xlA)
	nlB <- length(xlB)

	bad <- FALSE
	bad.h  <- (length( h)!=nh)
	bad.hA <- (length(hA)!=nhA)
	bad.hB <- (length(hB)!=nhB)
	bad.lA <- (length(lA)!=nlA)
	bad.lB <- (length(lB)!=nlB)

	bads <- c(bad.h,bad.hA,bad.hB,bad.lA,bad.lB)
	bad <- any(bads)
	msg <- c("h","hA","hB","lA","lB")[which(bads)]

	glue.2.string <- function (x, sep = "")
	{
	    str <- x[1]
	    n <- length(x)
	    if (n > 1)
	        for (i in 2:length(x)) str <- str %,% sep %,% x[i]
	    str
	}
	if(bad) stop("Check that lengths of time and hazard agree in "%,% glue.2.string(msg[bads], sep=" "))

        glegx <- glegx24
        glegw <- glegw24

	ngq <- length(glegx)
        
	ans <- .C("htilde",
	x = as.double(x),
	nx = as.integer(nx),
	gqxw = as.double(c(glegx,glegw)),
	ngq = as.integer(ngq),
	xh = as.double(xh),
	h = as.double(h),
	nh = as.integer(nh),
	xhA = as.double(xhA),
	hA = as.double(hA),
	nhA = as.integer(nhA),
	xhB = as.double(xhB),
	hB = as.double(hB),
	nhB = as.integer(nhB),
	xlA = as.double(xlA),
	lA = as.double(lA),
	nlA = as.integer(nlA),
	xlB = as.double(xlB),
	lB = as.double(lB),
	nlB = as.integer(nlB),
	gradual = as.integer(gradual),
	ftlde = as.double(rep(0,nx)),
	Stlde = as.double(rep(0,nx)),
	htlde = as.double(rep(0,nx)),
	PACKAGE = "PwrGSD")
	list(htilde = ans$htlde, ftilde = ans$ftlde, Stilde = ans$Stlde, call=.call.)
}

"trandfromh" <- function(N, tcut, h)
{
	.C("trandfromh",
	pn = as.integer(N),
	tcut = as.double(tcut),
	h = as.double(h),
	pncut = as.integer(length(tcut)),
	t = as.double(rep(0,N)),
	PACKAGE = "PwrGSD")
}

"trandhcdtl" <- function(N, tcut, h, tcutxA, hxA, tdA, tcutxB, hxB, tdB)
{
	.C("trandhcdtl",
	pn = as.integer(N),
	tcut = as.double(tcut),
	h = as.double(h),
	pncut = as.integer(length(tcut)),
	tcutxA = as.double(tcutxA),
	hxA = as.double(hxA),
	pncutxA = as.integer(length(tcutxA)),
	tdA = as.double(tdA),
	tcutxB = as.double(tcutxB),
	hxB = as.double(hxB),
	pncutXB = as.integer(length(tcutxB)),
	tdB = as.double(tdB),
	code = as.integer(rep(0,N)),
	t = as.double(rep(0,N)),
	PACKAGE = "PwrGSD")
}

"wtdlogrank" <- 
function(formula = formula(data), data = parent.frame(), WtFun = c("FH", "SFH", "Ramp"),
         param = c(0, 0), sided = c(2, 1), subset, na.action, w=FALSE) 
{
    if (missing(sided))
        sided <- 2
    m <- .call. <- match.call()
    m[[1]] <- as.name("model.frame")
    m$param <- m$sided <- m$WtFun <- m$w <- NULL
    m <- eval(m, parent.frame())
    mt <- attr(m, "terms")
    if (is.empty.model(mt)) 
        stop("No treatment indicator specified in model")
    Arm <- model.matrix(mt, m)[, -1]
    nlev <- length(unique(Arm))

    R <- model.extract(m, "response")
    if(nlev>2 || class(R)!="Surv")
      stop("Argument 'formula' is expected to have a response of class \"Surv\" and " %,%
           "a single predictor containing two levels")
    ind.too.small <- (R[,1]<1e-10)
    n.too.small <- sum(ind.too.small)
    if(n.too.small > 0) {
      warning("Deleting " %,% n.too.small %,% "observatioins that are less than 10^-10\n")
      R <- R[-which(ind.too.small),]
      Arm <- Arm[-which(ind.too.small)]
    }    
    TOS <- R[, 1]
    Event <- R[, 2]
    ntimes <- length(unique(TOS[Event!=0]))
    n <- length(Event)
    if(missing(WtFun)){
      param <- c(0,0)
      WtFun<-"FH"
    }
    wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1
    
    ans <- .C("WtdLogRank",
              TOS = as.double(TOS),
              Event = as.integer(Event), 
              Arm = as.integer(Arm),
              pn = as.integer(n),
              wttyp = as.integer(wttyp), 
              par = as.double(param),
              time = as.double(rep(0, ntimes)), 
              nrisk = as.integer(rep(0, 2*ntimes)),
              nevent = as.integer(rep(0, 2*ntimes)),
              wt = double(ntimes),              
              pntimes = as.integer(ntimes),              
              stat = double(1),
              var = double(1),
              m1 = double(1),
              UQt = double(ntimes),
              varQt = double(ntimes),
              var1t = double(ntimes),
              PACKAGE = "PwrGSD")

    out <- ans
    out$TOS <- out$Event <- out$Arm <- NULL
    out$time <- ans$time
    out$nrisk <- ans$nrisk[2*(1:ntimes) - 1] + ans$nrisk[2*(1:ntimes)]
    out$nevent <- ans$nevent[2*(1:ntimes)-1] + ans$nevent[2*(1:ntimes)]
    out$nrisk1 <- ans$nrisk[2*(1:ntimes)]
    out$nevent1 <- ans$nevent[2*(1:ntimes)]
    out$wt <- ans$wt
    out$pu0 <- sum(TOS[Arm==0])
    out$pu1 <- sum(TOS[Arm==1])
    sA <- support(Arm)
    out$n0 <- sum(Arm == sA[1])
    out$n1 <- sum(Arm == sA[2])
    out$n <- out$n0 + out$n1
    out$Z <- ans$stat/ans$var^0.5
    out$sided <- sided
    class(out) <- "survtest"
    out$call <- .call.
    out
}

"IntSurvDiff" <- 
function(formula = formula(data), data = parent.frame(), WtFun = c("FH", "SFH", "Ramp"),
         param = c(0, 0), sided = c(2, 1), subset, na.action, w=FALSE) 
{
    if (missing(sided))
        sided <- 2
    m <- .call. <- match.call()
    m[[1]] <- as.name("model.frame")
    m$param <- m$sided <- m$WtFun <- m$w <- NULL
    m <- eval(m, parent.frame())
    mt <- attr(m, "terms")
    if (is.empty.model(mt)) 
        stop("No treatment indicator specified in model")
    Arm <- model.matrix(mt, m)[, -1]
    nlev <- length(unique(Arm))

    R <- model.extract(m, "response")
    if(nlev>2 || class(R)!="Surv")
      stop("Argument 'formula' is expected to have a response of class \"Surv\" and " %,%
           "a single predictor containing two levels")
    ind.too.small <- (R[,1]<1e-10)
    n.too.small <- sum(ind.too.small)
    if(n.too.small > 0) {
      warning("Deleting " %,% n.too.small %,% "observatioins that are less than 10^-10\n")
      R <- R[-which(ind.too.small),]
      Arm <- Arm[-which(ind.too.small)]
    }    
    TOS <- R[, 1]
    Event <- R[, 2]
    ntimes <- length(unique(TOS[Event!=0]))
    n <- length(Event)
    if(missing(WtFun)){
      param <- c(0,0)
      WtFun<-"FH"
    }
    wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1

    ans <- .C("IntSurvDiff",
              TOS = as.double(TOS),
              Event = as.integer(Event), 
              Arm = as.integer(Arm),
              pn = as.integer(n),
              wttyp = as.integer(wttyp), 
              par = as.double(param),
              time = as.double(rep(0, ntimes)), 
              nrisk = as.integer(rep(0, 2*ntimes)),
              nevent = as.integer(rep(0, 2*ntimes)),
              pntimes = as.integer(ntimes),              
              stat = as.double(0),
              var = as.double(0),
              wt = as.double(rep(0, ntimes)),
              PACKAGE = "PwrGSD")
    out <- ans
    out$TOS <- out$Event <- out$Arm <- NULL
    out$time <- ans$time
    out$nrisk <- ans$nrisk[2*(1:ntimes) - 1] + ans$nrisk[2*(1:ntimes)]
    out$nevent <- ans$nevent[2*(1:ntimes)-1] + ans$nevent[2*(1:ntimes)]
    out$nrisk1 <- ans$nrisk[2*(1:ntimes)]
    out$nevent1 <- ans$nevent[2*(1:ntimes)]
    out$wt <- ans$wt
    out$pu0 <- sum(TOS[Arm==0])
    out$pu1 <- sum(TOS[Arm==1])
    sA <- support(Arm)
    out$n0 <- sum(Arm == sA[1])
    out$n1 <- sum(Arm == sA[2])
    out$n <- out$n0 + out$n1
    out$Z <- ans$stat/ans$var^0.5
    out$sided <- sided
    class(out) <- "survtest"
    out$call <- .call.
    out
  }

"print.survtest" <- function(x, ...) 
{
    sided <- x$sided
    Z <- x$stat/x$var^0.5
    pval <- (2^(sided == 2)) * (1 - pnorm(abs(Z)))
    cat("call:\n")
    print(x$call)
    z.nm <- as.character(x$call$formula[[3]])
    n0 <- x$n0
    n1 <- x$n1
    pu0 <- x$pu0
    pu1 <- x$pu1
    obs1 <- sum(x$nevent1)
    obs0 <- sum(x$nevent) - obs1
    exp0 <- sum(x$nevent * (1 - x$nrisk1/x$nrisk))
    exp1 <- sum(x$nevent * x$nrisk1/x$nrisk)
    tbl <- data.frame(n = c(n0, n1), pu = c(pu0, pu1), Obs = c(obs0, obs1), Exp = c(exp0, 
        exp1))
    dimnames(tbl) <- list(z.nm %,% c("=0:", "=1:"), names(tbl))
    print(tbl)
    ans <- c(Z, pval)
    names(ans) <- c("z", "p-val")
    print(ans)
    invisible(x)
}

"summary.survtest" <- function(object, ...) 
{
    sided <- object$sided
    Z <- object$stat/object$var^0.5
    pval <- (2^(sided == 2)) * (1 - pnorm(abs(Z)))
    cat("call:\n")
    print(object$call)
    z.nm <- as.character(object$call$formula[[3]])
    n0 <- object$n0
    n1 <- object$n1
    pu0 <- object$pu0
    pu1 <- object$pu1
    obs1 <- sum(object$nevent1)
    obs0 <- sum(object$nevent) - obs1
    exp0 <- sum(object$nevent * (1 - object$nrisk1/object$nrisk))
    exp1 <- sum(object$nevent * object$nrisk1/object$nrisk)
    tbl <- data.frame(n = c(n0, n1), pu = c(pu0, pu1), Obs = c(obs0, obs1), Exp = c(exp0, 
        exp1))
    dimnames(tbl) <- list(z.nm %,% c("=0:", "=1:"), names(tbl))
    object$tbl <- tbl
    ans <- c(Z, pval)
    names(ans) <- c("z", "p-val")
    object$ans
    object
}

"CondPower" <-
function(Z, frac, drift, drift.end, err.I, sided=1)
{
  mu.c.H0 <- Z * frac^0.5
  mu.c.HA <- mu.c.H0 + drift.end - drift
  rt <- drift.end >0
  z.c <- rt * qnorm(1 - err.I/sided) + (1 - rt) * qnorm(err.I/sided)
  Pr.cond.typeIIerr <- rt * pnorm((z.c - mu.c.HA)/(1 - frac)^0.5) +
    (1 - rt) * (1-pnorm((z.c - mu.c.HA)/(1 - frac)^0.5))
  Pr.cond.typeIerr <- rt * (1 - pnorm((z.c - mu.c.H0)/(1 -
        frac)^0.5)) + (1 - rt) * pnorm((z.c - mu.c.H0)/(1 - frac)^0.5)
  cbind(Pr.cond.typeIerr=Pr.cond.typeIerr, Pr.cond.typeIIerr=Pr.cond.typeIIerr)
}

"stpplt" <-
function (x, y, bw, stars, ...) 
{
    .call. <- match.call(expand.dots=FALSE)
    nms.dots <- names(.call.$...)   
    d.y <- dim(y)
    n <- d.y[1]
    d <- d.y[2]
    if(!bw)
    {
      cls <- list("magenta", "aquamarine")
      density <- NULL
      angle <- list(NULL, NULL)
    }
    if(bw==1)
    {
      cls <- list("grey90", NULL)
      density <- NULL
      angle <- list(NULL, NULL)
    }
    if(bw==2)
    {
      cls <- list(NULL, NULL)
      density <- diff(range(x))*10
      angle <- list(45, 135)
    }
    
    polygon(c(0, x, 1), c(0, y[, 1], 0), border = NA, col = cls[[1]],
            density=density, angle=angle[[1]])
    for (j in 2:d) polygon(c(x[n:1], x), c(y[n:1, j - 1], y[, 
        j]), border = NA, col = cls[[2 - (j%%2)]],
            density=density, angle=angle[[2 - (j%%2)]])
    polygon(c(x, 0), c(y[, d], 1), border = NA, col = cls[[2 - 
        ((d + 1)%%2)]], density=density, angle=angle[[2 - ((d + 1)%%2)]])
    lines.cmd <- as.call(expression(lines,x=x))
    ln.dots <- .call.$...

    if(!is.null(stars))
      ln.dots <- ln.dots[-which(nms.dots=="pch")]

    if(is.null(stars))
      stars <- FALSE

    lines.cmd$... <- ln.dots
    for(j in 1:d)
    {
      lines.cmd$type <- NULL
      lines.cmd$pch <- NULL
      if(stars && j==(2*stars+1))
      {
        lines.cmd$type <- "b"
        lines.cmd$pch <- 8
      }
    
      lines.cmd$y <- y[,j]
      eval(lines.cmd)
    }
}

"plot.cpd.PwrGSD" <-
function(x, formula, subset, na.action, ...)
{
  ow <- options("warn")
  options(warn = -1)
  .call. <- match.call(expand.dots=FALSE)
  dots <- .call.$...
  nms.dots <- names(dots)
  given.values <- rows <- columns <- show.given <- col <-
             pch <- bar.bg <- fac <- xlab <- ylab <- subscripts <- axlabels <-
             number <- overlap <- xlim <- ylim <- stars <- margin <- marFUN <- NULL
  bw <- FALSE
  ind.dots.dr <- NULL
  if("given.values" %in% nms.dots)
  {
    given.values <- dots[["given.values"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="given.values"))
  }
  if("rows" %in% nms.dots)
  {
    rows <- dots[["rows"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="rows"))
  }
  if("columns" %in% nms.dots)
  {
    columns <- dots[["columns"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="columns"))
  }
  if("show.given" %in% nms.dots)
  {
    show.given <- dots[["show.given"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="show.given"))
  }
  if("col" %in% nms.dots)
  {
    col <- dots[["col"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="col"))
  }
  if("pch" %in% nms.dots)
  {
    pch <- dots[["pch"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="pch"))
  }
  if("bar.bg" %in% nms.dots)
  {
    bar.bg <- dots[["bar.bg"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="bar.bg"))
  }
  if("fac" %in% nms.dots)
  {
    fac <- dots[["fac"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="fac"))
  }
  if("xlab" %in% nms.dots)
  {
    xlab <- dots[["xlab"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="xlab"))
  }
  if("ylab" %in% nms.dots)
  {
    ylab <- dots[["ylab"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="ylab"))
  }
  if("subscripts" %in% nms.dots)
  {
    subscripts <- dots[["subscripts"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="subscripts"))
  }
  if("axlabels" %in% nms.dots)
  {
    axlabels <- dots[["axlabels"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="axlabels"))
  }
  if("number" %in% nms.dots)
  {
    number <- dots[["number"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="number"))
  }
  if("overlap" %in% nms.dots)
  {
    overlap <- dots[["overlap"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="overlap"))
  }
  if("xlim" %in% nms.dots)
  {
    xlim <- dots[["xlim"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="xlim"))
  }
  if("ylim" %in% nms.dots)
  {
    ylim <- dots[["ylim"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="ylim"))
  }
  if("bw" %in% nms.dots)
  {
    bw <- dots[["bw"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="bw"))
  }
  if("stars" %in% nms.dots)
  {
    stars <- dots[["stars"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="stars"))
  }
  if("margin" %in% nms.dots)
  {
    margin <- dots[["margin"]]
    ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="margin"))
    if("marFUN" %in% nms.dots)
    {
      marFUN <- dots[["marFUN"]]
      ind.dots.dr <- c(ind.dots.dr, which(nms.dots=="marFUN"))
    }
    if(!("marFUN" %in% nms.dots))
    {
      marFUN <- mean
    }
  }
   
  if(!is.null(ind.dots.dr)) dots <- dots[-ind.dots.dr]
  .call.[[1]] <- as.name("costopplot")
  .call.$data <- as.call(expression(as.data.frame.cpd.PwrGSD))
  .call.$data$x <- .call.$x
  .call.$x <- NULL
  .call.$given.values <- given.values
  .call.$rows  <- rows
  .call.$columns <- columns
  .call.$show.given <- show.given
  .call.$col <- col
  .call.$pch <- pch
  .call.$bar.bg <- bar.bg
  .call.$fac <- fac
  .call.$xlab <- xlab
  .call.$ylab <- ylab
  .call.$subscripts <- subscripts
  .call.$axlabels <- axlabels
  .call.$number <- number
  .call.$overlap <- overlap
  .call.$xlim <- xlim 
  .call.$ylim <- ylim
  .call.$bw <- bw
  .call.$stars <- stars
  .call.$margin <- margin
  .call.$marFUN <- marFUN
  .call.$... <- dots
  eval(.call.)
  options(ow)
}

"costopplot" <- 
function (formula, data, given.values, rows, columns, show.given = TRUE, 
          col = par("fg"), pch = par("pch"), bar.bg = c(num = gray(0.8), fac = gray(0.95)),
          xlab = c(x.name, paste("Given :", a.name)),
          ylab = c("Type II Error Prob "%,%eIIlab%,%" & Power "%,%powlab, "Given : "%,% b.name),
          subscripts = FALSE, axlabels = function(f) abbreviate(levels(f)),
          number = 6, overlap = c(0,0), xlim, ylim, subset, na.action, bw=FALSE,
          stars=NULL, margin=NULL, marFUN=NULL,...) 
{
    .call. <- match.call(expand.dots=FALSE)
    if(missing(bw)) bw <- FALSE
    if(missing(stars)) stars <- NULL
    is.margin <- !missing(margin)
    if(!is.margin)
    {
      margin <- NULL
      marFUN <- NULL
    }
    if(!bw)
    {
      eIIlab <- "(Magenta)"
      powlab <- "(Aqua)"
    }
    if(bw==1)
    {
      eIIlab <- "(Grey)"
      powlab <- "(White)"
    }
    if(bw==2)
    {
      eIIlab <- "[//]"
      powlab <- "[\\\\]"
    }
    formula.new <- (y ~ x)
    formula.new[[3]] <- formula[[2]]
    formula.new[[2]] <- (I(cbind(dE, dP)) ~ x)[[2]]
    formula <- formula.new
    is.subset <- !missing(subset)
    is.na.action <- !missing(na.action)
    .call. <- match.call()
    .call.$formula <- formula
    deparen <- function(expr) {
        while (is.language(expr) && !is.name(expr) && deparse(expr[[1]]) == 
            "(") expr <- expr[[2]]
        expr
    }
    bad.formula <- function() stop("invalid conditioning formula")
    bad.lengths <- function() stop("incompatible variable lengths")
    formula <- deparen(formula)
    if (!inherits(formula, "formula")) 
        bad.formula()
    y <- deparen(formula[[2]])
    rhs <- deparen(formula[[3]])
    if (deparse(rhs[[1]]) != "|") 
        bad.formula()
    x <- deparen(rhs[[2]])
    rhs <- deparen(rhs[[3]])
    if (is.language(rhs) && !is.name(rhs) && (deparse(rhs[[1]]) == 
        "*" || deparse(rhs[[1]]) == "+")) {
        have.b <- TRUE
        a <- deparen(rhs[[2]])
        b <- deparen(rhs[[3]])
    }
    else {
        have.b <- FALSE
        a <- rhs
    }
    if (missing(data))
        data <- parent.frame()
    form <- (yv ~ xv + av)
    form[[2]] <- y
    form[[3]][[2]] <- x
    form[[3]][[3]] <- a
    if (have.b) {
        form <- (yv ~ xv + av * bv)
        form[[2]] <- y
        form[[3]][[2]] <- x
        form[[3]][[3]][[2]] <- a
        form[[3]][[3]][[3]] <- b
    }
    mdl <- as.call(expression(model.frame))
    mdl$formula <- form
    mdl$nlook <- as.name("nlook")
    mdl$data <- data
    if (is.subset) 
        mdl$subset <- .call.$subset
    if (is.na.action) 
        mdl$na.action <- .call.$na.action
    mdl <- eval(mdl, parent.frame())
    if(is.margin)
    {
      nms <- names(mdl)
      n.mdl <- nrow(mdl)
      margin.vars <- eval(.call.$margin)
      n.mv <- length(margin.vars)
      idx.vars <- setdiff(nms[-1], margin.vars)
      n.iv <- length(idx.vars)

      IND.vars <- list()
      cls <- NULL
      for(i in 1:n.iv)
      {
        iv.i <- idx.vars[i]
        X <- with(mdl, get(iv.i))
        cls <- c(cls, class(X))
        IND.vars[[iv.i]] <- as.factor(X)
      }
      names(cls) <- idx.vars
      
      tmp <- IND.vars[[1]]
      IND.vars[[1]] <- IND.vars[[n.iv]]
      IND.vars[[n.iv]] <- tmp

      tmp <- idx.vars[1]
      idx.vars[1] <- idx.vars[n.iv]
      idx.vars[n.iv] <- tmp

      tmp <- cls[1]
      cls[1] <- cls[n.iv]
      cls[n.iv] <- tmp
      names(cls) <- idx.vars
      
      tmp <- by(mdl[[1]], INDICES=IND.vars, FUN=function(x)apply(x, 2, FUN=marFUN))
      n.resp <- 2
      n.mgn <- n.mv
      d.tmp <- dim(tmp)
      dn.tmp <- dimnames(tmp)
      tmp <- unlist(tmp)
      n.tmp <- length(tmp)
      nr <- n.tmp/n.resp
      tmp <- t(matrix(tmp, n.resp, nr))
      idx <- 1:nr
      mdl <- data.frame(idx)
      nreps <- 1
      for(i in 1:n.iv)
      {
        lbls <- dn.tmp[[i]]
        if(cls[i]=="numeric") lbls <- as.numeric(lbls)
        mdl[[idx.vars[i]]]  <- lbls[gl(d.tmp[i], nreps, prod(d.tmp))]
        nreps <- nreps*d.tmp[i]
      }
      mdl <- cbind(mdl, tmp)
      nc <- ncol(mdl)
      names(mdl)[1] <- "index"
      names(mdl)[(nc-1):nc] <- c("dE","dP")
      mdl <- model.frame(I(cbind(dE,dP))~., data=mdl)
      mdl.nms <- names(mdl)
      mdl.nms[which(mdl.nms=="(nlook)")] <- "nlook"
      names(mdl) <- mdl.nms
    }
    y.name <- deparse(y)
    y <- model.extract(mdl, "response")
    nlook <- attr(data, "detail")["nlook"]
    d.y <- dim(y)
    ntot.orig <- d.y[1]
    n.conds <- ntot.orig/nlook
    y.new <- matrix(0, n.conds, 2 * nlook)
    for (j in 1:n.conds) {
        y.new[j, 2 * (1:nlook) - 1] <- y[nlook * (j - 1) + (1:nlook), 
            1]
        y.new[j, 2 * (1:nlook)] <- y[nlook * (j - 1) + (1:nlook), 
            2]
        y.new[j, ] <- cumsum(y.new[j, ])
    }
    y <- y.new
    mdl <- mdl[nlook * (0:((ntot.orig - 1)/nlook)) + 1, ]
    x.name <- deparse(x)
    x <- mdl[, x.name]
    nobs <- length(x)
    d.y <- c(nobs, length(y)/nobs)
    d <- d.y[2]
    if (length(y)%%nobs) 
        bad.lengths()
    a.name <- deparse(a)
    a <- as.factor(mdl[, a.name])
    if (length(a) != nobs) 
        bad.lengths()
    if (is.character(a)) 
        a <- as.factor(a)
    a.is.fac <- is.factor(a)
    if (have.b) {
        b.name <- deparse(b)
        b <- as.factor(mdl[, b.name])
        if (length(b) != nobs) 
            bad.lengths()
        if (is.character(b)) 
            b <- as.factor(b)
        b.is.fac <- is.factor(b)
        missingrows <- which(is.na(x) | is.na(y) | is.na(a) | 
            is.na(b))
    }
    else {
        missingrows <- which(is.na(x) | is.na(y) | is.na(a))
        b <- NULL
        b.name <- ""
    }
    number <- as.integer(number)
    if (length(number) == 0 || any(number < 1)) 
        stop("number must be integer >= 1")
    if (any(overlap >= 1)) 
        stop("overlap must be < 1 (and typically >= 0).")
    bad.givens <- function() stop("invalid given.values")
    if (missing(given.values)) {
        a.intervals <- if (a.is.fac) {
            i <- seq(along = a.levels <- levels(a))
            a <- as.numeric(a)
            cbind(i - 0.5, i + 0.5)
        }
        else co.intervals(a, number = number[1], overlap = overlap[1])
        b.intervals <- if (have.b) {
            if (b.is.fac) {
                i <- seq(along = b.levels <- levels(b))
                b <- as.numeric(b)
                cbind(i - 0.5, i + 0.5)
            }
            else {
                if (length(number) == 1) 
                  number <- rep.int(number, 2)
                if (length(overlap) == 1) 
                  overlap <- rep.int(overlap, 2)
                co.intervals(b, number = number[2], overlap = overlap[2])
            }
        }
    }
    else {
        if (!is.list(given.values)) 
            given.values <- list(given.values)
        if (length(given.values) != (if (have.b) 
            2
        else 1)) 
            bad.givens()
        a.intervals <- given.values[[1]]
        if (a.is.fac) {
            a.levels <- levels(a)
            if (is.character(a.intervals)) 
                a.intervals <- match(a.intervals, a.levels)
            a.intervals <- cbind(a.intervals - 0.5, a.intervals + 
                0.5)
            a <- as.numeric(a)
        }
        else if (is.numeric(a)) {
            if (!is.numeric(a.intervals)) 
                bad.givens()
            if (!is.matrix(a.intervals) || ncol(a.intervals) != 
                2) 
                a.intervals <- cbind(a.intervals - 0.5, a.intervals + 
                  0.5)
        }
        if (have.b) {
            b.intervals <- given.values[[2]]
            if (b.is.fac) {
                b.levels <- levels(b)
                if (is.character(b.intervals)) 
                  b.intervals <- match(b.intervals, b.levels)
                b.intervals <- cbind(b.intervals - 0.5, b.intervals + 
                  0.5)
                b <- as.numeric(b)
            }
            else if (is.numeric(b)) {
                if (!is.numeric(b.intervals)) 
                  bad.givens()
                if (!is.matrix(b.intervals) || ncol(b.intervals) != 
                  2) 
                  b.intervals <- cbind(b.intervals - 0.5, b.intervals + 
                    0.5)
            }
        }
    }
    if (any(is.na(a.intervals)) || (have.b && any(is.na(b.intervals)))) 
        bad.givens()
    if (have.b) {
        rows <- nrow(b.intervals)
        columns <- nrow(a.intervals)
        nplots <- rows * columns
        if (length(show.given) < 2) 
            show.given <- rep.int(show.given, 2)
    }
    else {
        nplots <- nrow(a.intervals)
        if (missing(rows)) {
            if (missing(columns)) {
                rows <- ceiling(round(sqrt(nplots)))
                columns <- ceiling(nplots/rows)
            }
            else rows <- ceiling(nplots/columns)
        }
        else if (missing(columns)) 
            columns <- ceiling(nplots/rows)
        if (rows * columns < nplots) 
            stop("rows * columns too small")
    }
    total.columns <- columns
    total.rows <- rows
    f.col <- f.row <- 1
    if (show.given[1]) {
        total.rows <- rows + 1
        f.row <- rows/total.rows
    }
    if (have.b && show.given[2]) {
        total.columns <- columns + 1
        f.col <- columns/total.columns
    }
    mar <- if (have.b) 
        rep.int(0, 4)
    else c(0.5, 0, 0.5, 0)
    oma <- c(5, 6, 5, 4)
    if (have.b) {
        oma[2] <- 5
        if (!b.is.fac) 
            oma[4] <- 5
    }
    if (a.is.fac && show.given[1]) 
        oma[3] <- oma[3] - 1
    opar <- par(mfrow = c(total.rows, total.columns), oma = oma, 
        mar = mar, xaxs = "r", yaxs = "r", new = FALSE)
    on.exit(par(opar))
    plot.new()
    if (missing(xlim)) 
        xlim <- range(as.numeric(x), finite = TRUE)
    if (missing(ylim)) 
        ylim <- range(as.numeric(y), finite = TRUE)
    pch <- rep(pch, length.out = nobs)
    col <- rep(col, length.out = nobs)
    do.panel <- function(index, subscripts = FALSE, id) {
        Paxis <- function(side, x) {
            if (nlevels(x)) {
                lab <- axlabels(x)
                axis(side, labels = lab, at = seq(lab), xpd = NA)
            }
            else axis(side, xpd = NA)
        }
        istart <- (total.rows - rows) + 1
        i <- total.rows - ((index - 1)%/%columns)
        j <- (index - 1)%%columns + 1
        par(mfg = c(i, j, total.rows, total.columns))
        plot.new()
        plot.window(xlim, ylim)
        if (any(is.na(id))) 
            id[is.na(id)] <- FALSE
        if (any(id)) {
            n.id <- sum(id)
            grid(lty = "solid")
            if (subscripts) 
                stpplt(x[id], y[id, ], subscripts = id, bw=bw, stars=stars, ...)
            else stpplt(x[id], y[id, ], bw=bw, stars=stars, ...)
        }
        if ((i == total.rows) && (j%%2 == 0)) 
            Paxis(1, x)
        else if ((i == istart || index + columns > nplots) && 
            (j%%2 == 1)) 
            Paxis(3, x)
        if ((j == 1) && ((total.rows - i)%%2 == 0)) 
            Paxis(2, y)
        else if ((j == columns || index == nplots) && ((total.rows - 
            i)%%2 == 1)) 
            Paxis(4, y)
        box()
    }
    if (have.b) {
        count <- 1
        for (i in 1:rows) {
            for (j in 1:columns) {
                id <- ((a.intervals[j, 1] <= a) & (a <= a.intervals[j, 
                  2]) & (b.intervals[i, 1] <= b) & (b <= b.intervals[i, 
                  2]))
                do.panel(count, subscripts, id)
                count <- count + 1
            }
        }
    }
    else {
        for (i in 1:nplots) {
            id <- ((a.intervals[i, 1] <= a) & (a <= a.intervals[i, 
                2]))
            do.panel(i, subscripts, id)
        }
    }
    mtext(xlab[1], side = 1, at = 0.5 * f.col, outer = TRUE, 
        line = 3.5, xpd = NA)
    mtext(ylab[1], side = 2, at = 0.5 * f.row, outer = TRUE, 
        line = 3.5, xpd = NA)
    if (length(xlab) == 1) 
        xlab <- c(xlab, paste("Given :", a.name))
    if (show.given[1]) {
        par(fig = c(0, f.col, f.row, 1), mar = mar + c(3 + (!a.is.fac), 
            0, 0, 0), new = TRUE)
        plot.new()
        nint <- nrow(a.intervals)
        a.range <- range(a.intervals, finite = TRUE)
        plot.window(a.range + c(0.03, -0.03) * diff(a.range), 
            0.5 + c(0, nint))
        rect(a.intervals[, 1], 1:nint - 0.3, a.intervals[, 2], 
            1:nint + 0.3, col = bar.bg[if (a.is.fac) 
                "fac"
            else "num"])
        if (a.is.fac) {
            text(apply(a.intervals, 1, mean), 1:nint, a.levels)
        }
        else {
            axis(3, xpd = NA)
            axis(1, labels = FALSE)
        }
        box()
        mtext(xlab[2], 3, line = 3 - a.is.fac, at = mean(par("usr")[1:2]), 
            xpd = NA)
    }
    else {
        mtext(xlab[2], 3, line = 3.25, outer = TRUE, at = 0.5 * 
            f.col, xpd = NA)
    }
    if (have.b) {
        if (length(ylab) == 1) 
            ylab <- c(ylab, paste("Given :", b.name))
        if (show.given[2]) {
            par(fig = c(f.col, 1, 0, f.row), mar = mar + c(0, 
                3 + (!b.is.fac), 0, 0), new = TRUE)
            plot.new()
            nint <- nrow(b.intervals)
            b.range <- range(b.intervals, finite = TRUE)
            plot.window(0.5 + c(0, nint), b.range + c(0.03, -0.03) * 
                diff(b.range))
            rect(1:nint - 0.3, b.intervals[, 1], 1:nint + 0.3, 
                b.intervals[, 2], col = bar.bg[if (b.is.fac) 
                  "fac"
                else "num"])
            if (b.is.fac) {
                text(1:nint, apply(b.intervals, 1, mean), b.levels, 
                  srt = 90)
            }
            else {
                axis(4, xpd = NA)
                axis(2, labels = FALSE)
            }
            box()
            mtext(ylab[2], 4, line = 3 - b.is.fac, at = mean(par("usr")[3:4]), 
                xpd = NA)
        }
        else {
            mtext(ylab[2], 4, line = 3.25, at = 0.5 * f.row, 
                outer = TRUE, xpd = NA)
        }
    }
    if (length(missingrows) > 0) {
        cat("\nMissing rows:", missingrows, "\n")
        invisible(missingrows)
    }
}

"agghaz" <- function(t.agg, time, nrisk, nevent)  
{
  n.t <- length(time)
  n.blocks <- ncol(nevent)
  result <-.C("agghaz",
              tagg=as.double(t.agg),
              time=as.double(time),
              nrisk=as.integer(nrisk),
              nevent=as.integer(nevent),
              pndth=as.integer(n.t),
              pnb=as.integer(n.blocks),
              timea=as.double(rep(0,n.t)),
              nriska=as.integer(rep(0,n.blocks*n.t)),
              neventa=as.integer(rep(0,n.blocks*n.t)),
              pnagg=as.integer(0),
              PACKAGE="PwrGSD")

  n.agg <- result$pnagg
  time.a <- result$timea[1:n.agg]
  nrisk.a <- matrix(result$nriska[1:(n.blocks*n.agg)], n.agg, n.blocks)
  nevent.a <- matrix(result$neventa[1:(n.blocks*n.agg)], n.agg, n.blocks)
  out <- as.data.frame(cbind(time.a, nrisk.a, nevent.a))
  names(out) <- c("time", "nrisk" %,% (1:n.blocks), "nevent" %,% (1:n.blocks))
  out
}

"mystack" <- 
function (object, fu.vars, create.idvar = FALSE) 
{
    d.object <- dim(object)
    n <- d.object[1]
    p <- d.object[2]
    for(k in 1:p)
      if(is.character(object[,k]))
        object[,k] <- as.factor(object[,k])

    if (create.idvar) 
        object$id <- 1:n
    out <- NULL
    nms <- names(object)
    n.fu.vars <- length(fu.vars)
    n.fus <- length(grep(fu.vars[1], nms))
    fu.var.pos <- matrix(0, n.fus, n.fu.vars)
    for (k in 1:n.fu.vars) fu.var.pos[, k] <- grep(fu.vars[k], 
        nms)
    base.vars <- object[, -c(fu.var.pos)]
    nbv <- dim(base.vars)[2]
    bvfactors <- NULL
    bvfactorlevs <- list()
    l <- 1
    for (k in 1:nbv) {
        isf <- is.factor(base.vars[, k])
        if (isf) {
            bvfactors <- c(bvfactors, k)
            bvfactorlevs[[l]] <- levels(base.vars[, k])
            base.vars[,k] <- as.numeric(base.vars[,k])
            l <- l + 1
        }
    }
    base.vars <- as.matrix(base.vars)
    nbvf <- l - 1
    fuvars <- object[, fu.var.pos]
    fufactors <- NULL
    fufactorlevs <- list()
    l <- 1
    for (k in 1:n.fu.vars) {
        isf <- is.factor(fuvars[, n.fus * (k - 1) + 1])
        if (isf) {
            fufactors <- c(fufactors, k)
            fufactorlevs[[l]] <- levels(fuvars[, n.fus * (k - 
                1) + 1])
            fuvars <- as.numeric(fuvars[,k])
        }
    }
    fuvars <- as.matrix(fuvars)
    nfuf <- l - 1
    out <- 
    .C("mystack", pn = as.integer(n), pnfus = as.integer(n.fus), 
        pnfuvars = as.integer(n.fu.vars), pnbasevars = as.integer(nbv), 
        basevars = as.double(base.vars), fuvars = as.double(fuvars), 
        out = double(n.fus * n * (nbv + 1 + n.fu.vars)), NAOK = TRUE, 
        PACKAGE = "PwrGSD")$out
    out <- as.data.frame(matrix(out, n.fus * n, nbv + 1 + n.fu.vars))
    names(out) <- c(nms[-fu.var.pos], "fu", fu.vars)
    if (nbvf > 0) 
        for (l in 1:nbvf) out[, bvfactors[l]] <- bvfactorlevs[[l]][out[,bvfactors[l]]]
    if (nfuf > 0) 
        for (l in 1:nfuf) out[, nbv + 1 + fufactors[l]] <- fufactorlevs[[l]][out[,nbv + 1 + fufactors[l]]]
    out
}

.onAttach <- function(libname, pkgname)
{
    options(stringsAsFactors=FALSE)
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
    packageStartupMessage(paste(pkgname, ver))
}
