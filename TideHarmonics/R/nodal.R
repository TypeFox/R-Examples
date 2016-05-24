### DEFINE GLOBAL VARIABLES ###

utils::globalVariables(c("harmonics", "hc60"))

### LAMBDA CREATION FUNCTION ###

lambdas <- function(dvec, astlon = c("task","cartwright"), ...) 
{
  astlon <- match.arg(astlon)
  dvec <- as.Date(dvec, ...)
  
  if(astlon == "task") {
    IY <- as.numeric(format(dvec,"%Y")) - 1900
    IL <- ifelse(IY > 0, trunc((IY-1)/4), trunc(IY/4))
    DL <- IL + as.POSIXlt(dvec)$yday
    coefs <- rbind(c(277.0247,129.38481,13.17639),
                 c(280.1895,-0.23872,0.98565),
                 c(334.3853,40.66249,0.11140),
                 c(259.1568,-19.32818,-0.05295),
                 c(281.2209,0.017192,0))
    out <- (coefs %*% rbind(1, IY, DL)) %% 360
  }
  if(astlon == "cartwright") {
    TT <- as.numeric(dvec - as.Date("1899-12-31"))/36525
    coefs <- 360 * rbind(c(0.751206,1336.855231,-0.000003),
                   c(0.776935,100.002136,0.000001),
                   c(0.928693,11.302872,-0.000029),
                   c(0.719954,-5.372617,0.000006),
                   c(0.781169,0.004775,0.000001))
    out <- (coefs %*% rbind(1, TT, TT*TT)) %% 360
  }
  out
}

### NODAL ADJUSTMENT FUNCTION ###

nodal_adj <- function(lambp, lambN, lambph, 
  indegree = TRUE, outdegree = TRUE)
{
  if(missing(lambp) || !is.numeric(lambp)) 
    stop("'lambp' must be a numeric vector")
  if(missing(lambN) || !is.numeric(lambN)) 
    stop("'lambN' must be a numeric vector")
  if(missing(lambph) || !is.numeric(lambph)) 
    stop("'lambph' must be a numeric vector")
  
  nt <- length(lambp)
  if(length(lambN) != nt)
    stop("'lambp' and 'lambN' must have the same length")
  if(length(lambph) != nt)
    stop("'lambp' and 'lambph' must have the same length")

  hc <- harmonics$name
  if(indegree) {
    lambN <- lambN * pi/180
    lambp <- lambp * pi/180
    lambph <- lambph * pi/180
  }
  
  fn <- un <- matrix(NA, nt, length(hc))
  colnames(fn) <- colnames(un) <- hc
  c1 <- cos(lambN); s1 <- sin(lambN)
  c2 <- cos(2*lambN); s2 <- sin(2*lambN)
  c3 <- cos(3*lambN); s3 <- sin(3*lambN)
  c1p <- cos(lambp); s1p <- sin(lambp)
  c2p <- cos(2*lambp); s2p <- sin(2*lambp)
  
  # ones and zeros as given by nodal z or nodal f
  
  ones <- c("aSa","Sa","Ssa","Sta","pi1","P1","xS1","S1","aS1","psi1",
    "NA2","aNA2","NB2","NA2*","MA2","MB2","MA2*","T2","S2","R2")
  fn[,ones] <- 1
  un[,ones] <- 0
  
  # ones and zeros from long-term constituents with nodal x
  
  ones <- c("MSm","Mnum","SM","KOo","MKo","Snu","SN","MStm","2SMN")
  fn[,ones] <- 1
  un[,ones] <- 0
  
  # fundamental constituents given by nodal y or Y
  
  fn[,"Mm"] <- 1.0000 - 0.1311*c1 + 0.0538*c2p + 0.0205*cos(2*lambp-lambN)
  un[,"Mm"] <- 0
  fn[,"Mf"] <- 1.084 + 0.415*c1 + 0.039*c2
  un[,"Mf"] <- -23.7*s1 + 2.7*s2 - 0.4*s3
  fn[,"O1"] <- 1.0176 + 0.1871*c1 - 0.0147*c2 + 0.0014*c3
  un[,"O1"] <- 10.80*s1 - 1.34*s2 + 0.19*s3
  fn[,c("K1","xK1")] <- 1.0060 + 0.1150*c1 - 0.0088*c2 + 0.0006*c3
  un[,c("K1","xK1")] <- -8.86*s1 + 0.68*s2 - 0.07*s3
  fn[,"J1"] <- 1.0129 + 0.1676*c1 - 0.0170*c2 + 0.0016*c3
  un[,"J1"] <- -12.94*s1 + 1.34*s2 - 0.19*s3
  fn[,"M2"] <- 1.0007 - 0.0373*c1 + 0.0002*c2
  un[,"M2"] <- -2.14*s1
  fn[,"K2"] <- 1.0246 + 0.2863*c1 + 0.0083*c2 - 0.0015*c3
  un[,"K2"] <- -17.74*s1 + 0.68*s2 - 0.04*s3
  
  # polar constituents given by nodal y or Y
  
  ss <- 2.783*s2p + 0.558*sin(2*lambp - lambN) + 0.184*s1
  cc <- 1 + 2.783*c2p + 0.558*cos(2*lambp - lambN) + 0.184*c1
  fn[,c("M1B","xM1B")] <- sqrt(cc^2 + ss^2)
  un[,c("M1B","xM1B")] <- atan2(ss, cc) * 180/pi
  
  ss <- s1p + 0.2*sin(lambp - lambN)
  cc <- 2.0*c1p + 0.4*cos(lambp - lambN)
  fn[,c("M1","xM1","aM1","M1C")] <- sqrt(cc^2 + ss^2)
  un[,c("M1","xM1","aM1","M1C")] <- atan2(ss, cc) * 180/pi
  
  ss <- -0.3593*s2p - 0.066*sin(2*lambp - lambN) - 0.2*s1
  cc <- 1 + 0.3593*c2p + 0.066*cos(2*lambp - lambN) + 0.2*c1
  fn[,"M1A"] <- sqrt(cc^2 + ss^2)
  un[,"M1A"] <- atan2(ss, cc) * 180/pi
  
  ss <- 0.147*sin(2*lambN - 2*lambp)
  cc <- 1 + 0.147*cos(2*lambN - 2*lambp)
  fn[,"gam2"] <- sqrt(cc^2 + ss^2)
  un[,"gam2"] <- atan2(ss, cc) * 180/pi
  
  ss <- -0.0446*sin(lambp - lambph)
  cc <- 1 - 0.0446*cos(lambp - lambph)
  fn[,"alp2"] <- sqrt(cc^2 + ss^2)
  un[,"alp2"] <- atan2(ss, cc) * 180/pi
  
  ss <- 0.477*s1
  cc <- 1 - 0.477*c1
  fn[,"del2"] <- sqrt(cc^2 + ss^2)
  un[,"del2"] <- atan2(ss, cc) * 180/pi
  
  ss <- -0.439*s1
  cc <- 1 + 0.439*c1
  fn[,c("xi2","eta2")] <- sqrt(cc^2 + ss^2)
  un[,c("xi2","eta2")] <- atan2(ss, cc) * 180/pi
  
  ss <- -0.2505*s2p - 0.1102*sin(2*lambp - lambN) - 0.0156*sin(2*lambp - 2*lambN) - 0.037*s1
  cc <- 1.0 - 0.2505*c2p - 0.1102*cos(2*lambp - lambN) - 0.0156*cos(2*lambp - 2*lambN) - 0.037*c1 
  fn[,"L2"] <- sqrt(cc^2 + ss^2)
  un[,"L2"] <- atan2(ss, cc) * 180/pi
  
  # single higher frequency terms 
  
  fn[,"MA4"] <- fn[,"M2"] 
  un[,"MA4"] <- un[,"M2"]
  fn[,c("M3","MB5")] <- fn[,"M2"]^1.5
  un[,c("M3","MB5")] <- 1.5*un[,"M2"]
  fn[,c("M4","N4","MA6")] <- fn[,"M2"]^2
  un[,c("M4","N4","MA6")] <- 2*un[,"M2"]
  fn[,c("M6","N6","MA8")] <- fn[,"M2"]^3
  un[,c("M6","N6","MA8")] <- 3*un[,"M2"]
  fn[,"M8"] <- fn[,"M2"]^4
  un[,"M8"] <- 4*un[,"M2"]
  fn[,c("M10","MA12")] <- fn[,"M2"]^5
  un[,c("M10","MA12")] <- 5*un[,"M2"]
  fn[,"M12"] <- fn[,"M2"]^6
  un[,"M12"] <- 6*un[,"M2"]
  fn[,c("M5","xM5","aM5")] <- (fn[,"M4"] + fn[,"M6"])/2
  un[,c("M5","xM5","aM5")] <- (un[,"M4"] + un[,"M6"])/2
  fn[,c("M7","aM7","MA9")] <- (fn[,"M6"] + fn[,"M8"])/2
  un[,c("M7","aM7","MA9")] <- (un[,"M6"] + un[,"M8"])/2
  
  fn[,c("S3","S4","S6","S8")] <- 1
  un[,c("S3","S4","S6","S8")] <- 0
  
  fn[,c("K3","xK3")] <- fn[,"K2"]^1.5
  un[,c("K3","xK3")] <- 1.5*un[,"K2"]
  fn[,"K4"] <- fn[,"K2"]^2
  un[,"K4"] <- 2*un[,"K2"]
  
  # specified terms based on simple (nodal abcjkmo)
  
  fn[,"Mfm"] <- fn[,"Mm"]
  un[,"Mfm"] <- un[,"Mm"]
  fn[,c("MSf","MSo","MSqm")] <- fn[,"M2"]
  un[,c("MSf","MSo","MSqm")] <- -un[,"M2"]
  fn[,"2SM"] <- fn[,"M2"]^2
  un[,"2SM"] <- -2*un[,"M2"]
  fn[,c("chi1","phi1","the1")] <- fn[,"J1"]
  un[,c("chi1","phi1","the1")] <- un[,"J1"]
  fn[,"tau1"] <- fn[,"K1"]
  un[,"tau1"] <- un[,"K1"]
  fn[,c("Mqm","eps2","2N2","mu2","N2","nu2","lam2")] <- fn[,"M2"]
  un[,c("Mqm","eps2","2N2","mu2","N2","nu2","lam2")] <- un[,"M2"]
  fn[,c("2Q1","sig1","Q1","rho1")] <- fn[,"O1"]
  un[,c("2Q1","sig1","Q1","rho1")] <- un[,"O1"]
  fn[,"O2"] <- fn[,"O1"]^2
  un[,"O2"] <- 2*un[,"O1"]
  
  # component terms
  fn2 <- fn[,c("M2","S2","T2","N2","K2","L2","R2","nu2","lam2","Q1","K1","J1","O1","P1","S1"),drop=FALSE]
  un2 <- un[,c("M2","S2","T2","N2","K2","L2","R2","nu2","lam2","Q1","K1","J1","O1","P1","S1"),drop=FALSE]
  colnames(fn2) <- colnames(un2) <- c("1M","1S","1T","1N","1K","1L","1R","1V","1D","1q","1k","1j","1o","1p","1s")
  fn2m <- fn2[,"1M"]; un2m <- un2[,"1M"]

  fn2 <- cbind(fn2, "2M" = 2*fn2m, "3M" = 3*fn2m, "4M" = 4*fn2m, 
               "5M" = 5*fn2m, "6M" = 6*fn2m, "7M" = 7*fn2m,
               "2S" = 2*fn2[,"1S"], "3S" = 3*fn2[,"1S"],
               "2N" = 2*fn2[,"1N"], "3N" = 3*fn2[,"1N"],
               "2k" = 2*fn2[,"1k"], "3k" = 3*fn2[,"1k"],
               "2K" = 2*fn2[,"1K"], "2p" = 2*fn2[,"1p"])
  un2 <- cbind(un2, "2M" = 2*un2m, "3M" = 3*un2m, "4M" = 4*un2m, 
               "5M" = 5*un2m, "6M" = 6*un2m, "7M" = 7*un2m,
               "2S" = 2*un2[,"1S"], "3S" = 3*un2[,"1S"],
               "2N" = 2*un2[,"1N"], "3N" = 3*un2[,"1N"],
               "2k" = 2*un2[,"1k"], "3k" = 3*un2[,"1k"],
               "2K" = 2*un2[,"1K"], "2p" = 2*un2[,"1p"])
  
  nms <- sub("x","",grep("\\.", hc, value=TRUE))
  pos <- substr(nms, 1, regexpr("\\.", nms)-1)
  neg <- substr(nms, regexpr("\\.", nms)+1, regexpr("\\d{1,2}$", nms)-1)
  pos <- strsplit(gsub("([[:alnum:]]{2})", "\\1 ", pos), " ")
  neg <- strsplit(gsub("([[:alnum:]]{2})", "\\1 ", neg), " ")
  
  fn2 <- log(fn2)
  tmpfn <- function(x) rowSums(fn2[,x,drop=FALSE])
  tmpun <- function(x) rowSums(un2[,x,drop=FALSE])
  fn[,grep("\\.",hc)] <- exp(sapply(pos, tmpfn) + sapply(neg, tmpfn))
  un[,grep("\\.",hc)] <- sapply(pos, tmpun) - sapply(neg, tmpun)
  
  # specified terms based on component (nodal dpq) 
  
  fn[,c("OO1","ups1")] <- fn[,"1K.1q1"]
  un[,c("OO1","ups1")] <- un[,"1K.1q1"]
  fn[,"L2A"] <- fn[,"2M.1N2"]
  un[,"L2A"] <- un[,"2M.1N2"]
  fn[,"L2B"] <- fn[,"1N1K.1M2"]
  un[,"L2B"] <- un[,"1N1K.1M2"]
  
  # return values
  
  if(!outdegree) un <- un * pi/180
  list(un = t(un), fn = t(fn))
}

