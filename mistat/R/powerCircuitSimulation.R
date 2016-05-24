powerCircuitSimulation <-function(
    rsA = 8200, 
           rsB = 220000,
           rsC = 1000,
           rsD = 33000,
           rsE = 56000,
           rsF = 5600,
           rsG = 3300,
           rsH = 58.5,
           rsI = 1000,
           rsJ = 120,
           trK = 130,
           trL = 100,
           trM = 130,
           tlA = 5,
           tlB = 10,
           tlC = 10,
           tlD = 5,
           tlE = 5,
           tlF = 5,
           tlG = 10,
           tlH = 5,
           tlI = 5,
           tlJ = 5,
           tlK = 5,
           tlL = 10,
           tlM = 5,
           each = 50, seed = NA){
    
    
    Lm <- length(rsA)
    if(any(c(length(rsB) != Lm, length(trM) != Lm, length(rsC) != Lm, 
             length(rsD) != Lm, length(rsE) != Lm, length(rsF) != Lm,
             length(rsG) != Lm, length(rsH) != Lm, length(rsI) != Lm,
             length(rsJ) != Lm, length(trK) != Lm, length(trL) != Lm,
             length(tlA) != Lm, length(tlB) != Lm, length(tlC) != Lm,
             length(tlD) != Lm, length(tlE) != Lm, length(tlF) != Lm,
             length(tlG) != Lm, length(tlH) != Lm, length(tlI) != Lm,
             length(tlJ) != Lm, length(tlK) != Lm, length(tlL) != Lm, length(tlM) != Lm)))
      stop("input values should be of same length")
    rm(Lm)
    
    Va <- 0
    if(any(c(rsA <= Va, rsB <= Va, trM <= Va, rsC <= Va, 
             rsD <= Va, rsE <= Va, rsF <= Va,
             rsG <= Va, rsH <= Va, rsI <= Va,
             rsJ <= Va, trK <= Va, trL <= Va, trM <= Va,
             tlA <= Va, tlB <= Va, tlC <= Va,
             tlD <= Va, tlE <= Va, tlF <= Va,
             tlG <= Va, tlH <= Va, tlI <= Va,
             tlJ <= Va, tlK <= Va, tlL <= Va, tlM <= Va)))
      stop("input values should be positive")
    rm(Va)
    
    Df <- data.frame(
      rsA =  rep(rsA,  each = each),
      rsB =  rep(rsB,  each = each),
      rsC =  rep(rsC,  each = each),
      rsD =  rep(rsD,  each = each),
      rsE =  rep(rsE,  each = each),
      rsF =  rep(rsF,  each = each),
      rsG =  rep(rsG,  each = each),
      rsH =  rep(rsH,  each = each),
      rsI =  rep(rsI,  each = each),
      rsJ =  rep(rsJ,  each = each),
      rsK =  rep(trK,  each = each),
      rsL =  rep(trL,  each = each),
      rsM =  rep(trM,  each = each),
      tlA =  rep(tlA,  each = each),
      tlB =  rep(tlB,  each = each),
      tlC =  rep(tlC,  each = each),
      tlD =  rep(tlD,  each = each),
      tlE =  rep(tlE,  each = each),
      tlF =  rep(tlF,  each = each),
      tlG =  rep(tlG,  each = each),
      tlH =  rep(tlH,  each = each),
      tlI =  rep(tlI,  each = each),
      tlJ =  rep(tlJ,  each = each),
      tlK =  rep(tlK,  each = each),
      tlL =  rep(tlL,  each = each),
      tlM =  rep(tlM,  each = each))
    
    n <- nrow(Df)
    
    X <- matrix(NA, 13*3, n)
    
    if(!is.na(seed))
      set.seed(seed)
    
    for(j in 1:13){
      X[j,] <- rnorm(n, mean = 0, sd = 1) 
    }
    
    X[1:13, ] <- X[1:13, ] * 
      # sd
      (2 * t(as.matrix(Df[,14:26])) * t(as.matrix(Df[,1:13])) / 600)+ 
      # mean
      t(as.matrix(Df[,1:13]))
    
    X[14:26,] <- t(as.matrix(Df[,1:13]))
    
    X[27:39,] <- t(as.matrix(Df[,14:26]))
    
    res <- as.data.frame(t(apply(X, 2, .vCircuit)))
    class(res) <- c(class(res), "mistatSimulation", "powerCircuitSimulation")
    return(res)
  }

.vCircuit <- function(x)
{
  a <- x[2] / (x[1] + x[2])
  b <- (x[1] * x[2] / (x[1] + x[2]) + x[3]) / (x[12] * x[13]) + x[9]
  c <- x[5] + 0.5 * x[7]
  d <- x[1] * x[2] * x[11] / (x[1] + x[2])
  e <- x[6] + 0.5 * x[7]
  f <- (c + e) * (1 + x[11]) * x[8] + c * e
  g <- x[8] + 0.6
  h <- 1.2
  
  aden <- 1 + d * e / f + b * (1 / x[10] + 0.006 * (1 + 13.67 / x[10])) + a * 0.1367 * 0.6
  anum <- (a + b / x[10]) * (138 - 1.33) + d * (c + e) * g / f - h
  
  res <- anum / aden
  res <- c(x[14:39], res)
  names(res) <- c("rsA", "rsB", "rsC", "rsD", "rsE", "rsF", 
                  "rsG", "rsH", "rsI", "rsJ", "trK", "trL", 
                  "trM", "tlA", "tlB", "tlC", "tlD", "tlE", 
                  "tlF", "tlG", "tlH", "tlI", "tlJ", "tlK", 
                  "tlL", "tlM", "volts") 
  return(res)
}
