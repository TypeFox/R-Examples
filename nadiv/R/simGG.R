simGG <- function(K, pairs, noff, g,
	nimm = 2, nimmG = seq(2, g-1, 1),
	VAf = 1, VAi = 1, VRf = 1, VRi = 1,
	mup = 20, muf = 0, mui = 0, murf = 0, muri = 0,
	d_bvf = 0, d_bvi = 0, d_rf = 0, d_ri = 0){
  if(pairs*2 > K) stop("pairs must be less than half of K")
  if(nimmG[1] == 1) stop("immigrants cannot arrive in the first generation")
  N <- pairs*noff
  da <- array(NA, dim = c(K, 10, g))
  dimnames(da) <- list(NULL,
     c("id", "dam", "sire", "parAvgU", "mendel", "u", "r", "p", "is", "gen"), seq(g))
  da[, "id", ] <- seq(K*g)
  # Assume last nimm rows in each generation are the immigrants: 
  ## 1=immigrant & 0=NOT immigrant
  da[, "is", ] <- 0
  da[(K-nimm+1):K, "is", nimmG] <- 1
  da[, "gen", ] <- rep(seq(g), each = K)
  # create standard normals for generation 1
  da[, "u", 1] <- rnorm(K, muf, sqrt(VAf))
  da[, "r", 1] <- rnorm(K, murf, sqrt(VRf))
 
  # Mating
  for(i in 2:g){
    if(i %in% nimmG) Knimm <- K-nimm else Knimm <- K
    # Parents ranked according to how close u is to the next generation's mean u
    ## Define function to do this
    prFun <- function(x, theta = 0){
      exp((-1 / (2*sd(x))) * (x - theta)^2)
    }
    ### End function definition
    if(d_bvf == 0) prin <- NULL else prin <- prFun(da[, "u", i-1],
	muf + d_bvf*sqrt(VAf)*(i-2)) # no selection in first generation hence i-2
    # Sampling WITH replacement to assign parents
    parPool <- sort(sample(x = da[, "id", i-1], size = pairs*2, replace = TRUE,
	prob = prin))
      pphalf <- length(parPool)/2
      pphalfId <- parPool[pphalf]
      while(pphalfId == parPool[pphalf + 1]){
        parPool[] <- parPool[c(which(parPool != pphalfId), which(parPool == pphalfId))]        
        pphalfId <- parPool[pphalf]
      }
      parPool[1:pphalf] <- parPool[1:pphalf][sample(length(parPool)/2)]
    iOff <- matrix(rep(parPool, each = noff), ncol = 2)
    da[1:Knimm, c("dam", "sire"), i] <- iOff[sort(sample(seq(nrow(iOff)),
	size = Knimm, replace = FALSE)), ]
    # Average of parent breeding values
    da[1:Knimm, "parAvgU", i] <- rowMeans(matrix(da[match(da[1:Knimm, c("dam", "sire"), i],
	da[, "id", i-1]), "u", i-1], ncol = 2, byrow = FALSE))
    # Assign Mendelian sampling variation
    # Within-family additive genetic variance
    ## p. 447, second eqn. in Verrier, Colleau, & Foulley. 1993. TAG
    if(i == 2){
      da[1:Knimm, "mendel", i] <- rnorm(Knimm, 0, sqrt(0.5 * VAf))
    } else{
      # Average of parent inbreeding coefficients
        tmpPed <- prunePed(data.frame(apply(da[, 1:3, 1:i-1], MARGIN = 2,
		FUN = function(x){x})), as.character(unique(c(da[1:Knimm, c("dam", "sire"), i]))))
        tmpParF <- makeAinv(tmpPed)$f
        da[1:Knimm, "mendel", i ] <- rnorm(Knimm, 0,
		sqrt(0.5 * VAf * (1 - rowMeans(matrix(tmpParF[match(as.character(c(da[1:Knimm,
		c("dam", "sire"), i])), tmpPed[, 1])], ncol = 2, byrow = FALSE)))))
      }
    # Total additive genetic effects
    da[1:Knimm, "u", i] <- rowSums(da[1:Knimm, c("parAvgU", "mendel"), i])
    # Residual deviations
    da[1:Knimm, "r", i] <- rnorm(Knimm, murf + d_rf*sqrt(VRf)*(i-1), sqrt(VRf))
    ############
    # Immigrants: assumed outbred
    if(i %in% nimmG){
      # no trend in first generation
      da[(Knimm+1):K, "u", i] <- rnorm(nimm, mui + d_bvi*sqrt(VAi)*(i-2), sqrt(VAi))
      # no trend in first generation
      da[(Knimm+1):K, "r", i] <- rnorm(nimm, muri + d_ri*sqrt(VRi)*(i-2), sqrt(VRi))
    }
  } 
  # create a data.frame out of the array
  df <- data.frame(apply(da, MARGIN = 2, FUN = function(x){x}))
  # Phenotypes
  df$p <- mup + rowSums(df[, c("u", "r")])
 df
} 
