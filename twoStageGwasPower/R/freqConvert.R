freqConvert <-
function(pD, pG, grr, inheritance="dominant") {
  #Take the prevalence, allele frequency, and genotypic relative risk
  # and compute the allele frequencies in cases and controls, and the
  # penetrances of the genotypes
  
  #pD  =  pr(D), prob of disease in population (prevalence)
  #pG  =  pr(G) in population
  #grr =  GRR for dominant disease, pr(D|GG)/pr(D|gg) = pr(D|Gg)/pr(D|gg)

  pGG <- pG^2
  pGg <- 2*pG*(1-pG)
  pgg <- (1-pG)^2

  # penetances (match CaTS!)
  if (inheritance=="dominant") {
    pD.given.gg <- pD/(grr*pGG + grr*pGg + pgg)
    pD.given.Gg <- pD.given.gg*grr
    pD.given.GG <- pD.given.gg*grr
    }
  else if (inheritance == "recessive") {
    pD.given.gg <- pD/(grr*pGG + pGg + pgg)
    pD.given.Gg <- pD.given.gg
    pD.given.GG <- pD.given.gg*grr
    }
  else if (inheritance == "multiplicative") {
    pD.given.gg <- pD/(grr^2*pGG + grr*pGg + pgg)
    pD.given.Gg <- pD.given.gg*grr
    pD.given.GG <- pD.given.gg*grr^2
    }
  else if (inheritance == "additive") {
    pD.given.gg <- pD/((2*grr-1)*pGG + grr*pGg + pgg)
    pD.given.Gg <- pD.given.gg*grr
    pD.given.GG <- pD.given.gg*(2*grr - 1)
    }
  else {
    stop("Invalid inheritance type")
     }
  pGG.given.D <- pD.given.GG*pGG/pD
  pGg.given.D <- pD.given.Gg*pGg/pD
  pgg.given.D <- pD.given.gg*pgg/pD

  pGG.given.notD <- (1 - pD.given.GG)*pGG/(1 - pD)
  pGg.given.notD <- (1 - pD.given.Gg)*pGg/(1 - pD)
  pgg.given.notD <- (1 - pD.given.gg)*pgg/(1 - pD)

  # allele frequency (match CaTS!)
  pG.given.D <- pGG.given.D + 0.5*pGg.given.D
  p1 <- pG.given.D
  pG.given.notD <- pGG.given.notD + 0.5*pGg.given.notD
  p0 <- pG.given.notD
  result <- list(pD.given.GG=pD.given.GG, pD.given.Gg=pD.given.Gg,
                 pD.given.gg=pD.given.gg, p1=p1, p0=p0)
  result
  }
