# CHNOSZ/util.protein.R
# calculate formulas and summarize properties of proteins
# MP90.cp: additive heat capacity from groups of Makhatadze and Privalov, 1990
# group.formulas: chemical makeup of the amino acid residues
# protein.formula: chemical makeup of the indicated proteins
# protein.length: lengths of the indicated proteins
# protein.info: summarize properties of proteins
# protein.basis: coefficients of basis species in formation reactions of [ionized] proteins [residues]
# protein.equil: step-by-step example of protein equilibrium calculation

MP90.cp <- function(protein, T) {
  # T (temperature, degrees C), protein (name of protein)
  # returns heat capacity of protein (kj/mol)
  # using algorithm of makhatadze and privalov, 1990.
  TMP <- c(5,25,50,75,100,125)
  A.cp <- splinefun(TMP,c(175.7,166.7,156.2,144.7,134.6,124.1))
  C.cp <- splinefun(TMP,c(225.4,237.6,250.8,260.7,268.2,276.1))
  D.cp <- splinefun(TMP,c( 72.8, 89.0,106.2,124.5,140.7,154.3))
  E.cp <- splinefun(TMP,c(168.3,179.0,192.0,203.7,211.4,217.8))
  F.cp <- splinefun(TMP,c(395.7,383.0,370.3,358.4,348.3,339.6))
  G.cp <- splinefun(TMP,c( 82.3, 78.0, 71.7, 66.4, 59.7, 53.9))
  H.cp <- splinefun(TMP,c(205.7,179.6,177.2,179.6,187.1,196.8))
  I.cp <- splinefun(TMP,c(406.8,402.3,397.1,390.8,386.0,380.8))
  K.cp <- splinefun(TMP,c(328.8,332.5,334.0,337.5,339.4,343.6))
  L.cp <- splinefun(TMP,c(385.9,381.7,377.8,372.9,369.4,365.5))
  M.cp <- splinefun(TMP,c(197.1,175.9,158.1,150.3,148.1,143.9))
  N.cp <- splinefun(TMP,c( 72.9, 88.8,109.8,125.2,140.5,154.2))
  P.cp <- splinefun(TMP,c(214.6,177.7,152.3,142.8,135.6,130.1))
  Q.cp <- splinefun(TMP,c(168.0,180.2,193.4,203.3,210.8,218.7))
  R.cp <- splinefun(TMP,c(204.6,273.4,305.8,315.1,318.7,318.5))
  S.cp <- splinefun(TMP,c( 75.6, 81.2, 85.7, 91.4, 97.3,102.1))
  T.cp <- splinefun(TMP,c(194.2,184.5,182.2,186.5,199.0,216.2))
  V.cp <- splinefun(TMP,c(324.6,314.4,305.0,294.7,285.7,269.6))
  W.cp <- splinefun(TMP,c(471.2,458.5,445.8,433.9,423.8,415.1))
  Y.cp <- splinefun(TMP,c(310.6,301.7,295.2,294.5,300.1,304.0))
  AA.cp <- splinefun(TMP,c(-158.3,-90.4,-21.5,-32.3,-92.4,-150.0))
  UPBB.cp <- splinefun(TMP,c(3.7,15.2,26.2,29.8,33.7,33.7))
  cnew <- numeric()
  for(i in 1:length(T)) {
    Ti <- T[i]
    cp <- c(A.cp(Ti),C.cp(Ti),D.cp(Ti),E.cp(Ti),F.cp(Ti),
            G.cp(Ti),H.cp(Ti),I.cp(Ti),K.cp(Ti),L.cp(Ti),
            M.cp(Ti),N.cp(Ti),P.cp(Ti),Q.cp(Ti),R.cp(Ti),
            S.cp(Ti),T.cp(Ti),V.cp(Ti),W.cp(Ti),Y.cp(Ti))
    # get the protein composition
    tt <- ip2aa(protein)[,6:25]
    cnew <- c(cnew, sum(cp * as.numeric(tt)) + sum(as.numeric(tt)) * UPBB.cp(Ti))
  }
  return(cnew)
}


group.formulas <- function() {
  # return a matrix with chemical formulas of residues
  # names of the sidechain groups
  groups <- paste("[", aminoacids(3), "]", sep="")
  # the indices of H2O, sidechain groups, and [UPBB]
  ig <- suppressMessages(info(c("H2O", groups, "[UPBB]")))
  # their formulas
  A <- i2A(ig)
  # add [UPBB] to the sidechain groups to get residues
  out <- A[1:21,]
  out[2:21,] <- t(t(A) + A[22,])[2:21,]
  # make "H2O" not "water"
  rownames(out)[1] <- "H2O"
  return(out)
}

protein.formula <- function(protein, organism=NULL, residue=FALSE) {
  # return a matrix with chemical formulas of proteins
  aa <- ip2aa(protein, organism)
  rf <- group.formulas()
  out <- as.matrix(aa[, 5:25]) %*% as.matrix(rf)
  if(residue) out <- out / rowSums(aa[, 6:25])
  row.names(out) <- make.unique(paste(aa$protein, aa$organism, sep="_"))
  return(out)
}

protein.length <- function(protein, organism=NULL) {
  # calculate the length(s) of proteins
  aa <- ip2aa(protein, organism)
  # use rowSums on the columns containing amino acid counts
  pl <- as.numeric(rowSums(aa[, 6:25]))
  return(pl)
}

protein.info <- function(protein, T=25, residue=FALSE, round.it=FALSE) {
  # make a table of selected properties for proteins
  # listed in protein
  aa <- ip2aa(protein)
  pname <- paste(aa$protein, aa$organism, sep="_")
  length <- protein.length(aa)
  pf <- protein.formula(aa)
  G <- unlist(subcrt(pname, T=T, property="G")$out)
  Z <- rep(NA, length(pname))
  G.Z <- rep(NA, length(pname))
  ZC <- ZC(pf)
  # run ionization calculations if we have H+
  thermo <- get("thermo")
  if(!is.null(thermo$basis)) {
    iHplus <- match("H+", rownames(thermo$basis))
    if(!is.na(iHplus)) {
      pH <- -thermo$basis$logact[iHplus]
      Z <- ionize.aa(aa, T=T, pH=pH)[1, ]
      G.ionization <- ionize.aa(aa, T=T, pH=pH, property="G")[1, ]
      G.Z <- G + G.ionization
      # add charge to the chemical formulas
      pf <- cbind(pf, Z)
      iH <- match("H", colnames(pf))
      pf[, iH] <- pf[, iH] + Z
    }
  }
  # take care of residue conversion
  if(residue) {
    pf <- pf / length
    G <- G / length
    Z <- Z / length
    G.Z <- G.Z / length
    length <- length / length
    # round the coefficients in the formulas
    if(round.it) pf <- round(pf, 3)
  }
  # convert each protein formula to a single line
  formula <- as.chemical.formula(round(pf, 3))
  if(round.it) {
    length <- round(length, 1)
    G <- round(G, 3)
    Z <- round(Z, 3)
    G.Z <- round(G.Z, 3)
    ZC <- round(ZC, 3)
  }
  out <- data.frame(protein=pname, length=length, formula=formula, G=G, Z=Z, G.Z=G.Z, ZC=ZC, stringsAsFactors=FALSE)
  rownames(out) <- NULL
  return(out)
}

protein.basis <- function(protein, T=25, normalize=FALSE) {
  # 20090902 calculate the coefficients of basis species in reactions
  # to form proteins (possibly per normalized by length) listed in protein
  # 20120528 renamed protein.basis from residue.info ...
  # what are the elemental compositions of the proteins
  aa <- ip2aa(protein)
  pf <- protein.formula(aa)
  # what are the coefficients of the basis species in the formation reactions
  sb <- species.basis(pf)
  # calculate ionization states if H+ is a basis species
  thermo <- get("thermo")
  iHplus <- match("H+", rownames(thermo$basis))
  if(!is.na(iHplus)) {
    pH <- -thermo$basis$logact[iHplus]
    Z <- ionize.aa(aa, T=T, pH=pH)[1, ]
    sb[, iHplus] <- sb[, iHplus] + Z
  }
  # compute per length-normalized coefficients if requested
  if(normalize) {
    # get lengths of proteins
    plen <- protein.length(aa)
    sb <- sb/plen
  }
  # return the result
  return(sb)
}

protein.equil <- function(protein, T=25, loga.protein=0, digits=4) {
  # show the individual steps in calculating metastable equilibrium among proteins
  msgout("protein.equil: temperature from argument is ", T, " degrees C\n")
  TK <- convert(T, "K")
  # get the amino acid compositions of the proteins
  aa <- ip2aa(protein)
  # get some general information about the proteins
  pname <- paste(aa$protein, aa$organism, sep="_")
  plength <- protein.length(aa)
  # use thermo$basis to decide whether to ionize the proteins
  thermo <- get("thermo")
  ionize.it <- FALSE
  iword <- "nonionized"
  bmat <- basis.elements()
  if("H+" %in% rownames(bmat)) {
    ionize.it <- TRUE
    iword <- "ionized"
    pH <- -thermo$basis$logact[match("H+", rownames(bmat))]
    msgout("protein.equil: pH from thermo$basis is ", pH, "\n")
  }
  # tell the user whose [Met] is in thermo$obigt
  info.Met <- info(info('[Met]', "aq"))
  msgout("protein.equil: [Met] is from reference ", info.Met$ref1, "\n")
  ## first set of output: show results of calculations for a single protein
  msgout("protein.equil [1]: first protein is ", pname[1], " with length ", plength[1], "\n")
  # standard Gibbs energies of basis species
  G0basis <- unlist(suppressMessages(subcrt(thermo$basis$ispecies, T=T, property="G")$out))
  # coefficients of basis species in formation reactions of proteins
  protbasis <- suppressMessages(protein.basis(aa, T=T))
  # sum of standard Gibbs energies of basis species in each reaction
  G0basissum <- colSums(t(protbasis) * G0basis)
  # standard Gibbs energies of nonionized proteins
  G0prot <- unlist(suppressMessages(subcrt(pname, T=T, property="G")$out))
  # standard Gibbs energy of formation reaction of nonionized protein, cal/mol
  G0protform <- G0prot - G0basissum
  msgout("protein.equil [1]: reaction to form nonionized protein from basis species has G0(cal/mol) of ", signif(G0protform[1], digits), "\n")
  if(ionize.it) {
    # standard Gibbs energy of ionization of protein, cal/mol
    G0ionization <- suppressMessages(ionize.aa(aa, property="G", T=T, pH=pH))[1, ]
    msgout("protein.equil [1]: ionization reaction of protein has G0(cal/mol) of ", signif(G0ionization[1], digits), "\n")
    # standard Gibbs energy of formation reaction of ionized protein, cal/mol
    G0protform <- G0protform + G0ionization
  }
  # standard Gibbs energy of formation reaction of non/ionized residue equivalents, dimensionless
  G0res.RT <- G0protform/thermo$opt$R/TK/plength
  msgout("protein.equil [1]: per residue, reaction to form ", iword, " protein from basis species has G0/RT of ", signif(G0res.RT[1], digits), "\n")
  # coefficients of basis species in formation reactions of residues
  resbasis <- suppressMessages(protein.basis(aa, T=T, normalize=TRUE))
  # logQstar and Astar/RT
  logQstar <- colSums(t(resbasis) * - thermo$basis$logact)
  msgout("protein.equil [1]: per residue, logQstar is ", signif(logQstar[1], digits), "\n")
  Astar.RT <- -G0res.RT - log(10)*logQstar
  msgout("protein.equil [1]: per residue, Astar/RT = -G0/RT - 2.303logQstar is ", signif(Astar.RT[1], digits), "\n")
  if(!is.numeric(protein)) msgout("protein.equil [1]: not comparing calculations with affinity() because 'protein' is not numeric\n")
  else {
    # for **Astar** we have to set the activities of the proteins to zero, not loga.protein!
    a <- suppressMessages(affinity(iprotein=protein, T=T, loga.protein=0))
    aAstar.RT <- log(10) * as.numeric(a$values) / plength
    msgout("check it!       per residue, Astar/RT calculated using affinity() is ", signif(aAstar.RT[1], digits), "\n")
    if(!isTRUE(all.equal(Astar.RT, aAstar.RT, check.attributes=FALSE)))
      stop("Bug alert! The same value for Astar/RT cannot be calculated manually as by using affinity()")
  }
  if(length(pname)==1) msgout("protein.equil [all]: all done... give me more than one protein for equilibrium calculations\n")
  else {
    ## next set of output: equilibrium calculations
    msgout("protein.equil [all]: lengths of all proteins are ", paste(plength, collapse=" "), "\n")
    msgout("protein.equil [all]: Astar/RT of all residue equivalents are ", paste(signif(Astar.RT, digits), collapse=" "), "\n")
    expAstar.RT <- exp(Astar.RT)
    sumexpAstar.RT <- sum(expAstar.RT)
    msgout("protein.equil [all]: sum of exp(Astar/RT) of all residue equivalents is ", signif(sumexpAstar.RT, digits), "\n")
    # boltzmann distribution
    alpha <- expAstar.RT / sumexpAstar.RT    
    msgout("protein.equil [all]: equilibrium degrees of formation (alphas) of residue equivalents are ", paste(signif(alpha, digits), collapse=" "), "\n")
    # check with equilibrate()
    if(is.numeric(protein)) {
      loga.equil.protein <- unlist(suppressMessages(equilibrate(a, normalize=TRUE))$loga.equil)
      # here we do have to convert from logarithms of activities of proteins to degrees of formation of residue equivalents
      a.equil.residue <- plength*10^loga.equil.protein
      ealpha <- a.equil.residue/sum(a.equil.residue)
      msgout("check it!     alphas of residue equivalents from equilibrate() are ", paste(signif(ealpha, digits), collapse=" "), "\n")
      if(!isTRUE(all.equal(alpha, ealpha, check.attributes=FALSE)))
        stop("Bug alert! The same value for alpha cannot be calculated manually as by using equilibrate()")
    }
    # total activity of residues
    loga.residue <- log10(sum(plength * 10^loga.protein))
    msgout("protein.equil [all]: for activity of proteins equal to 10^", signif(loga.protein, digits), ", total activity of residues is 10^", signif(loga.residue, digits), "\n")
    # equilibrium activities of residues
    loga.residue.equil <- log10(alpha*10^loga.residue)
    msgout("protein.equil [all]: log10 equilibrium activities of residue equivalents are ", paste(signif(loga.residue.equil, digits), collapse=" "), "\n")
    # equilibrium activities of proteins
    loga.protein.equil <- log10(10^loga.residue.equil/plength)
    msgout("protein.equil [all]: log10 equilibrium activities of proteins are ", paste(signif(loga.protein.equil, digits), collapse=" "), "\n")
    # check with equilibrate()
    if(is.numeric(protein)) {
      eloga.protein.equil <- unlist(suppressMessages(equilibrate(a, loga.balance=loga.residue, normalize=TRUE))$loga.equil)
      msgout("check it!    log10 eq'm activities of proteins from equilibrate() are ", paste(signif(eloga.protein.equil, digits), collapse=" "), "\n")
      if(!isTRUE(all.equal(loga.protein.equil, eloga.protein.equil, check.attributes=FALSE)))
        stop("Bug alert! The same value for log10 equilibrium activities of proteins cannot be calculated manually as by using equilibrate()")
    }
  }
}

