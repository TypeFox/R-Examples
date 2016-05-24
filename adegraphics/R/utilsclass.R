
#### revoir cette fonction ####
.util.ellipse <- function(mx, my, vx, vy, cxy, coeff) {
  if(!is.finite(mx) | !is.finite(my)) ## levels with no individuals
    return(NULL)
  lig <- 100
  epsi <- 1e-10
  x <- 0
  y <- 0
  if(vx < 0) 
    vx <- 0
  if(vy < 0) 
    vy <- 0
  if(vx == 0 && vy == 0) 
    return(NULL)
  delta <- (vx - vy) * (vx - vy) + 4 * cxy * cxy
  delta <- sqrt(delta)
  l1 <- (vx + vy + delta) / 2
  l2 <- vx + vy - l1
  if(l1 < 0) 
    l1 <- 0
  if(l2 < 0) 
    l2 <- 0
  l1 <- sqrt(l1)
  l2 <- sqrt(l2)
  test <- 0
  if(vx == 0) {
    a0 <- 0
    b0 <- 1
    test <- 1
  }
  if((vy == 0) && (test == 0)) {
    a0 <- 1
    b0 <- 0
    test <- 1
  }
  if(((abs(cxy)) < epsi) && (test == 0)) {
    a0 <- 1
    b0 <- 0
    test <- 1
  }
  if(test == 0) {
    a0 <- 1
    b0 <- (l1 * l1 - vx) / cxy
    norm <- sqrt(a0 * a0 + b0 * b0)
    a0 <- a0 / norm
    b0 <- b0 / norm
  }
  a1 <- 2 * pi / lig
  c11 <- coeff * a0 * l1
  c12 <- (-coeff) * b0 * l2
  c21 <- coeff * b0 * l1
  c22 <- coeff * a0 * l2
  angle <- 0
  for (i in 1:lig) {
    cosinus <- cos(angle)
    sinus <- sin(angle)
    x[i] <- mx + c11 * cosinus + c12 * sinus
    y[i] <- my + c21 * cosinus + c22 * sinus
    if(is.null(mx + c11 * cosinus + c12 * sinus) || is.null(y[i] <- my + c21 * cosinus + c22 * sinus))
      print("in util.ellipse x or y null")
    angle <- angle + a1
  }
  return(list(x = x, y = y, seg1 = c(mx + c11, my + c21, mx - c11, my - c21), seg2 = c(mx + c12, my + c22, mx - c12, my - c22)))
}

## Nouvelle version:
## principe:
## 1) calcul de distance entre les points appartenant a un groupe et le centroides du groupe
## 2) extraction du quantile correspondant a optchull (les % d les plus eloignes forment le polugfone
## x, y: points, mx, my: coordonnees des centroides, optchull: paramètre voulu pour lenvellope converxe, fac: facteur separeant les poitns
.util.chull <- function(x, y, mx, my, fac, chullSize) {
  ## pour chaque groupe calcul des distances
  chulls <- list()
  for(i in 1:nlevels(fac)) { ## attention fac est passe en facteur!
    index <- which(fac == levels(fac)[i])
    if(length(index) > 0) {
      x1 <- x[index]
      y1 <- y[index]
      dd <- sqrt((x1 - mx[i])^2 + (y1 - my[i])^2) ## distances chaque points a la moyenne
      tmp_quant <- list()
      for(quant in chullSize) { ## pour chaque envelope demandee
        selected <- which(dd <= quantile(dd, quant)) ## points en dessous du quant
        xin <- x1[selected]
        yin <- y1[selected]
        chullchoice <- chull(xin, yin) ## points formant la convex hull
        x2 <- xin[chullchoice]
        y2 <- yin[chullchoice]
        tmp_quant <- c(tmp_quant, list(cbind(x2, y2))) ## coord des points formant le polygone
      }
      names(tmp_quant) <- as.character(chullSize)
    } else
      tmp_quant <- NULL
       
    chulls <- c(chulls, list(tmp_quant))
  }
  names(chulls) <- as.character(levels(fac))
  return(chulls)
}

