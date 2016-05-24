#-----------------------------------------------------------------------------#
# Cooperative games - 3 players                                               #
# Imputation set and core                                                     #
#-----------------------------------------------------------------------------#

# point3Dto2D -----------------------------------------------------------------
# Transformar coordenadas de un punto de 3D a 2D
point3Dto2D <- function(x, alpha, beta) {
  
  # Calcular coordenadas en 2D
  x2d <- x[2]*cos(beta*(pi/180)) - x[1]*cos(alpha*(pi/180))
  y2d <- x[3] + x[1]*sin(alpha*(pi/180)) + x[2]*sin(beta*(pi/180))
  
  # Devolver punto en 2D
  return(c(x2d, y2d))
  
}
#-----------------------------------------------------------------------------#

# lineEq3D --------------------------------------------------------------------
# Ecuación de rectas en 3 dimensiones
lineEq3D <- function(A, B) {
  
  coefx <- B[1] - A[1]
  coefy <- B[2] - A[2]
  coefz <- B[3] - A[3]
  coefMat <- matrix(c(coefx, coefy, coefz))
  rhsMat <- matrix(c(A[1], A[2], A[3]))
  
  return(list(coef = coefMat, rhs = rhsMat))
  
}
#-----------------------------------------------------------------------------#

# intersec3D ------------------------------------------------------------------
intersec3D <- function(A, B) {
  
  # Construir matriz de coeficientes y de lados derechos
  coefMat <- cbind(-A$coef, B$coef); coefMat
  rhsMat <- A$rhs - B$rhs; rhsMat
  # Resolvemos el sistema para hallar el punto de corte
  a <- coefMat[1:2, ]; a
  b <- rhsMat[1:2, ]; b
  
  tparam <- solve(a, b); tparam
  
  # Puntode de intersección en cada recta
  pointA <- A$rhs + A$coef*tparam[1]
  pointB <- B$rhs + B$coef*tparam[2]
  
  # Si son el mismo las rectas se cortan en él
  if (all(pointA == pointB)) {
    intersecPoint <- c(pointA)
    return(intersecPoint)
  } else {
    cat("No se intersecan")
  }
  
}
#-----------------------------------------------------------------------------#

# getImputationSet ------------------------------------------------------------
# Conjunto de imputaciones

getImputationSet <- function(players, game) {
  
  # Vértices conjunto de imputaciones
  n <- players  # numero de jugadores
  v <- game  # valores del juego
  N <- length(v); N  # posición de la coalición total
  
  # El conjunto de imputaciones viene definido por los puntos:
  # 1 - (v(N)-x2-x3, x2, x3)
  v1 <- c(v[N] - v[2] - v[3], v[2], v[3])
  # 2 - (x1, v(N)-x1-x3, x3)
  v2 <- c(v[1], v[N] - v[1] - v[3], v[3])
  # 3 - (x1, x2, v(N)-x1-x2)
  v3 <- c(v[1], v[2], v[N] - v[1] - v[2])
  
  imputationSet <- matrix(c(v1, v2, v3), ncol = 3, byrow = TRUE)
  
  # Devolver vértices del conjunto de imputaciones
  return(imputationSet)
  
}
#-----------------------------------------------------------------------------#

# plotImputationSet -----------------------------------------------------------
# Dibujar conjunto de imputaciones
plotImputationSet <- function(imputationSet) {
  
  # Vértices a matriz de coordenadas en 2D
  coordMat <- matrix(ncol = 2)[-1, ]
  for (i in 1:nrow(imputationSet)) {
    coordMat <- rbind(coordMat, point3Dto2D(imputationSet[i, ], 30, 30))
  }
  coordMat
  
  # Axis limits
  minxlim <- min(coordMat[, 1]); minxlim
  maxxlim <- max(coordMat[, 1]); maxxlim
  extraxlim <- 0.005 * (maxxlim - minxlim); extraxlim
  minylim <- min(coordMat[, 2]); minylim
  maxylim <- max(coordMat[, 2]); maxylim
  extraylim <- 0.015 * (maxylim - minylim); extraylim
  
  # Representación
  par(mar = rep(1,4))
  plot(coordMat[, 1], coordMat[, 2], type = "n",
       axes = FALSE, xlab = NA, ylab = NA,
       xlim = c(minxlim - extraxlim, maxxlim + extraxlim), ylim = c(minylim - extraylim, maxylim + extraylim))
  box(col = "gray")
  polygon(coordMat[, 1], coordMat[, 2], col = "gray95")
  
  # Vertices labels
  for (i in 1:nrow(coordMat)) {
    vLabel <- paste("(", -imputationSet[i, 1], ",", -imputationSet[i, 2], ",", 
                    -imputationSet[i, 3], ")"); vLabel
    if (i == 1) {
      text(coordMat[i, 1], coordMat[i, 2], vLabel, cex = 0.6, pos = 1)
    } else if (i == 2) {
      text(coordMat[i, 1], coordMat[i, 2], vLabel, cex = 0.6, pos = 1)
    } else {
      text(coordMat[i, 1], coordMat[i, 2], vLabel, cex = 0.6, pos = 3)
    }
  }
  
}

# coreSet ---------------------------------------------------------------------
# Núcleo del juego cooperativo
getCoreSet <- function(players, game) {
  
  # Vértices conjunto de imputaciones
  n <- players  # numero de jugadores
  v <- game  # valores del juego
  N <- length(v); N  # posición de la coalición total
  
  # Puntos de partida para obtener el núcleo
  
  # Conjunto de imputaciones
  impSet <- getImputationSet(3, v); impSet
  
  # Puntos de corte de las rectas con el conjunto de imputaciones
  # x3
  x121 <- c(v[1], v[4] - v[1], v[N] - v[4]); x121
  x122 <- c(v[4] - v[2], v[2], v[N] - v[4]); x122
  # x2
  x131 <- c(v[1], v[N] - v[5], v[5] - v[1]); x131
  x132 <- c(v[5] - v[3], v[N] - v[5], v[3]); x132
  # x1
  x231 <- c(v[N] - v[6], v[2], v[6] - v[2]); x231
  x232 <- c(v[N] - v[6], v[6] - v[3], v[3]); x232
  # Juntarlos
  coreset <- matrix(c(x122, x121,
                      x131, x132,
                      x232, x231), ncol = 3, byrow = TRUE)
  csetTemp <- coreset; csetTemp
  
  # Puntos de intersección entre las rectas
  # Revisar tres puntos de corte de las rectas que limitan el núcleo
  # Recta x3 con x1
  A3D <- lineEq3D(coreset[1, ], coreset[2, ]); A3D
  B3D <- lineEq3D(coreset[5, ], coreset[6, ]); B3D
  if (!(all(A3D$coef == 0) || all(B3D$coef == 0))) {
    x3x1 <- intersec3D(A3D, B3D); x3x1
  } else {
    x3x1 <- NULL
  }
  # Recta x2 con x1
  A3D <- lineEq3D(coreset[3, ], coreset[4, ]); A3D
  B3D <- lineEq3D(coreset[5, ], coreset[6, ]); B3D
  if (!(all(A3D$coef == 0) || all(B3D$coef == 0))) {
    x2x1 <- intersec3D(A3D, B3D); x2x1
  } else {
    x2x1 <- NULL
  }
  # Recta x2 con x3
  A3D <- lineEq3D(coreset[3, ], coreset[4, ]); A3D
  B3D <- lineEq3D(coreset[1, ], coreset[2, ]); B3D
  if (!(all(A3D$coef == 0) || all(B3D$coef == 0))) {
    x2x3 <- intersec3D(A3D, B3D); x2x3
  } else {
    x2x3 <- NULL
  }
  # Unirlos
  interSet <- matrix(c(x3x1, x2x1, x2x3), ncol = 3, byrow = TRUE); interSet
  
  # Posibles núcleo
  coreTemp <- rbind(impSet, csetTemp, interSet); coreTemp
  
  # Borrar duplicados
  coreTemp <- unique(coreTemp); coreTemp
  
  # Selección de los puntos
  # Límites de x3
  x3min <- impSet[1, 3]; x3min
  x3max <- csetTemp[1, 3]; x3max
  # Límites de x2
  x2min <- impSet[1, 2]; x2min
  x2max <- csetTemp[3, 2]; x2max
  # Límites de x1
  x1min <- impSet[2, 1]; x1min
  x1max <- csetTemp[5, 1]; x1max
  
  # Construir nuevo coreSet
  coreSet <- matrix(ncol = 3)[-1, ]
  # Revisar punto por punto y eliminar el que no se encuentren entre los anteriores
  for (i in 1:nrow(coreTemp)) {
    revPoint <- coreTemp[i, ]; revPoint
    # Comprobar que queda entre los límites
    limx1 <- revPoint[1] >= x1min & revPoint[1] <= x1max
    limx2 <- revPoint[2] >= x2min & revPoint[2] <= x2max
    limx3 <- revPoint[3] >= x3min & revPoint[3] <= x3max
    # Si cumple todos los límites se añade al núcleo
    if (all(limx1, limx2, limx3)) {
      coreSet <- rbind(coreSet, revPoint)
    }
  }
  
  # Remove row names
  row.names(coreSet) <- NULL
  
  if (nrow(coreSet) == 0) {
    stop("Empty core")
  }
  
  return(coreSet)  
  
}
#-----------------------------------------------------------------------------#

# plotCoreSet -----------------------------------------------------------------
plotCoreSet <- function(coreSet) {
  
  coreset2d <- matrix(ncol = 2)[-1, ]
  for (i in 1:nrow(coreSet)) {
    coreset2d <- rbind(coreset2d, point3Dto2D(coreSet[i, ], 30, 30))
  }
  coreset2d
  
  # Estimate the center of the polygon with the mean of the x coords and
  # the mean of the y coords.
  centerx <- mean(coreset2d[, 1]); centerx
  centery <- mean(coreset2d[, 2]); centery
  center2d <- c(centerx, centery); center2d
  #points(center2d[1], center2d[2], pch = 19)
  
  # Calculate the angle of each point from that center point
  atanset <- matrix(ncol = 2)[-1, ]
  for (i in 1:nrow(coreset2d)) {
    atanset <- rbind(atanset, atan2(coreset2d[i, 2]-center2d[2], coreset2d[i, 1]-center2d[1]))
  }
  atanset
  atanset*180/pi
  
  # sort the data by the angle calculated
  coreset2d <- matrix(coreset2d[order(atanset[, 1]), ], ncol = 2); coreset2d
  
  points(coreset2d, pch = 19)
  polygon(coreset2d, lwd = 2, col = "gray75")
  
}
#-----------------------------------------------------------------------------#