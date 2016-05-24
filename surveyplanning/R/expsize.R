expsize <- function(Yh, H, s2h, poph, Rh = NULL, deffh = NULL, CVh,
                    dataset = NULL) {

  ### Checking

  if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Yh %in% names(dataset))!=1) stop("'Yh' does not exist in 'dataset'!")
      if (min(Yh %in% names(dataset))==1) Yh <- dataset[, Yh, with=FALSE]

      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE]}
      if(!is.null(s2h)) {
          if (min(s2h %in% names(dataset))!=1) stop("'s2h' does not exist in 'dataset'!")
          if (min(s2h %in% names(dataset))==1) s2h <- dataset[, s2h, with=FALSE] }
      if(!is.null(CVh)) {
          if (min(CVh %in% names(dataset))!=1) stop("'CVh' does not exist in 'dataset'!")
          if (min(CVh %in% names(dataset))==1) CVh <- dataset[, CVh, with=FALSE] }
      if(!is.null(poph)) {
          if (min(poph %in% names(dataset))!=1) stop("'poph' does not exist in 'dataset'!")
          if (min(poph %in% names(dataset))==1) poph <- dataset[, poph, with=FALSE] }
      if(!is.null(Rh)) {
          if (min(Rh %in% names(dataset))!=1) stop("'Rh' does not exist in 'dataset'!")
          if (min(Rh %in% names(dataset))==1) Rh <- dataset[, Rh, with=FALSE]}
      if(!is.null(deffh)) {
          if (min(deffh %in% names(dataset))!=1) stop("'deffh' does not exist in 'dataset'!")
          if (min(deffh %in% names(dataset))==1) deffh <- dataset[, deffh, with=FALSE] }
  }

  # Yh
  Yh <- data.table(Yh, check.names=TRUE)
  n <- nrow(Yh)
  m <- ncol(Yh)
  if (any(is.na(Yh))) stop("'Yh' has unknown values")
  if (!all(sapply(Yh, is.numeric))) stop("'Yh' must be all numeric values")
  if (is.null(names(Yh))) stop("'Yh' must be colnames")
  Yh <- data.table(sapply(Yh, as.numeric))

  s2h <- data.table(s2h, check.names=TRUE)
  if (nrow(s2h) != n) stop("'s2h' length must be equal with 'Yh' row count")
  if (ncol(s2h) != m) stop("'s2h' and 'Yh' must be equal column count")
  if (any(is.na(s2h))) stop("'s2h' has unknown values")
  if (!all(sapply(s2h, is.numeric))) stop("'s2h' must be numeric values")
  if (is.null(names(s2h))) stop("'s2h' must be colnames")

  # H
  H <- data.table(H)
  if (nrow(H) != n) stop("'H' length must be equal with 'Yh' row count")
  if (ncol(H) != 1) stop("'H' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(H))) stop("'H' has unknown values")
  if (is.null(names(H))) stop("'H' must be colnames")

  # CVh
  CVh <- data.frame(CVh)
  if (nrow(CVh) != n) stop("'CVh' must be equal with 'Yh' row count")
  if (ncol(CVh) != 1) stop("'CVh' must be 1 column data.frame, matrix, data.table")
  CVh <- CVh[,1]
  if (!is.numeric(CVh)) stop("'CVh' must be numerical")
  if (any(is.na(CVh))) stop("'CVh' has unknown values")

  # poph
  poph <- data.frame(poph)
  if (nrow(poph) != n) stop("'poph' must be equal with 'Yh' row count")
  if (ncol(poph) != 1) stop("'poph' must be vector or 1 column data.frame, matrix, data.table")
  poph <- poph[,1]
  if (!is.numeric(poph)) stop("'poph' must be numerical")
  if (any(is.na(poph))) stop("'poph' has unknown values")

  # Rh
  if (is.null(Rh)) Rh <- rep(1, n)
  Rh <- data.frame(Rh)
  if (nrow(Rh) != n) stop("'Rh' must be equal with 'Yh' row count")
  if (ncol(Rh) != 1) stop("'Rh' must be vector or 1 column data.frame, matrix, data.table")
  Rh <- Rh[, 1]
  if (!is.numeric(Rh)) stop("'Rh' must be numerical")
  if (any(is.na(Rh))) stop("'Rh' has unknown values")

  if (!is.null(deffh)) {
    deffh <- data.table(deffh, check.names=TRUE)
    if (nrow(deffh) != n) stop("'deffh' length must be equal with 'Yh' row count")
    if (ncol(deffh) != m) stop("'deffh' and 'Yh' must be equal column count")
    if (any(is.na(deffh))) stop("'deffh' has unknown values")
    if (!all(sapply(deffh, is.numeric))) stop("'deffh' must be numeric values")
    if (is.null(names(deffh))) stop("'deffh' must be colnames")
   }

  variable <- nh <- estim <- NULL

  CVh <- melt(data.table(H, CVh), id = c(names(H)))
  CVh[, variable := NULL]
  setnames(CVh, "value", "CVh")
  setkeyv(CVh, names(H))

  Rh <- melt(data.table(H, Rh), id = c(names(H)))
  Rh[, variable := NULL]
  setnames(Rh, "value", "Rh")
  setkeyv(Rh, names(H))
  resulth <- merge(CVh, Rh, all = TRUE)

  poph <- melt(data.table(H, poph), id = c(names(H)))
  poph[, variable := NULL]
  setnames(poph, "value", "poph")
  setkeyv(poph, names(H))
  resulth <- merge(resulth, poph, all = TRUE)

  Rh <- CVh <- poph <- NULL

  setnames(s2h, names(s2h), names(Yh))
  s2h <- melt(data.table(H, s2h), id=c(names(H)))
  setnames(s2h, "value", "s2h")
  setkeyv(s2h, c(names(H), "variable"))
  resulth <- merge(s2h, resulth, all = TRUE)
  setkeyv(resulth, c(names(H), "variable"))

  if (!is.null(deffh)) {
    setnames(deffh, names(deffh), names(Yh))
    deffh <- melt(data.table(H, deffh), id = c(names(H)))
    setnames(deffh, "value", "deffh")
    setkeyv(deffh, c(names(H), "variable"))
    resulth <- merge(deffh, resulth, all = TRUE)
  }

  if (is.null(deffh)) resulth[, deffh:=1]
  Yh <- melt(data.table(H, Yh), id=c(names(H)))
  setnames(Yh, "value", "estim")
  setkeyv(Yh, c(names(H), "variable"))
  resulth <- merge(Yh, resulth, all=TRUE)

  Yh <- deffh <- s2h <- NULL

  resulth[, nh := poph ^ 2 * s2h * deffh /
                    (Rh * ((estim * CVh/100) ^ 2 + poph * s2h  * deffh))]
  return(resulth)
}
