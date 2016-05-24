imputeBasic <- function(data, Init = "mean") {
  dat <- data
  dat1a <- dat[(sapply(dat, function(x) is.numeric(x)) == T)]
  if(ncol(dat1a) != 0) {
    if(Init == "mean") {
    dat1 <- as.data.frame(apply(dat1a, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
    } else {
    dat1 <- as.data.frame(apply(dat1a, 2, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)))
    }
    Initials <- dat1[sapply(dat1a, function(x) is.na(x))]
  } else {
    Initials <- "No Continous Variables in Dataset"
    dat1 <- dat1a
  }
  dat.factors.a <- dat[(sapply(dat, function(x) is.factor(x)) == T)]
  if(ncol(dat.factors.a) > 0) {
    dat.factors.b <- apply(dat.factors.a, 2, function(x) ifelse(x == "" | is.na(x), "Missing", x))
    dat.factors <- data.frame(sapply(1:ncol(dat.factors.b), function(x) as.factor(dat.factors.b[, x])))
    names(dat.factors) <- names(dat.factors.a)
    my.cols <- names(dat.factors)
    Imputed.Factor.Levels <- sapply(my.cols, function(this.col) {
      B <- table(dat.factors[, this.col])
      NL <- sort(B[names(B) != "Missing"], decreasing = T)
      NL2 <- names(NL)[1]
      # NL2 <- sample(names(NL), 1, prob = NL / sum(NL))
    })
    dat2 <- as.data.frame(sapply(1:length(Imputed.Factor.Levels), function(this.col) {
      # this.col <- 1
      NL2 <- levels(as.factor((as.character(Imputed.Factor.Levels))[this.col]))
      levels(dat.factors[, this.col])[levels(dat.factors[, this.col]) == "Missing"] <- NL2
      dat.factors[, this.col]
    })); names(dat2) <- names(dat.factors)
    dat[, names(dat2)] <- dat2
    dat.cat.imp <- dat
    dat[, names(dat1)] <- dat1
    output <- list(Imputed.DataFrame = dat, 
                   Imputed.Missing.Continous = Initials, 
                   Imputed.Missing.Factors = Imputed.Factor.Levels)
  } else {
    Imputed.Factor.Levels <- "No Categorical Variables in Dataset"
    dat <- dat1
    output <- list(Imputed.DataFrame = dat, 
                   Imputed.Missing.Continous = Initials, 
                   Imputed.Missing.Factors = NA)
  }
  class(output) <- "roughImputation"
  output
}

print.roughImputation <- function(x, ...) {
  print(x$Imputed.DataFrame)
}
