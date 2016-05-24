fitmm <-
function (data, marpars, homogeneous, restricted, add) 
{
    LORstr <- marpars$LORstr
    LORem <- marpars$LORem
    fmla <- marpars$fmla
    ncategories <- max(data$x)
    timepairs <- max(data$tp)
    if (any(data$counts == 0)) 
        data$counts <- data$counts + add
    LORterm <- matrix(0, timepairs, ncategories^2)
    if (LORstr == "uniform" | LORstr == "category.exch") {
        suppressWarnings(fitted.mod <- gnm(fmla, family = poisson, 
            data = data, verbose = FALSE, model = FALSE))
        if (is.null(fitted.mod)) 
            stop("gnm did not converge algorithm")
        coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
            "x:y")])
        if (LORem == "2way") 
            coefint <- mean(coefint)
        LORterm <- matrix(coefint, timepairs, ncategories^2)
        LORterm <- t(apply(LORterm, 1, function(x) exp(x * tcrossprod(1:ncategories))))
    }
    if (LORstr == "time.exch") {
           data$x <- factor(data$x)
           data$y <- factor(data$y)
        if (is.null(restricted)) {
            if (LORem == "3way") {
                if (homogeneous) {
                  suppressWarnings(fitted.mod <- gnm(fmla, family = poisson, 
                    data = data, verbose = FALSE, model = FALSE))
                  if (is.null(fitted.mod)) 
                    stop("gnm did not converge algorithm")
                  coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
                    "MultHomog")])
                  coefint <- c(tcrossprod(coefint))
                }
                else {
                  suppressWarnings(fitted.mod <- gnm(fmla, family = poisson, 
                    data = data, verbose = FALSE, model = FALSE))
                  if (is.null(fitted.mod)) 
                    stop("gnm did not converge algorithm")
                  coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
                    "Mult")])
                  coefint <- c(tcrossprod(coefint[-c(1:ncategories)], 
                    coefint[1:ncategories]))
                }
                LORterm <- exp(matrix(coefint, nrow = timepairs, 
                  ncol = ncategories^2, TRUE))
            }
            else {
                LORterm2 <- LORterm
                for (i in 1:timepairs) {
                  datamar <- data[data$tp == i, ]
                  suppressWarnings(fitted.mod <- gnm(fmla, family = poisson, 
                    data = datamar, verbose = FALSE, model = FALSE))
                  if (homogeneous) {
                    coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
                      "MultHomog")])
                    coefint <- c(tcrossprod(coefint))
                  }
                  else {
                    coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
                      "Mult")])
                    coefint <- c(tcrossprod(coefint[1:ncategories], 
                      coefint[-c(1:ncategories)]))
                  }
                  LORterm2[i, ] <- coefint
                }
                LORterm2 <- colMeans(LORterm2)
                LORterm <- exp(matrix(LORterm2, timepairs, 
                  ncategories^2, TRUE))
            }
        }
        else {
            if (LORem == "3way") {
                if (homogeneous) {
                  coefint <- RRChomog(fmla, data, ncategories)
                }
                else {
                  coefint <- RRCheter(fmla, data, ncategories)
                }
                LORterm <- exp(matrix(coefint, nrow = timepairs, 
                  ncol = ncategories^2, TRUE))
            }
            else {
                LORterm2 <- LORterm
                for (i in 1:timepairs) {
                  datamar <- data[data$tp == i, ]
                  if (homogeneous) {
                    coefint <- RRChomog(fmla, datamar, ncategories)
                  }
                  else {
                    coefint <- RRCheter(fmla, datamar, ncategories)
                  }
                  LORterm2[i, ] <- coefint
                }
                LORterm2 <- colMeans(LORterm2)
                LORterm <- exp(matrix(LORterm2, timepairs, 
                  ncategories^2, TRUE))
            }
        }
    }
    if (LORstr == "RC") {
        data$x <- factor(data$x)
        data$y <- factor(data$y)   
     for (i in 1:timepairs) {
            datamar <- data[data$tp == i, ]
            suppressWarnings(fitted.mod <- gnm(fmla, family = poisson, 
                data = datamar, verbose = FALSE, model = FALSE))
            if (is.null(restricted)) {
                if (homogeneous) {
                  coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
                    "MultHomog")])
                  coefint <- c(tcrossprod(coefint))
                }
                else {
                  coefint <- as.vector(coef(fitted.mod)[pickCoef(fitted.mod, 
                    "Mult")])
                  coefint <- c(tcrossprod(coefint[1:ncategories], 
                    coefint[-c(1:ncategories)]))
                }
            }
            else {
                coefint <- if (homogeneous) 
                  RRChomog(fmla, datamar, ncategories)
                else RRCheter(fmla, datamar, ncategories)
            }
            LORterm[i, ] <- exp(coefint)
        }
    }
    LORterm <- prop.table(LORterm, 1)
    LORterm
}

