L_hNV <- function(p, y = y, X = X, sc = sc) {
  b <- p[1:(length(p)-2)]
  sigmau2 <- p[(length(p)-1)]
  sigmav2 <- p[length(p)]
  epsilon <- y - X%*%b
  N <- length(y)
  ret = -(N * log(sqrt(2) / sqrt(pi)) + N * log(1 / (sqrt(sigmau2 + sigmav2)))
  + sum(log(pnorm(-sc*(epsilon * (sqrt(sigmau2) / sqrt(sigmav2))) / (sqrt(sigmau2 + sigmav2)))))
  - 1 / (2 * (sigmau2 + sigmav2)) * sum(epsilon^2))
  names(ret) <- "Log-Lik normal/half-normal distribution"
  return(ret)
  }
  
L_exp <- function(p, y = y, X = X, sc = sc) {
  b <- p[1:(length(p)-2)]
  sigmau2 <- p[(length(p)-1)]
  sigmav2 <- p[length(p)]
  epsilon <- y - X%*%b
  N <- length(y)
  ret = -(N * log(1/sqrt(sigmau2)) + N / 2 * (sigmav2 / sigmau2)
  + sum(log(pnorm((-sc*epsilon - (sigmav2 / (sqrt(sigmau2)))) / (sqrt(sigmav2)))))
  + sc / sqrt(sigmau2) * sum(epsilon))
  names(ret) <- "Log-Lik normal/exponential distribution"
  return(ret)
  }

L_trunc <- function(p, y = y, X = X, sc = sc) {
  b <- p[1:(length(p)-3)]
  sigmau2 <- p[length(p)-2]
  sigmav2 <- p[length(p)-1]
  mu <- p[length(p)]
  epsilon <- y - X%*%b
  N <- length(y)
  ret = -(N / 2* log(1 / (2 * pi)) + N * log(1 / (sqrt(sigmau2 + sigmav2)))
  - log(pnorm((mu/sqrt(sigmau2))))
  + sum(log(pnorm(((1 - (sigmau2/(sigmau2 + sigmav2)))
          - sc * (sigmau2 / (sigmau2 + sigmav2)) * epsilon)
        / (sqrt(sigmau2 * (1 - sigmau2 / (sigmau2 + sigmav2)))))))
  - 1 / (2 * (sigmau2 + sigmav2)) * sum((epsilon + sc * mu)^2))
  names(ret) <- "Log-Lik normal/truncated-normal distribution"
  return(ret)
  }

L_trunc_mufest <- function(p, mu = mu, y = y, X = X, sc = sc) {
  b <- p[1:(length(p)-2)]
  sigmau2 <- p[length(p)-1]
  sigmav2 <- p[length(p)]
  epsilon <- y - X%*%b
  N <- length(y)
  ret = -(N / 2* log(1 / (2 * pi)) + N * log(1 / (sqrt(sigmau2 + sigmav2)))
  - log(pnorm((mu/sqrt(sigmau2))))
  + sum(log(pnorm(((1 - (sigmau2/(sigmau2 + sigmav2)))
          - sc * (sigmau2 / (sigmau2 + sigmav2)) * epsilon)
        / (sqrt(sigmau2 * (1 - sigmau2 / (sigmau2 + sigmav2)))))))
  - 1 / (2 * (sigmau2 + sigmav2)) * sum((epsilon + sc * mu)^2))
  names(ret) <- "Log-Lik normal/truncated-normal distribution"
  return(ret)
  }
  
sfa.fit <- function(y, x, intercept = TRUE, fun = "hnormal",
              pars = NULL, par_mu = NULL, form = "cost", method = "BFGS", ...){
    # Anzahl Spalten in der Datenmatrix
    p <- ncol(x)
    X <- as.matrix(x)
    mu <- NULL
    if ((is.null(par_mu)) == FALSE && fun != "tnormal") {
        par_mu = NULL
        print("mu only expected for the truncatednormal case. Set par_mu = NULL")
        }
    if (fun == "hnormal") maxlik = L_hNV
    if (fun == "tnormal") {
        maxlik = L_trunc
        if(is.null(par_mu) == FALSE) {
            maxlik = L_trunc_mufest
            mu = par_mu
            }
        }
    if (fun == "exp") maxlik = L_exp
    if (form == "cost") sc = -1
    if (form == "production") sc = 1
    if (intercept) {
        X <- cbind(1, X)
        }
    daten <- data.frame(y, X)
    ols <- lm(y ~ ., data = daten)
       if (intercept) {
            ols <- lm(y ~ .- 1, data = daten)
            }
if (is.null(pars)) {
        pars <- coef(ols)
        if (intercept) {
            names(pars) <- c("Intercept", names(pars)[2:length(pars)])
            }
        res_ols <- resid(ols)
        su <- var(res_ols)*(1-2/pi)
        sv <- var(res_ols)
        pars <- c(pars, sigmau2 = su, sigmav2 = sv)
    if (fun == "tnormal" && is.null(par_mu)) {
        mu = 0
        pars <- c(pars, mu = mu)
        }
    }
    if (is.null(par_mu)) mod <- optim(pars, fn = maxlik, y = y, X = X, sc = sc, hessian = TRUE, method = method)
    if (is.null(par_mu) == FALSE) mod <- optim(pars, fn = maxlik, y = y, X = X, sc = sc, mu = par_mu, hessian = TRUE, method = method)
    
    # Koeffizienten aus mod herauslesen
    coef <- mod$par
    
    # Wert der Log-Liklihood zurückgeben
    if (is.null(par_mu)) LogLik <- - maxlik(coef, y = y, X = X, sc = sc)
    if (is.null(par_mu) == FALSE) LogLik <- - maxlik(coef, y = y, X = X, sc = sc, mu = par_mu)
    
    # return-Objekt vorbereiten - es werden Response, Designmatrix, Parameter, Wert der LogLiklihood,
    #         verwendete MaximumLiklihoodfunktion, name der Log-Lik Funkiton,
    #         Parameter für Art der Funktion (-1 = Kosten, 1 = Produktions), Die Hessematrix,
    #         Das Ergebnis des einfachen linearen Modells
    ret <- list(y = y, x = x, X = X, sigmau2 = coef[names(coef) == "sigmau2"],
                sigmav2 = coef[names(coef) == "sigmav2"], coef = coef, mu = mu,
                logLik = LogLik, maxlik = maxlik, fun = fun, sc = sc, hess = mod$hessian,
                ols = ols, par_mu = par_mu)
if(is.null(mu)){ret <- list(y = y, x = x, X = X, sigmau2 = coef[names(coef) == "sigmau2"],
				sigmav2 = coef[names(coef) == "sigmav2"], coef = coef,
				logLik = LogLik, maxlik = maxlik, fun = fun, sc = sc, hess = mod$hessian,
				ols = ols, par_mu = par_mu)
	}
    # Klasse von ret festlegen und returnieren
    class(ret) <- "sfa"
    return(ret)
}


# Hauptfunktion, die dem Anwender praesentiert wird
sfa <- function(formula, data = NULL, intercept = TRUE, fun = "hnormal", pars = NULL, par_mu = NULL, form = "cost", method = "BFGS", ...){

    # Extrahieren von Designmatrix, Response, Intercept (boolean) aus
    # formula-Argument
    mf <- model.frame(formula, data)
    y <- model.response(mf)
    x <- mf[-1]
    tm <- attr(mf, "terms")
    intercept <- attr(tm, "intercept") == 1

    # Fehlermeldungen abfangen
    if(!is.vector(y)) stop("y is not a vector")
    if(!is.data.frame(x)) stop("x is not a data frame")
    if(length(y) != nrow(x)) stop("x and y lengths differ")

    # Weitergabe der vorbereiteten Angaben an sfa.fit
    sfa.fit(y = y, x = x, intercept = intercept, fun = fun, pars = pars, par_mu = par_mu, form = form, method = method)
}


