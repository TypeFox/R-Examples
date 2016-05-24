### Gruppe 1: Frick, Kondofersky, Speidel

# Funktion zur Berechnung von Btilde (vgl. Paper)
# Uebergeben werden
# - ein Vektor x
# - Knoten
# - Hyperparameter lambda2
Btilde <- function(x, knots, lambda2){
    B <- bs(x, knots = knots[-c(1, length(knots))], intercept = FALSE)
    D <- diff(diag(ncol(B)), differences = 2)
    Omega <- crossprod(D)
    n <- length(x)
    M <- 1 / n * crossprod(B) + lambda2 * Omega
    R <- chol(M)
    R1 <- solve(R)
    Btilde <- B %*% R1
}



### Formel interface y ~ s(x1) + s(x2) ?
### lineare Terme?

# Funktion hgam.fit
hgam.fit <- function(y, x, weights = rep(1, nrow(x)), intercept = TRUE, 
                     nknots = 20, lambda1 = 2, lambda2 = 3, model = LinReg(), 
                     ...){

    # Anzahl Spalten in der Datenmatrix
    p <- ncol(x)

    # Berechnung der Knoten der einzelnen Kovariablen und Abspeichern in Liste
    knots <- lapply(x, quantile,
        prob = seq(from = 0, to = 1, length = nknots + 2))

    # Erweiterung von Btilde (durch lexical scoping)
    # myBtilde gibt Designmatrix aus, wobei Knoten und Hyperparameter
    # (und evtl. auch Daten) aus dem hoeheren environment (hgam.fit) benutzt werden
    # - xnew sind die (ggf. neuen) Daten
    # - which-Argument steuert, welche Kovariable gefittet werden soll und
    #   kann als charachter oder als Zahl angegeben werden.
    # - intercept (boolean) steuert, ob ein Intercept mitmodelliert werden soll
    myBtilde <- function(xnew = x, which = 1:ncol(x), 
                         intercept = TRUE, la2 = lambda2) {

        # Falls which als charachter angegeben wird,
        # wird es mittels match in einen Vektor ueberfuehrt, der angibt,
        # an welcher Stelle in den Daten die angegebenen Praediktoren stehen
        if (is.character(which))
        which <- match(which, names(x))

        # Falls mit which nicht korrekt angeben wurde, welche Praediktoren
        # verwendet werden sollen, und damit auch in den neuen Daten vorhanden
        # sein sollten, wird abgebrochen.
        stopifnot(all(names(x)[which] %in% names(xnew)))

        # Umsortieren von den (neuen) Daten - Reihenfolge kann vom Benutzer
        # festgelegt werden
        xnew <- xnew[, names(x)[which], drop = FALSE]

        # Berechnen der Designmatrix durch sukzessives anwenden von Btilde
        # auf einzelne Spalten von xnew. Dabei werden auch die jeweiligen
        # Knoten sukzessive verwendet. Ausgabe ist eine Liste, die je eine
        # Designmatrix pro Praediktor enthaelt
        X <- lapply(1:length(xnew), function(i)
               Btilde(x = xnew[, i], knots = knots[[which[i]]],
                      lambda2 = la2)
               )

        # Festlegen des index-Vektors
        index <- rep(which, sapply(X, ncol))

        # X aus Liste in eine grosse Designmatrix umwandeln
        X <- do.call("cbind", X)

        # Falls es einen Intercept gibt, die entsprechende Spalte in die
        # Design-Matrix einfuegen und den index-Vektor mit 0 erweitern
        if (intercept) {
            X <- cbind(1, X)
            index <- c(0, index)
        }

        # index als Attribut der Designmatrix mitgeben
        attr(X, "assign") <- index

        # X zurueckgeben
        return(X)
    }

    # Designmatrix fuer aktuelle Daten (also keine neuen) berechnen
    X <- myBtilde(intercept = intercept)

    # index Vektor aus X herausziehen
    index <- attr(X, "assign")

    # Falls intercept mitberechnet wurde, dann dessen index-Eintrag (0) auf
    # NA setzen, da grplasso dies so wuenscht
    index[index < 1] <- NA

    # grplasso aufrufen und Designmatrix, index und Hyperparameter uebergeben
    mod <- grplasso(x = X, y = y, weights = weights, index = index, 
                    model = model, lambda = lambda1, 
                    control = grpl.control(trace = 0), 
                    center = FALSE, standardize = FALSE)

    # Koeffizienten aus mod herauslesen
    coef <- coef(mod)

    # Koeffizienten benennen
    nm <- names(x)[index]
    nm[is.na(nm)] <- "Intercept"
    rownames(coef) <- nm

    # return-Objekt vorbereiten - es werden Response, Designmatrix,
    # Koeffizienten und Funktion myBtilde zurueckgegeben
    ret <- list(y = y, x = x, Btilde = X, coef = coef, 
                Btildenew = myBtilde, model = model, weights = weights)

    # Klasse von ret festlegen und returnieren
    class(ret) <- "hgam"
    return(ret)
}


# Hauptfunktion, die dem Anwender praesentiert wird
hgam <- function(formula, data = NULL, weights, model = LinReg(), 
                 nknots = 20, lambda1 = 2, lambda2 = 3, ...){
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
    if(nknots <= 0) stop("nknots should be positive")
    if(nknots %% 1 != 0) stop("nknots should be an integer")
    if(!is.numeric(lambda1)) stop("lambda1 should be numeric")
    if(lambda1 < 0) stop("lambda1 cannot be negative")
    if(!is.numeric(lambda2)) stop("lambda2 should be numeric")
    if(lambda2 < 0) stop("lambda2 cannot be negative")

    # Weitergabe der vorbereiteten Angaben an hgam.fit
    if (missing(weights)) weights <- rep(1, nrow(x))
    hgam.fit(y = y, x = x, weights = weights, model = model, 
             intercept = intercept, nknots = nknots, 
             lambda1 = lambda1, lambda2 = lambda2)
}



# Funktion zur Berechnung von Risiko
hrisk <- function(object, folds = 10,
                  type = c("cv", "bootstrap", "subsampling"),
                  nlambda1 = 10, lambda2 = 1:10, trace = TRUE,
                  papply = if (require("multicore")) mclapply else lapply){

    y <- model.response.hgam(object)
    n <- length(y)
    x <- object$x
    intercept <- rownames(coef(object))[1] == "Intercept"
    index <- attributes(object$Btilde)$assign
    index.lambda <- index                    
    if(intercept) index.lambda[1] <- NA
    model <- object$model
    w <- matrix(NA, ncol = folds, nrow = n)
    
    err <- matrix(0, nrow = nlambda1, ncol = length(lambda2))
    colnames(err) <- lambda2

    type <- match.arg(type)

    if (type == "cv") {
        ind <- split(sample(n), rep(1:folds, length = n))
        for (i in 1:folds) {
            w[,i] <- rep(1, n)
            w[ind[[i]], i] <- 0
        }
    }

    else if(type == "bootstrap"){
        for(i in 1:folds) w[, i] <- rmultinom(1, n, rep(1, n) / n)
    }

    else if(type == "subsampling"){
        for(i in 1:folds){
            ind <- split(sample(n), rep(1:2, length = n))
            w[,i] <- rep(1, n); w[ind[[1]], i] <- 0}
    }
    
    args <- expand.grid(fold = 1:folds, lambda2 = lambda2)

    risk <- papply(lambda2, function(la2) {
        l1.max <-  lambdamax(x = object$Btildenew(
            intercept = intercept, la2 = la2), 
            y = y, index = index.lambda, 
            center = FALSE, standardize = FALSE, model = model)
        lambda1 <- rev(exp(seq(log(1.2), log(l1.max + 1), 
                           length.out = nlambda1)) - 1)
        err <- 0
        for (i in 1:folds) {
            act.mod <- hgam.fit(y = y, x = x, weights = w[, i], 
                                nknots = 20,
                                model = object$model, intercept = intercept,
                                lambda1 = lambda1, 
                                lambda2 = la2)
                # nknots hier als 20 fest, evtl. muss es ubergeben werden
             f <- fitted(act.mod)
             oob <- w[, i] == 0
             err <- err + apply(f, 2, function(ff) 1/sum(w[oob,]) * act.mod$model@nloglik(
                          y[oob], ff[oob], rep(1, sum(oob))))
        }
        err <- err / folds
        list(lambda1 = lambda1, risk = err)
    })
    la2opt <- lambda2[mm <- which.min(sapply(risk, function(r) min(r$risk)))]
    la1opt <- risk[[mm]]$lambda1[which.min(risk[[mm]]$risk)]

    tmp <- c()
    for (i in 1:length(lambda2)) {
        tmp <- rbind(tmp, data.frame(lambda2 = lambda2[i],
                                     lambda1 = risk[[i]]$lambda1,
                                     risk = risk[[i]]$risk))
    }

    ret <- list(lambda1 = la1opt, lambda2 = la2opt, risk = tmp, folds = folds, type = type)
    class(ret) <- "hrisk"
    return(ret)
}


#set.seed(1234)
#test.d <- dgp(100)
#mod <- hgam(y ~ ., data = test.d)
#a <- hrisk(mod, folds = 3, nlambda1 = 10, lambda2 = 1:10)
