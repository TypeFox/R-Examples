fplsr = function (data, order = 6, type = c("simpls", "nipals"), unit.weights = TRUE, 
    weight = FALSE, beta = 0.1, interval = FALSE, method = c("delta", 
        "boota"), alpha = 0.05, B = 100, adjust = FALSE, backh = 10) 
{
    type = match.arg(type)
    rawdata = t(data$y)
    n = dim(rawdata)[1]
    Xtrain = rawdata[1:(n - 1), ]
    Ytrain = rawdata[2:n, ]
    Xtest = as.numeric(rawdata[n, ])
    if (interval == FALSE) {
        if (type == "simpls") {
            if (unit.weights == TRUE) {
                output = unitsimpls(Xtrain, Ytrain, Xtest, order, 
                  weight = weight, beta = beta)
                fitted = t(output$T %*% t(output$Q)) + colMeans(Ytrain)
                residuals = t(Ytrain) - fitted
                out = list(x1 = as.numeric(rownames(Xtrain)), 
                  y1 = as.numeric(colnames(Xtrain)), 
			ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
			y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
                  B = output$B, Ypred = fts(1:dim(Ytrain)[2], 
                    as.matrix(output$Ypred), xname = data$xname, 
                    yname = data$yname), P = output$P, Q = output$Q, 
                  T = output$T, R = output$R, fitted = fts(1:dim(Xtrain)[2], 
                    fitted, xname = data$xname, yname = "Fitted values"), 
                  residuals = fts(1:dim(Xtrain)[2], residuals, 
                    xname = data$xname, yname = "Residual"), 
                  meanX = fts(1:dim(Xtrain)[2], as.matrix(colMeans(Xtrain)), 
                    xname = data$xname, yname = data$yname), 
                  meanY = fts(1:dim(Ytrain)[2], as.matrix(colMeans(Ytrain)), 
                    xname = data$xname, yname = data$yname), 
                  call = match.call())
                return(structure(out, class = "fm"))
            }
            else {
                output = simpls(Xtrain, Ytrain, Xtest, order, 
                  weight = weight, beta = beta)
                fitted = t(output$T %*% t(output$Q)) + colMeans(Ytrain)
                residuals = t(Ytrain) - fitted
                out = list(x1 = as.numeric(rownames(Xtrain)), 
                  y1 = as.numeric(colnames(Xtrain)), 
			ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
			y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
                  B = output$B, Ypred = fts(1:dim(Ytrain)[2], 
                    as.matrix(output$Ypred), xname = data$xname, 
                    yname = data$yname), P = output$P, Q = output$Q, 
                  T = output$T, R = output$R, fitted = fts(1:dim(Xtrain)[2], 
                    fitted, xname = data$xname, yname = "Fitted values"), 
                  residuals = fts(1:dim(Xtrain)[2], residuals, 
                    xname = data$xname, yname = "Residual"), 
                  meanX = fts(1:dim(Xtrain)[2], as.matrix(colMeans(Xtrain)), 
                    xname = data$xname, yname = data$yname), 
                  meanY = fts(1:dim(Ytrain)[2], as.matrix(colMeans(Ytrain)), 
                    xname = data$xname, yname = data$yname), 
                  call = match.call())
                return(structure(out, class = "fm"))
            }
        }
        else {
            output = nipals(Xtrain, Ytrain, Xtest, order, weight = weight, 
                beta = beta)
            out = list(x1 = as.numeric(rownames(Xtrain)), y1 = as.numeric(colnames(Xtrain)), 
				ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
                y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, 
                  yname = data$yname), B = output$B, Ypred = fts(1:dim(Ytrain)[2], 
                  matrix(output$Ypred, dim(Ytrain)[2], ), xname = data$xname, 
                  yname = data$yname), P = output$P, Q = output$Q, 
                T = output$T, R = output$R, meanX = fts(1:dim(Xtrain)[2], 
                  as.matrix(colMeans(Xtrain)), xname = data$xname, 
                  yname = data$yname), meanY = fts(1:dim(Ytrain)[2], 
                  as.matrix(colMeans(Ytrain)), xname = data$xname, 
                  yname = data$yname), Yscores = output$Yscores, 
                projection = output$projection, fitted = fts(1:dim(Xtrain)[2], 
                  t(output$fitted.values[, , order]), xname = data$xname, 
                  yname = "Fitted values"), residuals = fts(1:dim(Xtrain)[2], 
                  t(output$residuals[, , order]), xname = data$xname, 
                  yname = "Residual"), Xvar = output$Xvar, Xtotvar = output$Xtotvar, 
                call = match.call())
            return(structure(out, class = "fm"))
        }
    }
    else {
        fplsrPI(t(Xtrain), t(Ytrain), Xtest, order, method = method, 
            alpha = alpha, B = B, weight = weight, beta = beta, 
            adjust = adjust, backh = backh)
    }
}
