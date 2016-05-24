stima <-
function (data, maxsplit, model="regtrunk", first = NULL, vfold = 10, CV = 1, Save = FALSE, 
    control = NULL, printoutput = TRUE)
{  #default model is a regression trunk; Classification trunk model is not available yet.    
    maxint <- maxsplit - 1
    n <- nrow(data)
    if (sum(apply(is.na(data[, ]), 2, sum)) != 0) {
        stop("Missing values in the data set are not allowed. Run the function rta again on the data set without missing values.")
    }
    type <- character((ncol(data)))
    ##short version
    type<-sapply(data,"class")
    ##end short version
    if (type[1] == "factor") {
        stop("The first column of the dataset is the response variable in the regression trunk model. This may not be a categorical variable.")
    }
    if (sum(type == "factor") != 0) {
        catmat <- data.frame(data[, type == "factor"])
        lev <- numeric(sum(type == "factor"))
        for (j in 1:ncol(catmat)) {
            lev[j] <- min(summary(catmat[, j]))
            if (lev[j] == 0) {
                catmat[, j] <- factor(catmat[, j])
                data[, type == "factor"][, j] <- factor(catmat[,
                  j])
                lev[j] <- min(summary(catmat[, j]))
            }
        }
        if (vfold != 0 & min(lev) <= 3) {
            stop("Some categories of the variable(s) '", names(catmat)[lev <=
                3], "' have small marginal frequencies. We recommend to increase the number of sets defined by vfold or merge categories.")
        }
    }
    if (is.null(control)) {
        control <- stima.control()
  
    }
    #print(control)
    ref <- control$ref
    if (is.null(first)) {
        first <- stima.pre(data, maxsplit, vfold = vfold, controlpre = control)
    }
    predtrunk <- control$predtrunk
    sel <- control$sel
    predsel <- control$predsel
    firstvar <- data.frame(data[, first])
    names(firstvar) <- names(data)[first]
    data0 <- data[, -first]
    typ <- type[-first] == "factor"
    type <- ifelse(type == "factor", 1, 2)
    class1 <- type[first]
     ##short version by Cor
     dataRTA<-cbind(as.data.frame(sapply(data0,"dummy")),firstvar)
	if(sum(typ) != 0){dataCAT<-data
	data[-first]<-sapply(data[-first],"as.numeric")}
	#end short version
    if (is.null(predtrunk)) {
        if (class1 == 1) {
            predtrunk <- c(1:ncol(data))[-c(1, first)]
        }
        else {
            predtrunk <- c(2:ncol(data))
        }
    }
    if (is.null(control$minbucket)) {
        minbucket <- floor(sqrt(n))
        if (vfold != 0) {
            minbucketCV <- floor(sqrt(n - (n/vfold)))
        }
        else {
            minbucketCV <- minbucket
        }
    }
    else {
        minbucket <- control$minbucket
        minbucketCV <- minbucket
    }
    minsplitCV <- 2 * minbucketCV
    minsplit <- 2 * minbucket
    
   
    if (class1 == 1) {
        CAT <- as.factor(dataRTA[, ncol(dataRTA)])
        indmat <- indicator(CAT)
        namesold <- names(dataRTA[, -ncol(dataRTA)])
        ncat <- ncol(indmat)
        namescat <- paste(names(firstvar), attributes(CAT)$levels[1:ncat],
            sep = "")
        data0 <- dataRTA[, -ncol(dataRTA)]
        data[, first] <- as.numeric(data[, first])
        dataRTA <- cbind(dataRTA[, -ncol(dataRTA)], indmat[,
            -ncat])
        names(dataRTA) <- c(namesold, namescat[-ncat])
    }
    if (class1 == 2) {
    	NUM <- as.double(dataRTA[, ncol(dataRTA)])
      	
	   split <- bos.f(dataRTA[,  - ncol(dataRTA)], rep(1, n), dataRTA[, ncol(dataRTA)], minbucket) 
		
		NUM <- ifelse(NUM <= split[1], 0, 1)
		indmat <- indicator(NUM)
		ncat <- ncol(indmat)
				dataRTA <- cbind(dataRTA,NUM)
     data0 <- dataRTA[, -ncol(dataRTA)]
    }
    y0rsq <- rsq.f(data0)
    predy0 <- fitted(lm(data0))
    SSdep <- ssq(dataRTA[, 1], mean(dataRTA[, 1]))
    relerr0<-1-y0rsq
    SE0<- sqrt(sum(((dataRTA[, 1] - predy0)^2 - ssqm(dataRTA[,
        1], predy0))^2))/SSdep
       
    y1rsq <- rsq.f(dataRTA)
    if (y1rsq >= 0.99) {
        stop("Linear main effects model explains all variance in y. Interaction effects are not needed.")
    }
    if (control$crit == "R2change") {
        critval <- y1rsq - y0rsq
    }
    if (control$crit == "F-value") {
        critval <- anova(lm(data0), lm(dataRTA))$"F Value"[2]
    }
    if (control$crit == "f2") {
        critval <- (y1rsq - y0rsq)/(1 - y1rsq)
    }
    relerr <- 1 - y1rsq
    predy<-fitted(lm(dataRTA))
    SE <- sqrt(sum(((dataRTA[, 1] - predy)^2 - ssqm(dataRTA[,
        1], predy))^2))/SSdep
    if (vfold == 0) {
        CV <- 1
    }
    if (CV == 0) {
        CV <- 1
    }
    if (vfold != 0) {
        if (is.null(control$cvvec)) {
            vecmat <- matrix(rep(0, n * CV), nrow = n)
            for (v in 1:CV) {
                vecmat[, v] <- vfoldvec(n, k = vfold, seed = control$seed +
                  v - 1)
            }
        }
        else {
            vecmat <- as.matrix(control$cvvec)
            if (ncol(vecmat) != CV) {
                stop("The size of the given cross-validation vector/matrix does not correspond to the value of CV")
            }
        }
    }
    PredCVarray <- array(0, dim = c(n, maxint + 1, CV))
    PredCVmat<-matrix(0,nrow=n,ncol=CV)
    CVm <- numeric(maxint + 1)
    SEm <- numeric(maxint + 1)
    relerrCV <- 0
    SECV <- 0
     relerrCV0 <- 0
    SECV0 <- 0
    ##crossvalidation of model in the rootnode
    if (vfold != 0) {
        for (v in 1:CV) {
            for (m in 1:vfold) {
                PredCVmat[vecmat[, v] == m, v] <- predict(lm(data0[vecmat[,
                  v] != m, ]), newdata = data0[vecmat[, v] ==
                  m, ])
            }
        }
        relerrCV0 <- ssq(data0[, 1], PredCVmat[,1])/SSdep
     SECV0 <- sqrt(sum(((data0[, 1] - PredCVmat[,1])^2 - ssqm(data0[,1], PredCVmat[,1]))^2))/SSdep
     
    }
    if(CV!=1){
    Resmat<-dataRTA[,1]-PredCVmat
    ssqres<-apply(Resmat^2,1,sum)
    CV0m<-(1/CV*sum(ssqres))/SSdep
    SE0m<-(1/CV*sqrt(sum((ssqres-mean(ssqres))^2)))/SSdep   }
   
  
    ##cross-validation of model after one split
    if (vfold != 0) {
        for (v in 1:CV) {
            for (m in 1:vfold) {
                PredCVarray[vecmat[, v] == m, 1, v] <- predict(lm(dataRTA[vecmat[,
                  v] != m, ]), newdata = dataRTA[vecmat[, v] ==
                  m, ])
            }
        }
    }
    result <- data.frame(matrix(0, ncat + maxint * 2, 7, byrow = TRUE),
        dup.row.names = FALSE)
    dimnames(result)[[1]] <- c(5e+05:(5e+05 + nrow(result) -
        1))
    gof <- data.frame(matrix(0, 1 + maxint, 6, byrow = TRUE))
    result[1:ncat, 1] <- names(firstvar)
    if (class1 == 1) {
        result[1:ncat, 2] <- c("=")
        result[1:ncat, 3] <- attributes(CAT)$levels[1:ncat]
    }
    else {
        result[1:ncat, 2] <- c("<=", ">")
        result[1:ncat, 3] <- split[1]
    }
    result[1:ncat, 4] <- c(apply(indmat, 2, sum))
    result[1:ncat, 6] <- rep(critval, ncat)
    for (i in 1:ncat) {
        result[i, 5] <- mean(dataRTA[indmat[, i] == 1, 1])
    }
    gof[1, 3:4] <- c(relerr, SE)
    dataCV <- dataRTA
    dataCV0 <- dataRTA
    lengthCV <- ncol(dataCV)
    k <- 1
    nodemat <- indmat
    splittingnode <- numeric(maxint)
    splittingvar <- numeric(maxint)
    rownumber <- (1:ncat) + 1
    index <- (1:ncat) + 1
    if (ncat > 2) {
        rownumber <- 1:ncat * 10000 + 1
        index <- 1:ncat * 10000 + 1
    }
    row.names(result)[1:length(index)] <- index
    indtermnodes <- c(paste("R", 1:length(rownumber), sep = ""))
    dataALL <- cbind(data0, nodemat)
    names(dataALL) <- c(names(data0), indtermnodes)
    if (sum(apply(nodemat, 2, sum) < minsplit) == ncol(nodemat)) {
        stop("The number of subjects in each category is lower than the minimum number of subjects before splitting: Decrease minbucket or merge categories.")
    }
    if(printoutput==TRUE){cat("splitting node", k, "is:  root node","\n")}
    #start splitting process:
    while (k <= maxint) {
        var.split <- splitc.f(dataRTA, data, nodemat, minsplit,
            minbucket, predtrunk, type, crit = control$crit)
        if (var.split[4] >= control$mincrit) {
            result[(ncat + k * 2 - 1):(ncat + k * 2), 1] <- names(data)[var.split[1]]
            if (type[var.split[1]] == 2) {
                Dnew <- ifelse(data[, var.split[1]] <= var.split[3] &
                  nodemat[, var.split[2]] == 1 & var.split[4] !=
                  0, 1, 0)
                result[(ncat + k * 2 - 1), 2] <- "<="
                result[(ncat + k * 2), 2] <- ">"
                result[(ncat + k * 2 - 1):(ncat + k * 2), 3] <- var.split[3]
            }
            else {
                bos <- bos.cat2(dataRTA, nodemat[, var.split[2]],
                  dataCAT[, var.split[1]])
                Dnew <- ifelse(bos$trx <= var.split[3] & nodemat[,
                  var.split[2]] == 1 & var.split[4] != 0, 1,
                  0)
                nam <- bos$ordernames
                result[(ncat + k * 2 - 1):(ncat + k * 2), 2] <- "="
                result[(ncat + k * 2 - 1), 3] <- paste(nam[c(1:length(nam)) <=
                  var.split[3]], collapse = ",")
                result[(ncat + k * 2), 3] <- paste(nam[c(1:length(nam)) >
                  var.split[3]], collapse = ",")
            }
            result[(ncat + k * 2 - 1), 4] <- sum(Dnew == 1)
            result[(ncat + k * 2), 4] <- sum(Dnew == 0 & nodemat[,
                var.split[2]] == 1)
            result[(ncat + k * 2 - 1):(ncat + k * 2), 5] <- c(mean(dataRTA[Dnew ==
                1, 1]), mean(dataRTA[Dnew == 0 & nodemat[, var.split[2]] ==
                1, 1]))
            result[(ncat + k * 2 - 1):(ncat + k * 2), 6] <- var.split[4]
            splittingnode[k] <- rownumber[var.split[2]]
            splittingvar[k] <- var.split[1]
            splitnode <- nodemat[, var.split[2]]
            if (printoutput == TRUE) {
                cat("splitting node",k+1,"is: ", paste(result[row.names(result) ==
                  splittingnode[k], 1], result[row.names(result) ==
                  splittingnode[k], 2], result[row.names(result) ==
                  splittingnode[k], 3]), "\n")
            }
            if (ncat == 2) {
                rownumber <- c(rownumber[-var.split[2]], splittingnode[k] *
                  2, splittingnode[k] * 2 + 1)
                index <- c(index, splittingnode[k] * 2, splittingnode[k] *
                  2 + 1)
            }
            else {
                rownumber <- c(rownumber[-var.split[2]], splittingnode[k] *
                  2 - floor(splittingnode[k]/10000) * 10000,
                  splittingnode[k] * 2 - floor(splittingnode[k]/10000) *
                    10000 + 1)
                index <- c(index, splittingnode[k] * 2 - floor(splittingnode[k]/10000) *
                  10000, splittingnode[k] * 2 - floor(splittingnode[k]/10000) *
                  10000 + 1)
            }
            nodemat <- cbind(nodemat[, -var.split[2]], ifelse(splitnode ==
                1 & Dnew == 1, 1, 0), ifelse(splitnode == 1 &
                Dnew == 0, 1, 0))
            dataRTA <- cbind(data0, nodemat[, -1])
            indtermnodes <- c(paste("R", 1:length(rownumber),
                sep = ""))
            names(dataRTA) <- c(names(data0), indtermnodes[-1])
            predy <- fitted(lm(dataRTA))
            relerr <- 1 - cor(dataRTA[, 1], predy)^2
            SE <- sqrt(sum(((dataRTA[, 1] - predy)^2 - ssqm(dataRTA[,
                1], predy))^2))/SSdep
            row.names(result)[1:length(index)] <- index
            gof[(1 + k), 3:4] <- c(relerr, SE)
            k <- k + 1
        }
        dataALL <- cbind(data0, nodemat)
        names(dataALL) <- c(names(data0), indtermnodes)
        if (sum(apply(nodemat, 2, sum) < minsplit) == ncol(nodemat) ||
            var.split[4] < control$mincrit) {
            k <- maxint + 1
        }
    }
    dataRTA <- cbind(data0, nodemat[, -ref])
    indtermnodes <- c(paste("R", 1:length(rownumber), sep = ""))
    names(dataRTA) <- c(names(data0), indtermnodes[-ref])
    if (!is.null(predsel)) {
        sel <- "manual"
    }
    if (sel != "none") {
        gofsel <- numeric(4)
        if (sel == "manual") {
            varsel <- c(1, predsel)
            datasel <- dataALL[, varsel]
            sel.lm <- lm(datasel)
        }
        if (sel == "backward") {
            sel1 <- lm(as.formula(lm(dataRTA)), data = dataRTA)
            stepsel <- step(sel1, k = control$ksel)
            predsel <- row.names(summary(stepsel)$coefficients)[-1]
            datasel <- cbind(dataRTA[, 1], dataRTA[, -1][, predsel])
            names(datasel)[1] <- names(dataRTA)[1]
            sel.lm <- lm(datasel)
        }
        resy2 <- (dataRTA[, 1] - fitted(sel.lm))^2
        gofsel[1] <- sum(resy2)/SSdep
        gofsel[2] <- sqrt(sum((resy2 - mean(resy2))^2))/SSdep
        predCVsel <- matrix(0, nrow = n, ncol = CV)
    }
    nodematCV <- indmat
    v <- 1
    m <- 1
    if (vfold != 0) {
        for (v in 1:CV) {
            for (m in 1:vfold) {
                k <- 1
                if (sum(apply(nodematCV[vecmat[, v] != m, ],
                  2, sum) < minsplitCV) == ncol(nodematCV[vecmat[,
                  v] != m, ])) {
                  stop("The number of subjects in each category is lower than the minimum number of subjects before splitting: Decrease minbucket or increase the number of sets defined by vfold or merge categories.")
                }
                while (k <= maxint) {
                  var.splitCV <- splitc.f(dataCV[vecmat[, v] !=
                    m, ], data[vecmat[, v] != m, ], nodematCV[vecmat[,
                    v] != m, ], minsplitCV, minbucketCV, predtrunk,
                    type, crit = control$crit)
                  if (var.splitCV[4] >= control$mincrit) {
                    if (type[var.splitCV[1]] == 2) {
                      dCV <- ifelse(data[, var.splitCV[1]] <=
                        var.splitCV[3] & nodematCV[, var.splitCV[2]] ==
                        1, 1, 0)
                    }
                    else {
                      bosCV <- bos.cat2(dataCV[vecmat[, v] !=
                        m, ], nodematCV[vecmat[, v] != m, var.splitCV[2]],
                        dataCAT[vecmat[, v] != m, var.splitCV[1]])
                      trx <- factor(dataCAT[, var.splitCV[1]],
                        levels = bosCV$ordernames)
                      trx <- as.numeric(trx)
                      trx[trx == "NA"] <- length(bosCV$ordernames) +
                        1
                      dCV <- ifelse(trx <= var.splitCV[3] & nodematCV[,
                        var.splitCV[2]] == 1, 1, 0)
                    }
                    splitnode <- nodematCV[, var.splitCV[2]]
                    nodematCV <- cbind(nodematCV[, -var.splitCV[2]],
                      ifelse(splitnode == 1 & dCV == 1, 1, 0),
                      ifelse(splitnode == 1 & dCV == 0, 1, 0))
                    dataCV <- cbind(data0, nodematCV[, -1])
                    indtermnodesCV <- c(paste("R", 1:ncol(nodematCV),
                      sep = ""))
                    names(dataCV) <- c(names(data0), indtermnodesCV[-1])
                    PredCVarray[vecmat[, v] == m, k + 1, v] <- predict(lm(dataCV[vecmat[,
                      v] != m, ]), newdata = dataCV[vecmat[,
                      v] == m, ])
                    k <- k + 1
                  }
                  if (sum(apply(nodematCV[vecmat[, v] != m, ],
                    2, sum) < minsplitCV) == ncol(nodematCV[vecmat[,
                    v] != m, ]) || var.splitCV[4] < control$mincrit) {
                    k <- maxint + 1
                  }
                  dataCVAll <- cbind(data0, nodematCV)
                  names(dataCVAll) <- c(names(data0), indtermnodesCV)
                }
                if (sel != "none") {
                  if (sel == "manual") {
                    dataCVsel <- dataCVAll[, varsel]
                    predCVsel[vecmat[, v] == m, v] <- predict(lm(dataCVsel[vecmat[,
                      v] != m, ]), newdata = dataCVsel[vecmat[,
                      v] == m, ])
                  }
                  if (sel == "backward") {
                    dataCVsel <- cbind(dataCVAll[, 1], dataCVAll[,
                      -1][, predsel])
                    names(dataCVsel)[1] <- names(dataCV)[1]
                    predCVsel[vecmat[, v] == m, v] <- predict(lm(dataCVsel[vecmat[,
                      v] != m, ]), newdata = dataCVsel[vecmat[,
                      v] == m, ])
                  }
                }
                dataCV <- dataCV0
                nodematCV <- indmat
            }
        }
    }
    predyCV <- PredCVarray[, , 1]
    
    if (vfold != 0 & maxsplit > 1) {
        for (k in 1:maxsplit) {
            relerrCV <- ssq(dataRTA[, 1], predyCV[, k])/SSdep
            SECV <- sqrt(sum(((dataRTA[, 1] - predyCV[, k])^2 -
                ssqm(dataRTA[, 1], predyCV[, k]))^2))/SSdep
            gof[k, 5:6] <- c(relerrCV, SECV)
        }
    }
    if (vfold != 0 & maxsplit == 1) {
        relerrCV <- ssq(dataRTA[, 1], predyCV)/SSdep
        SECV <- sqrt(sum(((dataRTA[, 1] - predyCV)^2 - ssqm(dataRTA[,
            1], predyCV))^2))/SSdep
        gof[1, 5:6] <- c(relerrCV, SECV)
    }
    names(result) <- c("Predictor", "Sign", "Splitpoint", "n",
        "MeanResponse", "Region", control$crit)
    names(gof) <- c("nsplit", control$crit, "RE", "SE", "REcv",
        "SEcv")
    Resarray <- dataRTA[, 1] - PredCVarray
    if (CV != 1) {
        for (k in 1:(maxint + 1)) {
            Rowres2 <- apply(Resarray[, k, ]^2, 1, sum)
            CVm[k] <- (1/CV * sum(Rowres2))/SSdep
            SEm[k] <- ((1/CV) * sqrt(sum((Rowres2 - mean(Rowres2))^2)))/SSdep
        }
        gof <- cbind(gof, REcvm = CVm, SEcvm = SEm)
    }
    if (sel != "none") {
        names(gofsel) <- c("RE", "SE", "REcv", "SEcv")
        if (vfold > 0) {
            resy2 <- (dataRTA[, 1] - predCVsel)^2
            gofsel[3] <- sum(resy2[, 1])/SSdep
            gofsel[4] <- sqrt(sum((resy2[, 1] - mean(resy2[,
                1]))^2))/SSdep
            if (CV != 1) {
                resy2 <- apply(resy2, 1, mean)
                CVselm <- sum(resy2)/SSdep
                SEselm <- sqrt(sum((resy2 - mean(resy2))^2))/SSdep
                gofsel <- c(gofsel, REcvm = CVselm, SEcvm = SEselm)
            }
        }
    }
    result <- result[result[, 4] != 0, ]
    for (i in 1:length(rownumber)) {
        result[row.names(result) == rownumber[i], 7] <- indtermnodes[i]
    }
    result[result[, 7] == 0, 7] <- ""
    result[, 6:7] <- result[, 7:6]
    rows <- seq(ncat, nrow(result), by = 2)
    nint <- 0:(length(rows) - 1)
    gof <- gof[1:length(nint), ]
    gof[, 1:2] <- cbind(nsplit = nint + 1, result[rows, 7])
    rootgof<-c(0,0,relerr0,SE0,relerrCV0,SECV0)
    
    if(CV!=1){ rootgof<-c(0,0,relerr0,SE0,relerrCV0,SECV0,CV0m,SE0m)  }
    gof<-rbind(rootgof,gof)
    rootrow <- c("root", "", "", n, mean(dataRTA[, 1]), "", "",
        "")
    result <- rbind(rootrow, result)
    fullcoef <- cbind(round(summary(lm(dataRTA))$coefficients[,
        1:2], digits = 4), round(summary(lm(stdata(dataRTA)))$coefficients[,
        1], digits = 4), round(summary(lm(dataRTA))$coefficients[,
        3:4], digits = 4))
    colnames(fullcoef) <- c("Coefficient", "Standard Error",
        " Std. Coef.", "  t value", "  Pr(>|t|)")
    if (Save == TRUE) {
        object <- list(call = match.call(), newdata = dataALL,
            trunk = result[, 1:6], splitsequence = splittingnode,
            goffull = gof, full = fullcoef)
    }
    else {
        object <- list(call = match.call(), trunk = result[,
            1:6], splitsequence = splittingnode, goffull = gof,
            full = fullcoef)
    }
    if (sel != "none") {
        selcoef <- cbind(round(summary(lm(datasel))$coefficients[,
            1:2], digits = 4), round(summary(lm(stdata(datasel)))$coefficients[,
            1], digits = 4), round(summary(lm(datasel))$coefficients[,
            3:4], digits = 4))
        colnames(selcoef) <- c("Coefficient", "Standard Error",
            " Std. Coef.", "  t value", "  Pr(>|t|)")
        if (Save == TRUE) {
            object <- list(call = match.call(), newdata = datasel,
                trunk = result[, 1:6], splitsequence = splittingnode,
                goffull = gof, full = fullcoef, gofsel = gofsel,
                selected = selcoef)
        }
        else {
            object <- list(call = match.call(), trunk = result[,
                1:6], splitsequence = splittingnode, goffull = gof,
                full = fullcoef, gofsel = gofsel, selected = selcoef)
        }
    }
    if (model=="regtrunk"){
    class(object) <- "rt"    }
        return(object)
}
