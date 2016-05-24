comparePredict <- function(valueName, pred, continue)
{
	stored <- c(
	  616,929,329,362,456,474,607,857,377,191,228,388,536,314,850,403,265,293,572,849,461,201,451,891,860,600,218,738,678,283,
    698,748,620,382,852,59,361,822,522,637,846,535,549,199,155,678,779,262,861,611,719,264,840,820,926,472,738,827,829,806,445,
    263,562,833,876,804,125,212,841,930,713,141,712,860,245,789,358,901,780,370,387,208,778,757,585,479,730,408,353,174,398,
    935,151,808,716,795,917,884,305,833,828,540,118,417,913,694,799,850,252,446,372,800,663,572,282,557,251,428,676,126,919,631,
    783,746,737,854,839,826,644,916,266,911,109,164,781,548,98,754,838,820,838,690,707,669,880,304,866,378,57,475,115,76,777,
    948,812,484,471,106,875,833,767,358,938,791,171,852,911,138,85,233,806,76,272,681,311,337,222,817,421,406,419,857,273,119,
    674,79,367,366,77,383,276,872,291,374,678,84,170,681,757,713)/1000

	#stored[3*(10:30)] <- 0.6 # simulate an error

  sdappr <- 0.095*(stored*(1 - stored))^0.61
	out <- pred$probabilities[, 2]
	err <- (out - stored)/sdappr
	curr <- mean(err^2)
	res <- curr < 2
    if (!res) {
		cat("Comparison FAILED for ", valueName, "\n", sep="")
        cat("mean(error^2) =", curr, "\n")
        ind <- order(err, decreasing=TRUE)[1:10]
        print(cbind(ind, stored=stored[ind], obtained=out[ind], error=err[ind]))
		if (continue) {
			cat("comparison FAILED\n")
		} else {
			stop("comparison FAILED")
		}
    }
	res
}

cmp.table <- function(a, b)
{
    aa <- unclass(a)
    bb <- unclass(b)
    if (identical(dim(aa), dim(bb))) {
        return(all(aa == bb))
    } else {
        return(FALSE)
    }
}

compareMEval <- function(valueName, pred, cl, mEval, continue)
{
    accuracy <- mean(pred$class == cl)
    res1 <- accuracy == mEval$accuracy
    aux.pred.mat <- table(cl, pred$class)
    res2 <- cmp.table(mEval$predictionMatrix, aux.pred.mat)
    res <- all(res1, res2)
    if (!res) {
		cat("Comparison FAILED for ", valueName, "\n", sep="")
        cat("accuracy =", accuracy, "\n")
        cat("mEval$accuracy =", mEval$accuracy, "\n")
        cat("aux.pred.mat\n")
        print(aux.pred.mat)
        cat("mEval$predictionMatrix\n")
        print(mEval$predictionMatrix)
		if (continue) {
			cat("comparison FAILED\n")
		} else {
			stop("comparison FAILED")
		}
    }
    res
}

compareClass <- function(valueName, value1, value2, continue)
{
    res <- identical(value1, value2)
	if (!res) {
		cat("Comparison FAILED for ", valueName, "\n", sep="")
        print(table(value1, value2))
		if (continue) {
			cat("comparison FAILED\n")
		} else {
			stop("Comparison FAILED")
		}
	}
	res
}

compareApprox <- function(valueName, value1, value2, tolerance, continue)
{
    res <- max(abs(value1 - value2)) <= tolerance
	if (!res) {
		cat("Comparison FAILED for ", valueName, "\n", sep="")
		cat(deparse(substitute(value1)), "\n")
		print(value1)
		cat("difference\n")
		print(value2 - value1)
		cat(deparse(substitute(value2)), "\n")
		print(value2)
		if (continue) {
			cat("comparison FAILED\n")
		} else {
			stop("Comparison FAILED")
		}
	}
	res
}

testCoreClass <- function(continue=TRUE)
{
	ncases <- 200
    RNGkind("Mersenne-Twister")
	set.seed(12345)
	train <- classDataGen(ncases)
	test <- classDataGen(ncases)
    model <- CoreModel(class ~ ., train, model="rf", minNodeWeightRF=5, minNodeWeightEst=2, rfNoTrees=50, maxThreads=1)
    pred <- predict(model, test, rfPredictClass=FALSE)
	destroyModels(model) 
    # consistency of predict output
    res1 <- compareClass("testCoreClass/pred$class", pred$class=="2", pred$probabilities[, 2] >= 0.5, continue)
    # compare with stored values
    res2 <- comparePredict("testCoreClass/pred", pred, continue)
    # modelEval test
    mEval <- modelEval(model, test$class, pred$class, pred$prob)
    res3 <- compareMEval("testCoreClass/modelEval", pred, test$class, mEval, continue)
    all(res1, res2, res3)
}

testCoreAttrEval <- function(continue=TRUE)
{
    ncases <- 200
    RNGkind("Mersenne-Twister")
    set.seed(0)
    train <- classDataGen(ncases)
    set.seed(0)
    estReliefF1 <- attrEval(class ~ ., train, estimator="ReliefFexpRank", maxThreads=1)
    set.seed(0)
    estReliefF0 <- attrEval(class ~ ., train, estimator="ReliefFexpRank", maxThreads=1) # makes sense to use maxThreads=2 once g++-5.0 is available
    resA1 <- compareApprox("testCoreAttrEval/estReliefF/threads", estReliefF0, estReliefF1, 1e-9, continue)
    stored <- c(0.07413109,0.0852054,0.05018575,0.02656779,0.0652197,0.03082657,-0.008773201,0.1002774,0.08263487,-0.00481844)
    resA2 <- compareApprox("testCoreAttrEval/estReliefF/stored", stored, estReliefF0, 1e-8, continue)
    set.seed(0)
    estMdl1 <- attrEval(class ~ ., train, estimator="MDL", minNodeWeightEst=5, maxThreads=1)
    set.seed(0)
    estMdl0 <- attrEval(class ~ ., train, estimator="MDL", minNodeWeightEst=5, maxThreads=1)
    resB1 <- compareApprox("testCoreAttrEval/estMdl/threads", estMdl0, estMdl1, 1e-8, continue)
    stored <- c(0.05068002,0.04873661,0.0236104,0.007925046,0.03501767,0.004991399,-0.02464427,0.1364728,0.1114686,0.00282483)
    resB2 <- compareApprox("testCoreAttrEval/estMdl/stored", stored, estMdl0, 1e-7, continue)
    all(resA1, resA2, resB1, resB2)
}

testCoreReg <- function(continue=TRUE)
{
    ncases <- 200
    RNGkind("Mersenne-Twister")
    set.seed(0)
    train <- regDataGen(ncases)
    test<- regDataGen(ncases)
    model <- CoreModel(response~., train, model="regTree", modelTypeReg=5, minNodeWeightEst=1, minNodeWeightTree=5)
    pred <- predict(model, test)
	destroyModels(model) 
    # Model evaluation
    mEval <- modelEval(model, test[["response"]], pred)
    stored <- c(0.7983146,0.7638713,0.5261048,0.7238808) # for ModelTypeReg=5, for modelTypeReg = 3: c(0.7326231, 0.7194654, 0.4705813, 0.6805264) 
    res1 <- compareApprox("testCoreReg/mEval", stored, c(mEval$MSE,mEval$RMSE,mEval$MAE,mEval$RMAE), 1e-6, continue)
    # Attribute evaluation with RReliefFexpRank
    estRReliefF <- attrEval(response~., train, estimator="RReliefFexpRank")
    stored <- c(0.03235394,0.02799911,0.004161262,-0.05577854,-0.04937824,0.05685567,-0.03946655,0.001630726,0.05570145,0.1200363)
    res2 <- compareApprox("testCoreReg/estReliefF", stored, estRReliefF, 1e-6, continue)
    # Attribute evaluation with MSEofMean
    stored <- c(-0.7192786,-0.686779,-0.7668486,-0.7625726,-0.7373601,-0.6837401,-0.753099,-0.7328526,-0.6412497,-0.6879387)
    estMSE <- attrEval(response~., train, estimator="MSEofMean")
    res3 <- compareApprox("testCoreReg/estMSE", stored, estMSE, 1e-6, continue)
    all(res1, res2, res3)
}


testCoreOrdEval <- function(continue=TRUE)
{
    ncases <- 200
    RNGkind("Mersenne-Twister")
    set.seed(0)
    train <-  ordDataGen(ncases)
    estOrdEval <- ordEval(class~., train, ordEvalNoRandomNormalizers=0 )
    stored <- c(0.5385996,0.6631206,0.344894,0.4327273,0.3623188,0.4223827,0.3440285,0.3432203)
    res1 <- compareApprox("testCoreOrdEval/reinfPosAttr", stored, estOrdEval$reinfPosAttr, 1e-7, continue)
    stored <- c(0.4653641,0.581854,0.3009524,0.3561151,0.2833333,0.3364486,0.2515213,0.2696177)
    res2 <- compareApprox("testCoreOrdEval/reinfNegAttr", stored, estOrdEval$reinfNegAttr, 1e-7, continue)
    stored <- c(0.4693182,0.525296,0.376569,0.4295302,0.3568282,0.4039517,0.3794926,0.4151309)
    res3 <- compareApprox("testCoreOrdEval/anchorAttr", stored, estOrdEval$anchorAttr, 1e-7, continue)
    all(res1, res2, res3)
}

#gener.Reg <- function(m,n)
#{
#    x <- matrix(runif(m*n),nrow=m,ncol=n)
#    data.frame(x,resp=rowSums(x))
#}

outputResult <- function(testName, status, failMessage, continue)
{
    if (!all(status)) {
        outMessage <- testName
		if (failMessage != "") {
			outMessage <- paste(outMessage, " (", failMessage, ")", sep="")
		}
		if (continue) {
			cat("Test FAILED: ", outMessage, "\n", sep="")
		} else {
			stop("Test FAILED: ", outMessage)
		}
	}
}

singleTestNA <- function(t, x)
{
    .C("testNA",
    as.integer(t),
    as.double(x),
    out=integer(2),
    NAOK=TRUE,
    PACKAGE="CORElearn")$out
}

testCoreNA <- function(continue=TRUE)
{
    a <- matrix(nrow=4, ncol=2)
    a[1,] <- singleTestNA(0, NA) # pass NA to CORElearn
    a[2,] <- singleTestNA(0, NaN) # pass NaN to CORElearn
    a[3,] <- singleTestNA(1, 0) # use internal NAcont
    a[4,] <- singleTestNA(2, 0) # generate NaN
    ok <- a == rbind(c(1,0), c(0,1), c(1,0), c(0,1))
	outputResult("testCoreNA", ok, "", continue)
	all(ok)
}

testCoreRPORT <- function(continue=TRUE)
{
    tmp <- .C("testRPORT", a=as.integer(2), PACKAGE="CORElearn")
    ok <- tmp$a == 1
	outputResult("testCoreRPORT", ok, tmp$a, continue)
	all(ok)
}

testCoreRand <- function(continue=TRUE)
{
    n <- 10
    runif(1)
    state <- .Random.seed
    x <- runif(n)
    .Random.seed <<- state
    y <- .C("testCoreRand", as.integer(n), a=double(n), PACKAGE="CORElearn")$a
	ok <- x == y
	outputResult("testCoreRand", ok, paste(x[1], x[2], y[1], y[2]), continue)
	all(ok)
}

asTxt <- function(ok)
{
	if (all(ok)) {
		"OK"
	} else {
		"FAIL"
	}
}

allTests <- function(continue=TRUE, timed=FALSE)
{
	
    t1 <- system.time(r1 <- testCoreClass(continue))
    cat("testCoreClass()    : ", asTxt(r1), "\n")
    if (timed) cat("Elapsed", t1["elapsed"],"sec\n")
    t2 <- system.time(r2 <- testCoreAttrEval(continue))
    cat("testCoreAttrEval() : ", asTxt(r2), "\n")
    if (timed) cat("Elapsed", t2["elapsed"],"sec\n")
    t3 <- system.time(r3 <- testCoreReg(continue))
    cat("testCoreReg()      : ", asTxt(r3), "\n")
    if (timed) cat("Elapsed", t3["elapsed"],"sec\n")
    t4 <- system.time(r4 <- testCoreOrdEval(continue))
    cat("testCoreOrdEval()  : ", asTxt(r4), "\n")
    if (timed) cat("Elapsed", t4["elapsed"],"sec\n")
    if (timed) {
        cat("system.time() summary\n")
        print(rbind(t1, t2, t3, t4)[, 1:3])
    }
    r5 <- testCoreNA(continue)
    cat("testCoreNA()       : ", asTxt(r5), "\n")
    r6 <- testCoreRPORT(continue)
    cat("testCoreRPORT()    : ", asTxt(r6), "\n")
    r7 <- testCoreRand(continue)
    cat("testCoreRand()     : ", asTxt(r7), "\n")
	result <- all(r1, r2, r3, r4, r5, r6, r7)
	outputResult("allTests", result, "", continue=FALSE)
	invisible(result)
}

testClassPseudoRandom <- function(s, k, m)
{
	n <- length(s)
	aux <- .C("testClassPseudoRandom",
		n = as.integer(n),
		s = as.integer(s),
		k = as.integer(k),
		m = as.integer(m),
		x = double(k*m),
		PACKAGE="CORElearn")
	matrix(aux$x, nrow=k, ncol=m)
}

testTime <- function()
{
	.C("testTime", x=double(1), PACKAGE="CORElearn")$x
}

