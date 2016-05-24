cv.customizedGlmnet <-
function(xTrain, yTrain, xTest, groupid = NULL, Gs = NULL,
    		dendrogram = NULL, dendrogramCV = NULL,
    		nfolds = 10, foldid = NULL,
    		family = c("gaussian", "binomial", "multinomial"))
{
    if (nrow(xTrain) != length(yTrain)) {
        stop(paste('num. of rows in xTrain (', nrow(xTrain),
        'does not match length of yTrain (', length(yTrain), ')', sep = ''))
    } else if (ncol(xTrain) != ncol(xTest)) {
        stop(paste('num. of cols of xTrain (', ncol(xTrain),
        'does not match num. of cols of xTest (', ncol(xTest), ')', sep = ''))
    } else if (!is.null(groupid) & nrow(xTest) != length(groupid)) {
        stop(paste('num. of rows of xTest (', nrow(xTest),
        'does not match length of groupid (', length(groupid), ')', sep = ''))
    }
    family = family[1]
    lambda = glmnet(xTrain, yTrain, family = family)$lambda*nrow(xTrain)

    if (family == "multinomial" | family == "binomial") {
    	yTrain = as.factor(yTrain)
    }

    if (!is.null(groupid)) {
        Gs = length(unique(groupid))
    } else {

        if (is.null(Gs)) {
            stop("Either groupid or Gs must be specified")}

        if (is.null(dendrogram)) {
    	    dendrogram = hclust(dist(rbind(xTrain, xTest)))}

        if (is.null(dendrogramCV)) {
    	    dendrogramCV = hclust(dist(xTrain))}
    }

    if (is.null(foldid)) {
    	foldid = sample(rep(1:nfolds, length.out = nrow(xTrain))) 
    }
    folds = sort(unique(foldid))

    error = matrix(NA, length(Gs), length(lambda))
    rownames(error) = Gs
    for (G in Gs) {
    	prediction = matrix(NA, length(yTrain), length(lambda))
    	for (fold in folds) {
    		xTrain_k = xTrain[foldid != fold, ]
    		yTrain_k = yTrain[foldid != fold]
    		xTest_k = xTrain[foldid == fold, ]
    		yTest_k = yTrain[foldid == fold]
            if (!is.null(groupid)) {
                groupid_k = foldid[foldid == fold]
            } else groupid_k = NULL
    		dendrogramTestIndices = array(FALSE, nrow(xTrain))
    		dendrogramTestIndices[foldid == fold] = TRUE
    		fit = customizedGlmnet(xTrain_k, yTrain_k, xTest_k, groupid_k, G,
    			family, dendrogramCV, dendrogramTestIndices)
            prediction[foldid == fold, ] = predict(fit, lambda)
    	}
        if (family == "gaussian") {
            error[as.character(G), ] =
                colMeans((prediction - yTrain)^2, na.rm = TRUE)
        } else if (family == "binomial") {
 		    error[as.character(G), ] =
                colMeans(1 + (prediction > 1/2) != yTrain, na.rm = TRUE)
        } else if (family == "multinomial") {
            error[as.character(G), ] =
                colMeans(prediction != yTrain, na.rm = TRUE)
        }
    }

    G.min = Gs[which.min(apply(error, 1, min))]
    lambda.min = lambda[which.min(error[as.character(G.min), ])]
    fit = customizedGlmnet(xTrain, yTrain, xTest, groupid, G.min,
    	dendrogram, family = family)
    selected = nonzero(fit, lambda = lambda.min)
    prediction = predict(fit, lambda = lambda.min)

    output = list(call = match.call(), prediction = prediction,
        selected = selected,
        fit = fit, G.min = G.min, lambda.min = lambda.min, error = error)
    class(output) = "cv.customizedGlmnet"
    output
}
