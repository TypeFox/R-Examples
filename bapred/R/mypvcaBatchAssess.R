mypvcaBatchAssess <-
function (X, factordata, threshold)
{
    ##require("lme4")

    X <- t(X)
    ###### X <- exprs(vsn2(abatch, verbose = FALSE))
    dataRowN <- nrow(X)
    dataColN <- ncol(X)

	########## Center the data (center rows) ##########
    XCentered <- matrix(data = 0, nrow = dataRowN, 
        ncol = dataColN)
    XCentered_transposed = apply(X, 1, 
        scale, center = TRUE, scale = FALSE)
    XCentered = t(XCentered_transposed)

	########## Compute correlation matrix &  Obtain eigenvalues ##########

    theDataCor <- cor(XCentered)
    eigenData <- eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n <- length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues/eigenValuesSum

	##===========================================
	##	Getting the experimental information
	##===========================================
    ###### expInfo <- pData(abatch)[, batch.factors]
    ###### factordata <- as.data.frame(expInfo)
    expDesignRowN <- nrow(factordata)
    expDesignColN <- ncol(factordata)

	########## Merge experimental file and eigenvectors for n components ##########

    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1) {
        my_sum_2 = my_sum_2 - percents_PCs[i]
        if ((my_sum_2) <= threshold) {
            my_counter_2 = my_counter_2 + 1
        }
    }
    if (my_counter_2 < 3) {
        pc_n = 3
    }
    else {
        pc_n = my_counter_2
    }

	## pc_n is the number of principal components to model

    pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
        pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n) {
        for (j in 1:expDesignRowN) {
            mycounter <- mycounter + 1
            pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                i]
        }
    }
    AAA <- factordata[rep(1:expDesignRowN, pc_n), ]
    Data <- cbind(AAA, pc_data_matrix)


	####### Edit these variables according to your factors #######

    variables <- c(colnames(factordata))
    for (i in 1:length(variables)) {
        Data$variables[i] <- as.factor(Data$variables[i])
    }


	########## Mixed linear model ##########
    op <- options(warn = (-1))
    effects_n = expDesignColN + choose(expDesignColN, 2) + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)

	##============================#
	##	Get model functions
	##============================#
    model.func <- c()
    index <- 1

	##	level-1
    for (i in 1:length(variables)) {
        mod = paste("(1|", variables[i], ")", sep = "")
        model.func[index] = mod
        index = index + 1
    }

	##	two-way interaction
    for (i in 1:(length(variables) - 1)) {
        for (j in (i + 1):length(variables)) {
            mod = paste("(1|", variables[i], ":", variables[j], 
                ")", sep = "")
            model.func[index] = mod
            index = index + 1
        }
    }
    function.mods <- paste(model.func, collapse = " + ")

	##============================#
	##	Get random effects	#
	##============================#

    for (i in 1:pc_n) {
        y = (((i - 1) * expDesignRowN) + 1)
        funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
        Rm1ML <- lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
            expDesignRowN), ], REML = TRUE, verbose = FALSE, 
            na.action = na.omit)
        randomEffects <- Rm1ML
		
		# From R-version 3.3.0 on lme4 no longer exports the function sigma().
		# Therefore we give a definition of the function here in the 
		# environment of mypvcaBatchAssess:
		sigmafromlme4 <- function (object, ...)
        {
          dc <- object@devcomp
          dd <- dc$dims
          if (dd[["useSc"]]) 
            dc$cmp[[if (dd[["REML"]]) 
              "sigmaREML"
            else "sigmaML"]]
          else 1
        }
		
        randomEffectsMatrix[i, ] <- c(unlist(lme4::VarCorr(Rm1ML)), 
            resid = sigmafromlme4(Rm1ML)^2)
    }
    effectsNames <- c(names(lme4::getME(Rm1ML, "cnms")), "resid")

	########## Standardize Variance ##########
    randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
        ncol = effects_n)
    for (i in 1:pc_n) {
        mySum = sum(randomEffectsMatrix[i, ])
        for (j in 1:effects_n) {
            randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                j]/mySum
        }
    }

	########## Compute Weighted Proportions ##########

    randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
        ncol = effects_n)
    for (i in 1:pc_n) {
        weight = eigenValues[i]/eigenValuesSum
        for (j in 1:effects_n) {
            randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                j] * weight
        }
    }

	########## Compute Weighted Ave Proportions ##########

    randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
    randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
    totalSum = sum(randomEffectsSums)
    randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
        ncol = effects_n)
    for (j in 1:effects_n) {
        randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
    }
    return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}
