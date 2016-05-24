brunch <- function(formula, data, phy, names.col, stand.contr = TRUE, robust=Inf, ref.var=NULL, 
	               node.depth=NULL, equal.branch.length=FALSE)
{


    # Program Flow:
    #   1) setup - check arguments, 
    #   2) use model functions to get design and response matrices, including all NA data
    #   3) feed the model matrices into a function to calculate nodal values and contrasts
    #   4) feed the returned contrast versions of the design and response matrices into lm.fit
    
    # TODO - return node age/height
    # TODO - allow caic to be used as a contrast calculator
    # TODO - explicit check for polytomy.brlen problems
    
    # CHECKS AND SETUP
    
	    # - test to see if there is a comparative data object and if not then
	    #   retrofit the remaining arguments into a comparative data object.
		if(! missing(data)){
			if(! inherits(data, 'comparative.data')){
				if(missing(names.col)) stop('names column is missing')
				names.col <- deparse(substitute(names.col))
				data <- caicStyleArgs(data=data, phy=phy, names.col=names.col, warn.dropped=TRUE)
			}
		}
	
		# extract the data and phylogeny
		cdata <- data # in case the original is needed later
		phy <- data$phy
		data <- data$data
	
	    # check node.depth is sensible
	    if(! is.null(node.depth)){
	        if(node.depth%%1 != 0 || node.depth < 1) stop("node.depth must be a positive integer greater than 1.")
	    }
    
		# set branch lengths doesn't get evaluated if FALSE or zero
        if(as.logical(equal.branch.length)) {
            phy$edge.length <- rep(2, nrow(phy$edge))
        } else {
            if(is.null(phy$edge.length)) stop("The phylogeny does not contain branch lengths and brunch has not been set to use equal branch lengths.")
            if(any(phy$edge.length <= 0)) stop("The phylogeny contains either negative or zero branch lengths and brunch has not been set to use equal branch lengths.")
        }
  		
	    # useful info...
	    root <- length(phy$tip.label) + 1
	    unionData <- nrow(data)
       
    # CALCULATE MODEL 
    # GET THE MODEL MATRIX and Model Response
        
        # ditch the intercept, if present
        formula <- update(formula, . ~ . - 1)
        
        # now we have the union of the phylogeny and data
        # get the model frame, matrix and response
        # these show the values at the tips for each term 
        mf <- model.frame(formula, data, na.action=na.pass)
        
        # TODO - think whether this check is always sufficient...
        mfComplete <- complete.cases(mf)
        if(sum(mfComplete) < 2 ) stop("Fewer than two taxa contain complete data for this analysis")


        # HANDLE CATEGORICAL VARIABLES:
         
        # currently, want to exclude the possibility of interactions in a factor
        # in Brunch because I don't have a clue how that should be handled
        # and also need a vector showing which terms are categorical and numeric
        # in order to allow correct standardization of contrasts
        # - do this step before turning the factors into numbers

        varClass <- attributes(attributes(mf)$terms)$dataClasses
        termFactors <- attributes(attributes(mf)$terms)$factors
       
        if(any(varClass %in% c("ordered","factor") & rowSums(termFactors) > 1)){
                stop("Interactions using categorical variables not supported in brunch analyses")}

        termClass <- apply(termFactors,2,function(X) unique(varClass[as.logical(X)]))

        # now modify the variables to be numeric for calculation 
        varLevels <- sapply(mf, function(x) length(levels(x)))
        varIsOrdered <- sapply(mf, is.ordered)
        
        # check for unordered multi states...
        if(any( varLevels > 2 & ! varIsOrdered )) stop("Unordered non-binary factors included in model formula.")
        
        # refit the model frame with numericized data
        data <- as.data.frame(lapply(mf, as.numeric))
        mf <- model.frame(formula, data, na.action=na.pass)

        # get the design matrix
        md <- model.matrix(formula, mf) 

        # sort out the reference variable and check for multiple or no factors...
        if(sum(varLevels > 0) == 0) stop("No factors specified in the model formula")        
        if(sum(varLevels > 0) > 1) warning("Multiple factors specified in the model formula")
        
        if( is.null(ref.var)){
            ref.var <- colnames(data)[which(varLevels > 0)[1]] # the (first) factor in the design matrix
        } else {
            ref.var <- deparse(substitute(ref.var))
            if(is.na(match(ref.var, colnames(md)))) stop("Reference variable not found in design matrix")
            if(! ref.var %in% colnames(data)[which(varLevels == 0)]) stop("The reference variable is not a factor")
        }

        # MODEL RESPONSE
        mr <- model.response(mf)
        # turn into a column matrix
        mr <- as.matrix(mr)
        colnames(mr) <- as.character(formula[2])
        # now that we have the model response for CAIC style contrasts we can substitute the reference variable
        # for empty models (i.e. models specified as X~1)
        if(is.empty.model(formula)) ref.var <- colnames(mr)

        # add to the design matrix - this strips the assign and contrast attributes so save...
		attrMD <- attributes(md)
        md <- cbind(mr, md)

    # NOW SETUP TO GET CONTRASTS AND NODAL VALUES
        # We know the tip values, the analysis tree         
        contr <- contrCalc(md, phy, ref.var, "brunch", 0) # brunch used an internal branch length of 0

    # GET RESPONSE MATRIX
        # first column of contrasts is response
        mrC <- contr$contr[,1,drop=FALSE]
        mdC <- contr$contr[,-1,drop=FALSE]

    # standardize the contrasts (but not the categorical) if required 
    
        if(stand.contr){
            notCateg <- ! termClass %in% c("factor","ordered")
            mdC[,notCateg] <- mdC[,notCateg, drop=FALSE]/sqrt(contr$var.contr)
            mrC <- mrC/sqrt(contr$var.contr)
        }
    
    # FEED THE RETURNED DESIGN AND RESPONSE MATRICES INTO THE MODELLING FUNCTIONS
        # assemble the data into a finished contrast object

        ContrObj <- list()
        ContrObj$contr$response <- mrC
        ContrObj$contr$explanatory <- mdC
        ContrObj$nodalVals$response <- contr$nodVal[,1,drop=FALSE]
        ContrObj$nodalVals$explanatory <- contr$nodVal[,-1,drop=FALSE]            
        ContrObj$contrVar <- contr$var.contr
        ContrObj$nChild <- contr$nChild
        ContrObj$nodeDepth <- contr$nodeDepth

		## need to keep the assign and contrasts attributes from the model 
		## matrix with the contrast object in order to get anova() methods to work
		## can't store assign permanently with explanatory contrasts because validNode subsetting strips attributes
		attr(ContrObj, 'assign') <- attrMD$assign
		if(! is.null(attrMD$contrasts)) attr(ContrObj, 'contrasts') <- attrMD$contrasts
		
        # gather the row ids of NA nodes to drop from the model
        validNodes <- with(ContrObj$contr, complete.cases(explanatory) & complete.cases(response))
       
        # enforce any node depth requirement
        if(! is.null(node.depth)){
            validNodes[ContrObj$nodeDepth > node.depth] <- FALSE
        }
                 
        # save for the user
        ContrObj$validNodes <- validNodes
        
        # feed the contr.model.response and contr.model.matrix
        # into lm.fit to get the model and then set up the lm object
        # - need to use lm.fit here rather than calling the model on 
        #   data=contrData because any functions in the formula are now
        #   set in the column names - don't want lm to try and reinterpret
        #   them in parsing the formula.
        # - the problem then becomes how to get the model to refer to the dataset

		## need to pass the assign and contrasts attributes over from the model 
		## matrix in order to get anova() methods to work
		contrMD <-  ContrObj$contr$explanatory[validNodes,,drop=FALSE]
		contrRS <-  ContrObj$contr$response[validNodes,,drop=FALSE]

		attr(contrMD, 'assign') <- attr(ContrObj, 'assign') ## replace attributes
		if(! is.null(attr(ContrObj, 'contrasts'))) attr(contrMD, 'contrasts') <- attr(ContrObj, 'contrasts')

       mod <- with(ContrObj$contr, lm.fit(contrMD, contrRS))
       class(mod) <-  "lm"
       
       # assemble the output
       # return fitted model and contrasts
       RET <- list(contrast.data=ContrObj, mod=mod,data=cdata)
       class(RET) <- c("caic")
       
        # convert the ContrObj into a data frame...
       contrData <- with(ContrObj$contr, as.data.frame(cbind(response,explanatory)))
       contrData <- contrData[validNodes, ,drop=FALSE]
       RET$mod$call <- substitute(lm(FORM, data=contrData), list(FORM=formula))
       RET$mod$terms <- attr(mf, "terms")
       
       # put the model.frame in to the lm object so that predict, etc. calls work
       RET$mod$model <- contrData
       attr(RET$mod$model, "terms") <- attr(mf, "terms")
 
	   ## Add studentized residuals: need to use matching in case of invalid nodes
       stRes <- rstudent(mod)
       SRallNodes <- rep(NA, length(RET$contrast.data$validNodes))
       names(SRallNodes) <- names(RET$contrast.data$contrVar)
	   SRallNodes[match(names(stRes), names(SRallNodes))] <- stRes
       RET$contrast.data$studentResid <- SRallNodes

       ## add some attributes
       attr(RET, "contr.method") <- "brunch"
       attr(RET, "macro.method") <- ""
       attr(RET, "stand.contr") <- stand.contr
       attr(RET, "robust") <- robust

	   # lastly, test for studentised outliers
	   if(any(stRes > robust)){
			RET <- caic.robust(RET, robust)
	   }

       return(RET)

}