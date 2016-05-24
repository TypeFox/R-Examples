crunch <- function(formula, data, phy, names.col, stand.contr = TRUE, robust=Inf, ref.var=NULL, node.depth=NULL,
                  polytomy.brlen=0, equal.branch.length=FALSE, factor.action="abort")
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
            if(is.null(phy$edge.length)) stop("The phylogeny does not contain branch lengths and crunch has not been set to use equal branch lengths.")
            if(any(phy$edge.length <= 0)) stop("The phylogeny contains either negative or zero branch lengths and crunch has not been set to use equal branch lengths.")
        }
            
        # check for factor.action
        factor.action <- match.arg(factor.action, c("abort", "warn", "allow"))

        # useful info...
        root <- length(phy$tip.label) + 1
        unionData <- nrow(data)

	# CALCULATE MODEL 
	# GET THE MODEL MATRIX and Model Response

	    # reduce to just the variables used in the formula so they can all be treated as numeric 
	    # but hang on to tip labels for subsetting the phylogeny down to complete tips
	    data <- subset(data, select=all.vars(formula))

    # HANDLE CATEGORICAL VARIABLES:
        # find the factors - (number of levels > 0)
        varLevels <- sapply(data, function(x) length(levels(x)))
        varIsOrdered <- sapply(data, is.ordered)
        
        if(any( varLevels > 0 )){
            
            # check for unordered multi states...
            if(any( varLevels > 2 & ! varIsOrdered )) stop("Unordered non-binary factors included in model formula.")
            
            # otherwise check for action on viable factors...
            
            if(factor.action == "abort"){
                stop("The formula includes factors. Change the factor.action argument to allow these to be fitted as numeric variables.")
            } else if(factor.action == "warn"){
                warning("The formula includes factors, which have been treated as continuous variables ")
            }
            
            data <- as.data.frame(lapply(data, as.numeric))
        }
                
        
        # ditch the intercept, if present
        formula <- update(formula, . ~ . - 1)
        
        # now we have the union of the phylogeny and data
        # get the model frame, matrix and response
        # these show the values at the tips for each term
        mf <- model.frame(formula, data, na.action=na.pass)

        # is there enough data in the model
        # TODO - think whether this check is always sufficient. 
        mfComplete <- complete.cases(mf)
        if(sum(mfComplete) < 2 ) stop("Fewer than two taxa contain complete data for this analysis")

        # get the design matrix
        md <- model.matrix(formula, mf) 

        # sort out the reference variable
        if( is.null(ref.var)){
            ref.var <- colnames(md)[1] # first column in the design matrix
        } else {
            ref.var <- deparse(substitute(ref.var))
            if(is.na(match(ref.var, colnames(md)))) stop("Reference variable not found in design matrix")
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
        # We know the tip values and have the tree         
        contr <- contrCalc(md, phy, ref.var, "crunch", polytomy.brlen)

    # GET RESPONSE MATRIX
        # standardize the contrasts if required     
        if(stand.contr)  contr$contr <- contr$contr/sqrt(contr$var.contr)
      
    # FEED THE RETURNED DESIGN AND RESPONSE MATRICES INTO THE MODELLING FUNCTIONS
        # assemble the data into a finished contrast object
        ContrObj <- list()
        ContrObj$contr$response <- contr$contr[,1,drop=FALSE]

		ContrObj$contr$explanatory <- contr$contr[,-1,drop=FALSE]

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

		contrMD <-  ContrObj$contr$explanatory[validNodes,,drop=FALSE]
		contrRS <-  ContrObj$contr$response[validNodes,,drop=FALSE]
		attr(contrMD, 'assign') <- attr(ContrObj, 'assign') ## replace attributes 
		if(! is.null(attr(ContrObj, 'contrasts'))) attr(contrMD, 'contrasts') <- attr(ContrObj, 'contrasts')
		
       mod <- lm.fit(contrMD, contrRS)
       class(mod) <-  "lm"
       
       # assemble the output
       # return fitted model and contrasts
       RET <- list(contrast.data=ContrObj, mod=mod, data=cdata)
       class(RET) <- c("caic")
       
	   ## THIS CAN BE DONE WITH MUCH MORE FINESSE!
       # convert the ContrObj into a data frame...
	   # - removed the OTT call to caic.table
	   # remove invalid Nodes to keep the residual and predict lengths the same
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
       attr(RET, "contr.method") <- "crunch"
       attr(RET, "macro.method") <- ""
       attr(RET, "stand.contr") <- stand.contr
       attr(RET, "robust") <- robust
       
	   # lastly, test for studentised outliers
	   if(any(stRes > robust)){
			RET <- caic.robust(RET, robust)
	   }

       return(RET)

}

