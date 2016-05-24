saemodel <-
function(formula, area, data, type="b", na.omit=FALSE){
   # get all variables in the formula object
   variables <- all.vars(formula)
   if (!is.language(area)) stop("area-specific random effect must be defined as formula \n")
   # Rao's model type "b" (basic unit-level model)
   if (type == "b" | type == "B"){
      type <- "B"
      # get the variable(s) that define the area-specifc random effect
      areavariable <- all.vars(area)
      if (!is.na(match(areavariable, variables))) stop("NOTE: '", areavariable, "' can be either in the formula or defining the area-specific random effect, but NOT both \n")
      # limit the number of area-defining variables to one
      if (length(areavariable) !=1 ) stop("area-specific random effect not properly specified \n")
      # generate the model frame (in doing so, we can specify the behavior in the case of NAs)
      na.action <- ifelse (na.omit == TRUE, "na.omit", "na.fail")
      mf <- model.frame(formula, data, na.action=na.action)
      # extract the terms
      mt <- attr(mf, "terms")
      # register if it has an intercept
      hasintercept <- attr(mt, "intercept")
      # extract the response
      y <- as.numeric(model.response(mf, "numeric"))
      # extract the design matrix
      X <- model.matrix(mt, mf)
      # dim of X
      p <- dim(X)[2]
      # check if X has full rank, else: stop
      qr <- qr(X)
      if (qr$rank < p) stop("Rank(design matrix) < p! Choose another model!\n")
      # generate area identifiers and sort the data 
      areaID <- data[, areavariable]
      # reshape areaID into a factor
      if (!is.factor(areaID)){
	 areaID <- as.factor(areaID)
      }
      # order the data along areaID so that within-area units form blocks
      ord <- order(areaID)
      mod <- data.frame(y, X, areaID)
      mod <- mod[ord, ]
      # get the area names
      areaNames <- levels(areaID)
      # vector of area size
      nsize <- unname(table(areaID))
      model <- list(X=mod[, 2:(p+1)], y=mod[, 1], nsize=nsize, areaID=mod[, (p+2)], g=length(nsize), p=p, n=length(y), intercept=hasintercept)
      attr(model, "areaNames") <- areaNames
      # additional stuff
      pf <- paste(formula)
      areadef <- paste(area[[2]])
      xnames <- colnames(X)
      yname <- pf[[2]]
   }
   # Rao's model type a (Fay-Herriot model, known variances must be given as argument of area)
   if (type == "a" | type == "A"){
      type <- "A"
      stop("not inplemented, yet! \n")
   }
   # build the model
   attr(model, "yname") <- yname 
   attr(model, "xnames") <- xnames
   attr(model, "areadef") <- areadef
   attr(model, "call") <- match.call()
   class(model) <- "saemodel"
   return(model)
}

