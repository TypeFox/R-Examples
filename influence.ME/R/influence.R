`influence` <-
  function(model, group=NULL, select=NULL, obs=FALSE, gf="single", count = FALSE, delete=TRUE, ...)
  {

	fixef <- NA
	rm(fixef)

    ## Checks, errors, and warnings
    # obs=TRUE cannot be used with delete=FALSE, group, gf,  parameters
    
    if(is.null(group) & !obs)
    {
      stop("Please specify either the 'group' parameter, or specify 'obs=TRUE'")
    }    
    
    if(!is.null(group) & obs)
    {
      stop("Either specify the 'group' parameter, or specify 'obs=TRUE', but not both.")
    }    
    
    
    # Defining Internal Variables
    data.adapted <- model.frame(model)
    original.no.estex <- which(substr(names(fixef(model)), 1,6) != "estex.")
    n.pred <- length(fixef(model)[original.no.estex])
    
    
    ####
	# Code kindly provided by Jennifer Bufford
	if("(weights)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(weights)"] <-
    as.character(model@call$weights)}
    if("(offset)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(offset)"] <-
    as.character(model@call$offset)}
    if(sum(grepl("offset", names(data.adapted)))>0) {
    names(data.adapted)[grep("offset", names(data.adapted))] <-
    gsub('offset\\(|\\)',"",names(data.adapted)[grep("offset", names(data.adapted))])}
	####
    
    
    
    if(!obs)
    {
      grouping.names <- grouping.levels(model, group)
      n.groups <- length(grouping.names)
    }
    
    if(obs)
    {
      n.obs <- nrow(data.adapted)
    }
    
    ###
    # Defining and naming the output elements for the original models
    ###
    # These apply to all
    ###
    
    # Fixed Estimates of the original model
    or.fixed <- matrix(ncol = n.pred , nrow = 1, data = fixef(model)[original.no.estex])
    dimnames(or.fixed) <- list(NULL, names(fixef(model))[original.no.estex])
    
    # Standard Error of the original model
    or.se <- matrix(ncol = n.pred , nrow = 1, data = se.fixef(model)[original.no.estex])
    dimnames(or.se) <- list(NULL, names(fixef(model))[original.no.estex])
    
    # Variance / Covariance Matrix of the original model
    or.vcov <- as.matrix(vcov(model)[original.no.estex, original.no.estex])
    dimnames(or.vcov) <- list(
      names(fixef(model)[original.no.estex]), 
      names(fixef(model)[original.no.estex]))
    
    # Test statistic of the original model
    or.test <- coef(summary(model))[original.no.estex,3]
    

    ###
    # Defining and naming the output elements for the adapted models
    ###
    # These procedures vary
    ###
    
    
    if(!obs)
    {
      
      if(is.null(select))
      {
        
        # Fixed Estimates of the modified model(s)
        alt.fixed <- matrix(ncol = n.pred, nrow = n.groups, data = NA)
        dimnames(alt.fixed) <- list(grouping.names, names(fixef(model))[original.no.estex])
        
        # Standard Error of the modified model(s)
        alt.se <- matrix(ncol = n.pred , nrow = n.groups, data = NA)
        dimnames(alt.se) <- list(grouping.names, names(fixef(model))[original.no.estex])  
        
        # Variance / Covariance Matrix of the modified model(s)
        alt.vcov <- list()  
        
        # Test statistic of the modified model(s)
        alt.test <- matrix(ncol = n.pred , nrow = n.groups, data = NA)
        dimnames(alt.test) <- list(grouping.names, names(fixef(model))[original.no.estex])  
        
        for (i in 1:n.groups)
        {
          
          if(count == TRUE) {print(n.groups + 1 - i)}

          model.updated <- exclude.influence(model=model, grouping=group, level=grouping.names[i], gf=gf, delete=delete)

          altered.no.estex <- which(substr(names(fixef(model.updated)), 1,6) != "estex.")
          
          alt.fixed[i,] 	<- as.matrix(fixef(model.updated)[altered.no.estex])
          alt.se[i,] 		<- as.matrix(se.fixef(model.updated)[altered.no.estex])
          alt.vcov[[i]] 	<- as.matrix(vcov(model.updated)[altered.no.estex, altered.no.estex])
          alt.test[i,] 		<- as.matrix(coef(summary(model.updated))[,3][altered.no.estex])
        }
        
       
       
      }
      
      if(!is.null(select))
      {  
        
        model.updated <- exclude.influence(model, group, select, gf=gf, delete=delete)
        altered.no.estex <- which(substr(names(fixef(model.updated)), 1,6) != "estex.")
        
        # Fixed Estimates of the modified model(s)
        alt.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model.updated)[altered.no.estex])
        dimnames(alt.fixed) <- list(
          "Altered model", 
          names(fixef(model.updated))[altered.no.estex])
        
        # Standard Error of the modified model(s)
        alt.se <- matrix(ncol = n.pred , nrow = 1, data = se.fixef(model.updated)[altered.no.estex])
        dimnames(alt.se) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
        
        # Variance / Covariance Matrix of the modified model(s)
        alt.vcov <- list()
        alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, altered.no.estex])
        dimnames(alt.vcov[[1]]) <- list(
          names(fixef(model.updated)[altered.no.estex]), 
          names(fixef(model.updated)[altered.no.estex]))
          
        # Test statistic of the modified model(s)
        alt.test <- matrix(ncol = n.pred , nrow = 1, data = coef(summary(model.updated))[,3][altered.no.estex])
        dimnames(alt.test) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
       
      }
    }
    
    
    if(obs)
    {
      
      if(is.null(select))
      {
        # Fixed Estimates of the modified model(s)
        alt.fixed <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
        dimnames(alt.fixed) <- list(1:n.obs, names(fixef(model))[original.no.estex])
        
        # Standard Error of the modified model(s)
        alt.se <- matrix(ncol = n.pred , nrow = n.obs, data = NA)
        dimnames(alt.se) <- list(1:n.obs, names(fixef(model))[original.no.estex])  
        
        # Variance / Covariance Matrix of the modified model(s)
        alt.vcov <- list()
        
        # Test statistic of the modified model(s)
        alt.test <- matrix(ncol = n.pred , nrow = n.obs, data = NA)
        dimnames(alt.test) <- list(1:n.obs, names(fixef(model))[original.no.estex])  
        
        for (i in 1:n.obs)
        {
          
          if(count == TRUE) {print(n.obs + 1 - i)}
          
          model.updated <- exclude.influence(model, obs=i)
          altered.no.estex <- which(substr(names(fixef(model.updated)), 1,6) != "estex.")
          
          alt.fixed[i,] 	<- as.matrix(fixef(model.updated)[altered.no.estex])
          alt.se[i,] 		<- as.matrix(se.fixef(model.updated)[altered.no.estex])
          alt.vcov[[i]] 	<- as.matrix(vcov(model.updated)[altered.no.estex, altered.no.estex])
          alt.test[i,]	 	<- as.matrix(coef(summary(model.updated))[,3][altered.no.estex])
        } 
        
        
      }
      
      if(!is.null(select))
      {
        model.updated <- exclude.influence(model, obs=select)
        altered.no.estex <- which(substr(names(fixef(model.updated)), 1,6) != "estex.")
        
        # Fixed Estimates of the modified model(s)
        alt.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model.updated)[altered.no.estex])
        dimnames(alt.fixed) <- list(
          "Altered model", 
          names(fixef(model.updated))[altered.no.estex])
        
        # Standard Error of the modified model(s)
        alt.se <- matrix(ncol = n.pred , nrow = 1, data = se.fixef(model.updated)[altered.no.estex])
        dimnames(alt.se) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
        
        # Variance / Covariance Matrix of the modified model(s)
        alt.vcov <- list()
        alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, altered.no.estex])
        dimnames(alt.vcov[[1]]) <- list(
          names(fixef(model.updated)[altered.no.estex]), 
          names(fixef(model.updated)[altered.no.estex]))
          
        # Test statistic of the modified model(s)
        alt.test <- matrix(ncol = n.pred , nrow = 1, data = coef(summary(model.updated))[,3][altered.no.estex])
        dimnames(alt.test) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
  
      }
      
      
      
    }
    
    
    estex <- list(
      or.fixed = or.fixed,
      or.se = or.se,
      or.vcov = or.vcov,
      or.test = or.test,
      alt.fixed = alt.fixed, 
      alt.se = alt.se,
      alt.vcov = alt.vcov,
      alt.test = alt.test)
    
    class(estex) <- "estex"
    return(estex)
    
  }



