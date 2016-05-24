
smnet<-function(formula, data.object, control = NULL, method = "AICC", netID = 1)
{  
  
  adjacency <- NULL
  
  # ERROR CHECK AND PROCESS THE control INPUT
  if(!is.null(control)){
    maxit <- control$maxit
    if(!is.null(maxit)) if(!is.wholenumber(maxit)) stop("maxit should be a positive integer or NULL", domain = NA, call. = FALSE)
    approx <- control$approx
    if(!is.null(approx)) if(!is.wholenumber(approx)) stop("approx should be a positive integer or NULL", domain = NA, call. = FALSE) 
    do.checks <- control$checks
    if(!is.null(do.checks)) if(!is.logical(do.checks)) stop("check should be logical or NULL", domain = NA, call. = FALSE) 
    if(is.null(do.checks)) do.checks <- TRUE
  }
  
  # ERROR CHECK AND PROCESS THE method INPUT
  if(!method %in% c("AICC", "AIC", "GCV")){
    err.msg <- "method should be either AIC, AICC or GCV"
    stop(err.msg, domain = NA, call. = FALSE)
  }
  
  # ERROR CHECK AND PROCESS THE DATA INPUT
  cls.data <- class(data.object)
  if(cls.data == "SpatialStreamNetwork"){
    data            <- getSSNdata.frame(data.object, Name = "Obs")
    data_netID      <- as.numeric(as.character(data$netID))
    data            <- data[data_netID == netID, ]
  } else if(cls.data == "data.frame"){
    data <- data.object
  } else stop("input data.object must be of class SpatialStreamNetwork or data.frame", domain = NA, call. = FALSE)
  
  # ERROR CHECK AND PROCESS THE FORMULA INPUT
  # if there is an error on the formula, print and error and show the list of available data columns
  formulaout <- try(formula_stuff <- get_formula_stuff(formula = formula, data = data), silent=T)
  if(class(formulaout) == "try-error"){
    var.names <- c(t(cbind(sort(as.character(names(data))), "\n")))
    err.msg   <- paste(c(var.names, "\nFormula error: check covariates are from the above"), collapse = " ")
    stop(err.msg, domain = NA, call. = FALSE)
  } 
  
  # construct the vector of weights, apply an ordering if appropriate
  cls.sm  <- lapply(formulaout$gp$smooth.spec, class)
  net     <- "network.spec" %in% cls.sm
  max.df  <- NULL
  if(net){
    networkSpecNo  <- which(cls.sm == "network.spec")
    networkSpec    <- formulaout$gp$smooth.spec[[networkSpecNo]]
    adjacency      <- networkSpec$adjacency
    if(!is.null(adjacency)){
      if(class(data.object) == "SpatialStreamNetwork"){
#         netID           <- networkSpec$netID
        # get the rids out of the SSN object associated with netID
        data            <- data[as.numeric(as.character(data$netID)) == netID,]
        rid_data        <- data.object@data
        rid_data        <- rid_data[as.numeric(as.character(rid_data$netID)) == netID,]
        rid_data$rid    <- re_map_rid(rid_data$rid, as.numeric(adjacency$rid_bid[,1]))
        if(nrow(rid_data) == 0) stop("No data associated with that netID", domain = NA, call. = FALSE)
        ord             <- order(rid_data$rid)
        response.locs   <- as.numeric(as.character(data$rid))
        response.locs   <- re_map_rid(response.locs, as.numeric(adjacency$rid_bid[,1]))
        weight          <- networkSpec$weight
        # set the maximum degrees of freedom for the network smooth (default of half number sites)
        if(is.null(networkSpec$max.df)) max.df <- length(unique(data$rid))/2 else max.df <- networkSpec$max.df
        # this happens if the user requests automatic stream weighting based on Shreve order
        if(weight == "autoShreve"){
          weight       <- adjacency_to_shreve(adjacency = adjacency)
          weight.type  <- check_weight(adjacency, weight)
          if(weight.type == "additive"){
            weight       <- get_shreve_weights(adjacency = adjacency$adjacency, shreve.order = as.matrix(weight))  
          } else if(weight.type == "unrecognised"){
            stop("supplied weight vector is neither additive or a network weighting", domain = NA, call. = FALSE)
          }
        } else {
          # this happens if the user requests the use of a ready made stream weighting or additive function
          weight       <- try(as.vector(as.matrix(rid_data[weight])), silent = T)
          if(class(weight) == "try-error"){
            # if the requested weight does not exist, show an error
            stop("Requested weight vector not found.  \n", domain = NA, call. = FALSE)
          }
          weight       <- weight[ord]
          weight.type  <- check_weight(adjacency, weight)
          if(weight.type == "additive"){
            weight     <- get_shreve_weights(adjacency = adjacency$adjacency, shreve.order = as.matrix(weight))  
          } else if(weight.type == "unrecognised"){
            stop("supplied weight vector is neither additive or a network weighting", domain = NA, call. = FALSE)
          }
        }
      }
      if(class(data.object) == "data.frame"){
        if(is.character(networkSpec$locs)){
          response.locs <- try(data[networkSpec$locs], silent = T)
          if(response.locs == "try-error"){
            stop("No data found for locs provided")
          }
        } else if(is.numeric(networkSpec$locs)){
          response.locs <- networkSpec$locs
        } else stop("locs should be numeric vector or a character string")
        weight <- networkSpec$weight
        if(!is.numeric(weight)) stop("If no data is not an SSN object, weight must be a vector", domain = NA, call. = FALSE)
        ord <- rid_data <- netID <- NULL
      }
      # check that weights vector and adjacency are the same dimension
      if(!nrow(adjacency$adjacency) == length(weight)) stop("weights and adjacency have different dimensions")
    } else {ord <- rid_data <- netID <- NULL; stop("network smooth requested, but adjacency matrix is missing", domain = NA, call. = FALSE)}
  } else ord <- rid_data <- weight <- NULL


  default.control = list(trace = 0, maxit = 500, start = NULL, approx = NULL, verbose = TRUE, tol = 10^-5)
  if(!is.null(control)) for(i in 1:length(control)) default.control[names(control)[i]]<-list(control[[i]])
  # create model objects
  model_objects <- get_model_objects(formula = formula, data = data, 
                                     adjacency = adjacency, 
                                     response.locs = response.locs, weight = weight,  
                                     rid_data = rid_data, netID = netID, ord = ord, 
                                     control = default.control, formulaout = formulaout)
  
  #  choose optimal smooth parameters using box constrained Nelder-Mead search
  opt <- with(model_objects, {
    if((n.smooth == 0) && (!net)){
      get_lm_fit(X.spam = X.spam, X.list = X.list, XTX.spam = XTX.spam, response = response, 
                 lin.means = lin.means, n.linear = n.linear)         
    } else {
      get_optimal_smooth(P.list = P.list, X.spam = X.spam, X.list=X.list, XTX.spam=XTX.spam, 
                         response=response, control = control, net=net, n.linear = n.linear,
                         lin.names = lin.names, sm.names = sm.names, lin.means = lin.means,
                         method = method, Pwee.list = Pwee.list, max.df = max.df)   
    }
    
#       P.list<-model_objects$P.list
#       Pwee.list<-model_objects$Pwee.list
#       X.spam<-model_objects$X.spam
#       X.list<-model_objects$X.list
#       XTX.spam<-model_objects$XTX.spam
#       response<-model_objects$response
#       control<-model_objects$control
#       n.linear<-model_objects$n.linear
#       sm.names<-model_objects$sm.names
#       smp<-model_objects$smp
#       crit<-model_objects$crit
#       verbose<-model_objects$verbose

  }
  )
  
  # calculate and print R2
  response      <- model_objects$response
  fitted.values <- opt$fit

  R2     <- 1 - sum((response - fitted.values)^2)/sum((response - mean(response))^2)
  np     <- ifelse(is.null(opt$ED), ncol(model_objects$X.spam), opt$ED)
  R2.adj <- R2 - ((1-R2)*np)/(nrow(model_objects$X.spam) - np - 1)
  cat(paste("   n = ", (nrow(model_objects$X.spam)), sep = ""))
  cat(paste("  R2.adj", " = ", round(R2.adj, 2), sep = ""))

  # create output list
  outputList               <- vector("list")
  outputList$response      <- response
  outputList$fitted.values <- fitted.values
  outputList$residuals     <- response - fitted.values
  outputList$coefficients  <- opt$beta_hat
  outputList$R2            <- R2
  outputList$R2.adj        <- R2.adj
  outputList$df.residual   <- length(response) - np
  outputList$ssn.object    <- data.object
  outputList$internals     <- c(opt, model_objects)
  class(outputList)        <- "smnet"
  invisible(outputList)
}

