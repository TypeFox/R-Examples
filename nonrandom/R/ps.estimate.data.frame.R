ps.estimate.data.frame <- function(object,
                                   resp,
                                   treat         = NULL,
                                   stratum.index = NULL,
                                   match.index   = NULL,
                                   adj           = NULL,
                                   weights       = "rr",
                                   family        = "gaussian",
                                   regr          = NULL,
                                   ...) 
{
  data <- object
  
  ## ##########
  ## find treat
  if (is.null(treat)){
    stop("Argument 'treat' is needed.")
  }else{
    if (is.character(treat) | is.numeric(treat)){
      A <- find.treat(data  = data,
                      treat = treat)
      treat      <- A[[1]]
      name.treat <- A[[2]]
    }else{
      stop("Argument 'treat' has to be either numeric or a string.") 
    }    
  }
  
  ## #########
  ## find resp
  if (missing(resp)){    
    stop("Argument 'resp' is needed.")
  }else{
    if (!is.character(resp) & !is.numeric(resp)){
      stop("Argument 'resp' has to be either numeric or a string.")
    }else{
      A <- find.resp(data = data,
                     resp = resp)
      resp      <- A[[1]]
      name.resp <- A[[2]]
    }   
  }

  if (name.treat == name.resp)
    stop("Arguments 'treat' and 'resp' are identical.")

  
  ## ########
  ## find adj
  if (!is.null(adj)){
    if (!inherits(adj, "formula")){ 
      adj <- find.sel(data,
                      adj,
                      "adj")
      name.adj <- names(adj)
      adj.form <- formula(paste(name.resp,
                                paste(c(name.treat, names(adj)),collapse="+"),
                                sep="~"))
    }else{
      adj.form <- adj
      adj <- model.frame(adj.form,data)[,-c(1:2)] ## w/o resp and treat
      name.adj <- names(model.frame(adj.form,data))[-c(1:2)]
    }
  }else{
    adj.form <- NULL
    name.adj <- NULL
  }
  
  if (any(name.adj == name.treat)){
    warning("Argmuent 'adj' contains argument 'treat'.")}
  if (any(name.adj == name.resp)){
    warning("Argmuent 'adj' contains argument 'resp'.")}
  

  
  ## #################################
  ## find stratum.index or match.index
  if ( is.null(stratum.index) & is.null(match.index) ){
    
    stop("Either argument 'stratum.index' or 'match.index' is needed.")
    
  }else{
    
    if ( !is.null(stratum.index) & !is.null(match.index) ){     
      warning("Both arguments 'stratum.index' and 'match.index' are given: 'match.index' is used.")
    }
    if ( !is.null(match.index) ){ ## matching
      
      stratum.index      <- NULL
      name.stratum.index <- NULL
      
      if ( is.character(match.index) | is.numeric(match.index) ){
        
        A <- find.strata(data    = data,
                         strata  = match.index,
                         match   = TRUE)    
        match.index      <- A[[1]]
        name.match.index <- A[[2]]
        
        match.strata <- c(rep(1,nrow(data)),
                          rep(2,nrow(data[match.index != 0,])))  
        data.new <- rbind(data, data[match.index != 0,])
        data.new$match.strata <- match.strata

        if ( !is.null(adj) ){         
          adj.match.index        <- as.data.frame(adj[match.index != 0,])
          names(adj.match.index) <- names(adj)
          adj1 <- rbind(adj, adj.match.index)
        }else{         
          adj1 <- NULL
        }
        
        treat1 <- c(treat, treat[match.index != 0])
        resp1 <- c(resp, resp[match.index != 0])
        
      }else{
        stop("Argument 'match.index' has to be either numeric or a string.")
      }
      
    }else{ ## stratification 
      
      name.match.index <- NULL
      
      if (is.character(stratum.index) | is.numeric(stratum.index)){
        
        A <- find.strata(data   = data,
                         strata = stratum.index,
                         match  = FALSE)    
        
        stratum.index        <- A[[1]]
        name.stratum.index   <- A[[2]]
        levels.stratum.index <- A[[3]]
        
      }else{
        
        stop("Argument 'stratum.index' has to be either numeric or a string.")

      }
      
      if (is.null(stratum.index)){
        if (any(names(adj) == name.treat) |
            any(names(adj) == name.match.index)){
          warning("Argmuent 'adj' contains argument 'treat' or 'match.index'.")
        }
      }else{
        if (any(names(adj) == name.treat) |
            any(names(adj) == name.stratum.index)){
          warning("Argmuent 'adj' contains argument 'treat' or 'stratum.index'.")
        }
      }
    }
  }
    
  if (is.null(stratum.index)){ ## matched

    data.lr  <- data
    data     <- data.new
    treat    <- treat1
    resp     <- resp1
    adj      <- adj1
    index    <- data.new$match.strata
    match.T  <- TRUE
    match.id <- match.index
    
  }else{ ## stratified

    data.lr  <- data
    index    <- stratum.index
    match.T  <- FALSE
    match.id <- NULL

  }
  

    
  ## #####################################################
  ## do stratified analysis depending on scale of response
  if (nlevels(as.factor(resp)) != 2){
    
    ps.estimation <- ps.estimate.gauss(data                 = data, 
                                       resp                 = resp,
                                       name.resp            = name.resp,
                                       treat                = treat,
                                       name.treat           = name.treat,
                                       stratum.index        = index,
                                       adj.form             = adj.form,
                                       adj                  = adj,
                                       name.adj             = name.adj,
                                       weights              = weights,
                                       levels.stratum.index = levels.stratum.index,
                                       family               = family,
                                       match.T              = match.T,
                                       match.id             = match.id,
                                       ...)
    
  }else{
    
    family <- "binomial"   
    ps.estimation <- ps.estimate.bin(data                 = data, 
                                     resp                 = resp,
                                     name.resp            = name.resp,
                                     treat                = treat,
                                     name.treat           = name.treat,
                                     stratum.index        = index,
                                     adj.form             = adj.form,
                                     adj                  = adj,
                                     name.adj             = name.adj,
                                     weights              = weights,
                                     levels.stratum.index = levels.stratum.index,
                                     family               = family,
                                     match.T              = match.T,
                                     match.id             = match.id,
                                     ...)
  }


  ## #########
  ## find regr
  if (is.null(regr)){

    lr.estimation <- "No regression"
    
  }else{

    if (!is.logical(regr)){
      if (!inherits(regr, "formula")){      
        regr <- find.sel(data.lr, ## data = original data
                         regr,
                         "regr")
        
        name.regr <- names(regr)
        
        lr.form <- formula(paste(name.resp,
                                 paste(c(name.treat, name.regr),collapse="+"), sep="~"))
        
        resp.regr  <- as.data.frame(model.frame(lr.form,data.lr)[,1]) ## resp
        treat.regr <- as.data.frame(model.frame(lr.form,data.lr)[,2]) ## treat
        regr       <- as.data.frame(model.matrix(lr.form,data.lr)[,-c(1:2)]) ## w/o resp, treat
        
        ## Falls nur eine Var in 'regr' ausgewaehlt, names(regr) nach
        ## model.frame funktioniert nicht, deswegen ==>
        if (length(names(regr)) != 1)
          name.regr <- names(regr)
        
        data.regr <- as.data.frame(cbind(resp.regr, treat.regr, regr))
        names(data.regr) <- c(name.resp, name.treat, name.regr)
        name.data.regr   <- names(data.regr)
        
        lr.form <- formula(paste(name.resp,
                                 paste(c(name.treat, name.regr),collapse="+"), sep="~"))
        
      }else{
        
        resp.regr  <- as.data.frame(model.frame(regr,data.lr)[,1]) ## resp
        treat.regr <- as.data.frame(model.frame(regr,data.lr)[,2]) ## treat
        name.regr <- names(as.data.frame(model.frame(regr,data.lr)))[-c(1:2)]

        if (any(names(as.data.frame(model.frame(regr,data.lr)))[c(1)] != name.resp))
          stop(paste("Argument 'regr' does not include ", name.resp,
                     " as response.", sep=""))
        
        if (any(names(as.data.frame(model.frame(regr,data.lr)))[c(2)] != name.treat))
          stop(paste("Argument 'regr' does not include ", name.treat,
                     " as treatment.", sep=""))
        
        regr       <- as.data.frame(model.matrix(regr,data.lr)[,-c(1:2)]) ## w/o resp
        data.regr <- as.data.frame(cbind(resp.regr, treat.regr, regr))
        
        ## Falls nur eine Var in 'regr' ausgewaehlt, names(regr) nach
        ## model.frame funktioniert nicht, deswegen ==>
        if (length(names(regr)) != 1)
          name.regr <- names(regr)
        
        names(data.regr) <- c(name.resp, name.treat, name.regr)
        name.data.regr   <- names(data.regr)
        
        lr.form <- formula(paste(name.resp,
                                 paste(c(name.treat, name.regr),collapse="+"), sep="~"))
        
      }
    }
    
    
    ## #############
    ## do regression
    lr.estimation <- ps.estimate.regr(data           = data.regr,
                                      name.data.regr = name.data.regr,
                                      regr           = regr,
                                      name.regr      = name.regr,
                                      treat          = treat,
                                      name.treat     = name.treat,
                                      resp           = resp,
                                      name.resp      = name.resp,
                                      lr.form        = lr.form,
                                      family         = family,
                                      ...)
  }

  ## Output
  object <- list(data               = data.lr,
                 resp               = resp,
                 name.resp          = name.resp,
                 treat              = treat,
                 name.treat         = name.treat,
                 name.stratum.index = name.stratum.index,
                 stratum.index      = stratum.index,
                 name.match.index   = name.match.index,
                 match.index        = match.index,
                 ps.estimation      = ps.estimation,
                 lr.estimation      = lr.estimation,
                 family             = family)

  class(object) <- c("est.data.frame",
                     class(object)[class(object)!="est.data.frame"])
  
  return(object)
  
}
