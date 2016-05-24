ps.estimate.stratified.data.frame <- function(object,
                                              resp,
                                              treat         = NULL,
                                              stratum.index = NULL,
                                              match.index   = NULL,
                                              adj           = NULL,
                                              weights       = "rr",
                                              family        = "gaussian",
                                              regr          = FALSE,
                                              ...)
{
  data <- object$data

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
      stop("Argument 'treat' must be either numeric or a string.")
    }
  }
    
  ## ###########
  ## find strata
  if (!is.character(stratum.index) & !is.numeric(stratum.index) & !is.null(stratum.index)){
    stop("Argument 'stratum.index' has to be either numeric or a string.")
  }else{
    if (is.null(stratum.index)){
      stratum.index        <- object$stratum.index
      name.stratum.index   <- object$name.stratum.index
      levels.stratum.index <- levels(as.factor(stratum.index))
    }else{
      A <- find.strata(data   = data,
                       strata = stratum.index,
                       match  = FALSE)    
      stratum.index        <- A[[1]]
      name.stratum.index   <- A[[2]]
      levels.stratum.index <- A[[3]]
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

  if (any(name.adj == name.stratum.index) |
      any(name.adj == name.treat)){
    stop("Argmuent 'adj' contains argument 'stratum.index' or 'treat'.")}

  if ((nlevels(as.factor(resp)) != 2 & family=="binomial") |
      (nlevels(as.factor(resp)) == 2 & family=="gaussian"))
    
    warning(paste("Discrepancy between argument 'family' = ",
                  family,
                  " and scale of argument 'resp' = ",
                  name.resp,
                  " ... internally change.", sep=""))
 

  ## #####################################################
  ## do stratified analysis depending on scale of response
  if (nlevels(as.factor(resp)) != 2){

    object$ps.estimation <- ps.estimate.gauss(data                 = data,
                                              resp                 = resp,
                                              name.resp            = name.resp,
                                              treat                = treat,
                                              name.treat           = name.treat,
                                              stratum.index        = stratum.index,
                                              adj.form             = adj.form,
                                              adj                  = adj,
                                              name.adj             = name.adj,
                                              weights              = weights,
                                              levels.stratum.index = levels.stratum.index,
                                              family               = family,
                                              match.T              = FALSE,
                                              match.id             = NULL,
                                              ...)  
  }else{
    
    family <- "binomial"    
    object$ps.estimation <- ps.estimate.bin(data                 = data,
                                            resp                 = resp,
                                            name.resp            = name.resp,
                                            treat                = treat,
                                            name.treat           = name.treat,
                                            stratum.index        = stratum.index,
                                            adj.form             = adj.form,
                                            adj                  = adj,
                                            name.adj             = name.adj,
                                            weights              = weights,
                                            levels.stratum.index = levels.stratum.index,
                                            family               = family,
                                            match.T              = FALSE,
                                            match.id             = NULL,
                                            ...)  
  }

  ## ####################################################
  ## find regr; regr/data.regr are reduced to data w/o NA
  if (!is.logical(regr)){
    if (!inherits(regr, "formula")){      
      regr <- find.sel(data, ## data = original data
                       regr,
                       "regr")

      name.regr <- names(regr)
      
      lr.form <- formula(paste(name.resp,
                            paste(c(name.treat, name.regr),collapse="+"), sep="~"))

      resp.regr  <- as.data.frame(model.frame(lr.form,data)[,1]) ## resp
      treat.regr <- as.data.frame(model.frame(lr.form,data)[,2]) ## treat
      regr       <- as.data.frame(model.matrix(lr.form,data)[,-c(1:2)]) ## w/o resp, treat

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

      resp.regr  <- as.data.frame(model.frame(regr,data)[,1]) ## resp
      treat.regr <- as.data.frame(model.frame(regr,data)[,2]) ## treat
      name.regr <- names(as.data.frame(model.frame(regr,data)))[-c(1:2)]

      if (any(names(as.data.frame(model.frame(regr,data)))[c(1)] != name.resp))
        stop(paste("Argument 'regr' does not include ", name.resp,
                   " as response.", sep=""))
      
      if (any(names(as.data.frame(model.frame(regr,data)))[c(2)] != name.treat))
        stop(paste("Argument 'regr' does not include ", name.treat,
                   " as treatment.", sep=""))

      regr       <- as.data.frame(model.matrix(regr,data)[,-c(1:2)]) ## w/o resp
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
  if (identical(regr,FALSE)){    
    object$lr.estimation <- "No regression"
  }else{

    object$lr.estimation <- ps.estimate.regr(data           = data.regr,
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
  object$resp <- resp
  object$name.resp <- name.resp
  
  object$treat <- treat
  object$name.treat <- name.treat
  
  object$stratum.index <- stratum.index
  object$name.stratum.index <- name.stratum.index

  object$family <- family

  class(object) <- c("est.stratified.data.frame",
                     class(object)[class(object)!="est.stratified.data.frame"])

  return(object)
  
}
