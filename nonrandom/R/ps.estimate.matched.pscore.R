ps.estimate.matched.pscore <- function(object,
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
  ## #################
  ## Check match.index
  if (is.null(match.index)){
    d1 <- object$data
    d2 <- object$data.matched
 
    d1[,"match.strata"] <- rep(1,nrow(d1))  ## unmatched
    d2[,"match.strata"] <- rep(2,nrow(d2))  ## matched
  
    data <- (rbind(d1,d2))
    data$match.strata <- as.factor(data$match.strata)

    match.index <- data$match.strata
    name.match.index <- "match.strata"

    match.id <- d1[,object$name.match.index]
     
  }else{
    if (is.character(match.index) | is.numeric(match.index)){
      if (is.character(match.index)){

        match.id <- object$data[, match.index]

        A <- find.strata(data   = object$data,
                         strata = match.index,
                         match  = TRUE)

        d1 <- object$data             ## A[[1]] contains match index
        d2 <- object$data[A[[1]]>0,]  ## A[[1]]>0 : matched obs

        d1[,"match.strata"] <- rep(1,nrow(d1))  ## unmatched
        d2[,"match.strata"] <- rep(2,nrow(d2))  ## matched
        
        data <- (rbind(d1,d2))
        data$match.strata <- as.factor(data$match.strata)

        match.index <- data$match.strata
        name.match.index <- A[[2]]
        
      }else{

        match.id <- object$data[, match.index]

        A <- find.strata(data   = object$data,
                         strata = match.index,
                         match  = TRUE)    

        d1 <- object$data
        d2 <- object$data[A[[1]]>0,]  ## A[[1]]>0 : matched obs

        d1[,"match.strata"] <- rep(1,nrow(d1))  ## unmatched
        d2[,"match.strata"] <- rep(2,nrow(d2))  ## matched
        
        data <- (rbind(d1,d2))
        data$match.strata <- as.factor(data$match.strata)

        match.index <- data$match.strata
        name.match.index <- A[[2]]
        
      }
    }else{
      stop("Argument 'match.index' must be numeric or a string.")
    } 
  }



  ## ##########
  ## find treat
  if (is.null(treat)){   
    name.treat <- object$name.treat
    treat      <- data[,name.treat]
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
                                              stratum.index        = data$match.strata,
                                              adj.form             = adj.form,
                                              adj                  = adj,
                                              name.adj             = name.adj,
                                              levels.stratum.index = 1,
                                              family               = family,
                                              match.T              = TRUE,
                                              match.id             = match.id,
                                              ...)

  }else{
    
    family <- "binomial"
  
    object$ps.estimation <- ps.estimate.bin(data                 = data, 
                                            resp                 = resp,
                                            name.resp            = name.resp,
                                            treat                = treat,
                                            name.treat           = name.treat,
                                            stratum.index        = data$match.strata,
                                            adj.form             = adj.form,
                                            adj                  = adj,
                                            name.adj             = name.adj,
                                            levels.stratum.index = 1,
                                            family               = family,
                                            match.T              = TRUE,
                                            match.id             = match.id,
                                            ...)
  }

  ## ####################################################
  ## find regr; regr/data.regr are reduced to data w/o NA
  if (!is.logical(regr)){
    if (!inherits(regr, "formula")){      
      regr <- find.sel(d1, ## d1 = original data
                       regr,
                       "regr")

      name.regr <- names(regr)
      
      lr.form <- formula(paste(name.resp,
                            paste(c(name.treat, name.regr),collapse="+"), sep="~"))

      resp.regr  <- as.data.frame(model.frame(lr.form,d1)[,1]) ## resp
      treat.regr <- as.data.frame(model.frame(lr.form,d1)[,2]) ## treat
      regr       <- as.data.frame(model.matrix(lr.form,d1)[,-c(1:2)]) ## w/o resp, treat

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

      resp.regr  <- as.data.frame(model.frame(regr,d1)[,1]) ## resp
      treat.regr <- as.data.frame(model.frame(regr,d1)[,2]) ## treat
      name.regr <- names(as.data.frame(model.frame(regr,d1)))[-c(1:2)]

      if (any(names(as.data.frame(model.frame(regr,data)))[c(1)] != name.resp))
        stop(paste("Argument 'regr' does not include ", name.resp,
                   " as response.", sep=""))
      
      if (any(names(as.data.frame(model.frame(regr,data)))[c(2)] != name.treat))
        stop(paste("Argument 'regr' does not include ", name.treat,
                   " as treatment.", sep=""))
      
      regr       <- as.data.frame(model.matrix(regr,d1)[,-c(1:2)]) ## w/o resp
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

  object$match.index  <- data$match.strata

  object$family <- family

  class(object) <- c("est.matched.pscore",
                     class(object)[class(object)!="est.matched.pscore"])


  return(object)
  

}
