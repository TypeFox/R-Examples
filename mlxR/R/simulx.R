#' Simulation of mixed effects models and longitudinal data
#'
#' Compute predictions and sample data from \code{Mlxtran}, \code{R} and \code{PharmML} models
#' 
#' simulx takes advantage of the modularity of hierarchical models for simulating 
#' different components of a model: models for population parameters, individual 
#' covariates, individual parameters and longitudinal data.
#' 
#' Furthermore, \code{simulx} allows to draw different types of longitudinal data, 
#' including continuous, count, categorical, and time-to-event data.
#' 
#' The models are encoded using either the model coding language \samp{Mlxtran}, \samp{R} or the 
#' markup language \samp{PharmML}. \samp{Mlxtran} models are automatically converted into C++ codes, 
#' compiled on the fly and linked to R using the \samp{Rcpp} package. 
#' That allows one to implement very easily complex models and to take advantage 
#' of the numerical sovers used by the C++ \samp{mlxLibrary}.
#' 
#' See http://simulx.webpopix.org for more details.      
#' @param model a \code{Mlxtran}, \code{R} or \code{PharmML} model used for the simulation
#' @param parameter a vector of parameters with their names and values
#' @param output a list (or list of lists) with fields: 
#' \itemize{
#'   \item \code{name}: a vector of output names
#'   \item \code{time}: a vector of times (only for the longitudinal outputs)
#' }
#' @param treatment a list with fields
#' \itemize{
#'   \item \code{time} : a vector of input times,
#'   \item \code{amount} : a scalar or a vector of amounts,
#'   \item \code{rate} : a scalar or a vector of infusion rates (default=\code{Inf}),
#'   \item \code{tinf} : a scalar or a vector of infusion times (default=0),
#'   \item \code{type} : the type of input (default=1),
#'   \item \code{target} : the target compartment (default=NULL). 
#' }
#' @param regressor a list, or a list of lists, with fields
#' \itemize{
#'   \item \code{name} : a vector of regressor names,
#'   \item \code{time} : a vector of times,
#'   \item \code{value} : a vector of values.
#' }
#' @param varlevel (not supported by mlxLibrary 2016R1)
#' @param group a list, or a list of lists, with fields: 
#' \itemize{
#'   \item \code{size} : size of the group (default=1),
#'   \item \code{level} : level(s) of randomization,
#'   \item \code{parameter} : if different parameters per group are defined,
#'   \item \code{output} : if different outputs per group are defined,
#'   \item \code{treatment} : if different treatements per group are defined,
#'   \item \code{regressor} : if different regression variables per group are defined.
#' }
#' @param data a list (output of simulx when settings$data.in==TRUE)
#' @param project the name of a Monolix project
#' @param nrep number of replicates
#' @param npop number of population parameters to draw randomly 
#' @param fim a string with the Fisher Information Matrix to be used 
#' @param result.folder the name of the folder where the outputs of simulx should be stored
#' @param result.file the name of the single file where the outputs of simulx should be saved
#' @param stat.f a R function for computing some summary (mean, quantiles, survival,...) of the simulated data. Default = "statmlx".
#' @param settings a list of optional settings
#' \itemize{
#'   \item \code{seed} : initialization of the random number generator (integer),
#'   \item \code{load.design} : TRUE/FALSE (if load.design is not defined, a test is automatically performed to check if a new design has been defined),
#'   \item \code{data.in} : TRUE/FALSE (default=FALSE)
#'   \item \code{id.out}  : add columns id (when N=1) and group (when #group=1), TRUE/FALSE (default=FALSE)
#'   \item \code{kw.max} : maximum number of trials for generating a positive definite covariance matrix (default = 100) 
#'   \item \code{sep} : the field separator character (default = ",") 
#'   \item \code{digits} : number of decimal digits in output files (default = 5) 
#'   \item \code{disp.iter} : TRUE/FALSE (default = FALSE) display replicate and population numbers
#'   \item \code{replacement} : TRUE/FALSE (default = FALSE) sample id's with/without replacement
#'   \item \code{out.trt} : TRUE/FALSE (default = TRUE) output of simulx includes treatment
#' }       
#' 
#' @return A list of data frames. Each data frame is an output of simulx
#' 
#' @examples
#' \dontrun{
#' myModel <- inlineModel("
#' [LONGITUDINAL]
#' input = {A, k, c, a}
#' EQUATION:
#' t0    = 0 
#' f_0   = A
#' ddt_f = -k*f/(c+f)
#' DEFINITION:
#' y = {distribution=normal, prediction=f, sd=a}
#' [INDIVIDUAL]
#' input = {k_pop, omega}
#' DEFINITION:
#' k = {distribution=lognormal, prediction=k_pop, sd=omega}
#' ")
#' f <- list(name='f', time=seq(0, 30, by=0.1))
#' y <- list(name='y', time=seq(0, 30, by=2))
#' res <- simulx(model     = 'model/home.txt', 
#'               parameter = c(A=100, k_pop=6, omega=0.3, c=10, a=2), 
#'               output    = list(f,y,"k"),
#'               group     = list(size=4, level='individual'))
#' 
#' plot(ggplotmlx() + geom_line(data=res$f, aes(x=time, y=f, colour=id)) +
#'      geom_point(data=res$y, aes(x=time, y=y, colour=id)))
#' print(res$parameter)
#' }
#' 
#' @importFrom stats runif
#' @export

simulx <- function(model=NULL, parameter=NULL, output=NULL,treatment=NULL, 
                   regressor=NULL, varlevel=NULL, group=NULL, 
                   data=NULL, project=NULL, nrep=1, npop=NULL, fim=NULL, 
                   result.folder=NULL, result.file=NULL, stat.f="statmlx",
                   settings=NULL)
{ 
  #--------------------------------------------------
  #  simulx.R is governed by the CeCILL-B license. 
  #  You can  use, modify and/ or redistribute the software under the terms of 
  #  the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
  #  http://www.cecill.info/index.en.html
  #
  #  simulx.R was developed by Marc Lavielle and the Inria popix team for the DDMoRe project. 
  #--------------------------------------------------
  
  myOldENVPATH = Sys.getenv('PATH');
  initMlxLibrary()
  session=Sys.getenv("session.simulx")
  if (!is.null(varlevel) && grepl('MonolixSuite2016R1',session))
    # if (!is.null(varlevel) && grepl('mlxsuite-release',session))
  {
    cat("\nvarlevel is not supported with this version of mlxLibrary\n")
    return()
  }
  Sys.setenv(LIXOFT_HOME=session)
  
  if (is.null(settings$seed))
    settings$seed <- round(runif(1)*100000)
  set.seed(settings$seed)
  disp.iter <- ifelse((!is.null(settings$disp.iter) && settings$disp.iter==TRUE), TRUE, FALSE)
  sep <- settings$sep
  if (is.null(settings$sep)) sep <- ","
  digits <- settings$digits
  if (is.null(settings$digits)) digits <- 5
  kw.max <- settings$kw.max
  if (is.null(settings$kw.max)) kw.max <- 500
  replacement <- settings$replacement
  if (is.null(settings$replacement)) replacement <- FALSE
  out.trt <- settings$out.trt
  if (is.null(settings$out.trt))  out.trt <- T
  
  if (!is.null(data))
  {
    r <- simulxunit(data=data,settings=settings)
    Sys.setenv(LIXOFT_HOME="")
    Sys.setenv('PATH'=myOldENVPATH);
    return(r)
  }
  
  if (isfield(settings,"record.file"))  
    warning("\n\n 'record.file' is a deprecated option. Use 'result.file' instead.")
  #   if (!is.null(result.file) &&  !is.null(result.folder))
  #     warning("\n\n You can define either 'result.file' or 'result.folder', not both...")
  
  #--------------------------------------------------
  #    R MODEL
  #--------------------------------------------------
  Rmodel <- FALSE
  if (identical(file_ext(model),"R")) 
    Rmodel <- TRUE
  else if ( !is.null(model) && exists(model, mode="function") )
    Rmodel <- TRUE
  
  #--------------------------------------------------
  #     reshape inputs
  #--------------------------------------------------
  group     <- mklist(group)
  parameter <- mklist(parameter)
  treatment <- mklist(treatment)
  regressor <- mklist(regressor)
  varlevel  <- mklist(varlevel)
  output    <- mklist(output)
  
  #--------------------------------------------------
  #     stat output
  #--------------------------------------------------
  if (!is.null(result.folder) || !is.null(result.file))
    write.simul=TRUE
  else
    write.simul=FALSE
  
  stat.n <- NULL
  stat.0 <- NULL
  stat.a <- list()
  for (k in (1:length(output))){
    outk <- output[[k]]
    outk$type <- NULL
    if (is.null(outk$time))
      outk.n <- "parameter"
    else
      outk.n <- outk$name
    outk$name <- NULL
    outk$time <- NULL
    if (length(outk)>0)
    {
      stat.n <- c(stat.n, outk.n)
      stat.a <- c(stat.a, rep(list(outk),length(outk.n)))
    } else if (write.simul==T) {
      stat.0 <- c(stat.0, outk.n)
    }
  }
  if (write.simul==T)
    stat.0 <- c(stat.0, "treatment")
  names(stat.a)=stat.n
  
  #--------------------------------------------------
  #     Monolix project
  #--------------------------------------------------
  if (!(is.null(project)))  
  {
    
    if (!is.null(npop))
    {
      iproj.pop <- T
      if (is.null(fim))  fim <- "needed"
    }
    else
      iproj.pop <- F
    
    if (is.list(parameter[[1]]) && !is.null(parameter[[1]]$pop))
      p.pop <- parameter[[1]]
    else
      p.pop <- NULL
    
    ans <- processing_monolix( project=project,
                               treatment=treatment,
                               parameter=parameter,
                               output=output,
                               group=group,
                               fim=fim
    )
    model     <- ans$model
    treatment <- ans$treatment
    parameter <- ans$param
    output    <- ans$output
    group     <- ans$group
    regressor <- ans$regressor
    varlevel  <- ans$occasion
    fim       <- ans$fim
    infoParam <- ans$infoParam
    id        <- as.factor(ans$id$oriId)
    N         <- nlevels(id)
    if (iproj.pop==T)
    {
      if(is.null(fim))
        stop('The covariance matrix of the population parameters is requested for simulating several replicates 
           of the population parameters')
      else 
        parameter[[1]] <- sim.pop(npop,parameter[[1]],infoParam,fim,kw.max=kw.max)
    } 
    else if (!is.null(p.pop))
    {
      parameter[[1]] <- p.pop
      iproj.pop <- T
    }
    test.project <- T
    test.N <- T
    test.pop <- F
  }
  else
  {
    test.project <- F
    #--------------------------------------
    
    l.input <- c('parameter', 'treatment', 'regressor', 'varlevel', 'output')
    popid <- list()
    N <- NULL
    id <- NULL
    test.N <- F
    test.pop <- F
    for (k in (1:length(l.input)))
    {
      lk <- l.input[k]
      eval(parse(text=paste0('pk <- ',lk))) 
      pk<- dpopid(pk,lk)
      if (!is.null(pk$N))
      {
        test.N <- T
        if (!is.null(pk$id)) 
          id <- unique(c(id, pk$id))
      }
      if (!is.null(pk$npop))
      {
        test.pop <- T
        npop <- unique(c(npop,pk$npop))
        if (length(npop)>1)
          stop('\n\nDifferent numbers of populations are defined in different inputs')
      }
      popid[[k]] <- pk
    }
    if (test.N==T)
      id <- as.factor(id)
    N <- length(id)
  }
  #--------------------------------------------------
  #     Pop parameters
  #--------------------------------------------------
  if (test.pop == T)
  {
    for (k in (1: length(parameter)))
    {
      paramk <- parameter[[k]]
      if (isfield(paramk,"pop"))
      {
        if (isfield(paramk,"id"))
          stop("\n\n Both 'id' and 'pop' cannot be defined in the same data frame")
        npop <- nrow(paramk)
        paramk$pop <- NULL
        parameter[[k]] <- paramk[1,]
        if (npop>1){
          test.pop <- T
          k.pop <- k
          pop.mat <- paramk
        }
      }
    } 
  }
  
  lv <- list(treatment=treatment,
             parameter=parameter,
             output=output,
             regressor=regressor,
             varlevel=varlevel,
             id=id)
  
  if (test.N==T && !is.null(group))
  {
    if (any(sapply(group, function(x) is.null(x$size))))
      stop("'size' is missing in group")
    g.size <- sapply(group, function(x) x$size)
    if (any(sapply(group, function(x) !is.null(x$level))))
      warning("\n\n'level' in group is ignored when id's are defined in the inputs of simulx")
    ng <- length(group)
    if (ng==1)
    {
      if (!identical(names(group[[1]]),'size'))
        stop("\n\nOnly 'size' can be defined in group when a single group is created and when id's are defined in the inputs of simulx")
    } else {
      u.name <- unique(unlist(sapply(group, function(x) names(x))))
      if (!all(u.name %in% c("size","treatment")))
        stop("\n\nOnly 'size' and 'treatment' can be defined in group when several groups are created and when id's are defined in the inputs of simulx")
      if ("treatment" %in% u.name)
      {
        tr <- NULL
        for (k in (1:ng))
        {
          tk <- as.data.frame(group[[k]]$treatment)
          tk <- tk[rep(seq.int(1,nrow(tk)), group[[k]]$size), ]
          #           tk$group <- k
          if (is.null(tr))
            tr <- tk
          else
            tr <- merge(tr,tk, all=T)
        }
        tr[is.na(tr)]='.'
        tr <- cbind(data.frame(id=as.factor(seq(1:sum(g.size)))),tr)
        lv$treatment <- tr
      }
      gr.ori <- NULL
      for (k in (1:ng))
        gr.ori <- c(gr.ori, rep(k, group[[k]]$size))
      lv$gr.ori <- as.factor(gr.ori)
    }
  }
  
  if (is.null(N)) N<-1
  if (is.null(npop)) npop<-1
  
  test.rep <- FALSE
  if (nrep>1){
    if (test.N == F)
      test.rep <- TRUE
    else if (is.null(group))
      test.rep <- TRUE
  }
  
  R.complete <- list()
  rs <- NULL
  
  
  for (ipop in (1:npop))
  {
    if (disp.iter==TRUE) 
    {
      if (nrep>1) 
        cat("\n")
      if (npop>1)
        cat("population: ",ipop,"\n")
    }
    irw <- 0
    if (test.pop == T)  parameter[[k.pop]] <- pop.mat[ipop,]
    if (test.rep == T)
    {
      if (test.N==F)  
        lv$group <- group
      dataIn <- simulxunit(model=model,lv=lv,settings=c(settings, data.in=T))
      settings$data.in=F
      settings$load.design=F
    }
    for (irep in (1:nrep))
    {
      irw <- irw + 1
      settings$seed <- settings$seed +12345
      
      
      if (disp.iter==TRUE && nrep>1)  
        cat("replicate: ",irw,"\n")
      
      if (test.rep == T)
      {
        r <- simulxunit(data = dataIn,settings=settings, out.trt=out.trt)
        # r$treatment <- NULL
      } 
      else 
      {
        if (test.N==T && !is.null(group))  
          lv <- resample.data(lv,id,sum(g.size),replacement)
        if (test.N==F)  
          lv$group <- group
        r <- simulxunit(model=model,lv=lv,settings=settings, out.trt=out.trt)
      }
      
      rs <- r
      rs[stat.0] <- NULL
      if (length(stat.n) > 0)
      {
        for (k in (1:length(stat.n)))
        {
          rnk <- stat.n[k]
          if (!is.null(rnk))
          {
            resak <- stat.a[[rnk]]
            rs[[rnk]] <- do.call(stat.f, c(list(rs[[rnk]]),resak))
          }
        }
      }
      
      r.attr <- sapply(r,attr,"type")
      if (nrep>1)
      {
        for (k in (1:length(rs)))
          if (is.data.frame(rs[[k]]))
            rs[[k]] <- cbind(list(rep=as.factor(irw)), rs[[k]])
      }
      
      
      if (write.simul==TRUE)
      {
        for (k in (1:length(r)))
        {
          r[[k]] <- data.frame(lapply(r[[k]], function(y) if(is.numeric(y)) signif(y, 5) else y)) 
          if (npop>1)  
            r[[k]] <- cbind(list(pop=as.factor(ipop)),r[[k]])
          if (nrep>1)  
            r[[k]] <- cbind(list(rep=as.factor(irw)), r[[k]])
          attr(r[[k]],"type") <- r.attr[k]
        }
        if (ipop==1 & irep==1)
          app <- F
        else
          app <- T
        writeDatamlx(r,result.folder=result.folder,result.file=result.file,
                     sep=sep,digits=digits,app.dir=app,app.file=app)
        
      } 
      if (irep==1) 
      {
        res <- rs
      } else 
      {
        for(k in (1:length(rs)))
          if (is.data.frame(rs[[k]]))
            res[[k]] <- rbind(res[[k]],rs[[k]])
          else
            res[[k]] <- rs[[k]]
      }  
    } # irep
    if (length(res)>0)
    {
      for (k in (1:length(res)))
      {
        Rk <- res[[k]]
        if (is.data.frame(Rk))
        {
          if (npop>1)
            Rk <- cbind(list(pop=as.factor(ipop)),Rk)
          if (ipop==1)
            R.complete[[k]] <- Rk
          else
            R.complete[[k]] <- rbind(R.complete[[k]],Rk)
        } else
          R.complete[[k]] <- Rk
      }
    }
    
  } # ipop
  
  if (length(res)>0)
  {
    names(R.complete) <- names(res)
    for (k in (1:length(res)))
    {
      attrk <- attr(r[[names(res)[k]]],'type')
      if (!is.null(attrk))
        attr(R.complete[[k]],"type") <- attrk
    } 
  }
  
  pop <- NULL
  if (test.pop == T)
  {
    pop <- as.data.frame(pop.mat)
    pop <- format(pop, digits = 5, justify = "left")
    pop <- cbind(pop=(1:npop),pop)
  }
  if (!(is.null(project))) 
    pop <- parameter[[1]]
  if (!is.null(pop))
  {
    if (write.simul==TRUE)
    {
      r <- list(population=pop)
      writeDatamlx(r,result.folder=result.folder,sep=sep,digits=digits,app.dir=T)
    } 
    R.complete$population <- pop
  }
  
  Sys.setenv(LIXOFT_HOME="")
  Sys.setenv('PATH'=myOldENVPATH);
  return(R.complete)
}


#--------------------------------------------------
#--------------------------------------------------
#       Simulxunit
#--------------------------------------------------
#--------------------------------------------------

simulxunit <- function(model=NULL, lv=NULL, data=NULL, settings=NULL, out.trt=T)
{ 
  #--------------------------------------------------
  # Manage settings
  #--------------------------------------------------
  cc  <-  processing_setting(settings)
  s <- cc[[1]]
  data.in <- cc[[2]]
  id.out  <- cc[[3]]
  id.ori  <- NULL
  
  if(is.null(data)){
    
    #--------------------------------------------------
    #    MODEL
    #--------------------------------------------------
    if (!(is.null(model))) {
      if(model=="pkmodel"){
        model = generateModelFromPkModel(lv$parameter[[1]],lv$output[[1]]) 
      } else {
        model_ext <- file_ext(model)
        if(model_ext=="xml"){
          model = pharmml2mlxtran(model)
        }
      }
    }
    
    
    dataIn <- list()    
    
    iop.group <- 1
    if (is.null(lv$group))
      iop.group <- 0
    #--------------------------------------------------
    
    #     doseRegimen <- lv$treatment
    lv$model <- model
    id.ori <- lv$id
    gr.ori <- lv$gr.ori
    lv  <- hformat(lv)
    
    dataIn  <-  hgdata(lv)
    dataIn$model <- model
    dataIn$trt <- lv$treatment
    
    if(data.in==T) {
      dataIn$iop.group <- iop.group
      dataIn$id.ori <- id.ori
      s$loadDesign <- TRUE
      #       return(dataIn)
    }    
  }else{
    dataIn <- data
    iop.group <- data$iop.group
    dataIn$iop.group <- NULL
    id.ori <- data$id.ori
    dataIn$id.ori <- NULL
  }
  
  
  
  if (out.trt==T)
    trt <- dataIn$trt
  else
    trt <- NULL
  
  # dataIn$trt <- NULL
  if (length(s)==0){
    argList <- list(DATA=dataIn) 
  } else {
    argList <- list(DATA=dataIn, SETTINGS=s)       
  }
  
  if (identical(file_ext(model),"R")) {Rfile <- TRUE} else {Rfile <- FALSE}
  if ( !is.null(model) && exists(model, mode="function") ){Rsource <- TRUE} else {Rsource <- FALSE}
  if (Rfile || Rsource){
    dataOut <- simulR(argList)
    return(dataOut)
  }else{
    dot_call <- .Call
    dataOut  <- dot_call( "mlxComputeR", argList, PACKAGE = "mlxComputeR" )
    if(data.in==T)
      return(dataIn)
    if (!exists("gr.ori"))
      gr.ori <- NULL
    dataOut  <- convertmlx(dataOut,dataIn,trt,iop.group,id.out,id.ori,gr.ori)
    if(!is.null(id.ori) && is.data.frame(id.ori))
      dataOut$originalId <- id.ori
    return(dataOut)
  }
}

mergeres <- function(r,s,m=NULL,N=NULL){
  K <- length(s)
  if (is.null(r)){
    if (!is.null(m)){
      for (k in (1:K))
        if ("id" %in% names(s[[k]])){
          s[[k]]$group <- factor(m) 
          i1 <- which(names(s[[k]])=="id")
          i2 <- which(names(s[[k]])=="group")
          s[[k]] <- s[[k]][,c((1:i1),i2,((i1+1):(i2-1)))]
          if (!is.null(N))
            s[[k]] <- s[[k]][(1:N),]
        }
    }
    u <- s
  }else{
    u <- r
    for (k in (1:K)){
      if (is.data.frame(s[[k]])){
        if (!is.null(m))
          s[[k]]$group <- factor(m)      
        
        if ("id" %in% names(s[[k]])){
          skid <- as.numeric(as.character(s[[k]]$id))
          if (!is.null(N))
            ikN <-which(skid <= N)
          ir <- as.numeric(tail(r[[k]]$id,1))
          skid <- skid + ir
          s[[k]]$id <- factor(skid) 
          if (!is.null(N))
            s[[k]] <- s[[k]][ikN,]
        }
        u[[k]] <- rbind(r[[k]],s[[k]])
        attr(u[[k]],"name")=attr(s[[k]],"name")
      }
    }
  } 
  return(u)
}


sim.pop <- function(n,mu,infop,fim,kw.max)
{
  p1.name <- names(fim$se)
  np <- length(p1.name)
  p1.trans=rep("N",np)
  p2.name <- sub("_pop","",p1.name)
  i.pop <- match(infop$name,p2.name)
  i1 <- which(!is.na(i.pop))
  p1.trans[i.pop[i1]] <- infop$trans[i1]
  i.omega <- grep("omega_",p1.name)
  p1.trans[i.omega] <- "L"
  i.omega2 <- grep("omega2_",p1.name)
  p1.trans[i.omega2] <- "L"
  param <- data.frame(pop.param=mu,sd=fim$se,trans=p1.trans)
  if (is.null(kw.max))
    x <- simpopmlx(n=n,parameter=param,corr=fim$mat)
  else
    x <- simpopmlx(n=n,parameter=param,corr=fim$mat,kw.max=kw.max)
  return(x)
}

