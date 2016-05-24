plot.trtsel <-
function(x, bootstraps = 500,  
            plot.type = c("risk", "treatment effect", "cdf", "selection impact"), 
            ci = "default", 
            alpha = .05, 
            fixed.values = NULL,  
            offset = 0.01, 
            conf.bands = TRUE, 
            conf.bandsN = 100, 
            trt.names = c("Treatment", "    No \nTreatment"),
            xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL, show.marker.axis = TRUE, 
            ...)

{
 
  if(!is.trtsel(x)) stop("x must be an object of class 'trtsel' created by using the function 'trtsel' see ?trtsel for more help")
  
  plot.type <- match.arg(plot.type) 
  
  if(!is.element(plot.type, c("risk", "treatment effect", "cdf", "selection impact"))){ 
    stop("plot.type must be one of \"risk\", \"treatment effect\", \"cdf\" or \"selection impact\"")
  }
  stopifnot(length(plot.type) ==1)
  
  if(!is.element(ci, c("default", "horizontal", "vertical", "none"))){ 
    stop("ci must be one of \"default\",  \"horizontal\", or  \"vertical\" or \"none\" ")
  }
  

  if(alpha<0 | alpha > 1) stop("Error: alpha should be between 0 and 1")
  if(bootstraps < 2) warning("Number of bootstraps must be greater than 1, bootstrap confidence intervals will not be computed") 
  if(ci == "none") bootstraps = 0; 
  # theta curves can only be plotted for cohort data and with vertical ci's
  if(substr(plot.type, 1, 3)=="sel"){
   # if(substr(x$model.fit$study.design, 1, 3) != "ran") stop("theta curves cannot be created for subcohort designs")
    if(substr(ci, 1, 1) =="h") {
      
      warning("horizontal ci's are not available for selection impact curves, vertical ci's will be calculated")
      ci = "vertical"
    }
  }
  #set default ci's 
  
  if(ci =="default"){
    #continuous marker
    if(is.null(x$model.fit$disc.marker.neg)){
      if(substr(plot.type, 1, 3) =="ris") ci = "horizontal"
      if(substr(plot.type, 1, 3) =="tre") ci = "horizontal"
      if(substr(plot.type, 1, 3) =="cdf") ci = "vertical"
      if(substr(plot.type, 1, 3) =="sel") ci = "vertical"
    }else{
      
      if(substr(plot.type, 1, 3) =="ris") ci = "vertical"
      if(substr(plot.type, 1, 3) =="tre") ci = "vertical"
      if(substr(plot.type, 1, 3) =="cdf") ci = "horizontal"
    }

  
  }
  #no cdf plots for binary marker
  if(!is.null(x$model.fit$disc.marker.neg)){
  if(is.element(substr(plot.type, 1, 3), c("cdf", "sel"))) stop("cdf or selection impact plots cannot be created for a binary marker. Please choose plot.type to be \"risk\" or \"treatment effect\" ")
  }
  #save the current plot parameters
  #old.par <- par(no.readonly = TRUE)
 

#extract the needed data from x, which is our TrtSel object

  trt <- x$derived.data$trt
  event <- x$derived.data$event

  rho  <- x$model.fit$cohort.attributes
  if(length(rho) == 4) rho = c(rho, -9999, -9999, -9999) #accomodate the nested case-control design
  
  study.design <- x$model.fit$study.design
  link <- x$model.fit$link
  boot.sample <- x$functions$boot.sample
  get.F <- x$functions$get.F
  delta <- x$derived.data$trt.effect  

  if(link == "risks_provided") {
    marker = NULL
    show.marker.axis = FALSE
  }else{
    marker <- x$derived.data$marker 
  }

 if(is.null(x$model.fit$disc.marker.neg)){  #continuous marker 
    plot.functions <- list(  predcurvePLOT_gg, trteffectPLOT_gg, CDFdeltaPLOT_gg, SelectionImpactPLOT_gg)
    
  if(length(fixed.values)!=0) conf.bands = FALSE
  if(conf.bands & length(fixed.values)==0 ){
   #create fixed values

   if(substr(ci, 1,1 )=="v"){
       
      if(is.element(substr(plot.type, 1,3), c("tre", "ris", "sel"))) fixed.values = seq(from = 1, to = 100, length.out = conf.bandsN)
      else if(substr(plot.type, 1,3)=="cdf") fixed.values = seq(from = min(delta), to = max(delta), length.out = conf.bandsN)
      
   }else{
      
      if(is.element(substr(plot.type, 1,3), c("cdf"))) fixed.values = seq(from = 1, to = 100, length.out = conf.bandsN)
      else if(substr(plot.type, 1,3)=="tre") fixed.values = seq(from = min(delta), to = max(delta), length.out = conf.bandsN)
      else if(substr(plot.type, 1,3)=="ris"){
       allrisks <- c(x$derived.data$fittedrisk.t0, 
                     x$derived.data$fittedrisk.t1)
       fixed.values = seq(from = min(allrisks), to = max(allrisks), length.out = conf.bandsN)

      }
      ## add case for theta curve

   }
   offset = 0

   }
    ##Bootstrap for confidence intervals...
    #bootstrapping done by fixing response and strapping marker 
    
    if((conf.bands & bootstraps>1) | (length(fixed.values)>0 & bootstraps > 1)){ 

      ci.bounds <- get.plot.ci(x, 
                               plot.type, 
                               ci, bootstraps, 
                               fixed.values =fixed.values,  
                               alpha = alpha)
    }else{ 
      ci.bounds <- NULL
      conf.bands = FALSE
    }
    ## end Bootstrap
    
    
    
    
  }else{
    #discrete marker...ignore fixed values and ci bands, we only need ci's around the observed points
    
    plot.functions <- list(  predcurvePLOT_gg_disc, trteffectPLOT_gg_disc, CDFdeltaPLOT_gg_disc)

    ##Bootstrap for confidence intervals...
    #much simpler for a discrete marker

    if(( conf.bands & bootstraps>1) ){ 
      
      ci.bounds <- get.plot.ci_disc(x, plot.type, 
                                    ci, bootstraps, 
                                    alpha)
      ci = ci.bounds$newci
      ci.bounds = ci.bounds$myconf.ints
      }else{ 
      ci.bounds <- NULL
      conf.bands = FALSE
    }
    ## end Bootstrap
    
    
  }

  tmp.plotfun <- plot.functions[[match(plot.type, c("risk", "treatment effect", "cdf", "selection impact"))]]
   

if(is.null(x$model.fit$disc.marker.neg)){  # continous marker 
  if(substring(plot.type, 1, 4) == "risk"){

  curves <- tmp.plotfun(x, ci, ci.bounds, get.F, fixed.values, conf.bands,  rho, trt.names, xlab, ylab, xlim, ylim, main, show.marker.axis, offset = offset)

    if(!is.null(ci.bounds)){
      ci.bounds <- data.frame(t(ci.bounds))
      ci.bounds <- cbind(fixed.values, ci.bounds)
      names(ci.bounds) <- c("fixed.values", "trt0.lower", "trt0.upper", "trt1.lower", "trt1.upper")
      }
    }else{

      curves <- tmp.plotfun(x, ci, ci.bounds, get.F, fixed.values, conf.bands, rho, xlab, ylab, xlim, ylim, main)
      if(!is.null(ci.bounds)){
        ci.bounds <- data.frame(t(ci.bounds))
        ci.bounds <- cbind(fixed.values, ci.bounds)
        names(ci.bounds) <- c("fixed.values", "lower", "upper")
      }
    }
  
}else{  # discrete marker
  

    
    curves <- tmp.plotfun(x=x, ci = ci, ci.bounds = ci.bounds, get.F = get.F, 
                          xlab = xlab, ylab=ylab, 
                          xlim = xlim, ylim=ylim, 
                          main=main, trt.names = trt.names)

    if(!is.null(ci.bounds)){
      ci.bounds <- data.frame(curves[[2]])

    }
   
}    

  #par(old.par)

invisible(list("plot" = curves[[1]], "ci.bounds" = ci.bounds))
}
