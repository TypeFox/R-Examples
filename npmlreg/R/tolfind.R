"tolfind" <-
function(formula,
                    random=~1,
                    family = gaussian(),
                    data ,
                    k = 4,
                    random.distribution="np",
                    offset,
                    weights,
                    na.action,
                    EMmaxit=500,
                    EMdev.change=0.001,
                    lambda=0,
                    damp=TRUE,
                    damp.power=1,
                    spike.protect=1,
                    sdev,
                    shape,
                    plot.opt=1,
                    steps=15,
                    find.in.range=c(0.05,0.8),
                    verbose=FALSE,
                    noformat=FALSE,
                    ...)

{
  # JE/JH, 2005/06.

  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }
  data <- as.data.frame(data)
  ddim <- dim(data)
  mf   <- match.call(expand.dots = FALSE)
  
  # Test for k>1 and for inadmissibly removed intercept term
  if (k==1){
      stop("Please choose k > 1.")
  } else if (k>1 && random.distribution=='np' && max(length(grep('- 1', deparse(formula(mf)))),length(grep('-1', deparse(formula(mf))))) >0 ){
      stop(" term '-1'  in model formula not supported for k>1 & random.distribution='np'. ")
  }
  
  # Test for incorrect offset specification in formula object
  testoffset <- try(is.null(attr(terms(formula(mf)),"offset")),silent=TRUE)
  if (!(class(testoffset)=="try-error" || testoffset)){
      stop("Please specify offset as separate argument outside the model formula.")
  }
  
  # Extract variables from call and set up initial (fixed effect) model
  m    <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0)
  mf   <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- try(eval(mf,  parent.frame()), silent=TRUE)
  if (class(mf)=="try-error"){
      if (!missing(offset) && length(offset) != ddim[1]) {
      stop("Number of offsets is ", length(offset), ", should equal ", ddim[1], " (number of observations)")
      }
      if (!missing(weights) && length(weights) != ddim[1]){
      stop("Number of weights is ", length(weights), ", should equal ", ddim[1], " (number of observations)")
      }
      stop(geterrmessage())
  }
  
  offset  <- model.offset(mf) 
  weights <- model.weights(mf) 
  Y  <- model.response(mf, "numeric") # response
  Ym <- is.matrix(Y)
  N  <- NROW(Y)  # corresponds to ddim[1] if there are no missing values
  
  if (is.null(offset)){offset  <-rep(0,N) }
  if (is.null(weights)){weights<-rep(1,N)}   
  data$poffset <- numeric(1); data$pweights<-numeric(1)
  data  <- if (is.matrix(Y)) data[dimnames(Y)[[1]],] else  data[names(Y),] # omit missing values   
  poffset<-data$poffset<-offset;  pweights<-data$pweights<-weights # global binding only to avoid R CMD Note.
  
  all.Disparities <- all.converged <- rep(0,steps)
  min.Disp        <- min.Disp.conv <- 10^8
  step.min        <- step.min.conv <- 1
  mform           <- strsplit(as.character(random)[2],'\\|')[[1]]
  mform           <- gsub(' ', '',mform)
  
  # Configure graphics window
  if (!noformat){
      if (steps >8) par(mfrow=c(4,4), cex=0.5) else par(mfrow=c(3,3), cex=0.5,cex.axis=1.1 )
  }

  # Run alldist/allvc for a grid of tol values
  for (t in 1: steps){
     tol  <- find.in.range[1]+ (find.in.range[2]-find.in.range[1])*t/steps
     tol0 <- tol  #for main title of graphical output
     if(length(mform)==1){ 
         npfit<- try(alldist(formula=formula,random=formula(random),family=family, data=data,k=k,random.distribution=random.distribution,  offset=poffset,  weights=pweights, na.action=na.action,  tol=tol, EMmaxit=EMmaxit, EMdev.change=EMdev.change, lambda=lambda, damp=damp, damp.power=damp.power, spike.protect=spike.protect, sdev=sdev, shape=shape, plot.opt=plot.opt, verbose=verbose))         
         if(class(npfit)=="try-error"){
                cat("tolfind failed using tol=", tol, ". Hint:  specify another range of tol values and try again. "); return()
         }  
         all.Disparities[t]<-  npfit$disparity
         all.converged[t]<-  npfit$EMconverged
     } else { 
         vcfit<- try(allvc(formula=formula,random=formula(random),family=family, data=data,k=k,random.distribution=random.distribution,  offset=poffset, weights=pweights, na.action=na.action, tol=tol, EMmaxit=EMmaxit,EMdev.change=EMdev.change, lambda=lambda, damp=damp, damp.power=damp.power, spike.protect=spike.protect, shape=shape, sdev=sdev, plot.opt=plot.opt, verbose=verbose))
         if(class(vcfit)=="try-error"){
                cat("tolfind failed using tol=", tol,  ". Hint: specify another range of tol values and try again. "); return()
         } 
         all.Disparities[t]<-  vcfit$disparity
         all.converged[t]<-vcfit$EMconverged
     }
     
     if (all.Disparities[t]< min.Disp){min.Disp<- all.Disparities[t]; step.min<-t }
     if (all.Disparities[t]< min.Disp.conv && all.converged[t]){min.Disp.conv<- all.Disparities[t]; step.min.conv<-t }
  }
  
  tol.min       <- find.in.range[1] + (find.in.range[2]-find.in.range[1])*step.min/steps
  tol.min.conv  <- find.in.range[1] + (find.in.range[2]-find.in.range[1])*step.min.conv/steps
  npcolors      <- 2 + all.converged

  plot(find.in.range[1]+(find.in.range[2]-find.in.range[1])*(1:steps)/steps,all.Disparities,type="o",xlab="tol",ylab="Disparity",col=npcolors)
  segments(tol.min, min.Disp, tol.min, 1.1*min.Disp,col=4)
  cat("Minimal Disparity:", min.Disp, "at tol=", tol.min, "\n")

  # Print on screen
  if(max(all.converged)==1){
      cat("Minimal Disparity with EM converged:", min.Disp.conv, "at tol=", tol.min.conv, "\n")
  } else {
      cat(" No convergence achieved for any choice of tol.", "\n")
  }

  list("MinDisparity"=min.Disp.conv, "Mintol"=tol.min.conv, "AllDisparities"= all.Disparities , "Alltol"=  find.in.range[1]+ (find.in.range[2]-find.in.range[1])* (1:steps)/steps, "AllEMconverged"=all.converged==TRUE)
}

