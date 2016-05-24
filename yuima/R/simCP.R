simCP<-function(xinit,yuima,env){

	
	sdeModel<-yuima@model
	
	modelstate <- sdeModel@solve.variable
	modeltime <- sdeModel@time.variable
	Terminal <- yuima@sampling@Terminal[1]
    Initial <- yuima@sampling@Initial[1]
	dimension <- yuima@model@dimension
    dummy.val <- numeric(dimension)
    if(length(xinit) !=  dimension)
    xinit <- rep(xinit,  dimension)[1:dimension]
    if(length(unique(as.character(xinit)))==1 &&
       is.numeric(tryCatch(eval(xinit[1],envir=env),error=function(...) FALSE))){
           dX_dummy<-xinit[1]
           dummy.val<-eval(dX_dummy, envir=env)
           if(length(dummy.val)==1){
               dummy.val<-rep(dummy.val,dimension)
           }
           for(i in 1:length(modelstate)){
                assign(modelstate[i],dummy.val[i] ,envir=env)
           }

       } else {
           for(i in 1:dimension){
               dummy.val[i] <- eval(xinit[i], envir=env)
           }
           
       }

###    Simulation of CP using Lewis' method

##:: Levy
    JP <- eval(sdeModel@jump.coeff[[1]], envir=env)
    mu.size <- length(JP)
    #  print(str(JP))
   
   #assign(sdeModel@measure$intensity, env) ## intensity param
   .CPintensity <- function(.t) {
       assign(modeltime, .t, envir=env)
       eval(sdeModel@measure$intensity, envir=env)
   }
   
   
      dummyList<-as.list(env)
      
      lgth.meas<-length(yuima@model@parameter@measure)
      if(lgth.meas>1){
        for(i in c(2:lgth.meas)){
          idx.dummy<-yuima@model@parameter@measure[i]
          assign(idx.dummy,as.numeric(dummyList[idx.dummy]))
        }
      }
      
      # we use Lewis' acceptance/rejection method
      
      #if(grep("^[dexp|dnorm|dgamma|dconst]", sdeModel@measure$df$expr)){
          ##:: e.g. dnorm(z,1,1) -> rnorm(mu.size*N_sharp,1,1)
          F <- suppressWarnings(parse(text=gsub("^d(.+?)\\(.+?,", "r\\1(mu.size*N_sharp,", sdeModel@measure$df$expr, perl=TRUE)))
          #} else{
          #stop("Sorry. CP only supports dconst, dexp, dnorm and dgamma yet.")
          #}
      
      ell <- optimize(f=.CPintensity, interval=c(Initial, Terminal), maximum = TRUE)$objective
      ellMax <- ell * 1.01


      time <- Initial
      E <- Initial
      
      while(time < Terminal) {
          time <- time - 1/ellMax * log(runif(1))
          if(runif(1) < .CPintensity(time)/ellMax)
          E <- c(E, time)
      }
      N_sharp <- length(E)-1
      
      F.env <- new.env(parent=env)
      assign("mu.size", mu.size, envir=F.env)
      assign("N_sharp", N_sharp, envir=F.env)
      randJ <- eval(F, envir=F.env)  ## this expression is evaluated in the F.env
   
      randJ <- JP[1]*randJ
      randJ <- as.matrix(randJ, ncol=yuima@dimension)
  
      randJ <- rbind(dummy.val, randJ)
      CP <- apply(randJ,2,cumsum)
      tsX <- zoo(x=CP, order.by=E)
      yuimaData <- setYuima(data=setData(tsX))
      yuimaData <- subsampling(yuimaData, sampling=yuima@sampling)
      return(yuimaData@data)
}
