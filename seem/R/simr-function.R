simr <- function(model, file, param=NULL, pdfout=F, lab.out="") {

 # multiple runs
 # required: model as a list
 #           path/name of input file
 # optional: param set for multiple runs 
 #           switch to produce pdf files, default F
 #           label to identify output files, deafult none

 # read input, file format is csv 
 dat <- c("character","numeric","character","character")
 input <- read.csv(file,header=T,colClasses=dat)
 # input values and labels
 v <- t(input$Val); vl <- t(input$Lab)
 # t0,tf,dt,tw are in v[1:4]
 t <- seq(v[1],v[2],v[4]); dt <- v[3]
 # position of first initial condition
 iNS <- max(which((vl=="X1.0")),which((vl=="X0")))
 last <- length(vl); NS <- last-iNS # number of states
 mp <- c(5:(iNS-1)); mX0 <- c(iNS:(iNS+NS-1)); digX <- v[last]
 # units ind and dep vars
 if(NS>1){
 unit <- c(input$Uni[1], rep(input$Uni[iNS],NS))
 tXlab <- paste(c("t",paste("X",seq(1:NS),sep="")),"[",unit,"]", sep="")
 } else {
 unit <- c(input$Uni[1], input$Uni[iNS],
          paste("(",input$Uni[iNS],")/",input$Uni[1],sep=""))
 tXlab <- paste(c("t","X","dX/dt"),"[",unit,"]", sep="")
}
# number of time stamps and runs  
  nt <- length(t)
  if(is.null(param)) np <- 1 else {
  if(length(param$plab)==1) param$pval<- matrix(param$pval)
  np <- dim(param$pval)[1]
  }
# arrays to store results
  X <-  structure(1:(np*nt*NS),dim=c(nt,np,NS))
  if(NS==1) dX.dt <- matrix(nrow=nt,ncol=np)
# loop for parameter p that varies
  for(i in 1:np) { 
     # pick a value of X0 or p from the param set
     if(is.null(param)==F){
      for(k in 1:length(param$plab)){
      for(j in mp[1]:mX0[length(mX0)])
        if(param$plab[k]==vl[j]) v[j]<- param$pval[i,k]
      }
     }
     # parameters
     p <- v[mp[1]:mp[length(mp)]]
     # X0 initial condition 
     X0 <- v[mX0[1]:mX0[length(mX0)]]
     # integration (edit RK4 for ramos)
     # out <- RK4(X0, t, model$f, p, dt)
     out <- ramos(X0, t, model$f, p, dt)
     # rename x part of out for simpler use
     if(NS>1) for(k in 1:NS) X[,i,k] <- out[,k]
     # calculate rate for one-state systems
     else{
      X[,i,1] <- out
      for(j in 1:nt) dX.dt[j,i] <- model$f(t[j],p,X[j,i,1])[1]
     }
  } # end of parameter p loop 

 # prep and organize output as data frame
  if(NS==1) X <- cbind(X[,,1],dX.dt)
  X <- signif(X,digX)
  output <- data.frame(t,X)

# forms names for output vars and runs 
  nv <- length(tXlab)
  var.run <- c(tXlab[1], paste(tXlab[2],paste(".Run",c(1:np),sep=""),sep=""))
  for(i in 3:nv) var.run <- c(var.run, paste(tXlab[i],paste(".Run",c(1:np),sep=""),sep="")) 
  names(output) <- var.run

 # call output function to generate output files
  vars.out(file, lab.out, input, param, output, tXlab, pdfout)

  return(list(input=input, param=param, output=output))
} # end of function

