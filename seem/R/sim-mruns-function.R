sim.mruns <- function(model, file, param, pdfout=F, lab.out="") {

 # multiple runs
 # required: model as a list
 #           path/name of input file
 #           param set for multiple runs 
 # optional: switch to produce pdf files, default F
 #           label to identify output files, default none
 
 # read input, file format is csv 
 dat <- c("character","numeric","character","character")
 input <- read.csv(file,header=T,colClasses=dat)
 # input values and labels
 v <- t(input$Val); vl <- t(input$Lab)

 # t0,tf,dt,tw are in v[1:4]
 t <- seq(v[1],v[2],v[4]); dt <- v[3]
 digX <- v[7]
 # units ind and dep vars
 unit <- c(input$Uni[1], input$Uni[6])

 # number of time stamps and param values 
  nt <- length(t); np <- length(param$pval)
 # array to store results
  X <- matrix(ncol=np,nrow=nt)
 # loop for parameter p
  for(i in 1:np) { 
   # pick a value of X0 or p from param 
    for(j in 5:6) if(param$plab==vl[j]) v[j]<- param$pval[i]
   # p,X0, in v[5:6]
    p <- v[5]; X0 <- v[6]
   # integration
     X[,i] <- RK4(X0, t, model$f, p, dt)
  } # end of parameter p loop

 # prep and organize output
  X <- signif(X,digX); output <- data.frame(t,X)
  tXlab <- paste(c("t","X"),"[",unit,"]", sep="")
  var.run <- c(tXlab[1], paste(tXlab[2],paste(".Run",c(1:np),sep=""),sep=""))
  names(output) <- var.run

 # call output function to generate output files
 mruns.out(prefix=file, lab.out=param$plab, input, param, output, tXlab, pdfout)

 return(list(input=input, param=param, output=output))

} # end of function

