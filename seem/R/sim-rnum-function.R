sim.rnum <- function(model, file, pdfout=F,lab.out=""){

 # required: model as a list
 #           path/name of input file
 # optional: switch to produce pdf files, default F
 #           label to identify output files, deafult none
 
 # read input, file format is csv 
 dat <- c("character","numeric","character","character")
 input <- read.csv(file,header=T,colClasses=dat)
 # input values and labels
 v <- t(input$Val); vl <- t(input$Lab)

 # t0,tf,dt,tw are in v[1:4]
 t <- seq(v[1],v[2],v[4]); dt <- v[3]
 # mu,sd,n in v[5:7]
 p <- c(v[5],v[6]); n <- v[7]; X0 <- v[8]; digX<-v[9]
 unit <- c(input$Uni[1], rep(input$Uni[8],2))

 # call integration
 X <- rnum(X0, t, model$f, p, dt, n)
 # organize output
 X <- signif(X,digX); output <- data.frame(t,X)  
 names(output) <- paste(names(output),"[",unit,"]", sep="")
 tXlab <- paste(c("t","X"),"[",unit,"]", sep="")
 # generate output files
 rnum.out(prefix=file, lab.out, input, output, tXlab, pdfout)

 return(list(input=input,output=output))

} # end of function


