sim.comp <- function(model, file, method= "RK4", pdfout=F, lab.out=""){

 # required: model as a list
 #           path/name of input file
 # optional: name of function for num method, default RK4 
 #           switch to produce pdf files, default F
 #           label to identify output files, default none
       
 # read input, file format is csv 
 dat <- c("character","numeric","character","character")
 input <- read.csv(file,header=T,colClasses=dat)
 # input values and labels
 v <- t(input$Val); vl <- t(input$Lab)

 # t0,tf,dt,tw are in v[1:4]
 t <- seq(v[1],v[2],v[4]); dt <- v[3]
 # p,X0, digX in v[5:7]
 p <- v[5]; X0 <- v[6]; digX <- v[7]
 # units ind and dep vars
 unit <- c(input$Uni[1], input$Uni[6])

 # call integration
 if(method=="euler")  X <- euler(X0, t, model$f, p, dt)
 else  X <- RK4(X0, t, model$f, p, dt)
 # organize output
 X <- signif(X,digX); output <- data.frame(t,X)  
 names(output) <- paste(names(output),"[",unit,"]", sep="")
 # generate output files
 onerun.out(prefix=file, lab.out, input, output, pdfout)

 return(list(input=input,output=output))
} # end of function


