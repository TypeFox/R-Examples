example.as.ESS.object <- function(){
path.input <- system.file("Input", package="R2GUESS")
path.output <- system.file("Output", package="R2GUESS")
path.par <- system.file("extdata", package="R2GUESS")
file.par <- "Par_file_example_Hopx.xml"
root.file.output <- "Example-GUESS-Y-Hopx"
label.Y <- c("ADR","Fat","Heart","Kidney")
my.env <- new.env()
data(MAP.file,envir=my.env)
MAP.file <- my.env$MAP.file
modelY_Hopx <-as.ESS.object(dataY="data-Y-ALL-C-CODE.txt",dataX="data-X-C-CODE.txt",path.input=path.input,
                            path.output=path.output,root.file.output=root.file.output,label.X=NULL,
                            label.Y=label.Y,path.par=path.par,file.par=file.par,MAP.file=MAP.file)
return(modelY_Hopx)
}
