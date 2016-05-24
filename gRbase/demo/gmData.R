## ###################################################################
## examples:

## add new type of variable
oldtypes <- validVarTypes()
validVarTypes <- function() {
  c(oldtypes,"MyVarType")
}


newgmData(1:3,"MyVarType")

mydata <- newgmData(paste("x",1:6,sep=""),"Discrete",2,NA,list(x3=c("lab1","lab2")))

valueLabels(mydata)

data(iris)
mydata <- as.gmData(iris)
mydata

observations(mydata) <- observations(mydata)[1:10,]

mydata$mimName <- letters[1:nrow(mydata)]

varNames(mydata)
varTypes(mydata)
nLevels(mydata)
valueLabels(mydata)

varTypes(mydata)[3]     <- "Discrete"
nLevels(mydata)[3] <- 2


data(rats)
gmd.rats <- as.gmData(rats)
data(HairEyeColor)
gmd.hec  <- as.gmData(HairEyeColor)

gmd.rats.nodata  <-  newgmData(
                         varNames=c("Sex","Drug","W1","W2"),
                         varTypes=c("Discrete","Discrete","Continuous","Continuous"),
                         nLevels=c(2,3,NA,NA), 
                         valueLabels=list(Sex=c("M","F"), Drug=c("D1","D2","D3")))
observations(gmd.rats.nodata) <- rats
valueLabels(gmd.rats.nodata) <- list(Sex=c("Male","Female"))

nVar <- nrow(  gmd.rats.nodata ) 
gmd.rats.nodata$mimName <-  letters[1:nVar]


