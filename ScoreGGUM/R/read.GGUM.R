read.GGUM <-
function(file,
           numItems,
           numCats,
           model = 8) {

### Import data from GGUM item parameter output file
  
ipar.list <- scan(file, skip=4, what=list( NAME=character(0),VALUE=numeric(0)))
  
ipar <- matrix(0,numItems,max(numCats)+3)
ipar.new <- matrix(0,numItems,max(numCats)+3)
pos <- 0

# Stop program and report error if input values are not correct dimensions

if (length(numCats)==1) {numCats <- rep(numCats,numItems)}
if (length(numCats)!=numItems) {
  stop("ERROR: length of numCats not equal to numItems") }
model <- as.integer(model)
if ((length(model)!=1) | (model < 1) | (model > 8)) {
  stop("ERROR: model must be an integer between 1 and 8") }

### Model 8

if (model==8){
ref <- matrix(c(1,3,5,10+0:9*4),ncol=1)
for (i in 1:numItems) {
  ipar[i,1:(numCats[i]+3)] <- pos+ref[1:(numCats[i]+3),1]
  ipar.new[i,1:(numCats[i]+3)] <- ipar.list$VALUE[ipar[i,1:(numCats[i]+3)]]
  pos <- pos + ref[(numCats[i]+3),1]+1
  }
}

### Model 7 

if (model==7){
ref <- matrix(c(1,3,5,7*numItems+0:9*4+3),ncol=1)
for (i in 1:numItems) {
  ipar[i,1:3] <- pos+ref[1:3,1]
  ipar[i,4:(numCats[i]+3)] <- ref[4:(numCats[i]+3),1]
  ipar.new[i,1:(numCats[i]+3)] <- ipar.list$VALUE[ipar[i,1:(numCats[i]+3)]]
  pos <- pos + 7
  }
}

### Models 6 and 5

if (model==6 | model==5){
ref <- matrix(c(1,3,5,rep(9,10)),ncol=1)
for (i in 1:numItems) {
  ipar[i,1:(numCats[i]+3)] <- pos+ref[1:(numCats[i]+3),1]
  ipar.new[i,1:(numCats[i]+3)] <- ipar.list$VALUE[ipar[i,1:(numCats[i]+3)]]
  pos <- pos + 10
  thresh <- numCats[i]
  ipar.new[i,4] <- 0
  for (j in 4:(numCats[i]+3)) {
    ipar.new[i,j] <- ipar.new[i,j]*(-2)*thresh
    thresh = thresh - 1
    }
  }
}

### Model 4

if (model==4){ 
ref <- matrix(c(1,3,7+0:9*4),ncol=1)
for (i in 1:numItems) {
  ipar[i,c(1:2,4:(numCats[i]+3))] <- pos+ref[1:(numCats[i]+2),1]
  ipar.new[i,c(1:2,4:(numCats[i]+3))] <- ipar.list$VALUE[ipar[i,1:(numCats[i]+3)]]
  pos <- pos + ref[(numCats[i]+2),1]+1
  }
ipar.new[,3] <- 1
}

### Model 3 

if (model==3){ 
ref <- matrix(c(1,3,4*numItems+0:9*4+3),ncol=1)
for (i in 1:numItems) {
  ipar[i,1:2] <- pos+ref[1:2,1]
  ipar[i,4:(numCats[i]+3)] <- ref[3:(numCats[i]+2),1]
  ipar.new[i,c(1:2,4:(numCats[i]+3))] <- ipar.list$VALUE[ipar[i,1:(numCats[i]+3)]]
  pos <- pos + 4
  }
ipar.new[,3] <- 1
}

### Models 2 and 1

if (model==2 | model==1){ 
ref <- matrix(c(1,3,rep(6,10)),ncol=1)
for (i in 1:numItems) {
  ipar[i,c(1:2,4:(numCats[i]+3))] <- pos+ref[1:(numCats[i]+2),1]
  ipar.new[i,c(1:2,4:(numCats[i]+3))] <- ipar.list$VALUE[ipar[i,1:(numCats[i]+3)]]
  pos <- pos + 7 
  thresh <- numCats[i]
  ipar.new[i,4] <- 0
  for (j in 4:(numCats[i]+3)) {
    ipar.new[i,j] <- ipar.new[i,j]*(-2)*thresh
    thresh = thresh - 1
    }
  }
ipar.new[,3] <- 1
}

### Name columns of item parameter file

ipar.names <- c("ITEM","DELTA","ALPHA","TAU0","TAU1","TAU2","TAU3","TAU4","TAU5","TAU6","TAU7","TAU8","TAU9","TAU10")
colnames(ipar.new)[1:(max(numCats)+3)] <- ipar.names[1:(numCats[i]+3)]

### Return I x (C+3) item parameter matrix, where I is the number of items and C is the number of response categories 

return(ipar.new)

}
