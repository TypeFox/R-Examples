##load the data
data(mesa.model)

##create a data matrix
M1 <- createDataMatrix(mesa.model)
dim(M1)
head(M1)

##create data matrix for only a few locations
M2 <- createDataMatrix(mesa.model, subset =
                         c("60370002","60370016","60370113","60371002",
                           "60371103","60371201","L001","L002"))
dim(M2)
head(M2)
\dontshow{
  if( (dim(M1)[1]!=dim(mesa.model$trend)[1]) ||
     (dim(M1)[2]!=dim(mesa.model$locations)[1]) ){
    stop("createDataMatrix: dimension missmatch - M1")
  }
  if( (dim(M2)[1]!=dim(mesa.model$trend)[1]) || (dim(M2)[2]!=8) ){
    stop("createDataMatrix: dimension missmatch - M2")
  }
#  if( isTRUE( all.equal(M1[,colnames(M2)],M2) ) ){
#    stop("createDataMatrix: M1!=M2")
#  }
}

