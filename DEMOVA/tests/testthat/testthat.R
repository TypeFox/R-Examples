library(testthat)
library(leaps)
library(DEMOVA)


test_that("Testing DEMOVA", {

property<-"NameProperty"
nom<-"trainingSet.csv"
data<-read.csv(nom , header = TRUE , sep=" ") 

dim<-dim(data)
desc<-data[,3:dim[2]]
id<-data[,1]
y<-data[,2] 

d<-preselection(desc)
expect_is(d,"data.frame")

desc<-select_variables(id,y,d,0.99,auto=TRUE)
expect_is(desc,"data.frame")
expect_true(file.exists("variables_selected.csv"))

MLR<-select_MLR(y,desc,5)
expect_is(MLR,"data.frame")
expect_true(file.exists("_MLR"))

mydata<-cbind(y,MLR)

fit<-fitting(mydata,dim(MLR)[2],property)
expect_is(fit,"lm")
file<-paste('prediction_TrainSet_',property,'.csv',sep="")
expect_true(file.exists(file))
file<-paste(property,'_TrainingSet','.tiff',sep="")
expect_true(file.exists(file))

q<-LOO(mydata,dim(MLR)[2])
expect_is(q,"numeric")
expect_less_than(q, 1)

q<-LMO(mydata,5,dim(MLR)[2])
expect_is(q,"numeric")
expect_less_than(q, 1)

q<-LMO(mydata,10,dim(MLR)[2])
expect_is(q,"numeric")
expect_less_than(q, 1)

YS<-scramb(mydata,1000,dim(MLR)[2],cercle=TRUE)
expect_true(file.exists("Scramb.csv"))
expect_true(file.exists("Scramb.tiff"))
expect_is(YS,"numeric")
expect_true(is.vector(YS))
expect_equal(length(YS),2)
expect_less_than(YS[1], 1)

new_nom<-"testSet.csv"
newdata<-read.csv(new_nom,header=TRUE , sep=" ")
mynewdata<-newdata[,2:dim(newdata)[2]]

Rext<-prediction(fit,mydata,mynewdata,dim(MLR)[2])
expect_true(file.exists("predictions_TestSet.csv"))
expect_true(file.exists("Exp.vs.Pred.tiff"))
expect_is(Rext,"numeric")
expect_less_than(Rext, 1)

new_nom2<-"ExternalSet.csv"	
newdata2<-read.csv(new_nom2,header=TRUE , sep=" ")
mynewdata2<-newdata2[,2:dim(newdata2)[2]]
coeff_val<-graphe_3Sets(fit,mydata,mynewdata,mynewdata2,dim(MLR)[2])
expect_is(coeff_val,"numeric")
expect_true(is.vector(coeff_val))
expect_equal(length(coeff_val),2)
expect_true(file.exists("Graphe_3sets.tiff"))
expect_less_than(coeff_val[1], 1)
expect_less_than(coeff_val[2], 1)

})



