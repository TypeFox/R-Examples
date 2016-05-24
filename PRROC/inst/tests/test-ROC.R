context("PRROC-ROC")

test_that("AUC-ROC - perfect, hard-labeled",{
	scores0<-runif(100);
	scores1<-(-runif(100));
	auc<-roc.curve(scores.class0 = scores0, scores.class1 = scores1)
	expect_that(auc$auc,equals(1))
})


test_that("AUC-ROC - class and fields",{
	scores0<-c(1);
	scores1<-c(-1);
	auc<-roc.curve(scores.class0 = scores0, scores.class1 = scores1)
	expect_that(auc,is_a("PRROC"))
	expect_that(auc$type,matches("ROC"))
	expect_that(auc$curve,equals(NULL))
	expect_that(auc$auc,is_a("numeric"))
	
	auc.curve<-roc.curve(scores.class0 = scores0, scores.class1 = scores1, curve = TRUE)
	expect_that(auc.curve,is_a("PRROC"))
	expect_that(auc.curve$type,matches("ROC"))
	expect_that(auc.curve$curve,is_a("matrix"))
	expect_that(dim(unique( auc.curve$curve ) ),equals(c(3,3)))
	expect_that(auc.curve$auc,is_a("numeric"))
	
	auc.curve<-roc.curve(scores.class0 = scores0, scores.class1 = scores1, curve = T, max.compute = T, min.compute = T, rand.compute = T);
	expect_that(auc.curve,is_a("PRROC"))
	expect_that(auc.curve$type,matches("ROC"))
	expect_that(auc.curve$curve,is_a("matrix"))
	expect_that(dim( unique(auc.curve$curve) ),equals(c(3,3)))
	expect_that(auc.curve$auc,is_a("numeric"))
	expect_that(auc.curve$max,is_a("PRROC"))
	expect_that(auc.curve$max$type,matches("ROC"))
	expect_that(auc.curve$min,is_a("PRROC"))
	expect_that(auc.curve$min$type,matches("ROC"))
	expect_that(auc.curve$rand,is_a("PRROC"))
	expect_that(auc.curve$rand$type,matches("ROC"))
})

test_that("AUC-ROC - worst, hard-labeled",{
	scores0<-(-runif(100));
	scores1<-(runif(100));
	auc<-roc.curve(scores.class0 = scores0, scores.class1 = scores1)
	expect_that(auc$auc,equals(0))
	
})


test_that("AUC-ROC - max, min, rand",{
	scores0<-(-runif(100));
	scores1<-(runif(100));
	auc<-roc.curve(scores.class0 = scores0, scores.class1 = scores1, min.compute=T, max.compute=T, rand.compute=T)
	expect_that(auc$min$auc,equals(0))
	expect_that(auc$max$auc,equals(1))
	expect_that(auc$rand$auc,equals(0.5))
})

test_that("AUC-ROC - identical scores",{
	scores0<-runif(100);
	scores1<-scores0;
	auc<-roc.curve(scores.class0 = scores0, scores.class1 = scores1)
	expect_that(auc$auc,equals(0.5))
})


test_that("ROC curve - perfect, hard-labeled",{
	scores0<-c(1)
	scores1<-c(-1)
	curve<-roc.curve(scores.class0 = scores0, scores.class1 = scores1,curve=T)$curve
	curve<-unique(curve);
	expect_that(curve,equals(matrix(c(1,0,0,1,1,0,-1,-1,1),nrow=3)))
})


test_that("AUC-ROC - ROCR, if available",{
	
	if( suppressWarnings( require("ROCR",character.only = T,quietly = T) ) ){
		
		scores0<-rnorm(n = 100, mean = 2);
		scores1<-rnorm(n = 100, mean = 0);
		scores<-c(scores0,scores1);
		lab<-c(rep(1,length(scores0)),rep(0,length(scores1)))
		pred<-prediction(predictions = scores,labels = lab)
		auc.rocr<-performance(pred,measure = "auc")@y.values[[1]]
		auc.prroc<-roc.curve(scores0,scores1)$auc
		expect_that(auc.prroc,equals(auc.rocr))
	}
	
})


test_that("AUC-ROC - value for random data (with specified seed)",{
	set.seed(101);
	scores0<-rnorm(n = 100, mean = 2);
	scores1<-rnorm(n = 100, mean = 0);
	auc.prroc<-roc.curve(scores0,scores1)$auc
	expect_that(auc.prroc,equals(0.9236))# value checked against ROCR
})


test_that("AUC-ROC - value for weighted, random data (with specified seed)",{
	set.seed(113);
	scores<-c(rnorm(n = 100, mean = 2), rnorm(n = 100, mean = 0));
	wfg<-c(runif(100,min=0.5,max=1),runif(100,min=0,max=0.5))
	auc.prroc<-roc.curve(scores.class0 = scores, weights.class0 = wfg)$auc
	expect_that(auc.prroc,equals(0.7287332599))
})

test_that("AUC-ROC - value of min, max, rand for weighted data",{
	scores<-c(rnorm(n = 100, mean = 2), rnorm(n = 100, mean = 0));
	wfg<-c(runif(100,min=0.5,max=1),runif(100,min=0,max=0.5))
	auc.prroc<-roc.curve(scores.class0 = scores, weights.class0 = wfg, max.compute = T, min.compute = T, rand.compute = T);
	expect_that(auc.prroc$rand$auc,equals(0.5))
	expect_that(auc.prroc$max$auc,equals(1-auc.prroc$min$auc))
})


test_that("AUC-ROC - consistency between weighted and hard-labeled interface",{
	scores0<-rnorm(n = 100, mean = 2);
	scores1<-rnorm(n = 100, mean = 0);
	scores<-c(scores0,scores1);
	lab<-c(rep(1,length(scores0)),rep(0,length(scores1)))
	auc.roc1<-roc.curve(scores0,scores1)$auc;
	auc.roc2<-roc.curve(scores,weights.class0 = lab)$auc
	expect_that(auc.roc1,equals(auc.roc2))
})



