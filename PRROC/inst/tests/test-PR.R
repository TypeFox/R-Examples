context("PRROC-PR")

test_that("AUC-PR - perfect, hard-labeled",{
	scores0<-runif(100);
	scores1<-(-runif(100));
	auc<-pr.curve(scores.class0 = scores0, scores.class1 = scores1)
	auc.in<-auc$auc.integral
	auc.dg<-auc$auc.davis.goadrich
	expect_that(auc.in,equals(1))
	expect_that(auc.dg,equals(1))
})


test_that("AUC-PR - class and fields",{
	scores0<-c(1);
	scores1<-c(-1);
	auc<-pr.curve(scores.class0 = scores0, scores.class1 = scores1)
	expect_that(auc,is_a("PRROC"))
	expect_that(auc$type,matches("PR"))
	expect_that(auc$curve,equals(NULL))
	expect_that(auc$auc.integral,is_a("numeric"))
	expect_that(auc$auc.davis.goadrich,is_a("numeric"))
	
	auc.curve<-pr.curve(scores.class0 = scores0, scores.class1 = scores1, curve = TRUE)
	expect_that(auc.curve,is_a("PRROC"))
	expect_that(auc.curve$type,matches("PR"))
	expect_that(auc.curve$curve,is_a("matrix"))
	expect_that(dim(auc.curve$curve),equals(c(3,3)))
	expect_that(auc.curve$auc.integral,is_a("numeric"))
	expect_that(auc.curve$auc.davis.goadrich,is_a("numeric"))
	
	auc.curve<-pr.curve(scores.class0 = scores0, scores.class1 = scores1, curve = T, max.compute = T, min.compute = T, rand.compute = T);
	expect_that(auc.curve,is_a("PRROC"))
	expect_that(auc.curve$type,matches("PR"))
	expect_that(auc.curve$curve,is_a("matrix"))
	expect_that(dim(auc.curve$curve),equals(c(3,3)))
	expect_that(auc.curve$auc.integral,is_a("numeric"))
	expect_that(auc.curve$auc.davis.goadrich,is_a("numeric"))
	expect_that(auc.curve$max,is_a("PRROC"))
	expect_that(auc.curve$max$type,matches("PR"))
	expect_that(auc.curve$min,is_a("PRROC"))
	expect_that(auc.curve$min$type,matches("PR"))
	expect_that(auc.curve$rand,is_a("PRROC"))
	expect_that(auc.curve$rand$type,matches("PR"))
})

test_that("AUC-PR - worst, hard-labeled",{
	scores0<-(-runif(1000));
	scores1<-(runif(1000));
	auc<-pr.curve(scores.class0 = scores0, scores.class1 = scores1)
	expect_that(auc$auc.integral,equals(0.3068528194))
})


test_that("AUC-PR - max, min, rand",{
	scores0<-(-runif(100));
	scores1<-(runif(100));
	auc<-pr.curve(scores.class0 = scores0, scores.class1 = scores1, min.compute=T, max.compute=T, rand.compute=T)
	expect_that(auc$min$auc.integral,equals(0.3068528194))
	expect_that(auc$auc.integral,equals(auc$min$auc.integral))
	expect_that(auc$max$auc.integral,equals(1))
	expect_that(auc$rand$auc.integral,equals(0.5))
	
	scores0<-(-runif(100));
	scores1<-(runif(900));
	auc<-pr.curve(scores.class0 = scores0, scores.class1 = scores1, rand.compute=T)
	expect_that(auc$rand$auc.integral,equals(0.1))
})


test_that("AUC-PR - identical scores",{
	scores0<-runif(100);
	scores1<-scores0;
	auc<-pr.curve(scores.class0 = scores0, scores.class1 = scores1)$auc.integral
	expect_that(auc,equals(0.5))
})


test_that("PR curve - perfect, hard-labeled",{
	scores0<-c(1)
	scores1<-c(-1)
	curve<-pr.curve(scores.class0 = scores0, scores.class1 = scores1,curve=T)$curve
	curve<-unique(curve);
	expect_that(curve,equals(matrix(c(1,1,0,0.5,1,1,-1,1,1),nrow=3)))
})


test_that("AUC-PR - value for random data (with specified seed)",{
	set.seed(101);
	scores0<-runif(n = 100, min = 0.3, max = 1.0)# we need scores between 0 and 1 for auc.jar
	scores1<-runif(n = 100, min = 0.0, max = 0.6)
	auc.prroc<-pr.curve(scores0,scores1)
	expect_that(auc.prroc$auc.davis.goadrich,equals(0.9018680216532841))# value from Davis & Goadrich (auc.jar)
	expect_that(auc.prroc$auc.integral,equals(0.901869384113))
})


test_that("AUC-PR - value for weighted, random data (with specified seed)",{
	set.seed(113);
	scores<-c(rnorm(n = 100, mean = 2), rnorm(n = 100, mean = 0));
	wfg<-c(runif(100,min=0.5,max=1),runif(100,min=0,max=0.5))
	auc.prroc<-pr.curve(scores.class0 = scores, weights.class0 = wfg)$auc.integral
	expect_that(auc.prroc,equals(0.7053830324))
})


test_that("AUC-PR - consistency between weighted and hard-labeled interface",{
	scores0<-rnorm(n = 100, mean = 2);
	scores1<-rnorm(n = 100, mean = 0);
	scores<-c(scores0,scores1);
	lab<-c(rep(1,length(scores0)),rep(0,length(scores1)))
	auc.pr1<-pr.curve(scores0,scores1)$auc.integral;
	auc.pr2<-pr.curve(scores,weights.class0 = lab)$auc.integral
	expect_that(auc.pr1,equals(auc.pr2))
})



