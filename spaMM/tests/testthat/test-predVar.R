cat("\ntest of prediction variance:")
# examples from Booth & Hobert 1998 JASA

npos <- c(11,16,14,2,6,1,1,4,10,22,7,1,0,0,1,6)
ntot <- c(36,20,19,16,17,11,5,6,37,32,19,17,12,10,9,7)
treatment <- c(rep(1,8),rep(0,8))
clinic <-c(seq(8),seq(8))
clinics <- data.frame(npos=npos,nneg=ntot-npos,treatment=treatment,clinic=clinic)
fitobject <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(),data=clinics,HLmethod="ML")
res <- sqrt(get_predVar(fitobject))
# res: cf row 2, table 6; some discrepancies (surely stemming from point estimates).
expect_equal(res[6],c("6"=0.7316108),tolerance=2e-4)

if(require("rsae", quietly = TRUE)) {
  data(landsat)
  fitobject <- HLfit(HACorn ~ PixelsCorn + PixelsSoybeans + (1|CountyName),data=landsat[-33,],HLmethod="ML")
  newXandZ <- unique(data.frame(PixelsCorn=landsat$MeanPixelsCorn,PixelsSoybeans=landsat$MeanPixelsSoybeans,CountyName=landsat$CountyName))
  res <- get_predVar(fitobject,newdata=newXandZ)
  expect_equal(res[12],c("32"=27.80684),tolerance=2e-4)
  res <- get_predVar(fitobject,newdata=newXandZ,variances=list(disp=FALSE)) ## sum of cols nu and beta in Table 7
  expect_equal(res[12],c("32"=26.56215),tolerance=2e-4)
} else {
  cat( "package 'rsae' not available, cannot run prediction variance test.\n" )
}

# simple LM
set.seed(123)
n <- 100
k <- rep(1:(n/10), each=10)
data.test <- data.frame(y=rnorm(n, mean=k), x=factor(k))
hlfit <- HLfit(y ~ x, data=data.test)
p <- predict(hlfit, newdata=data.frame(x=unique(data.test$x)), variances=list(respVar=TRUE))
expect_equal(length(attr(p, "residVar")),10L)
expect_equal(attr(p, "residVar")[1],c(`1`=0.8347247),tolerance=1e-6)
expect_equal(attr(p, "predVar")[1],c(`1`=0.08347247),tolerance=1e-6)
expect_equal(attr(p, "respVar")[1],c(`1`=0.918197),tolerance=1e-6)
## also tests that a single level of the factor in the newdata is not a problem
pp <- predict(hlfit, newdata=data.test[1,], variances=list(respVar=TRUE))
expect_equal(pp[1],c(`1`=1.074626),tolerance=1e-6)

## multiple tests with two ranefs 
set.seed(123)
data(Loaloa)
ll <- cbind(Loaloa,idx=sample(2,size=nrow(Loaloa),replace=TRUE))
hl <- HLCor(cbind(npos,ntot-npos)~1+Matern(1|longitude+latitude)+(1|idx),data=ll,
            family=binomial(),ranPars=list(nu=0.5,rho=1/0.7),method="PQL/L")
# verif single point input
p1 <- predict(hl,variances=list(respVar=TRUE))
p2 <- predict(hl,newdata=ll[1,],variances=list(respVar=TRUE))
expect_equal(attr(p1,"respVar")[1],attr(p2,"respVar")[1])
expect_equal(p1[1],p2[1])
# verif 'slice' mechanism including Evar (hence new levels of ranef) 
lll <- ll
lll$idx <- lll$idx+1
lll$latitude <- lll$latitude+1
p4 <- predict(hl,newdata=lll[1:101,],variances=list(respVar=TRUE))
p5 <- predict(hl,newdata=lll[101:102,],variances=list(respVar=TRUE))
expect_equal(attr(p4,"respVar")[101],attr(p5,"respVar")[1])
expect_equivalent(p4[101],p5[1]) ## _equivalent does not check names and other attributes



