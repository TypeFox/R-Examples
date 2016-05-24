
library(survey)
data(api)
options(survey.replicates.mse=TRUE)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

a<-svyby(~api00+api99, ~comp.imp+sch.wide,design=rclus1,svymean,
         covmat=TRUE,drop.empty.groups=FALSE)
b<-svyby(~api00+api99, ~comp.imp+sch.wide,design=rclus1,svymean,
         covmat=TRUE,drop.empty.groups=TRUE)

stopifnot(all(na.omit(
              as.vector(as.matrix(SE(a)))==sqrt(diag(vcov(a)))
)))
stopifnot(all(
              as.vector(as.matrix(SE(b)))==sqrt(diag(vcov(b)))
              ))

rat <- svyratio(~ell+mobility, ~mobility+meals, dclus1,covmat=TRUE)
all <- svytotal(~ell+mobility+meals, dclus1)

stopifnot(all(abs(vcov(svycontrast(all,
                                   list(quote(ell/mobility),
                                        quote(mobility/mobility),
                                        quote(ell/meals),quote(mobility/meals))))
                  -vcov(rat))<1e-10))

stopifnot(all(abs(SE(rat)-sqrt(diag(vcov(rat))))<1e-10))

rat <- svyratio(~ell+mobility, ~mobility+meals, rclus1,covmat=TRUE)
all <- svytotal(~ell+mobility+meals, rclus1, return.replicates=TRUE)

con<-svycontrast(all,
                 list(quote(ell/mobility),
                      quote(mobility/mobility),
                      quote(ell/meals),quote(mobility/meals)))

stopifnot(all(abs(survey:::svrVar(con$replicates, rclus1$scale,rclus1$rscales,mse=rclus1$mse, coef=coef(con))-vcov(rat))<1e-10))

options(survey.replicates.mse=FALSE)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

a<-svyby(~api00+api99, ~comp.imp+sch.wide,design=rclus1,svymean,
         covmat=TRUE,drop.empty.groups=FALSE)
b<-svyby(~api00+api99, ~comp.imp+sch.wide,design=rclus1,svymean,
         covmat=TRUE,drop.empty.groups=TRUE)

stopifnot(all(na.omit(
              as.vector(as.matrix(SE(a)))==sqrt(diag(vcov(a)))
)))
stopifnot(all(
              as.vector(as.matrix(SE(b)))==sqrt(diag(vcov(b)))
              ))

rat <- svyratio(~ell+mobility, ~mobility+meals, dclus1,covmat=TRUE)
all <- svytotal(~ell+mobility+meals, dclus1)

stopifnot(all(abs(vcov(svycontrast(all,
                                   list(quote(ell/mobility),
                                        quote(mobility/mobility),
                                        quote(ell/meals),quote(mobility/meals))))
                  -vcov(rat))<1e-10))

stopifnot(all(abs(SE(rat)-sqrt(diag(vcov(rat))))<1e-10))

rat <- svyratio(~ell+mobility, ~mobility+meals, rclus1,covmat=TRUE)
all <- svytotal(~ell+mobility+meals, rclus1, return.replicates=TRUE)

con<-svycontrast(all,
                 list(quote(ell/mobility),
                      quote(mobility/mobility),
                      quote(ell/meals),quote(mobility/meals)))

stopifnot(all(abs(survey:::svrVar(con$replicates, rclus1$scale,rclus1$rscales,mse=rclus1$mse, coef=coef(con))-vcov(rat))<1e-10))


