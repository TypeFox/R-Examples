### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: vignette.Rnw:20-31
###################################################
library(ezsim)
ezsim_basic<-ezsim(
    m             = 50,
    run           = TRUE,
    display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
    parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
    dgp           = function() rnorm(n,mu,sigma),
    estimator     = function(x) c(mean_hat = mean(x), 
                                 sd_mean_hat=sd(x)/sqrt(length(x)-1)),
    true_value    = function() c(mu, sigma / sqrt(n-1))
)


###################################################
### code chunk number 2: vignette.Rnw:32-34
###################################################
summary_ezsim_basic<-summary(ezsim_basic)
head(summary_ezsim_basic,16)


###################################################
### code chunk number 3: vignette.Rnw:35-37 (eval = FALSE)
###################################################
## plot(ezsim_basic)
## plot(ezsim_basic,"density")


###################################################
### code chunk number 4: vignette.Rnw:40-41
###################################################
print(plot(ezsim_basic,return_print=TRUE)[[1]])


###################################################
### code chunk number 5: vignette.Rnw:44-45
###################################################
print(plot(ezsim_basic,return_print=TRUE)[[2]])


###################################################
### code chunk number 6: vignette.Rnw:48-49
###################################################
print(plot(ezsim_basic,"density",return_print=TRUE)[[1]])


###################################################
### code chunk number 7: vignette.Rnw:52-53
###################################################
print(plot(ezsim_basic,"density",return_print=TRUE)[[2]])


###################################################
### code chunk number 8: vignette.Rnw:82-84
###################################################
par_def<-createParDef(selection=list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5)))
par_def


###################################################
### code chunk number 9: vignette.Rnw:87-88
###################################################
generate(par_def)[1:3]


###################################################
### code chunk number 10: vignette.Rnw:92-95
###################################################
par_def2<-createParDef(selection=list(mu1=5,mu2=3,n=c(10,20)), 
    banker=list(Sigma=matrix(c(1,.4,.4,1),nrow=2)))
generate(par_def2)


###################################################
### code chunk number 11: vignette.Rnw:101-104
###################################################
dgp<-function(){
	rnorm(n,mu,sigma)
}


###################################################
### code chunk number 12: vignette.Rnw:106-109
###################################################
evalFunctionOnParameterDef(par_def,dgp,index=1)
evalFunctionOnParameterDef(par_def,dgp,index=2)



###################################################
### code chunk number 13: vignette.Rnw:111-119
###################################################
dgp_2<-function(){
    z1<-rnorm(n)
    z2<-rnorm(n)
    cbind(x1=mu1+z1*Sigma[1,1], 
          x2=mu2+ Sigma[2,2]*(Sigma[1,2]*z1+ sqrt(1-Sigma[1,2]^2)*z2 ))
}
evalFunctionOnParameterDef(par_def2,dgp_2)



###################################################
### code chunk number 14: vignette.Rnw:123-127
###################################################
estimator<-function(x){
    c(mean_hat = mean(x), sd_mean_hat=sd(x)/sqrt(length(x)-1))
}
estimator(evalFunctionOnParameterDef(par_def,dgp,index=1))


###################################################
### code chunk number 15: vignette.Rnw:132-136
###################################################
true<-function(){
    c(mu, sigma / sqrt(n-1))
}
evalFunctionOnParameterDef(par_def,true)


###################################################
### code chunk number 16: vignette.Rnw:141-142
###################################################
display_name<-c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])")


###################################################
### code chunk number 17: vignette.Rnw:147-153 (eval = FALSE)
###################################################
## estimator_lm <- function(x) {
##     lm(y~x1+x2, data=x)
## }
## estimator_parser_lm <- function(x){
##     coef(x)
## }


###################################################
### code chunk number 18: vignette.Rnw:161-173 (eval = FALSE)
###################################################
## library(ezsim)
## ezsim_basic<-ezsim(
##     m             = 50,
##     run           = FALSE,
##     display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
##     parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
##     dgp           = function() rnorm(n,mu,sigma),
##     estimator     = function(x) c(mean_hat = mean(x), 
##                                  sd_mean_hat=sd(x)/sqrt(length(x)-1)),
##     true_value    = function() c(mu, sigma / sqrt(n-1))
## )
## ezsim_basic = run(ezsim_basic)


###################################################
### code chunk number 19: vignette.Rnw:187-202 (eval = FALSE)
###################################################
## library(ezsim)
## ezsim_basic<-ezsim(
##     m             = 50,
##     run           = TRUE,
##     use_core      = 4, 
##     use_seed      = 123,
##     cluster_packages = "Jmisc",
##     display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
##     parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
##     dgp           = function() rnorm(n,mu,sigma),
##     estimator     = function(x) c(mean_hat = mean(x), 
##                                  sd_mean_hat=sd(x)/sqrt(length(x)-1)),
##     true_value    = function() c(mu, sigma / sqrt(n-1))
## )
## summary(ezsim_basic)


###################################################
### code chunk number 20: vignette.Rnw:211-230 (eval = FALSE)
###################################################
## library(parallel)
## my_cluster = makeCluster(4)
## clusterSetRNGStream(my_cluster, 123)
## 
## library(ezsim)
## ezsim_basic<-ezsim(
##     m             = 50,
##     run           = TRUE,
##     use_core      = 4, 
##     cluster       = my_cluster,
##     display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
##     parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
##     dgp           = function() rnorm(n,mu,sigma),
##     estimator     = function(x) c(mean_hat = mean(x), 
##                                  sd_mean_hat=sd(x)/sqrt(length(x)-1)),
##     true_value    = function() c(mu, sigma / sqrt(n-1))
## )
## stopCluster(my_cluster)
## summary(ezsim_basic)


###################################################
### code chunk number 21: vignette.Rnw:247-248
###################################################
summary(ezsim_basic,subset=list(estimator="mean_hat",n=c(20,40),sigma=c(1,3)))


###################################################
### code chunk number 22: vignette.Rnw:252-254
###################################################
summary(ezsim_basic,simple=FALSE,
        subset=list(estimator="mean_hat",n=c(20,40),sigma=c(1,3)))


###################################################
### code chunk number 23: vignette.Rnw:258-262
###################################################
summary(ezsim_basic,stat=c("q25","median","q75"),
        Q025=quantile(value_of_estimator,0.025),
        Q975=quantile(value_of_estimator,0.975),
        subset=list(estimator="mean_hat",n=c(20,40),sigma=c(1,3)))


###################################################
### code chunk number 24: vignette.Rnw:272-273 (eval = FALSE)
###################################################
## plot(ezsim_basic,subset=list(estimator="sd_mean_hat",mu=3))


###################################################
### code chunk number 25: vignette.Rnw:274-276
###################################################
print(plot(ezsim_basic,subset=list(estimator="sd_mean_hat",mu=0),return_print=TRUE)[[1]])



###################################################
### code chunk number 26: vignette.Rnw:277-278 (eval = FALSE)
###################################################
## plot(ezsim_basic,subset=list(estimator="mean_hat",sigma=3))


###################################################
### code chunk number 27: vignette.Rnw:279-280
###################################################
print(plot(ezsim_basic,subset=list(estimator="mean_hat",sigma=3),return_print=TRUE)[[1]])


###################################################
### code chunk number 28: vignette.Rnw:286-288 (eval = FALSE)
###################################################
## plot(ezsim_basic,subset=list(estimator="sd_mean_hat",mu=0),
## parameters_priority=c("sigma","n"))


###################################################
### code chunk number 29: vignette.Rnw:289-290
###################################################
print(plot(ezsim_basic,subset=list(estimator="sd_mean_hat",mu=0),parameters_priority="sigma",return_print=TRUE)[[1]])


###################################################
### code chunk number 30: vignette.Rnw:291-292 (eval = FALSE)
###################################################
## plot(ezsim_basic,subset=list(estimator="mean_hat",sigma=c(1,3)),parameters_priority="mu")


###################################################
### code chunk number 31: vignette.Rnw:293-294
###################################################
print(plot(ezsim_basic,subset=list(estimator="mean_hat",sigma=c(1,3)),parameters_priority="mu",return_print=TRUE)[[1]])


###################################################
### code chunk number 32: vignette.Rnw:300-303 (eval = FALSE)
###################################################
## plot(ezsim_basic,"density",
##      subset=list(estimator="mean_hat",sigma=3),
##      parameters_priority="n",benchmark=dnorm)


###################################################
### code chunk number 33: vignette.Rnw:304-305
###################################################
print(plot(ezsim_basic,"density",benchmark=dnorm,subset=list(estimator="mean_hat",sigma=3),parameters_priority="n",return_print=TRUE)[[1]])


###################################################
### code chunk number 34: vignette.Rnw:306-309 (eval = FALSE)
###################################################
## plot(ezsim_basic,"density",
##      subset=list(estimator="mean_hat",mu=0),
##      parameters_priority="n" ,benchmark=dnorm)


###################################################
### code chunk number 35: vignette.Rnw:310-311
###################################################
print(plot(ezsim_basic,"density",benchmark=dnorm,subset=list(estimator="mean_hat",mu=0),parameters_priority="n",return_print=TRUE)[[1]])


###################################################
### code chunk number 36: vignette.Rnw:316-317 (eval = FALSE)
###################################################
## plot(summary(ezsim_basic,c("q25","q75")))


###################################################
### code chunk number 37: vignette.Rnw:318-319
###################################################
print(plot(summary(ezsim_basic,c("q25","q75")),return_print=TRUE))


###################################################
### code chunk number 38: vignette.Rnw:320-321 (eval = FALSE)
###################################################
## plot(summary(ezsim_basic,c("q25","q75"),subset=list(estimator="mean_hat")))


###################################################
### code chunk number 39: vignette.Rnw:322-323
###################################################
print(plot(summary(ezsim_basic,c("q25","q75"),subset=list(estimator="mean_hat")),return_print=TRUE))


###################################################
### code chunk number 40: vignette.Rnw:324-325 (eval = FALSE)
###################################################
## plot(summary(ezsim_basic,c("median"),subset=list(estimator="sd_mean_hat")))


###################################################
### code chunk number 41: vignette.Rnw:326-327
###################################################
print(plot(summary(ezsim_basic,c("median"),subset=list(estimator="sd_mean_hat")),return_print=TRUE))


###################################################
### code chunk number 42: vignette.Rnw:334-356
###################################################
ez_powerfun<-ezsim(
    m             = 100,
    run           = TRUE,
    display_name  = c(b="beta",es="sigma[e]^2",xs="sigma[x]^2"),
    parameter_def = createParDef(selection=list(xs=1,n=50,es=5,b=seq(-1,1,0.1))),
    dgp           = function(){
                        x<-rnorm(n,0,xs)
                        e<-rnorm(n,0,es)
                        y<-b * x + e
                        data.frame(y,x)
                    },
    estimator     = function(d){
                        r<-summary(lm(y~x-1,data=d))
                        stat<-r$coef[,1]/r$coef[,2]

                        # test whether b > 0
                        # level of significance : 5%
                        out <- stat > c(qnorm(.95), qt(0.95,df=r$df[2]))
                        names(out)<-c("z-test","t-test")
                        out
                    }
)


###################################################
### code chunk number 43: vignette.Rnw:359-360 (eval = FALSE)
###################################################
## plot(ez_powerfun,"powerfun",null_hypothesis=0)


###################################################
### code chunk number 44: vignette.Rnw:361-362
###################################################
print(plot(ez_powerfun,"powerfun",null_hypothesis=0,return_print=TRUE))


