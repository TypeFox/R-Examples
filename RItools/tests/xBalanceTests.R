require("RItools")

data(nuclearplants)

##################################################
### Basic uses
##################################################
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data=nuclearplants)

xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test'))

xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=data.frame(unstrat=factor(character(32)),
           pt=factor(nuclearplants$pt)),
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test'))

(xb0 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=list(unstrat=NULL, pt=~pt),
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test'))
)
##########################################################################################
### Oddness on LHS of formula
##########################################################################################Q
xBalance(I(pr==1) ~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=list(unstrat=NULL, pt=~pt),
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test')
         )


###(b0 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                ~factor(pt), nuclearplants))
#####################################################
###
###
#####################################################
###### Test balance on v-bles one at a time.
#####################################################
### xBalance(pr~ date, ~factor(pt), nuclearplants)
###
#####################################################
###### Factor handling
#####################################################
###a0 <- xBalance(pr~ date + t1 + t2 + cap + CUM.N,
###               ~factor(pt),
###               data.frame(nuclearplants,
###                          CUM.N=cut(nuclearplants$cum.n, c(0,2,5,21), include=TRUE)),
###               chisq=FALSE)
###xBalance(pr~ date + t1 + t2 + cap + CUM.N,
###         ~factor(pt),
###         data.frame(nuclearplants,
###                    CUM.N=cut(nuclearplants$cum.n, c(0,2,5,21), include=TRUE)),
###         chisq=TRUE)
###a1 <- xBalance(pr~ date + t1 + t2 + cap +
###               cut(cum.n, c(0,2,5,21), include=TRUE),
###               ~factor(pt),
###               nuclearplants)
###attributes(a0) <- attributes(a1) <- NULL
###all.equal(a0,a1)
###
###Also, handling bad factors in naImpute (i.e. only one level).
#####################################################
###### When there are NAs in the grouping variable and we're doing
###### na.omit=TRUE, we don't want to omit rows with NA on
###### the grouping before evaluating fmla.  For me this is of concern
###### when fmla has a call to ns() in it, but it's easier to test the issue
###### using rank() (-BH)
#####################################################
###nuclearplants$pt1 <- factor(nuclearplants$pt)
###nuclearplants$pt1[1:3] <- NA
###nuclearplants$rank.t1 <- rank(nuclearplants$t1)
###all.equal(nuclearplants$rank.t1[-(1:3)], rank(nuclearplants$t1[-(1:3)]))
###### note that the above is not TRUE -- omitting changes the ranks
###b1 <- xBalance(pr~rank.t1+t2, ~pt1, data=nuclearplants, na.rm=TRUE)[1,]
###b2 <- xBalance(pr~rank(t1)+t2, ~pt1, data=nuclearplants, na.rm=TRUE)[1,]
###attributes(b1) <- attributes(b2) <- NULL
###all.equal(b2,b1)
###
#####################################################
###### finding stratum variable that's in calling
###### frame but not in data
#####################################################
###mypt <- factor(nuclearplants$pt)
###all.equal(b0, xBalance(pr~date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                       ~mypt, nuclearplants))
#####################################################
###### * Handling stratum weights                 ###
#####################################################
###### ** Default weightings should evaluate before NAs are removed,
###### in my view.  (Note that this will require some re-coding, even before
###### nonstandard weighting functionality is added.  However, that
###### same re-coding ought to streamline integration of the new functionality.)
###b3 <- xBalance(pr~ifelse(pt,rank.t1,0)+t2, ~factor(pt), data=nuclearplants, na.rm=TRUE)[1,c(2,4)]
###b4 <- xBalance(pr~ifelse(pt,rank.t1,0)+pt1, ~factor(pt), data=nuclearplants, na.rm=TRUE)[1,c(2,4)]
###attributes(b3) <- attributes(b4) <- NULL
###all.equal(b3, b4) # Should be TRUE
###
###### ** Provided they're all nonnegative, they shouldn't
###### have to sum to 1.  Ie, the function should handle the re-normalizing.
###(my.hmnic.wts <- tapply(nuclearplants$pr, nuclearplants$pt,
###                        function(x) {1/(1/sum(x) + 1/sum(!x))}) )
###identical(b0,
###          xBalance(pr~date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                   ~factor(pt), nuclearplants,
###                   stratum.weights=my.hmnic.wts) )
###
###### **... unless user specifies "normalize.weights=FALSE"!  (For group-level
###### assmt, users may want to work with group-aggregated covariates while folding
###### de-aggregation into the weights.)
#####  NEED TEST HERE
###
###### ** If stratum.weights is a function, it should
###### find v-bles within its data argument.
###pr <- nuclearplants$pr ; pr[which.max(pr)] <- 0
###(not.hmnic.wts <- tapply(pr, nuclearplants$pt,
###                        function(x) {1/(1/sum(x) + 1/sum(!x))}) )
###identical(b0, # This should be FALSE
###          xBalance(pr~date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                   ~factor(pt), nuclearplants,
###                   stratum.weights=not.hmnic.wts) )
###identical(b0, # Whereas this should be TRUE, if evaluation happens properly
###          xBalance(pr~date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                   ~factor(pt), nuclearplants,
###                   stratum.weights=function(data){tapply(data$pr, data$pt,
###                     function(x) {1/(1/sum(x) + 1/sum(!x))})} ) )
###
###### ** If stratum.weights is a function, its
###### evaluation should within xBalance's
###### parent frame for objects.
###PT <- nuclearplants$pt
###identical(b0, # Should be TRUE
###          xBalance(pr~date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                   ~factor(pt), nuclearplants,
###                   stratum.weights=function(data){tapply(data$pr, PT,
###                     function(x) {1/(1/sum(x) + 1/sum(!x))})} ) )
###
###### ** If stratum.weights is function to be evaluated within
###### data, then it needs to evaluate before na's are removed. (IMHO. -BH)
###
###### ** If stratum.weights is not a function, require it to specify a nonnegative
###### real number weight for each stratum in the reduced levels
###### of groups. The following tries should fail.
###try(xBalance(pr~.-cost-pt, ~factor(pt), nuclearplants,
###             stratum.weights=my.hmnic.wts[-1]) )
###try(xBalance(pr~.-cost-pt, ~factor(pt), nuclearplants,
###             stratum.weights=c(NA,my.hmnic.wts[-1]) ))
###try(xBalance(pr~.-cost-pt, ~factor(pt), nuclearplants,
###             stratum.weights=-my.hmnic.wts) )
###try(xBalance(pr~.-cost-pt, ~factor(pt), nuclearplants,
###             stratum.weights=unlist(list('0'='a', '1'='b')) ) )
###
###### Non-binary "treatment" variable.
###xBalance(cost~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants)
###
###xBalance(cost~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants)
###
###xBalance(cost~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants,
###         stratum.weights=my.hmnic.wts)
###
###xBalance(cost~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants,
###         stratum.weights=function(data){tapply(data$pr, data$pt,
###                     function(x) {1/(1/sum(x) + 1/sum(!x))})} )
###
#####################################################
######        covariate.scaling                   ###
#####################################################
###xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants,
###         covariate.scaling=1)
###b0[1:2,]
###
#####################################################
######   gracefully handle case of no strata w/   ###
######   variation in Tx group                    ###
#####################################################
###xBalance(pr~ date + t1,
###         ~factor(pt), nuclearplants[c(1,32),])
###
###
#####################################################
######      PRINT JUST THE CHISQUARE TESTS        ###
#####################################################
###b0[numeric(0),]
###
#####################################################
######               CHISQ="only"                 ###
#####################################################
###xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###                ~factor(pt), nuclearplants, chisq='only')
###
#####################################################
######               xtable method                ###
#####################################################
if (require('xtable'))
  {
  xtablea <- xtable(xb0)
  xtableb <- xtable(xb0, caption="Caption!", label="thetable", digits=1,
       align=rep('l', prod(dim(xb0$result)[2:3])+1),
       display=c('s', rep(c(rep('fg', dim(xb0$result)[2]-1), 's'),
         dim(xb0$result)[3]) ) #,col.labels= do this one later
       )
  }


#####################################################
######               include.means                ###
#####################################################
###xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants,
###         covariate.scaling=1, include.means=TRUE)
###xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants, include.means=TRUE)
###
###
#####################################################
######  na.rm=FALSE with missing covariates       ###
#####################################################
### Should create a new variable (0=not missing,1=missing) and impute missing values with the mean (median is new default)

set.seed(123)
testdata<-nuclearplants
testdata$date[sample(1:32,10)]<-NA

xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data = testdata,
    na.rm = FALSE,impfn=mean.default) ##first using the mean to match up with previous versions

#####################################################
######  handling factor with no obs for a level in a strata  ###
#####################################################

testdata<-nuclearplants
testdata$cum.nF<-factor(testdata$cum.n)

##Notice that for most levels of cum.n, there are no obs in one of the two strata
table(testdata$pt,testdata$cum.n)

##    1 2 3 5 6 7 8 11 12 14 15 16 17 18 19 20 21
##  0 4 3 4 2 1 1 1  0  2  1  1  1  1  1  1  1  1
##  1 0 0 0 0 0 1 2  3  0  0  0  0  0  0  0  0  0

##First no missing levels, same in both strata --- looks ok
xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.nF, strata=list(nostrata=NULL,thept=~pt), data = testdata,impfn=mean.default)

##Second two missing levels, same in both strata
##This doesn't look as good --- we'd prefer to drop levels that don't exist.

testdata$cum.nF[testdata$cum.n>16]<-NA
testdata$cum.nF[testdata$cum.n==7]<-NA
table(testdata$pt,testdata$cum.nF,exclude=c()) ##Notice that the levels don't disappear by default.
xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.nF, strata=list(nostrata=NULL,thept=~pt), data = testdata,na.rm=FALSE,impfn=mean.default)

##This isn't right either.
xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.nF, strata=list(nostrata=NULL,thept=~pt), data = testdata,na.rm=TRUE,impfn=mean.default)


#####################################################
######  handling factor with no levels strata=argument  ###
#####################################################
try(xBalance(pr ~ date,
             strata=data.frame(nastrat=factor(rep(NA,nrow(testdata))) ),
             data=testdata), FALSE)

xBalance(pr ~ date,
         strata=data.frame(nostrata=factor(rep('a',nrow(testdata))),
           nastrat=factor(rep(NA,nrow(testdata))) ),
         data=testdata)
#####################################################
######             WISHLIST                       ###
#####################################################
###
###
###
###
