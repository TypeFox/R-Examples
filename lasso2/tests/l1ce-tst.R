####------ l1ce() tests ------------
library(lasso2)

###-------- +/- the same as  example(l1ce) :
data(Iowa)

l1c.I <- l1ce(Yield ~ ., Iowa, bound = 10, trace = TRUE, absolute.t=TRUE)

## the next ones give a l1ce LIST (one for each bound)
l1c.liI <- l1ce(Yield ~ ., Iowa, bound = seq(1e-6,1, len= 17))
print(plot(l1c.liI))

data(Prostate)
## '0' fails on 64b: l1c.P <- l1ce(lpsa ~ ., Prostate, bound= c(0, 1, len= 17))
l1c.P <- l1ce(lpsa ~ ., Prostate, bound= c(1e-6,seq(1/16, 1, by=1/16)))
##                    complicated bound (for the comparison test below)
print(plot(l1c.P))
## test  multi-/single- bound problem:
l1c.P.25 <- l1ce(lpsa ~ ., Prostate, bound= 0.25)# 0.25 is nr. [5] above
l1c.P.1 <-  l1ce(lpsa ~ ., Prostate, bound= 1)
lm.P     <- lm  (lpsa ~ ., Prostate)
stopifnot(all.equal(coef(l1c.P.25), coef(l1c.P)[ 5,], tol= 1e-13),
          all.equal(coef(l1c.P.1 ), coef(l1c.P)[17,], tol= 1e-13),
          all.equal(coef(lm.P),     coef(l1c.P.1),    tol= 1e-13)
          )

###-------- Try a case where p > n :
n <- 100
p <- 120

RNGversion("1.6.0")
set.seed(n+p)
x <- matrix(runif(n*p), n,p, dimnames= list(NULL, paste("x",1:p,sep="")))
x <- as.data.frame(apply(x, 2, sort))
with(x,
     y <<- 5 + 4*x1 -3*x2 + 10*x3  + (eps <<- rnorm(n, sd = 1/200)))
d.ex <- cbind(y = y, x)
dim(d.ex)# 100 x 121
if(FALSE)
summary(lm(y ~ ., data = d.ex))# almost complete nonsense

## both these give something, but not at all the true model ...
l20 <- l1ce(y ~ ., data = d.ex, bound = 20, absolute.t = TRUE)
coef(l20)[coef(l20) > 0]

l15 <- l1ce(y ~ ., data = d.ex, bound = 15, absolute.t = TRUE)
coef(l15)[coef(l15) > 0]

sum(eps^2)
sum(resid(l20)^2) / sum(eps^2)
sum(resid(l15)^2) / sum(eps^2)

## Lower the bounds dramatically now:
l1.lis <- l1ce(y ~ ., data = d.ex, bound = seq(1e-5, 0.1, len=21))

pl1lis <- plot(l1.lis)#ylim = c(-10,10))
round(pl1lis$mat,3)
bnds <- pl1lis$bounds[,"rel.bound"]
round(1000 * t(pl1lis$mat[0.01 < bnds & bnds < 0.06 ,]))

### ---------- Check that l1ce(.) can be used inside functions -------

## {has not worked for a long time, till 2005-06-10}:
myLasso <- function(formula, data, sweep.out = ~ 1, boundset) {
    ## Of course this is silly: calling l1ce() for each bound separately:
    sapply(boundset, function(B)
           l1ce(formula, data = data, sweep.out = sweep.out, bound = B),
           simplify = FALSE)
}
my.lis <- myLasso(y ~ ., data = d.ex, boundset = seq(1e-5, 0.1, len=6))
## should have a subset from the models in  l1.lis above !


### ------------------------------------------------------------------
value = c(-1.2229125582, 0.0291875346, -1.1354055994, -1.1586843969,
0.0265310929, -1.2761698858, -0.9567990729, 0.0001296789,
-0.9174081609, -0.8564132893 , 0.0825918391, -0.8348783996,
-0.8967906196, -0.1102660185, -0.9392872627, -0.6826339823,
0.0521741319, -0.7737865885)

P3 =  c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1)
P4 =  c(0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1)

P1 = c(1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0)

df = data.frame(value,P1,P3,P4)

try( l1ce(value~-1+(P3+P4+P1)^2,data=df,na.action=na.omit,absolute.t=TRUE,standardize=FALSE,bound=0.4,sweep.out=~-1+ (P3+P4+P1) ))

try( l1ce(value~-1+(P3+P4+P1)^2,data=df,na.action=na.omit,absolute.t=TRUE,standardize=FALSE,bound=0.5,sweep.out=~-1+ (P3+P4+P1) ))

try(l1ce(value~-1+(P3+P4+P1)^2,data=df,na.action=na.omit,absolute.t=TRUE,standardize=FALSE,bound=0.6,sweep.out=~-1+ (P3+P4+P1) ))
