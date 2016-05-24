## The names of the components of trellis objects
## are different in R than in S-Plus.

## These data reprinted in \cite{hand:1994} are originally from Pearce,
## S.C., 1983, The Agricultural Field Experiment, Wiley.
##
## The response is crop yield in pounds and the covariable is yield
## in bushels in a prior period under the same growing conditions.
## The treatments are growing conditions, where level 6 is a control.
## There are 4 blocks.  Hand implies that treat is significant iff
## the covariable is taken into account.
##

## if.R(r=data(apple),
##      s={
##        apple <- read.table(hh("datasets/apple.dat"), header=TRUE)
##        apple$treat <- factor(apple$treat)
##        contrasts(apple$treat) <- contr.treatment(6)
##        apple$block <- factor(apple$block)
##      })
data(apple)

apple.ancova.1 <- aov(yield ~ block + pre*treat, data=apple)
anova(apple.ancova.1)

apple.ancova.2 <- aov(yield ~ block + pre + treat, data=apple)
anova(apple.ancova.2)

apple.ancova.2b <- aov(yield ~ block + treat + pre, data=apple)
anova(apple.ancova.2b)

apple.ancova.2 <- update(apple.ancova.2, x=TRUE)
apple.ancova.2$x
coef(apple.ancova.2)
predict(apple.ancova.2)


## find and remove block effect from response variable and covariable
yield.block.effect <- fitted(lm(yield ~ block, data=apple))-mean(apple$yield)
pre.block.effect   <- fitted(lm(pre   ~ block, data=apple))-mean(apple$pre)
yield.block        <- apple$yield-yield.block.effect
pre.block          <- apple$pre-pre.block.effect
apple <- cbind(apple, yield.block=yield.block, pre.block=pre.block)


## Same sums of squares as apple.ancova.1 and apple.ancova.2
## for pre and treat adjusted for block
## The sum of the pre:treat and residual sum of squares is correct.
## The residual Df includes the block df and is therefore wrong.
## Therefore we suppress the residual Means Square and the F tests
apple.ancova.3 <- ancova(yield.block ~ pre.block*treat, data=apple,
                         blocks=apple$block)
tmp3 <- anova(apple.ancova.3)[,1:3]
tmp3[4,3] <- NA
tmp3
apple.ancova.3b <- ancova(yield.block ~ treat*pre.block, data=apple,
                         blocks=apple$block)
tmp3b <- anova(apple.ancova.3b)[,1:3]
tmp3b[4,3] <- NA
tmp3b

## Same sums of squares as apple.ancova.1 and apple.ancova.2
## for pre and treat adjusted for block
## The residual sum of squares is correct.
## The residual Df includes the block df and is therefore wrong.
## Therefore we suppress the residual Means Square and the F tests
apple.ancova.4 <- ancova(yield.block ~ pre.block + treat, data=apple)
tmp4 <- anova(apple.ancova.4)[,1:3]
tmp4[3,3] <- NA
tmp4
apple.ancova.4b <- ancova(yield.block ~ treat + pre.block, data=apple)
tmp4b <- anova(apple.ancova.4b)[,1:3]
tmp4b[3,3] <- NA
tmp4b

apple.ancova.6 <- ancova(yield.block ~ treat, x=pre.block, data=apple)
tmp6 <- anova(apple.ancova.6)[,1:3]
tmp6[2,3] <- NA
tmp6

predict.lm(apple.ancova.4, type="terms")
yield.block.pre <-
  apple$yield.block -
  predict.lm(apple.ancova.4, type="terms", terms="pre.block")

apple <- cbind(apple, yield.block.pre=as.vector(yield.block.pre))
apple.ancova.5 <- ancova(yield.block.pre ~ treat, x=pre.block, data=apple)
tmp5 <- anova(apple.ancova.5)[,1:2]
tmp5
if.R(r=
     attr(apple.ancova.5, "trellis")$y.limits <- attr(apple.ancova.3, "trellis")$y.limits
     , s=
     attr(apple.ancova.5, "trellis")$ylim <- attr(apple.ancova.3, "trellis")$ylim
     )
  
apple.ancova.7 <- ancova(yield.block ~ pre.block, groups=treat, data=apple)
tmp7 <- anova(apple.ancova.7)[,1:3]
tmp7[2,3] <- NA
tmp7

apple.ancova.8 <- ancova(yield ~ pre * treat, data=apple,
                         blocks=apple$block)
tmp8 <- anova(apple.ancova.8)[,1:2]
tmp8

## first step at printing all 6 panels together
simplify.legend <- function(x) {
  if.R(r={
    x$legend <- NULL
    x$sub <- NULL
  }, s={
    x$key <- NULL
    x$sub <- NULL
  })
  x
}

if.R(r={
  bot <- 3
  mid <- 2
  top <- 1
}, s={
  bot <- 1
  mid <- 2
  top <- 3
})

aa5 <- attr(apple.ancova.5, "trellis")
aa4 <- attr(apple.ancova.4, "trellis")
aa6 <- attr(apple.ancova.6, "trellis")
aa7 <- attr(apple.ancova.7, "trellis")
aa3 <- attr(apple.ancova.3, "trellis")
aa8 <- attr(apple.ancova.8, "trellis")

print(simplify.legend(aa5), split = c(1,bot,1,3), more = TRUE) # bottom of 6
print(simplify.legend(aa4), split = c(1,mid,1,3), more = TRUE) # middle
print(simplify.legend(aa6), split = c(1,top,1,3), more = FALSE)# middle
print(simplify.legend(aa7), split = c(1,bot,1,3), more = TRUE) # middle
print(simplify.legend(aa3), split = c(1,mid,1,3), more = TRUE) # middle
print(simplify.legend(aa8), split = c(1,top,1,3), more = FALSE)# top of 6
## export.eps(hh("dsgntwo/figure/apple.ancova0.eps"))



## modify trellis parameters
## second step at printing all 6 panels together
simplify.legend.labels <- function(x) {
  if.R(r={
    x$legend <- NULL
    x$sub <- NULL
    x$par.strip.text <- list(cex=.5)
    tmp.scales <- list(alternating=FALSE,
                       x=list(cex=.5),
                       y=list(at=seq(200,375,25), cex=.6))
    x <- update(x, scales=tmp.scales)
    x$main$cex <- 1
    x$xlab <- NULL
    x$ylab <- NULL
  }, s={
    x$key <- NULL
    x$sub <- NULL
    x$par.strip.text <- list(cex=1)
    tmp.scales <- list(alternating=FALSE,
                       x=list(cex=.8),
                       y=list(at=seq(200,375,25), cex=.9))
    x$scales <- tmp.scales
    x$main$cex <- 1.8
    x$xlab <- NULL
    x$ylab <- NULL
  })
  x
}



print(simplify.legend.labels(aa5), position = c(0, 0.00/1.35, 1, 0.40/1.35), more = TRUE)  # bottom of 6
print(simplify.legend.labels(aa4), position = c(0, 0.40/1.35, 1, 0.80/1.35), more = TRUE)  # middle
print(simplify.legend.labels(aa6), position = c(0, 0.80/1.35, 1, 1.20/1.35), more = FALSE) # middle
print(simplify.legend.labels(aa7), position = c(0, 0.00/1.35, 1, 0.40/1.35), more = TRUE)  # middle
print(simplify.legend.labels(aa3), position = c(0, 0.40/1.35, 1, 0.80/1.35), more = TRUE)  # middle
print(simplify.legend.labels(aa8), position = c(0, 0.80/1.35, 1, 1.35/1.35), more = FALSE) # top of 6
## export.eps(hh("dsgntwo/figure/apple.ancova.eps"))


a.xlim <- range(apple$pre, pre.block)
a.ylim <- range(apple$yield, yield.block)


a.y <- if.R(s=
            t(bwplot(block ~ yield, data=apple,
                     main="yield --- observed by block",
                     par.strip.text=list(cex=1),
                     xlim=a.ylim,
                     strip=function(...)
                     strip.default(..., strip.names = c(TRUE, TRUE))))
            ,r=
            bwplot(yield ~ block, data=apple,
                     main="yield --- observed by block",
                     par.strip.text=list(cex=1),
                     ylim=a.ylim,
                     strip=function(...)
                     strip.default(..., strip.names = c(TRUE, TRUE)))
            )

a.p <- if.R(s=
            t(bwplot(block ~ pre, data=apple,
                     main="pre --- observed by block",
                     par.strip.text=list(cex=1),
                     xlim=a.xlim,
                     strip=function(...)
                     strip.default(..., strip.names = c(TRUE, TRUE))))
            ,r=
            bwplot(pre ~ block, data=apple,
                   main="pre --- observed by block",
                   par.strip.text=list(cex=1),
                   ylim=a.xlim,
                   strip=function(...)
                   strip.default(..., strip.names = c(TRUE, TRUE)))
            )


a.y.b <- if.R(s=
              t(bwplot(block ~ yield.block, data=apple,
                       main="yield --- adjusted for block",
                       par.strip.text=list(cex=1),
                       xlim=a.ylim,
                       strip=function(...)
                       strip.default(..., strip.names = c(TRUE, TRUE))))
              ,r=
              bwplot(yield.block ~ block, data=apple,
                       main="yield --- adjusted for block",
                       par.strip.text=list(cex=1),
                       ylim=a.ylim,
                       strip=function(...)
                       strip.default(..., strip.names = c(TRUE, TRUE)))
              )

a.p.b <- if.R(s=
              t(bwplot(block ~ pre.block, data=apple,
                       main="pre --- adjusted for block",
                       par.strip.text=list(cex=1),
                  xlim=a.xlim,
                       strip=function(...)
                       strip.default(..., strip.names = c(TRUE, TRUE))))
              ,r=
              bwplot(pre.block ~ block, data=apple,
                     main="pre --- adjusted for block",
                     par.strip.text=list(cex=1),
                     ylim=a.xlim,
                     strip=function(...)
                     strip.default(..., strip.names = c(TRUE, TRUE)))
              )

if.R(r={
  bot <- 2
  top <- 1
}, s={
  bot <- 1
  top <- 2
})
print(a.y,   split = c(1,top,2,2), more = TRUE)  # left  top
print(a.p,   split = c(2,top,2,2), more = TRUE)  # right top
print(a.y.b, split = c(1,bot,2,2), more = TRUE)  # left  bottom
print(a.p.b, split = c(2,bot,2,2), more = FALSE) # right bottom
## export.eps(hh("dsgntwo/figure/apple.y.p.eps"))






## apple.ancova.2 and apple.ancova.4 have the same Sums of Squares in
## the anova table and the same regression coefficients.
summary.lm(apple.ancova.2, corr=FALSE)
summary.lm(apple.ancova.4, corr=FALSE)
## apple.ancova.2 has the correct residual df, hence Mean Squares and F tests.
## apple.ancova.4 has the wrong   residual df, hence Mean Square and F tests.

if.R(r={
  ## glht must be done with apple.ancova.2
  
  tmp <-
    glht(apple.ancova.2, linfct=mcp(treat=contrMat(table(apple$treat), type="Dunnett", base=6)))
  confint(tmp)
  plot(tmp)
  
  apple.mmc <-
    mmc(apple.ancova.2, linfct=mcp(treat=contrMat(table(apple$treat), type="Dunnett", base=6)))
  plot(apple.mmc, x.offset=12, ry=c(245,310))
  ## export.eps(hh("dsgntwo/figure/apple.mmc.eps"))
  print(apple.mmc)
  plotMatchMMC(apple.mmc$mca)
  ## export.eps(hh("dsgntwo/figure/apple.multicomp.mca.eps"))
}, s={
## multicomp must be done with apple.ancova.2

tmp <-
multicomp(apple.ancova.2, comparisons="mcc", method="dunnett", valid.check=FALSE,
          focus="treat")
tmp
plot(tmp)
## export.eps(hh("dsgntwo/figure/apple.multicomp.eps"))

## find out which rows of lmat we need
zapsmall(tmp$lmat)
## keep just the treatment rows
apple.mmc <-
multicomp.mmc(apple.ancova.2,
              comparisons="mcc", method="dunnett", valid.check=FALSE,
              focus="treat", lmat.rows=7:12, x.offset=10, plot=FALSE)
plot(apple.mmc, col.iso=16, x.offset=10)
## export.eps(hh("dsgntwo/figure/apple.mmc.eps"))
print(apple.mmc)
plotMatchMMC(apple.mmc$mca)
})
## export.eps(hh("dsgntwo/figure/apple.multicomp.mca.eps"))
