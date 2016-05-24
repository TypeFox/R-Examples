data(hotdog, package="HH")

CpT <- ancovaplot(Sodium ~ Calories + Type, data=hotdog)
## CpT
anova(aov(Sodium ~ Calories + Type, data=hotdog))

CsT <- ancovaplot(Sodium ~ Calories * Type, data=hotdog)
## CsT
anova(aov(Sodium ~ Calories * Type, data=hotdog))

CgT <- ancovaplot(Sodium ~ Calories, groups=Type, data=hotdog)
## CgT
anova(aov(Sodium ~ Calories, data=hotdog))

TxC <- ancovaplot(Sodium ~ Type, x=Calories, data=hotdog)
## TxC
anova(aov(Sodium ~ Type, data=hotdog))

## pdf(width=9, height=7)
removeAnnotation <-
  function(x) {
    update(x,
           main=list(x$main, cex=1.1),
           ylab=NULL,
           xlab=NULL,
           scales=list(alternating=0, tck=0),
           par.strip.text=list(cex=.9, lines=1.1))
  }

## 2 x 3, with empty spots
print(position=c(.03, .31,  .53, .62), removeAnnotation(CgT), more=TRUE )
print(position=c(.50, .00, 1.00, .31), removeAnnotation(TxC), more=TRUE )
print(position=c(.50, .31, 1.00, .62), removeAnnotation(CpT), more=TRUE )
print(position=c(.50, .62, 1.00, .93), removeAnnotation(CsT), more=FALSE)

## grid writes on the current page, even though lattice has been told more=FALSE
## column labeling
grid.text(x=c(.29, .75), y=.02,
          c(expression("constant intercept" ~~ alpha),
            expression("variable intercept" ~~ alpha)))

## row labeling
grid.text(x=.02, y=c(.15, .45, .75), rot=90,
          c(expression("zero slope" ~~ beta==0),
            expression("constant slope" ~~ beta),
            expression("variable slope" ~~ beta)))

## main title
grid.text(x=.5, y=.98, gp=gpar(fontsize=15),
          "Composite graph illustrating four models with a factor and a covariate")
