## library(gmodels)
## library(lme4)
## set.seed(314159)

## sleepstudy$dayGroup <- cut(sleepstudy$Days, seq(-1,9,by=2), include=T)

## # ci example
## fm2 <- lmer(Reaction ~ dayGroup + (1|Subject) + (0+Days|Subject), sleepstudy)
## ci(fm2)


## # estimable examples
## estimable(fm2, c( 0, -1, 1, 0,  0 )  ) # list all terms
## estimable(fm2, c("dayGroup(1,3]"=-1, "dayGroup(3,5]"=1)) # just the nonzero terms
## estimable(fm2, c("dayGroup(1,3]"=-1, "dayGroup(3,5]"=1), n.sim=5000 ) # more simulations...


## # fit.contrast example
## fit.contrast( fm2, "dayGroup",
##               rbind("0-1 vs 3-4"=c(-1,0,1,0,0),
##                     "3-4 vs 5-6"=c(0,0,-1,1,0)
##                   ),
##             conf=0.95 )

## # Example from Ariel.Muldoon@oregonstate.edu
## homerange=c(
##         "male","1","fall","0.1",
## 	"male","1","winter","0.3",
## 	"male","1","spring","5.2",
## 	"male","1","summer","3.1",
## 	"male","2","fall","3.4",
## 	"male","2","winter","1.3",
## 	"male","2","spring","4.8",
## 	"male","2","summer","4.3",
## 	"male","3","fall","3.9",
## 	"male","3","winter","3.8",
## 	"male","3","spring","5.7",
## 	"male","3","summer","2.0",
## 	"male","4","fall","3.7",
## 	"male","4","winter","4.3",
## 	"male","4","spring","6.0",
## 	"male","4","summer","1.8",
## 	"female","5","fall","4.3",
## 	"female","5","winter","1.9",
## 	"female","5","spring","7.2",
## 	"female","5","summer","6.9",
## 	"female","6","fall","5.3",
## 	"female","6","winter","4.3",
## 	"female","6","spring","6.2",
## 	"female","6","summer","4.8",
## 	"female","7","fall","7.1",
## 	"female","7","winter","4.9",
## 	"female","7","spring","8.3",
## 	"female","7","summer","7.7"
## 	)

## homerange <- data.frame(matrix(homerange,ncol=4, byrow=T))
## names(homerange) <- c("sex", "animal", "season", "area")
## homerange$area = as.numeric(as.character(homerange$area))

## fit1 <- lmer(area ~ sex*season + (1|animal), data=homerange)
## summary(fit1)
## anova(fit1)


## #matrix to give estimable for making estimates

## spr <- rbind(c(1,0,1,0,0,0,0,0),
##              c(1,1,0,0,0,1,0,0),
##              c(1,0,0,1,0,0,0,0),
##              c(1,1,0,0,0,0,1,0),
##              c(1,0,0,0,1,0,0,0),
##              c(1,1,0,0,0,0,0,1))


## estimable(fit1, spr)
