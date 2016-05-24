## stage.classify.R

## copied from demo graphics
#if(dev.cur() <= 1) get(getOption("device"))()

#opar <- par(ask = interactive() &&
#            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))


options(digits=3)

# ---------------------------------------------------------------- #
# Data sets

data(aq.trans)
data(aq.census)

head(aq.trans)
#  ADD logical vectors for survival and flowering

aq.trans$survived<-aq.trans$fate!="dead"

aq.trans$flowered<-aq.trans$fate=="flower"

head(aq.trans,1)


## pause AFTER tables
caption<-function (x) 
{
    cat("Press <Enter> to continue...")
    readline()
    invisible()
}


# ---------------------------------------------------------------- #
# TRANSITION FREQUENCY TABLES
surv <-subset(aq.trans, stage %in% c("flower", "large", "small") & !is.na(rose) )

## show how cut() works
table(surv$rose)
table(cut(surv$rose,c(-1:7,45)))


# DEAD 
z<-table(cut(surv$rose,c(0:7,45), labels=c(1:7, "8+")), surv$stage=="flower", surv$survived, dnn=c("State", "Flower", "Alive"))

ftable(z, row.vars=c(2,1), col.vars=c(3))

caption("Transition frequency table with initial state (flowering and rosette number) and fate (dead or alive)")

y<-table(cut(surv$rose,c(0:7,45), labels=c(1:7, "8+")),
        cut(surv$rose2,c(0:7,48), labels=c(1:7, "8+")),
        surv$stage=="flower", surv$flowered, dnn=c("State", "Fate", "Flower", "Flower"))

ftable(y, row.vars=c(3,1), col.vars=c(4,2))

caption("Transition frequency table with initial state and fate of live plants")





# ---------------------------------------------------------------- #
#  Logistic regression survival vs. leaf



a<-subset(aq.trans, leaf<50 & stage!="recruit", c(leaf,survived))

logi.hist.plot(a$leaf,  a$survived, 
type="hist", boxp=FALSE, counts=TRUE, las.h=0, int=10, 
ylabel="Survival probability", ylabel2="Number of plants", xlab="Number of leaves", main="Logistic regression leaf vs. survival")


# ---------------------------------------------------------------- #
# SPLIT into FLOWERING AND VEGETATIVE states
## rosettes not counted before 1998!

flwr<-subset(surv, rose>0 & rose<=7 & stage=="flower" & year>1997)
non.flwr<-subset(surv, rose>0 & rose<=7 & stage!="flower" & year>1997)


# counts
rbind(flower=table(flwr$rose), vegetative=table(non.flwr$rose))

# ---------------------------------------------------------------- #
# Check FLOWERING sample sizes of potential stage class boundaries using cut



a<-table( flwr$year, cut(flwr$rose,c(0,1,3,10), labels=c("rose1", "rose2-3", "rose3+")) )


a

# or just  two flowering classes

table( flwr$year, cut(flwr$rose,c(0,2,10), labels=c("rose1-2", "rose3+")))

caption("Potential stage class boundaries and sample sizes of flowering plants using cut")




# ---------------------------------------------------------------- #
# SURVIVAL rates by rosette number and initial flowering state  USING proportions

 
tapply(flwr$survived, flwr$rose, mean)

tapply(non.flwr$survived, non.flwr$rose, mean)
 
##  TWO columns for binomial regression (counts of dead and live plants)

x.flwr<-table(flwr$rose, flwr$survived)
x.flwr
## non flower
x.non.flwr<-table(non.flwr$rose, non.flwr$survived)


rose<-1:7

glm.flwr<-glm(x.flwr~rose, binomial)   
glm.non.flwr<-glm(x.non.flwr~rose, binomial) 

summary(glm.flwr)

xrange<-c(.5,7.5)

## plot
plot(rose, prop.table(x.flwr,1)[,2], ylim=c(0.45,1), xlim=xrange, pch=1,
      col="royalblue", las=1, 
          ylab=expression(paste("Proportion surviving (year ", italic("t+1"), ")")),
     xlab=expression(paste("Rosette number (year ", italic("t"), ")")),
     main="Survival rates by flowering state" )


points(rose, prop.table(x.non.flwr,1)[,2], pch=3, col="darkgreen")


xp<-seq(1,max(rose),.1)
lines(xp, 1-predict(glm.flwr, list(rose=xp), type="response"), col="royalblue", lwd=2)
lines(xp, 1-predict(glm.non.flwr, list(rose=xp), type="response"), col="darkgreen", lwd=2)


## add space to increase box size
legend(.5,1, c("Flowering", "Vegetative "), pch=c(1,3),  col=c("royalblue", "darkgreen"))


# ---------------------------------------------------------------- #
#  GROWTH rates by rosette number and initial flowering state
#  1. probability of FLOWERING (ie, GROW into reproductive stage) the next year 


x.flwr<-table(flwr$rose, flwr$flowered)
x.non.flwr<-table(non.flwr$rose, non.flwr$flowered)


glm.flwr<-glm(x.flwr~rose, binomial)        
glm.non.flwr<-glm(x.non.flwr~rose, binomial) 

plot(rose, prop.table(x.flwr,1)[,2], ylim=c(0,1), xlim=xrange, pch=1, col="royalblue", las=1,
     ylab=expression(paste("Proportion flowering (year ", italic("t+1"), ")")),
     xlab=expression(paste("Rosette number (year ", italic("t"), ")")),
     main="Growth into reproductive stage")


points(rose, prop.table(x.non.flwr,1)[,2], pch=3, col="darkgreen")


lines(xp, 1-predict(glm.flwr, list(rose=xp), type="response"), col="royalblue", lwd=2)
lines(xp, 1-predict(glm.non.flwr, list(rose=xp), type="response"), col="darkgreen", lwd=2)


legend(.5,1, c("Flowering", "Vegetative "), pch=c(1,3),  col=c("royalblue", "darkgreen"))


# ---------------------------------------------------------------- #
#  2. probability of growing new rosettes (growth), losing rosettes (regression),
# or staying the same (stasis)) the next year
     
# ADD column growth to dataframe 

flwr$growth<-flwr$rose2-flwr$rose
non.flwr$growth<-non.flwr$rose2-non.flwr$rose

table(flwr$growth)

#  CUT into groups corresponding to  0=statis, <0= regression, >0 = growth??

x.flwr<-table(flwr$rose, cut(flwr$growth,c(-10,-.5,.5,10), labels=c("regress", "stasis", "growth")))
#non-flower
x.non.flwr<-table(non.flwr$rose, cut(non.flwr$growth,c(-15,-.5,.5,15), labels=c("regress", "stasis", "growth")))

x.flwr

caption("Growth rates of flowering plants by rosette number" )


## TWO columns for binomial regression (successes and failures)

a1.f<-cbind(regress=x.flwr[,1],    b=x.flwr[,2]+x.flwr[,3])
a2.f<-cbind(stasis= x.flwr[,2],    b=x.flwr[,1]+x.flwr[,3])
a3.f<-cbind(grow=   x.flwr[,3],    b=x.flwr[,1]+x.flwr[,2])

a1.nf<-cbind(regress=x.non.flwr[,1],  b=x.non.flwr[,2]+x.non.flwr[,3])
a2.nf<-cbind(stasis =x.non.flwr[,2],  b=x.non.flwr[,1]+x.non.flwr[,3])
a3.nf<-cbind(grow   =x.non.flwr[,3],  b=x.non.flwr[,1]+x.non.flwr[,2])

# INITIAL flowering
m1.f<-glm(a1.f~rose, binomial)
m2.f<-glm(a2.f~rose, binomial)
m3.f<-glm(a3.f~rose, binomial)

## proportions for easier plotting

p.flwr<-prop.table(x.flwr,1)
p.non.flwr<-prop.table(x.non.flwr,1)
 
plot(rose, p.flwr[,1], type="n", ylim=c(0,1),xlim=xrange, las=1,
     ylab=expression(paste("Proportion (year ", italic("t+1"), ")")),
     xlab=expression(paste("Rosette number (year ", italic("t"), ")")),
     main=paste("Growth, stasis, and regression\nof flowering plants"))

points(rose, p.flwr[,1], pch=1, col="darkgreen")
lines(xp, predict(m1.f, list(rose=xp), type="response"), col="darkgreen", lwd=2)

points(rose, p.flwr[,2],  pch=3, col="red")
lines(xp, predict(m2.f, list(rose=xp), type="response"), col="red", lwd=2)

points(rose, p.flwr[,3], pch=6, col="royalblue" )
lines(xp, predict(m3.f, list(rose=xp), type="response"), col="royalblue" , lwd=2)

legend(2.5,1,c("Growth", "Stasis", "Regress  "), pch=c(6,3,1), col=c("royalblue", "red","darkgreen")  )



# ---------------------------------------------------------------- #
#  NON-flowering growth rates


m1.nf<-glm(a1.nf~rose, binomial)
m2.nf<-glm(a2.nf~rose, binomial)
m3.nf<-glm(a3.nf~rose, binomial)

plot(rose, p.non.flwr[,1], type="n", ylim=c(0,1),xlim=xrange, las=1,
     ylab=expression(paste("Proportion (year ", italic("t+1"), ")")),
     xlab=expression(paste("Rosette number (year ", italic("t"), ")")),
         main=paste("Growth, stasis, and regression\nof vegetative plants"))

points(rose, p.non.flwr[,1], pch=1, col="darkgreen")
lines(xp, predict(m1.nf, list(rose=xp), type="response"), col="darkgreen", lwd=2)

points(rose, p.non.flwr[,2],  pch=3, col="red")
lines(xp, predict(m2.nf, list(rose=xp), type="response"), col="red", lwd=2)

points(rose, p.non.flwr[,3], pch=6, col="royalblue" )
lines(xp, predict(m3.nf, list(rose=xp), type="response"), col="royalblue" , lwd=2)


legend(2,1,c("Growth", "Stasis", "Regress  "), pch=c(6,3,1), col=c("royalblue", "red","darkgreen")  )



# ---------------------------------------------------------------- #
# REPRODUCTION - number of fruits by rosette (complete fruit surveys at FILLMORE Canyon only)

fruit<- subset(aq.trans, stage=="flower" & fruits>=0 & rose<=7 & rose>0, select=c(plant, year, rose,leaf, fruits) )

table(fruit$rose)  

fruit1<- glm(fruits~rose, poisson, data=fruit)
summary(fruit1)

caption("glm using poisson family")

## correct for overdispersion

fruit1<- glm(fruits~rose, quasipoisson, data=fruit)
summary(fruit1)

caption("glm using quasipoisson to correct for overdispersion")

## PLOT methods?  side-by-side boxplots  or try plotCI in gplots

boxplot(fruit$fruits ~ fruit$rose, xlab="Rosette number", ylab="Mature fruits", col="green" , ylim=c(0,12), main=paste("Side-by-side boxplots of\nmature fruits by plant size"))



lines(xp, exp(predict(fruit1,list(rose=xp))), col="red", lwd=2)




#par(opar)







