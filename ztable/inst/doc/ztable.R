## ----,results='asis'-----------------------------------------------------
require(ztable)
options(ztable.type="html")
z=ztable(head(iris))
z

## ----,results='asis'-----------------------------------------------------
z=ztable(head(iris),align="cccccc")
z

## ----,results='asis'-----------------------------------------------------
cgroup=c("Sepal","Petal","Species")
n.cgroup=c(2,2,1)
z=addcgroup(z,cgroup=cgroup,n.cgroup=n.cgroup)
z

## ----,results='asis'-----------------------------------------------------
rgroup=c("OneToThree","Four","FiveToSix")
n.rgroup=c(3,1,2)

z=addrgroup(z,rgroup=rgroup,n.rgroup=n.rgroup,cspan.rgroup=1)
z
print(z,type="latex")

## ----,results='asis'-----------------------------------------------------
ncount=c(123,120,123,124)
sub=paste("(N=",ncount,")",sep="")
z=addSubColNames(z,c(sub,NA))
z

## ----,results='asis'-----------------------------------------------------
z=spanRow(z,col=2,from=4,to=7,"orange")
z=spanRow(z,col=3,from=5,to=7,"platinum")
z=spanRow(z,col=4,from=6,to=7,"cyan")
z=spanRow(z,col=5,from=5,to=7,"yellow")
z=spanRow(z,col=6,from=3,to=5,"yellow")
z

z=spanCol(z,row=2,from=3,to=4,"yellow")
z=spanCol(z,row=3,from=4,to=5,"lightblue")
z

## ----,results='asis'-----------------------------------------------------
vlines(z,type="all")       # type=1 gets same result
z=vlines(z,type="none")      # type=0 gets same result
z
z=vlines(z,add=c(1,2,5))
z

## ----,results='asis'-----------------------------------------------------
t1=head(iris,10)[,c(1,3,5)]
t2=tail(iris,10)[,c(1,3,5)]
t=cbind(t1,t2)
z=ztable(t,caption="Table 1. Top 10 and Last 10 Data from iris",align="ccccccc")
z

## ----,results='asis'-----------------------------------------------------
cgroup=c("Top 10","Last 10")
n.cgroup=c(3,3)
z=addcgroup(z,cgroup=cgroup,n.cgroup=n.cgroup)
z 
rgroup=c("Top 1-3","Top 4-5",NA," Top 7-10")
n.rgroup=c(3,2,1,4)
z=addrgroup(z,rgroup=rgroup,n.rgroup=n.rgroup,cspan.rgroup=1)
z
z=addRowColor(z,c(5,10),"pink")
z=addColColor(z,4,"amber")
z=addCellColor(z,rows=c(5,10),cols=4,"orange")
z
z=spanCol(z,row=2,from=2,to=3,color="lightcyan")
z=spanRow(z,col=7,from=7,to=8,color="cyan")
z
hlines(z,type=1)

## ----,results='asis'-----------------------------------------------------
vlines(z,type=0)  # No vertical lines
vlines(z,type=1)  # Vertical lines for all column

## ----, eval=FALSE--------------------------------------------------------
#  options(ztable.type="html")

## ----,results="asis",message=FALSE---------------------------------------
require(ztable)
options(ztable.type="html")
options(ztable.zebra=1)
options(ztable.zebra.color="platinum")
options(ztable.colnames.bold=TRUE)
ztable(head(mtcars))

## ----,results='asis'-----------------------------------------------------
ztable(head(mtcars),zebra=NULL,size=3,
       caption="Table 1. Non-zebra Table with small size")

## ----,results='asis'-----------------------------------------------------
ztable(head(mtcars[c(1:7)]),zebra=2,zebra.color="lightcyan",size=7,
       caption="Table 2. Left-sided caption at botom with large font",
       caption.placement="bottom",caption.position="l") 

## ----,results="asis"-----------------------------------------------------
out <- aov(mpg ~ ., data=mtcars)
ztable(out)

## ----,results='asis'-----------------------------------------------------
fit <- lm(mpg ~ cyl + disp + wt + drat + am, data=mtcars)
ztable(fit)

## ----,results='asis'-----------------------------------------------------
a=anova(fit)
ztable(a)

## ----,results='asis'-----------------------------------------------------
fit2 <- lm(mpg ~ cyl+wt, data=mtcars)
b=anova(fit2,fit)
ztable(b)
ztable(b,show.heading=FALSE)

## ----,results='asis',warning=FALSE---------------------------------------
require(survival)
data(colon)
attach(colon)
out <- glm(status ~ rx+obstruct+adhere+nodes+extent, data=colon, family=binomial)
ztable(out)

## ----,results='asis'-----------------------------------------------------
ztable(anova(out))

## ----,results='asis'-----------------------------------------------------
op <- options(contrasts = c("contr.helmert", "contr.poly"))
npk.aov <- aov(yield ~ block + N*P*K, npk) 
ztable(npk.aov,zebra=1)

## ----,results='asis'-----------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
ztable(lm.D9)
ztable(anova(lm.D9),align="|c|rrrr|r|")

## ----,results='asis'-----------------------------------------------------
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
d.AD <- data.frame(treatment, outcome, counts)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
ztable(glm.D93)

## ----,results='asis',message=FALSE---------------------------------------
data(USArrests)
pr1 <- prcomp(USArrests) 
ztable(pr1)
ztable(summary(pr1))

## ----,results='asis',message=FALSE---------------------------------------
colon$TS = Surv(time,status==1) 
out=coxph(TS~rx+obstruct+adhere+differ+extent+surg+node4,data=colon)
ztable(out)

## ----,comment=NA---------------------------------------------------------
require(graphics)

DNase1 <- subset(DNase, Run == 1)

## using a selfStart model
fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal),DNase1)
summary(fm1DNase1)

## ----,results='asis',message=FALSE---------------------------------------
ztable(fm1DNase1)

## ----,results='asis'-----------------------------------------------------
require(MASS)
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
a=fitdistr(x, "gamma")
ztable(a)
x3 <- rweibull(100, shape = 4, scale = 100)
b=fitdistr(x3, "weibull")
ztable(b)

## ----,results='asis',message=FALSE---------------------------------------
ztable(head(mtcars,15),zebra=0,zebra.color=NULL) 

## ----,results='asis'-----------------------------------------------------
z1=ztable(head(iris),zebra=2)
z1
print(z1,zebra.type=2)
print(z1,zebra=1,zebra.type=2,zebra.colnames=TRUE)

## ----,results='asis'-----------------------------------------------------
options(ztable.zebra.color=NULL)
(z1=ztable(head(iris),zebra=0,zebra.type=2))

## ----,results='asis'-----------------------------------------------------
update_ztable(z1,colnames.bold=TRUE,zebra.colnames=TRUE)

## ----,results='asis'-----------------------------------------------------
print(z1,zebra.color=c(rep("white",5),"peach"),zebra.colnames=TRUE)

## ----,results='asis'-----------------------------------------------------
ztable(head(iris),zebra=0,zebra.type=0)
ztable(head(iris),zebra=0,zebra.type=0,zebra.color=zcolors$name,zebra.colnames=TRUE)

## ----,results='asis'-----------------------------------------------------
ztable(head(iris),zebra=0,zebra.type=0,zebra.color=1:7,zebra.colnames=TRUE)
ztable(head(mtcars[,1:9]),zebra=0,zebra.type=0,zebra.color=1:9,zebra.colnames=TRUE)

## ----,results='asis'-----------------------------------------------------
mycolor=rep("white",6)
for(i in 1:149){
    mycolor=c(mycolor,"white",zcolors$name[((i-1)*5+1):((i-1)*5+5)])
}
mycolor=c(mycolor,"white",zcolors$name[c(746:749,1)])
a=c(zcolors$name[1:5])
for(i in 2:149){
    a=rbind(a,zcolors$name[((i-1)*5+1):((i-1)*5+5)])
}
a=rbind(a,zcolors$name[c(746:749,1)])
a=data.frame(a,stringsAsFactors=FALSE,row.names=NULL)
ztable(a,zebra=0,zebra.type=0,zebra.color=mycolor,include.rownames=FALSE,
       include.colnames=FALSE,longtable=TRUE)

## ----,results='asis'-----------------------------------------------------
z=ztable(head(mtcars[1:3]),tabular=TRUE,zebra.color="peach-orange")
z1=ztable(head(iris[1:3]),tabular=TRUE,zebra=2)

parallelTables(width=c(0.5,0.5),list(z,z1),type="html")
parallelTables(width=c(0.5,0.5),list(z,"figures/ztable3.png"),type="html")

## ----,results='asis'-----------------------------------------------------
require(moonBook)
res=mytable(Dx~.,data=acs)
options(ztable.zebra=NULL)
z=ztable(res)
z
vlines(z,type="all")


## ----,results='asis'-----------------------------------------------------
res1=mytable(sex+DM~.,data=acs)
z=ztable(res1)
z
vlines(z,type="all")

## ----,results='asis'-----------------------------------------------------
z=addRowColor(z,c(13,16),"platinum")
z=addColColor(z,c(5,8),"pink")
z=addCellColor(z,rows=16,cols=c(5,8),color="orange")
z=addCellColor(z,rows=13,cols=5,color="orange")
z

