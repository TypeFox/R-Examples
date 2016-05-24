### R code from vignette source 'Quick_Start.Rnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
options(prompt=" ", continue=" ")


###################################################
### code chunk number 2: model.spec (eval = FALSE)
###################################################
## B1=matrix(list("b",0,0,"b"),2,2)
## U1=matrix(0,2,1)
## Q1=matrix(c("q11","q12","q12","q22"),2,2)
## Z1=matrix(c(1,0,1,1,1,0),3,2)
## A1=matrix(list("a1",0,0),3,1)
## R1=matrix(list("r11",0,0,0,"r",0,0,0,"r"),3,3)
## pi1=matrix(0,2,1); V1=diag(1,2)
## model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)


###################################################
### code chunk number 3: marss.call (eval = FALSE)
###################################################
## fit=MARSS(data, model=model.list)


###################################################
### code chunk number 4: model.spec2 (eval = FALSE)
###################################################
## TT=dim(data)[2]
## B1=array(list(),dim=c(2,2,TT))
## B1[,,1:20]=matrix(list("b",0,0,"b_1"),2,2)
## B1[,,21:TT]=   matrix(list("b",0,0,"b_2"),2,2)


###################################################
### code chunk number 5: model.spec (eval = FALSE)
###################################################
## C1=matrix(c("temp1","temp2"),2,1)
## model.list=list(B=B1,U=U1,C=C1,c=temp,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)


###################################################
### code chunk number 6: Reset
###################################################
options(prompt="> ", continue="+ ")


