reg.poly <-
function(resp, treat, DFerror, SSerror, DFtreat, SStreat) { 

MSerror<-SSerror/DFerror

cat('\nAdjustment of polynomial models of regression\n------------------------------------------------------------------------\n')

X<-matrix(1,length(treat),4)
X[,2]<-treat
X[,3]<-treat^2
X[,4]<-treat^3
mean.table<-tapply.stat(resp,treat,mean)
colnames(mean.table)<-c('  Levels','   Observed Means')

###############################################################################
#reta
###############################################################################
b=ginv(t(X[,1:2])%*%X[,1:2], tol=.Machine$double.eps)%*%t(X[,1:2])%*%resp
ep = sqrt(diag(ginv(t(X[,1:2])%*%X[,1:2], tol=.Machine$double.eps)*MSerror))
tc = b/ep
pv = 2*pt(abs(tc),DFerror,lower.tail=FALSE)     
tm1<-data.frame('Estimate' = round(b,8),'Standard Error' = round(ep,5),'tc'=round(tc,5),'p-value' = round(pv,5))
rownames(tm1)<-c('b0','b1')

aov.m1<-anova(lm(resp~treat))
if(dim(mean.table)[1]==2){r2m1<-1}
if(dim(mean.table)[1]>2) {r2m1<-aov.m1[1,2]/SStreat}

#ANAVA da regressao linear
nomes1<-c("Linear Effect","Lack of fit","Residuals")
anava1<-data.frame('DF'=c(1,(gld=c(DFtreat-1)),(glr=DFerror)),
                  'SS'=c(round(c(aov.m1[[2]][1],(sqd=c(SStreat-aov.m1[[2]][1])),SSerror),5)),
                  'MS'=c(round(c(aov.m1[[3]][1],(if(gld==0){qmd=0} else{qmd=(sqd/gld)}), MSerror),5)),
                  'Fc'=c(round(c((fcl=aov.m1[[3]][1]/MSerror),(fc=qmd/MSerror)),2),''),
                  'p-value'=c(round(c(pf(fcl,1,glr,lower.tail=FALSE),(if(gld==0){pv=1}else{pv=1-pf(fc,gld,glr)})),5),''))
rownames(anava1)<-nomes1

output1<-list('Linear Model
------------------------------------------------------------------------' = tm1,
              'R2 of linear model' = r2m1,
              'Analysis of Variance of linear model' = anava1)
print(output1,right=TRUE)
cat('------------------------------------------------------------------------\n')
###############################################################################
  if(dim(mean.table)[1]>2) {
###############################################################################
#parabola
###############################################################################
b2=ginv(t(X[,1:3])%*%X[,1:3], tol=.Machine$double.eps)%*%t(X[,1:3])%*%resp
ep2 = sqrt(diag(ginv(t(X[,1:3])%*%X[,1:3], tol=.Machine$double.eps)*MSerror))
tc2 = b2/ep2
pv2 = 2*pt(abs(tc2),DFerror,lower.tail=FALSE)  
tm2<-data.frame('Estimate' = round(b2,8),'Standard Error' = round(ep2,5),'tc'=round(tc2,5),'p-value' = round(pv2,5))
rownames(tm2)<-c('b0','b1','b2')

t2<-treat^2
aov.m2<-anova(lm(resp~treat+t2))
if(dim(mean.table)[1]==3){r2m2<-1}
if(dim(mean.table)[1]>3) {r2m2<-(aov.m2[1,2]+aov.m2[2,2])/SStreat}

#ANAVA da regressao quadratica
nomes2<-c("Linear Effect","Quadratic Effect","Lack of fit","Residuals")
anava2<-data.frame('DF'=c(aov.m2[[1]][1:2],(gld=c(DFtreat-2)),(glr=DFerror)),
                  'SS'=c(round(c(aov.m2[[2]][1:2],(sqd=c(SStreat-sum(aov.m2[[2]][1:2]))),SSerror),5)),
                  'MS'=c(round(c(aov.m2[[3]][1:2],(if(gld==0){qmd=0} else{qmd=(sqd/gld)}), MSerror),5)),
                  'Fc'=c(round(c((fcl=aov.m2[[3]][1:2]/MSerror),(fc=qmd/MSerror)),2),''),
                  'p-value'=c(round(c(pf(fcl,1,glr,lower.tail=FALSE),(if(gld==0){pv=1}else{pv=1-pf(fc,gld,glr)})),5),''))
rownames(anava2)<-nomes2

output2<-list('Quadratic Model
------------------------------------------------------------------------' = tm2,
              'R2 of quadratic model' = r2m2,
              'Analysis of Variance of quadratic model' = anava2)
print(output2,right=TRUE)
cat('------------------------------------------------------------------------\n')

                            }
###############################################################################
  if(dim(mean.table)[1]>3) {
###############################################################################
#cubica
###############################################################################
b3=ginv(t(X[,1:4])%*%X[,1:4], tol=.Machine$double.eps)%*%t(X[,1:4])%*%resp
ep3 = sqrt(diag(ginv(t(X[,1:4])%*%X[,1:4], tol=.Machine$double.eps)*MSerror))
tc3 = b3/ep3
pv3 = 2*pt(abs(tc3),DFerror,lower.tail=FALSE)  
tm3<-data.frame('Estimate' = round(b3,8),'Standard Error' = round(ep3,5),'tc'=round(tc3,5),'p-value' = round(pv3,5))
rownames(tm3)<-c('b0','b1','b2','b3')

t3<-treat^3
aov.m3<-anova(lm(resp~treat+t2+t3))
if(dim(mean.table)[1]==4){r2m3<-1}
if(dim(mean.table)[1]>4) {r2m3<-(aov.m3[1,2]+aov.m3[2,2]+aov.m3[3,2])/SStreat}

#ANAVA da regressao cubica
nomes3<-c("Linear Effect","Quadratic Effect","Cubic Effect","Lack of fit","Residuals")
anava3<-data.frame('DF'=c(aov.m3[[1]][1:3],(gld=c(DFtreat-3)),(glr=DFerror)),
                  'SS'=c(round(c(aov.m3[[2]][1:3],(sqd=c(SStreat-sum(aov.m3[[2]][1:3]))),SSerror),5)),
                  'MS'=c(round(c(aov.m3[[3]][1:3],(if(gld==0){qmd=0} else{qmd=(sqd/gld)}), MSerror),5)),
                  'Fc'=c(round(c((fcl=aov.m3[[3]][1:3]/MSerror),(fc=qmd/MSerror)),2),''),
                  'p-value'=c(round(c(pf(fcl,1,glr,lower.tail=FALSE),(if(gld==0){pv=1} else {pv=1-pf(fc,gld,glr)})),5),''))
rownames(anava3)<-nomes3

output3<-list('Cubic Model
------------------------------------------------------------------------' = tm3,
             'R2 of cubic model' = r2m3, 
             'Analysis of Variance of cubic model' = anava3
              )
print(output3,right=TRUE)
cat('------------------------------------------------------------------------\n')
                            }
###############################################################################
print(mean.table)

cat('------------------------------------------------------------------------\n\n')

}
