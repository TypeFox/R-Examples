RE.Johnson <-
function(x){

sort.x=sort(x)
#######VARIABLES######
z<-seq(0.25,1.25,0.01) # z values.
QR<-matrix(0,1,101)
q<-matrix(0,101,4) #quartile
j<-matrix(0,101,4) # element of x relative to q
y<-matrix(0,101,4)

xl=matrix(0,101,1)
xm=matrix(0,101,1)
xu=matrix(0,101,1)

b.eta=matrix(0,101,1)
b.gamma=matrix(0,101,1)
b.lambda=matrix(0,101,1)
b.epsilon=matrix(0,101,1)

l.eta=matrix(0,101,1)
l.gamma=matrix(0,101,1)
l.lambda=matrix(0,101,1)
l.epsilon=matrix(0,101,1)

u.eta=matrix(0,101,1)
u.gamma=matrix(0,101,1)
u.lambda=matrix(0,101,1)
u.epsilon=matrix(0,101,1)

xsb=matrix(0,length(x),101)
xsl=matrix(0,length(x),101)
xsu=matrix(0,length(x),101)

xsb.valida=matrix(0,1,101)
xsl.valida=matrix(0,1,101)
xsu.valida=matrix(0,1,101)

xsb.adtest=matrix(0,1,101)
xsl.adtest=matrix(0,1,101)
xsu.adtest=matrix(0,1,101)

f.gamma<-0
f.lambda<-0
f.epsilon<-0
f.eta<-0

#################################

for(i in 1:101) {
  q[i,1]<-pnorm(-3*z[i])
  q[i,2]<-pnorm(-1*z[i])
  q[i,3]<-pnorm(1*z[i])
  q[i,4]<-pnorm(3*z[i])

j[i,1]<-length(x)*q[i,1]+0.5 
j[i,2]<-length(x)*q[i,2]+0.5 
j[i,3]<-length(x)*q[i,3]+0.5 
j[i,4]<-length(x)*q[i,4]+0.5 
  
ifelse(j[i,1]<1,(y[i,1]<-min(sort.x)),(y[i,1]<-((sort.x[ceiling(j[i,1])]-sort.x[floor(j[i,1])])/(ceiling(j[i,1])-floor(j[i,1])))*(j[i,1]-floor(j[i,1]))+sort.x[floor(j[i,1])]))
ifelse(j[i,2]>length(x),(y[i,2]<-max(sort.x)),(y[i,2]<-((sort.x[ceiling(j[i,2])]-sort.x[floor(j[i,2])])/(ceiling(j[i,2])-floor(j[i,2])))*(j[i,2]-floor(j[i,2]))+sort.x[floor(j[i,2])]))
ifelse(j[i,3]>length(x),(y[i,3]<-max(sort.x)),(y[i,3]<-((sort.x[ceiling(j[i,3])]-sort.x[floor(j[i,3])])/(ceiling(j[i,3])-floor(j[i,3])))*(j[i,3]-floor(j[i,3]))+sort.x[floor(j[i,3])]))
ifelse(j[i,4]>length(x),(y[i,4]<-max(sort.x)),(y[i,4]<-((sort.x[ceiling(j[i,4])]-sort.x[floor(j[i,4])])/(ceiling(j[i,4])-floor(j[i,4])))*(j[i,4]-floor(j[i,4]))+sort.x[floor(j[i,4])]))


QR[i]<-((y[i,4]-y[i,3])*(y[i,2]-y[i,1]))/((y[i,3]-y[i,2])^2)

xl[i]<-y[i,2]-y[i,1]
xm[i]<-y[i,3]-y[i,2]
xu[i]<-y[i,4]-y[i,3]

}

######### SB,SL,SU
for(i in 1:101) {

ifelse((.5*(((1+xm[i]/xu[i])*(1+xm[i]/xl[i]))^.5))<1,(b.eta[i,1]<--1000),(b.eta[i,1]<-z[i]/(acosh(.5*(((1+xm[i]/xu[i])*(1+xm[i]/xl[i]))^.5))))) 
ifelse((.5*(((1+xm[i]/xu[i])*(1+xm[i]/xl[i]))^.5))<1,(b.gamma[i,1]<--1000),(b.gamma[i,1]<-b.eta[i,1]*asinh(((xm[i]/xl[i]-xm[i]/xu[i])*(((1+xm[i]/xu[i])*(1+xm[i]/xl[i])-4)^.5))/(2*(((xm[i]^2)/(xl[i]*xu[i]))-1))))) 
ifelse((((((1+xm[i]/xu[i])*(1+xm[i]/xl[i])-2)^2)-4))<0,b.lambda[i,1]<-1000,b.lambda[i,1]<-(xm[i]*(((((1+xm[i]/xu[i])*(1+xm[i]/xl[i])-2)^2)-4)^.5))/(((xm[i]^2)/(xl[i]*xu[i]))-1))
ifelse((((((1+xm[i]/xu[i])*(1+xm[i]/xl[i])-2)^2)-4))<0,b.epsilon[i,1]<-1000,b.epsilon[i,1]<-.5*(y[i,2]+y[i,3]-b.lambda[i]+((xm[i]*(xm[i]/xl[i]-xm[i]/xu[i]))/(((xm[i]^2)/(xl[i]*xu[i]))-1))))


l.eta[i]<- 2*z[i]/(log(xu[i]/xm[i]))
ifelse((xu[i]/xm[i]-1)/((xu[i]*xm[i])^.5)<=0,l.gamma[i]<-1000,l.gamma[i]<- l.eta[i]*log((xu[i]/xm[i]-1)/((xu[i]*xm[i])^.5)))
l.epsilon[i]<- .5*(y[i,2]+y[i,3]-xm[i]*((xu[i]/xm[i]+1)/(xu[i]/xm[i]-1)))


ifelse((.5*(xu[i]/xm[i]+xl[i]/xm[i]))<1,(u.eta[i,1]<--1000),u.eta[i,1]<-2*z[i]/(acosh(.5*(xu[i]/xm[i]+xl[i]/xm[i]))))
ifelse( ((xu[i]*xl[i]/(xm[i]^2))-1)<0,(u.gamma[i,1]<--1000),(u.gamma[i,1]<-u.eta[i]*asinh((xl[i]/xm[i]-xu[i]/xm[i])/(2*(((xu[i]*xl[i]/(xm[i]^2))-1)^.5))))) 
ifelse(((xu[i]*xl[i]/(xm[i]^2))-1)<0,u.lambda[i,1]<-1000,u.lambda[i]<- (2*xm[i]*(((xu[i]*xl[i]/(xm[i]^2))-1)^.5))/((xu[i]/xm[i]+xl[i]/xm[i]-2)*((xu[i]/xm[i]+xl[i]/xm[i]+2)^.5)))
u.epsilon[i]<- .5*(y[i,2]+y[i,3]+((xm[i]*(xl[i]/xm[i]-xu[i]/xm[i]))/((xu[i]/xm[i]+xl[i]/xm[i])-2)))

for(o in 1:length(x)) {
 ifelse(((x[o]-b.epsilon[i])/(b.lambda[i]+b.epsilon[i]-x[o]))<=0,xsb.valida[1,i]<- xsb.valida[1,i]+1,xsb[o,i]<-b.gamma[i]+b.eta[i]*log((x[o]-b.epsilon[i])/(b.lambda[i]+b.epsilon[i]-x[o])))
 ifelse((x[o]-l.epsilon[i])<=0,xsl.valida[1,i]<- xsl.valida[1,i]+1,xsl[o,i]<- l.gamma[i]+l.eta[i]*log(x[o]-l.epsilon[i]))
 ifelse(((xu[i]*xl[i]/(xm[i]^2))-1)<0,xsu.valida[1,i]<- xsu.valida[1,i]+1,xsu[o,i]<- u.gamma[i]+u.eta[i]*asinh((x[o]-u.epsilon[i])/u.lambda[i]))
 
 }
if(xsb.valida[1,i]==0) xsb.adtest[1,i]<-(RE.ADT(xsb[,i])$p)
if(xsl.valida[1,i]==0) xsl.adtest[1,i]<-(RE.ADT(xsl[,i])$p)
if(xsu.valida[1,i]==0) xsu.adtest[1,i]<-(RE.ADT(xsu[,i])$p)
}
#####
	xsb.adtest[which(is.na(xsb.adtest))] <- 0 # insertion by Boda Martin to handle with NA/NaN occurences
	xsl.adtest[which(is.na(xsl.adtest))] <- 0 # insertion by Boda Martin to handle with NA/NaN occurences
	xsu.adtest[which(is.na(xsu.adtest))] <- 0 # insertion by Boda Martin to handle with NA/NaN occurences



ifelse((max(xsb.adtest)>max(xsl.adtest)& max(xsb.adtest)>max(xsu.adtest)),{p<-max(xsb.adtest);fun<-"SB";transformed<-xsb[,max.col(xsb.adtest)];f.gamma<-b.gamma[max.col(xsb.adtest)];f.lambda<-b.lambda[max.col(xsb.adtest)];f.epsilon<-b.epsilon[max.col(xsb.adtest)];f.eta<-b.eta[max.col(xsb.adtest)]} ,ifelse(max(xsl.adtest)>max(xsu.adtest),{p<-max(xsl.adtest);fun<-"SL";transformed<-xsl[,max.col(xsl.adtest)];f.gamma<-l.gamma[max.col(xsl.adtest)];f.lambda<-l.lambda[max.col(xsl.adtest)];f.epsilon<-l.epsilon[max.col(xsl.adtest)];f.eta<-l.eta[max.col(xsl.adtest)]},{p<-max(xsu.adtest); fun<-"SU";transformed<-xsu[,max.col(xsu.adtest)];f.gamma<-u.gamma[max.col(xsu.adtest)];f.lambda<-u.lambda[max.col(xsu.adtest)];f.epsilon<-u.epsilon[max.col(xsu.adtest)];f.eta<-u.eta[max.col(xsu.adtest)]}))

outList = list("Johnson Transformation","function"=fun,p=p,transformed=transformed,f.gamma=f.gamma,f.lambda=f.lambda,f.epsilon=f.epsilon,f.eta=f.eta)
    invisible(outList)

}




