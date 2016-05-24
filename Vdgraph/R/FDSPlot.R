FDSPlot <- function(des, mod=2) 
{
#####################################################################
# This function creates a fraction of design space plot for a design
# stored in the matrix or data frame des
#  Inputs;
#    des - a matrix or data frame containing a experimental design
#          in coded or uncoded units. There should be one column
#          for each column in the matrix and one row for each run 
#          in the design. The maximum number of columns allowed is 
#          7.
#
#    mod - the model to be represented:
#          0 = linear model
#          1 = linear main effects plus linear by linear 2-factor
#               interactions
#          2 = full quadratic response surface model (default)
####################################################################

#First check to see if des is a matrix, and if so convert it to a data frame
     if(is.matrix(des)) {des <- as.data.frame(des)}
#Next check to see if des is centered, and if not center it
   chkm <-apply(des, 2, f)
  if (mean(chkm) !=0) des<-scale(des,scale=FALSE)
#Check to see if scaled, and if not scale it
   chkm <-apply(des, 2, mx)
   if(mean(chkm) !=1) des<-scale(des, center=FALSE, scale=chkm)
   rc <-dim(des)
   ndpts<-rc[1]
   kvar1<-rc[2]
  if (kvar1<2) {
    stop("The minimum number of independent variables must be 2","\n")
               }
  if (kvar1>7) {
    stop("The maximum number of independent variables allowed is 7","\n")
               }
des<-as.data.frame(des)
# get the model
y<-runif(nrow(des),0,1)
### Case for two variables #############
 if(kvar1==2) {
names<-list("x1","x2")
names(des)<-names
lmod<-lm(y~x1+x2,data=des)
imod<-lm(y~x1+x2+x1:x2,data=des)
qmod<-lm(y~x1+x2+I(x1^2)+I(x2^2)+x1:x2,data=des)

 if(mod==0){
X<-model.matrix(lmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
fX<-cbind(Int,X1,X2)
           }

 if(mod==1){
X<-model.matrix(imod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X1X2<-X1*X2
fX<-cbind(Int,X1,X2,X1X2)
           }
 if(mod==2){
X<-model.matrix(qmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X1sq<-X1*X1
X2sq<-X2*X2
X1X2<-X1*X2
fX<-cbind(Int,X1,X2,X1sq,X2sq,X1X2)
            }
              }
# Case for 3 variables ########################
 if(kvar1==3)    {
names<-list("x1","x2","x3")
names(des)<-names
lmod<-lm(y~x1+x2+x3,data=des)
imod<-lm(y~x1+x2+x3+x1:x2+x1:x3+x2:x3,data=des)
qmod<-lm(y~x1+x2+x3+I(x1^2)+I(x2^2)+I(x3^2)+x1:x2+x1:x3+x2:x3,data=des)
             
#get the design matrix and fX
 if(mod==0) {
X<-model.matrix(lmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
fX<-cbind(Int,X1,X2,X3)
             }

 if(mod==1)  {
X<-model.matrix(imod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X1X2<-X1*X2
X1X3<-X1*X3
X2X3<-X2*X3
fX<-cbind(Int,X1,X2,X3,X1X2,X1X3,X2X3)
            }

 if(mod==2) {
X<-model.matrix(qmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X1sq<-X1*X1
X2sq<-X2*X2
X3sq<-X3*X3
X1X2<-X1*X2
X1X3<-X1*X3
X2X3<-X2*X3
fX<-cbind(Int,X1,X2,X3,X1sq,X2sq,X3sq,X1X2,X1X3,X2X3)
            }
               }

# Case for 4 variables ########################
 if(kvar1==4)    {
names<-list("x1","x2","x3","x4")
names(des)<-names
lmod<-lm(y~x1+x2+x3+x4,data=des)
imod<-lm(y~x1+x2+x3+x4+x1:x2+x1:x3+x1:x4+x2:x3+x2:x4+x3:x4,data=des)
qmod<-lm(y~x1+x2+x3+x4+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+x1:x2+x1:x3+x1:x4+x2:x3+x2:x4+x3:x4,data=des)
             
#get the design matrix and fX
 if(mod==0) {
X<-model.matrix(lmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
fX<-cbind(Int,X1,X2,X3,X4)
             }

 if(mod==1)  {
X<-model.matrix(imod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X2X3<-X2*X3
X2X4<-X2*X4
X3X4<-X3*X4
fX<-cbind(Int,X1,X2,X3,X4,X1X2,X1X3,X1X4,X2X3,X2X4,X3X4)
            }
 if(mod==2) {
X<-model.matrix(qmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X1sq<-X1*X1
X2sq<-X2*X2
X3sq<-X3*X3
X4sq<-X4*X4
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X2X3<-X2*X3
X2X4<-X2*X4
X3X4<-X3*X4
fX<-cbind(Int,X1,X2,X3,X4,X1sq,X2sq,X3sq,X4sq,X1X2,X1X3,X1X4,X2X3,X2X4,X3X4)
            }
               }
# Case for 5 variables ########################
 if(kvar1==5)    {
names<-list("x1","x2","x3","x4","x5")
names(des)<-names
lmod<-lm(y~x1+x2+x3+x4+x5,data=des)
imod<-lm(y~x1+x2+x3+x4+x5+x1:x2+x1:x3+x1:x4+x1:x5+x2:x3+x2:x4+x2:x5+x3:x4+x3:x5+x4:x5,data=des)
qmod<-lm(y~x1+x2+x3+x4+x5+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+x1:x2+x1:x3+x1:x4+x1:x5+x2:x3+x2:x4+x2:x5+x3:x4+x3:x5+x4:x5,data=des)
             
#get the design matrix and fX
 if(mod==0) {
X<-model.matrix(lmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
fX<-cbind(Int,X1,X2,X3,X4,X5)
             }

 if(mod==1)  {
X<-model.matrix(imod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X1X5<-X1*X5
X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X3X4<-X3*X4
X3X5<-X3*X5
X4X5<-X4*X5
fX<-cbind(Int,X1,X2,X3,X4,+X5+X1X2,X1X3,X1X4,X1X5,X2X3,X2X4,X2X5,X3X4,X3X5,X4X5)
            }
 if(mod==2) {
X<-model.matrix(qmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X1sq<-X1*X1
X2sq<-X2*X2
X3sq<-X3*X3
X4sq<-X4*X4
X5sq<-X5*X5
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X1X5<-X1*X5
X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X3X4<-X3*X4
X3X5<-X3*X5
X4X5<-X4*X5
fX<-cbind(Int,X1,X2,X3,X4,X5,X1sq,X2sq,X3sq,X4sq,X5sq,X1X2,X1X3,X1X4,X1X5,X2X3,X2X4,X2X5,X3X4,X3X5,X4X5)
            }
               }
# Case for 6 variables ########################
 if(kvar1==6)    {
names<-list("x1","x2","x3","x4","x5","x6")
names(des)<-names
lmod<-lm(y~x1+x2+x3+x4+x5+x6,data=des)
imod<-lm(y~x1+x2+x3+x4+x5+x6+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x2:x3+x2:x4+x2:x5+x2:x6+x3:x4+x3:x5+x3:x6+x4:x5+x4:x6+x5:x6,data=des)
qmod<-lm(y~x1+x2+x3+x4+x5+x6+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+I(x6^2)+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x2:x3+x2:x4+x2:x5+x2:x6+x3:x4+x3:x5+x3:x6+x4:x5+x4:x6+x5:x6,data=des)
             
#get the design matrix and fX
 if(mod==0) {
X<-model.matrix(lmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X6<-runif(5000,-1,1)
fX<-cbind(Int,X1,X2,X3,X4,X5,X6)
             }

 if(mod==1)  {
X<-model.matrix(imod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X6<-runif(5000,-1,1)
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X1X5<-X1*X5
X1X6<-X1*X6
X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X2X6<-X2*X6
X3X4<-X3*X4
X3X5<-X3*X5
X3X6<-X3*X6
X4X5<-X4*X5
X4X6<-X4*X6
X5X6<-X5*X6
fX<-cbind(Int,X1,X2,X3,X4,X5,X6,X1X2,X1X3,X1X4,X1X5,X1X6,X2X3,X2X4,X2X5,X2X6,X3X4,X3X5,X3X6,X4X5,X4X6,X5X6)
            }
 if(mod==2) {
X<-model.matrix(qmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X6<-runif(5000,-1,1)
X1sq<-X1*X1
X2sq<-X2*X2
X3sq<-X3*X3
X4sq<-X4*X4
X5sq<-X5*X5
X6sq<-X6*X6
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X1X5<-X1*X5
X1X6<-X1*X6
X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X2X6<-X2*X6
X3X4<-X3*X4
X3X5<-X3*X5
X3X6<-X3*X6
X4X5<-X4*X5
X4X6<-X4*X6
X5X6<-X5*X6
fX<-cbind(Int,X1,X2,X3,X4,X5,X6,X1sq,X2sq,X3sq,X4sq,X5sq,X6sq,X1X2,X1X3,X1X4,X1X5,X1X6,X2X3,X2X4,X2X5,X2X6,X3X4,X3X5,X3X6,X4X5,X4X6,X5X6)
            }
               }
# Case for 7 variables ########################
 if(kvar1==7)    {
names<-list("x1","x2","x3","x4","x5","x6","x7")
names(des)<-names
lmod<-lm(y~x1+x2+x3+x4+x5+x6+x7,data=des)
imod<-lm(y~x1+x2+x3+x4+x5+x6+x7+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x1:x7+x2:x3+x2:x4+x2:x5+x2:x6+x2:x7+x3:x4+x3:x5+x3:x6+x3:x7+x4:x5+x4:x6+x4:x7+x5:x6+x5:x7+x6:x7,data=des)
qmod<-lm(y~x1+x2+x3+x4+x5+x6+x7+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+I(x6^2)+I(x7^2)+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x1:x7+x2:x3+x2:x4+x2:x5+x2:x6+x2:x7+x3:x4+x3:x5+x3:x6+x3:x7+x4:x5+x4:x6+x4:x7+x5:x6+x5:x7+x6:x7,data=des)
             
#get the design matrix and fX
 if(mod==0) {
X<-model.matrix(lmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X6<-runif(5000,-1,1)
X7<-runif(5000,-1,1)
fX<-cbind(Int,X1,X2,X3,X4,X5,X6,X7)
             }

 if(mod==1)  {
X<-model.matrix(imod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X6<-runif(5000,-1,1)
X7<-runif(5000,-1,1)
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X1X5<-X1*X5
X1X6<-X1*X6
X1X7<-X1*X7
X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X2X6<-X2*X6
X2X7<-X2*X7
X3X4<-X3*X4
X3X5<-X3*X5
X3X6<-X3*X6
X3X7<-X3*X7
X4X5<-X4*X5
X4X6<-X4*X6
X4X7<-X4*X7
X5X6<-X5*X6
X5X7<-X5*X7
X6X7<-X6*X7
fX<-cbind(Int,X1,X2,X3,X4,X5,X6,X7,X1X2,X1X3,X1X4,X1X5,X1X6,X1X7,X2X3,X2X4,X2X5,X2X6,X2X7,X3X4,X3X5,X3X6,X3X7,X4X5,X4X6,X4X7,X5X6,X5X7,X6X7)
            }
 if(mod==2) {
X<-model.matrix(qmod)
XtXI<-solve(t(X)%*%X)
Int<-c(rep(1,5000))
X1<-runif(5000,-1,1)
X2<-runif(5000,-1,1)
X3<-runif(5000,-1,1)
X4<-runif(5000,-1,1)
X5<-runif(5000,-1,1)
X6<-runif(5000,-1,1)
X7<-runif(5000,-1,1)
X1sq<-X1*X1
X2sq<-X2*X2
X3sq<-X3*X3
X4sq<-X4*X4
X5sq<-X5*X5
X6sq<-X6*X6
X7sq<-X7*X7
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X1X5<-X1*X5
X1X6<-X1*X6
X1X7<-X1*X7
X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X2X6<-X2*X6
X2X7<-X2*X7
X3X4<-X3*X4
X3X5<-X3*X5
X3X6<-X3*X6
X3X7<-X3*X7
X4X5<-X4*X5
X4X6<-X4*X6
X4X7<-X4*X7
X5X6<-X5*X6
X5X7<-X5*X7
X6X7<-X6*X7
fX<-cbind(Int,X1,X2,X3,X4,X5,X6,X7,X1sq,X2sq,X3sq,X4sq,X5sq,X6sq,X7sq,X1X2,X1X3,X1X4,X1X5,X1X6,X1X7,X2X3,X2X4,X2X5,X2X6,X2X7,X3X4,X3X5,X3X6,X3X7,X4X5,X4X6,X4X7,X5X6,X5X7,X6X7)
            }
               }

# get v and vi
v<-diag(fX%*%XtXI%*%t(fX))
vi<-sort(v)
# calculate fraction of design space             
fs<-(1:length(vi))/length(vi)
# make FDS Plot
plot(fs,vi,type='l',ylim=c(0.0,(max(vi)+.15)),xlab="Fraction of Space",ylab="Relative Prediction Variance")
abline(v=.5,lty=3)
abline(h=(vi[2500]+vi[2501])/2,lty=3)
}

