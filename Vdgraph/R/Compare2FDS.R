Compare2FDS <- function(des1, des2, name1, name2, mod=2) 
{
#####################################################################
# This function creates a fraction of design space plot for a design
# stored in the matrix or data frame des
#  Inputs;
#   des1 - a matrix or data frame containing a experimental design
#          in coded or uncoded units. There should be one column
#          for each column in the matrix and one row for each run 
#          in the design. The maximum number of columns allowed is 
#          7.
#   des2 - a matrix or data frame containing a experimental design
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
#
# First check calculate the coordinates for des1
####################################################################
#First check to see if des1 is a matrix, and if so convert it to a data frame
   if(is.matrix(des1)) des1 <- as.data.frame(des1)
#Next check to see if des1 is centered, and if not center it
   chkm <-apply(des1, 2, f)
  if (mean(chkm) !=0) des1<-scale(des1, scale=FALSE)
#Check to see if scaled, and if not scale it
   chkm <-apply(des1, 2, mx)
   if(mean(chkm) !=1) des1<-scale(des1, center=FALSE, scale=chkm)
   rc <-dim(des1)
   ndpts<-rc[1]
   kvar1<-rc[2]
  if (kvar1<2) {
    stop("The minimum number of independent variables must be 2","\n")
               }
  if (kvar1>7) {
    stop("The maximum number of independent variables allowed is 7","\n")
               }
des1<-as.data.frame(des1)
# get the model
y<-runif(nrow(des1),0,1)
### Case for two variables #############
 if(kvar1==2) {
names<-list("x1","x2")
names(des1)<-names
lmod<-lm(y~x1+x2,data=des1)
imod<-lm(y~x1+x2+x1:x2,data=des1)
qmod<-lm(y~x1+x2+I(x1^2)+I(x2^2)+x1:x2,data=des1)

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
names(des1)<-names
lmod<-lm(y~x1+x2+x3,data=des1)
imod<-lm(y~x1+x2+x3+x1:x2+x1:x3+x2:x3,data=des1)
qmod<-lm(y~x1+x2+x3+I(x1^2)+I(x2^2)+I(x3^2)+x1:x2+x1:x3+x2:x3,data=des1)
             
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
names(des1)<-names
lmod<-lm(y~x1+x2+x3+x4,data=des1)
imod<-lm(y~x1+x2+x3+x4+x1:x2+x1:x3+x1:x4+x2:x3+x2:x4+x3:x4,data=des1)
qmod<-lm(y~x1+x2+x3+x4+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+x1:x2+x1:x3+x1:x4+x2:x3+x2:x4+x3:x4,data=des1)
             
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
names(des1)<-names
lmod<-lm(y~x1+x2+x3+x4+x5,data=des1)
imod<-lm(y~x1+x2+x3+x4+x5+x1:x2+x1:x3+x1:x4+x1:x5+x2:x3+x2:x4+x2:x5+x3:x4+x3:x5+x4:x5,data=des1)
qmod<-lm(y~x1+x2+x3+x4+x5+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+x1:x2+x1:x3+x1:x4+x1:x5+x2:x3+x2:x4+x2:x5+x3:x4+x3:x5+x4:x5,data=des1)
             
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
names(des1)<-names
lmod<-lm(y~x1+x2+x3+x4+x5+x6,data=des1)
imod<-lm(y~x1+x2+x3+x4+x5+x6+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x2:x3+x2:x4+x2:x5+x2:x6+x3:x4+x3:x5+x3:x6+x4:x5+x4:x6+x5:x6,data=des1)
qmod<-lm(y~x1+x2+x3+x4+x5+x6+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+I(x6^2)+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x2:x3+x2:x4+x2:x5+x2:x6+x3:x4+x3:x5+x3:x6+x4:x5+x4:x6+x5:x6,data=des1)
             
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
names(des1)<-names
lmod<-lm(y~x1+x2+x3+x4+x5+x6+x7,data=des1)
imod<-lm(y~x1+x2+x3+x4+x5+x6+x7+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x1:x7+x2:x3+x2:x4+x2:x5+x2:x6+x2:x7+x3:x4+x3:x5+x3:x6+x3:x7+x4:x5+x4:x6+x4:x7+x5:x6+x5:x7+x6:x7,data=des1)
qmod<-lm(y~x1+x2+x3+x4+x5+x6+x7+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+I(x6^2)+I(x7^2)+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x1:x7+x2:x3+x2:x4+x2:x5+x2:x6+x2:x7+x3:x4+x3:x5+x3:x6+x3:x7+x4:x5+x4:x6+x4:x7+x5:x6+x5:x7+x6:x7,data=des1)
             
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

# get v1 and v1i
v1<-diag(fX%*%XtXI%*%t(fX))
v1i<-sort(v1)
# calculate fraction of design space             
fs1<-(1:length(v1i))/length(v1i)
####################################################################
#
# Next check calculate the coordinates for des2
####################################################################
#First check to see if des2 is a matrix, and if so convert it to a data frame
   if(is.matrix(des2)) des2 <- as.data.frame(des2)
#Next check to see if des2 is centered, and if not center it
   chkm <-apply(des2, 2, f)
  if (mean(chkm) !=0) des2<-scale(des2, scale=FALSE)
#Check to see if scaled, and if not scale it
   chkm <-apply(des2, 2, mx)
   if(mean(chkm) !=1) des2<-scale(des2, center=FALSE, scale=chkm)
   rc <-dim(des2)
   ndpts<-rc[1]
   kvar1<-rc[2]
  if (kvar1<2) {
    stop("The minimum number of independent variables must be 2","\n")
               }
  if (kvar1>7) {
    stop("The maximum number of independent variables allowed is 7","\n")
               }
des2<-as.data.frame(des2)
# get the model
y<-runif(nrow(des2),0,1)
### Case for two variables #############
 if(kvar1==2) {
names<-list("x1","x2")
names(des2)<-names
lmod<-lm(y~x1+x2,data=des2)
imod<-lm(y~x1+x2+x1:x2,data=des2)
qmod<-lm(y~x1+x2+I(x1^2)+I(x2^2)+x1:x2,data=des2)

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
names(des2)<-names
lmod<-lm(y~x1+x2+x3,data=des2)
imod<-lm(y~x1+x2+x3+x1:x2+x1:x3+x2:x3,data=des2)
qmod<-lm(y~x1+x2+x3+I(x1^2)+I(x2^2)+I(x3^2)+x1:x2+x1:x3+x2:x3,data=des2)
             
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
names(des2)<-names
lmod<-lm(y~x1+x2+x3+x4,data=des2)
imod<-lm(y~x1+x2+x3+x4+x1:x2+x1:x3+x1:x4+x2:x3+x2:x4+x3:x4,data=des2)
qmod<-lm(y~x1+x2+x3+x4+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+x1:x2+x1:x3+x1:x4+x2:x3+x2:x4+x3:x4,data=des2)
             
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
names(des2)<-names
lmod<-lm(y~x1+x2+x3+x4+x5,data=des2)
imod<-lm(y~x1+x2+x3+x4+x5+x1:x2+x1:x3+x1:x4+x1:x5+x2:x3+x2:x4+x2:x5+x3:x4+x3:x5+x4:x5,data=des2)
qmod<-lm(y~x1+x2+x3+x4+x5+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+x1:x2+x1:x3+x1:x4+x1:x5+x2:x3+x2:x4+x2:x5+x3:x4+x3:x5+x4:x5,data=des2)
             
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
names(des2)<-names
lmod<-lm(y~x1+x2+x3+x4+x5+x6,data=des2)
imod<-lm(y~x1+x2+x3+x4+x5+x6+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x2:x3+x2:x4+x2:x5+x2:x6+x3:x4+x3:x5+x3:x6+x4:x5+x4:x6+x5:x6,data=des2)
qmod<-lm(y~x1+x2+x3+x4+x5+x6+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+I(x6^2)+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x2:x3+x2:x4+x2:x5+x2:x6+x3:x4+x3:x5+x3:x6+x4:x5+x4:x6+x5:x6,data=des2)
             
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
names(des2)<-names
lmod<-lm(y~x1+x2+x3+x4+x5+x6+x7,data=des2)
imod<-lm(y~x1+x2+x3+x4+x5+x6+x7+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x1:x7+x2:x3+x2:x4+x2:x5+x2:x6+x2:x7+x3:x4+x3:x5+x3:x6+x3:x7+x4:x5+x4:x6+x4:x7+x5:x6+x5:x7+x6:x7,data=des2)
qmod<-lm(y~x1+x2+x3+x4+x5+x6+x7+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)+I(x6^2)+I(x7^2)+x1:x2+x1:x3+x1:x4+x1:x5+x1:x6+x1:x7+x2:x3+x2:x4+x2:x5+x2:x6+x2:x7+x3:x4+x3:x5+x3:x6+x3:x7+x4:x5+x4:x6+x4:x7+x5:x6+x5:x7+x6:x7,data=des2)
             
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

# get v2 and v2i
v2<-diag(fX%*%XtXI%*%t(fX))
v2i<-sort(v2)
# calculate fraction of design space             
fs2<-(1:length(v2i))/length(v2i)

# make Comparison FDS Plots
maxv<-max(max(v1i),max(v2i))
plot(fs1,v1i,type='l',ylim=c(0.0,(maxv+.15)),xlab="Fraction of Space",ylab="Relative Prediction Variance")
abline(v=.5,lty=3)
abline(h=(v1i[2500]+v1i[2501])/2,lty=3)
lines(fs2,v2i,col="red",lty=2)
abline(h=(v2i[2500]+v2i[2501])/2,lty=3)
lines(c(0,.10),c(maxv,maxv),lty=1,col="black")
text(.11,maxv,name1,adj=-.1)
lines(c(0,.10),c(maxv-maxv/20,maxv-maxv/20),lty=3,col="red")
text(.11,maxv-maxv/20,name2,col="red",adj=-.1)
}
