#
#   code used to adjust parameter lambda in awsimage and awspimage  
#   values are currently adjusted for default settings of lkern and plateau   
#   and additive independent Gaussian noise 
#   
#   deviations from these settings may require to in-/decrease lambda by a factor specified in ladjust
#
library(adimpro)
set.seed(1);timgc <- make.image(array(as.integer(rnorm(512^2*3,32000,6000)),c(512,512,3)),gammatype="None")
set.seed(1);timgg <- make.image(array(as.integer(rnorm(512^2,32000,6000)),c(512,512)),gammatype="None")

# local constant color  alpha=0.09
unix.time(t0imghatInf <- awsimage(timgc,hmax=10,aws=FALSE,earlystop=FALSE,homogen=FALSE))
r0inf <- sqrt(var(as.vector(extract.image(t0imghatInf))))
unix.time(t0imghat1 <- awsimage(timgc,hmax=10,earlystop=FALSE,homogen=FALSE,ladjust=1))
r01 <- sqrt(var(as.vector(extract.image(t0imghat1))))
r01/r0inf-1

# local constant greyvalue   alpha=0.1
unix.time(t0imghatInf <- awsimage(timgg,hmax=10,aws=FALSE,earlystop=FALSE,homogen=FALSE))
r0inf <- sqrt(var(as.vector(extract.image(t0imghatInf))))
unix.time(t0imghat1 <- awsimage(timgg,hmax=10,earlystop=FALSE,homogen=FALSE,ladjust=1))
r01 <- sqrt(var(as.vector(extract.image(t0imghat1))))
r01/r0inf-1

# local polynomial (degree=2) color  alpha=0.3

unix.time(timghatInf <- awspimage(timgc,hmax=20,degree=2,aws=FALSE,earlystop=FALSE,homogen=FALSE))
rinf <- sqrt(var(as.vector(extract.image(timghatInf))))
unix.time(timghat1 <- awspimage(timgc,hmax=20,degree=2,earlystop=FALSE,homogen=FALSE,ladjust=1))
r1 <- sqrt(var(as.vector(extract.image(timghat1))))
r1/rinf-1

# local polynomial (degree=2) greyvalue    alpha=0.2

unix.time(timghatInf <- awspimage(timgg,hmax=20,degree=2,aws=FALSE,earlystop=FALSE,homogen=FALSE))
rinf <- sqrt(var(as.vector(extract.image(timghatInf))))
unix.time(timghat1 <- awspimage(timgg,hmax=20,degree=2,earlystop=FALSE,homogen=FALSE,ladjust=1))
r1 <- sqrt(var(as.vector(extract.image(timghat1))))
r1/rinf-1

# local polynomial (degree=1) color  alpha=0.2

unix.time(timghatInf <- awspimage(timgc,hmax=15,degree=1,aws=FALSE,earlystop=FALSE,homogen=FALSE))
rinf <- sqrt(var(as.vector(extract.image(timghatInf))))
unix.time(timghat1 <- awspimage(timgc,hmax=15,degree=1,earlystop=FALSE,homogen=FALSE,ladjust=1))
r1 <- sqrt(var(as.vector(extract.image(timghat1))))
r1/rinf-1

# local polynomial (degree=1) greyvalue  alpha=0.2

unix.time(timghatInf <- awspimage(timgg,hmax=15,degree=1,aws=FALSE,earlystop=FALSE,homogen=FALSE))
rinf <- sqrt(var(as.vector(extract.image(timghatInf))))
unix.time(timghat1 <- awspimage(timgg,hmax=15,degree=1,earlystop=FALSE,homogen=FALSE,ladjust=1))
r1 <- sqrt(var(as.vector(extract.image(timghat1))))
r1/rinf-1

# local constant color  alpha=0.2
unix.time(t0imghatInf <- awsaniso(timgc,hmax=10,aws=FALSE))
r0inf <- sqrt(var(as.vector(extract.image(t0imghatInf))))
unix.time(t0imghat1 <- awsaniso(timgc,hmax=10,ladjust=1.25))
r01 <- sqrt(var(as.vector(extract.image(t0imghat1))))
r01/r0inf-1

# local constant greyvalue   alpha=0.19
unix.time(t0imghatInf <- awsaniso(timgg,hmax=10,aws=FALSE))
r0inf <- sqrt(var(as.vector(extract.image(t0imghatInf))))
unix.time(t0imghat1 <- awsaniso(timgg,hmax=10,ladjust=1.5))
r01 <- sqrt(var(as.vector(extract.image(t0imghat1))))
r01/r0inf-1
