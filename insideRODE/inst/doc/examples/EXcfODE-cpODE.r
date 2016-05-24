require(lattice)
require(deSolve)
require(compiler)
require(insideRODE)
data(Theoph)
TheophODE <- Theoph
TheophODE$Dose[TheophODE$Time!=0] <- 0
TheophODE$Cmt <- rep(1,dim(TheophODE)[1])

# model files
OneComp <- list(DiffEq=list(
                            dy1dt = ~ -ka*y1 ,
                            dy2dt = ~ ka*y1-ke*y2),
                ObsEq=list(
                            c1 = ~ 0,
                            c2 = ~ y2/CL*ke),
                Parms=c("ka","ke","CL"),
                States=c("y1","y2"),
                Init=list(0,0))
                
sink("mymod.c")
cat("
/* file mymod.c */
#include <R.h>
#include <math.h>
static double parms[3];
#define ka parms[0]
#define ke parms[1]
#define CL parms[2]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
   int N=2;
   odeparms(&N, parms);
}

/* names for states and derivatives */
#define y1 y[0]
#define y2 y[1]
#define dy1 ydot[0]
#define dy2 ydot[1]
#define c1 yout[0]
#define c2 yout[1]

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    dy1 = -exp(ka)*y1;
    dy2 = exp(ka)*y1-exp(ke)*y2;
    c1 = 0.0;
    c2 = y2/exp(CL)*exp(ke);
}

/* END file mymod1.c */
",fill=TRUE)
sink()
system("RCMD SHLIB mymod.c")
dllname<-dyn.load("mymod.dll")[[1]]


TheophModelc <- cfODE(OneComp,TheophODE,dllname=dllname)

Theoph.nlme <- nlme(conc ~ TheophModelc(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)

plot(augPred(Theoph.nlme,level=0:1))
dyn.unload("mymod.dll")
unlink("mymod.dll")
unlink("mymod.o")



#cfLSODA,cfLSODE, cfODE, cfVODE SOLVER
sink("mymodff.f")
cat("
c file mymodf.f
        subroutine initmod(odeparms)
        external odeparms
        double precision parms(3)
        common /myparms/parms
        call odeparms(2, parms)
        return
        end
        subroutine derivs (neq, t, y, ydot, yout, ip)
        double precision t, y, ydot, ka, ke, CL
        integer neq, ip(*)
        dimension y(2), ydot(2), yout(2)
        common /myparms/ka,ke,CL
        ydot(1) = -exp(ka)*y(1)
        ydot(2) = exp(ka)*y(1)-exp(ke)*y(2)
        yout(1) = 0
        yout(2) = y(2)/exp(CL)*exp(ke)
        return
        end
",fill=TRUE)
sink(file = NULL, append = FALSE, type = c("output"),split = FALSE)

system("RCMD SHLIB mymodff.f")
dllname<-dyn.load("mymodff.dll")[[1]]

TheophModel <- cfODE(OneComp,TheophODE,dllname=dllname)# cpLSODA

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)

plot(augPred(Theoph.nlme,level=0:1))
dyn.unload("mymodff.dll")
unlink("mymodff.o")
unlink("mymodff.dll")
