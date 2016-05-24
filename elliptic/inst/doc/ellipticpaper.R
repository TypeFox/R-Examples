### R code from vignette source 'ellipticpaper.Rnw'

###################################################
### code chunk number 1: requirepackage
###################################################
require(elliptic,quietly=TRUE)


###################################################
### code chunk number 2: setOverallImageQuality
###################################################
n <- 400
n_BACCO <- 40


###################################################
### code chunk number 3: require_packages
###################################################



###################################################
### code chunk number 4: ellipticpaper.Rnw:234-237
###################################################
require(elliptic)
require(emulator)
require(calibrator)


###################################################
### code chunk number 5: simple_usage_of_P
###################################################
z <- 1.9+1.8i
P(z,g=c(1,1i))
P(z,Omega=c(1,1i))


###################################################
### code chunk number 6: define_maxdiff
###################################################
maxdiff <- function(x,y){max(abs(x-y))}


###################################################
### code chunk number 7: laurent
###################################################
g <- c(3,2+4i)
z <- seq(from=1,to=0.4+1i,len=34)


###################################################
### code chunk number 8: maxdiff_laurent
###################################################
maxdiff(P(z,g), P.laurent(z,g))


###################################################
### code chunk number 9: abs_e18.10.9
###################################################
abs( e18.10.9(parameters(g=g)))


###################################################
### code chunk number 10: lattice_figure
###################################################
jj <- parameters(g=c(1+1i,2-3i))$Omega
latplot(jj,xlim=c(-4,4),ylim=c(-4,4),xlab="Re(z)",
     ylab="Im(z)")
polygon(Re(c(jj[1],sum(jj),jj[2],0)),
        Im(c(jj[1],sum(jj),jj[2],0)),lwd=2,col="gray90",pch=16,cex=3)

kk <- -c(3*jj[1] + 2*jj[2] , jj[1] + jj[2])  #det(matrix(c(3,2,1,1),2,2,T)==1

polygon(Re(c(kk[1],sum(kk),kk[2],0)),
        Im(c(kk[1],sum(kk),kk[2],0)),lwd=2,col="gray30",pch=16,cex=3)


###################################################
### code chunk number 11: congruence
###################################################
M <- congruence(c(4,9))


###################################################
### code chunk number 12: define_o
###################################################
o <- c(1,1i)


###################################################
### code chunk number 13: maxdiff_o
###################################################
maxdiff(g.fun(o), g.fun(M %*% o,maxiter=800))


###################################################
### code chunk number 14: u_udash
###################################################
u <- function(x){exp(pi*2i*x)}
udash <- function(x){pi*2i*exp(pi*2i*x)}
Zeta <- function(z){zeta(z,g)}
Sigma <- function(z){sigma(z,g)}
WeierstrassP <- function(z){P(z,g)}


###################################################
### code chunk number 15: integrate
###################################################
jj <- integrate.contour(Zeta,u,udash)


###################################################
### code chunk number 16: maxdiff_integrate
###################################################
maxdiff(jj, 2*pi*1i)


###################################################
### code chunk number 17: abs_integrate
###################################################
abs(integrate.contour(WeierstrassP,u,udash))


###################################################
### code chunk number 18: jj_omega
###################################################
jj.omega <- half.periods(g=c(1+1i,2-3i))


###################################################
### code chunk number 19: calculate_wp_figure
###################################################
x <- seq(from=-4, to=4, len=n)
y <- x
z <- outer(x,1i*x, "+")
f <- P(z, c(1+1i,2-3i))


###################################################
### code chunk number 20: wp_figure_file
###################################################
png("wp_figure.png",width=800,height=800)


###################################################
### code chunk number 21: wp_figure_plot
###################################################
persp(x, y, limit(Re(f)), xlab="Re(z)",ylab="Im(z)",zlab="Re(P(z))",
theta=30, phi=30, r=1e9, border=NA, shade=0.8, expand=0.6)


###################################################
### code chunk number 22: wp_figure_close
###################################################
null <- dev.off()


###################################################
### code chunk number 23: thallerfig_file
###################################################
png("thallerfig.png",width=800,height=800)


###################################################
### code chunk number 24: thallerfig_plot
###################################################
par(pty="s")
view(x,y,f,code=0,real.contour=FALSE, imag.contour=FALSE,drawlabel=FALSE,col="red",axes=FALSE,xlab="Re(z)",ylab="Im(z)")
axis(1,pos = -4)
axis(2,pos = -4)
lines(x=c(-4,4),y=c(4,4))
lines(y=c(-4,4),x=c(4,4))


###################################################
### code chunk number 25: thallerfig_close
###################################################
null <- dev.off()


###################################################
### code chunk number 26: sigma_green_calc
###################################################
x <- seq(from= -12, to=12, len=n)
y <- x
z <- outer(x, 1i*y, "+")
f <- sigma(z, c(1+1i,2-3i))


###################################################
### code chunk number 27: sigma_green_file
###################################################
png("sigma_green.png",width=800,height=800)


###################################################
### code chunk number 28: sigma_green_plot
###################################################
par(pty="s")
view(x,y,f,scheme=4,real.contour=FALSE,drawlabels=FALSE,axes=FALSE, xlab="Re(z)",ylab="Im(z)")
axis(1,pos= -12)
axis(2,pos= -12)
lines(x=c(-12,12),y=c(12,12))
lines(y=c(-12,12),x=c(12,12))
lines(x=c(-12,12),y=-c(12,12))
lines(y=c(-12,12),x=-c(12,12))


###################################################
### code chunk number 29: sigma_green_close
###################################################
null <- dev.off()


###################################################
### code chunk number 30: calculate_zeta
###################################################
zeta.z <- zeta(z, c(1+1i,2-3i))


###################################################
### code chunk number 31: zetafig_file
###################################################
png("zetafig.png",width=800,height=800)


###################################################
### code chunk number 32: zetafig_plot
###################################################
par(pty="s")
view(x,y,zeta.z,scheme=0,real.contour=FALSE,drawlabels=FALSE,xlab="Re(z)",ylab="Im(z)")


###################################################
### code chunk number 33: zetafig_close
###################################################
null <- dev.off()


###################################################
### code chunk number 34: calculate_sn
###################################################
jj <- seq(from=-40,to=40,len=n)
m <- outer(jj,1i*jj,"+")
f <- sn(u=5-2i,m=m,maxiter=1000)


###################################################
### code chunk number 35: sn_figure_file
###################################################
png("sn_figure.png",width=800,height=800)


###################################################
### code chunk number 36: sn_figure_plot
###################################################
par(pty="s")
 view(jj,jj,f,scheme=0,r0=1/5,real=T,imag=F,levels=c(0,-0.4,-1),drawlabels=F,xlab="Re(m)",ylab="Im(m)")


###################################################
### code chunk number 37: sn_figure_close
###################################################
null <- dev.off()


###################################################
### code chunk number 38: stag_calc
###################################################
     f <- function(z){1i*z^2}
     x <- seq(from=-6,to=6,len=n)
     y <- seq(from=-6,to=6,len=n)
     z <- outer(x,1i*y,"+")


###################################################
### code chunk number 39: stag_point_file
###################################################
png("stag_point.png",width=800,height=800)


###################################################
### code chunk number 40: stag_point_plot
###################################################
par(pty="s")
view(x,y,f(z),nlevels=14,imag.contour=TRUE,real.cont=TRUE,scheme=-1, 
     drawlabels=FALSE,
     axes=FALSE,xlab="Re(z)",ylab="Im(z)")
axis(1,pos=-6)
axis(2,pos=-6)
lines(x=c(-6,6),y=c(6,6))
lines(y=c(-6,6),x=c(6,6))
d1 <- c(-0.1,0,0.1)
d2 <- c( 0.1,0,0.1)
lines(x=d1,y=1+d2)
lines(x=d1,y=-1-d2)
lines(x=1-d2,y=d1)
lines(x=-1+d2,y=d1)


###################################################
### code chunk number 41: stag_point_close
###################################################
null <- dev.off()


###################################################
### code chunk number 42: two_calc
###################################################
     f <- function(z){1i*log((z-1.7+3i)*(z-1.7-3i)/(z+1-0.6i)/(z+1+0.6i))}
     x <- seq(from=-6,to=6,len=n)
     y <- seq(from=-6,to=6,len=n)
     z <- outer(x,1i*y,"+")


###################################################
### code chunk number 43: two_sources_two_sinks_file
###################################################
png("two_sources_two_sinks.png",width=800,height=800)


###################################################
### code chunk number 44: two_sources_two_sinks_plot
###################################################
par(pty="s")
view(x,y,f(z),nlevels=24,imag.contour=TRUE,real.cont=TRUE,scheme=17,power=0.1,drawlabels=FALSE,axes=FALSE,xlab="Re(z)",ylab="Im(z)")
axis(1,pos=-6)
axis(2,pos=-6)
lines(x=c(-6,6),y=c(6,6))
lines(y=c(-6,6),x=c(6,6))


###################################################
### code chunk number 45: two_sources_two_sinks_close
###################################################
null <- dev.off()


###################################################
### code chunk number 46: rect_calc3
###################################################
m <- 0.5
K <- K.fun(m)
iK <- K.fun(1-m)

#b <- sn(1.8 + 0.8i, m=m) # 1.8 to the right and 0.8 up.
#b <- 0 # middle bottom
b <- sn(0 + 1i*iK/2,m=m)  #dead centre of the rectangle.
#b <- -1 # lower left
#b <- 1/sqrt(m) # top right
#b <- -1/sqrt(m) # top left
#b <- 1e9*1i # top centre


a <- 1   #bottom right hand side corner


f <- function(z){1i*log((z-a)*(z-Conj(a))/(z-b)/(z-Conj(b)))}

     x <- seq(from=-K,to=K,len=n)
     y <- seq(from=0,to=iK,len=n)
     z <- outer(x,1i*y,"+")
    fsn <- f(sn(z,m=m))


###################################################
### code chunk number 47: rectangle_pot_flow_file
###################################################
png("rectangle_pot_flow.png",width=800,height=800)


###################################################
### code chunk number 48: rectangle_pot_flow_plot
###################################################
view(x,y,fsn,nlevels=44,imag.contour=FALSE,real.cont=TRUE,scheme=17,power=0.1,drawlabels=FALSE,axes=FALSE,xlab="",ylab="")
rect(-K,0,K,iK,lwd=3)


###################################################
### code chunk number 49: rectangle_pot_flow_close
###################################################
null <- dev.off()


###################################################
### code chunk number 50: bacco_flow
###################################################
# Choose the size of the computational mesh:
n <- n_BACCO

# Choose the number of code observations for D1:
n.code.obs <- 30

# And the number of field observations for D2:
n.field.obs <- 31

# First, up the D1 design matrix.  Recall that D1 is the set of code
# observations, which here means the observations of fluid speed when
# the sink is at a known, specified, position.

set.seed(0)

latin.hypercube <- function (n, d){
  sapply(seq_len(d), function(...) { (sample(1:n) - 0.5)/n })
}


D1.elliptic <- latin.hypercube(n.code.obs , 4)
colnames(D1.elliptic) <- c("x","y","x.sink","y.sink")
D1.elliptic[,c(1,3)] <- (D1.elliptic[,c(1,3)] -0.5)*2
#D1.elliptic[,c(2,4)] <- D1.elliptic[,c(2,4)] *iK

# now a D2 design matrix.  This is field observations: observations of
# fluid speed when the sink is at the true, unknown, specified,
# position.
D2.elliptic <- latin.hypercube(n.field.obs , 2)
colnames(D2.elliptic) <- c("x","y")
D2.elliptic[,1] <- (D2.elliptic[,1] -0.5)*2


# Now a function that, given x and y and x.sink and y.sink, returns
# the log of the fluid speed at x,y:

fluid.speed <- function(x.pos, y.pos, x.sink, y.sink){
  
  a <- 1   #bottom right hand side corner
  b <- sn(x.pos/K + 1i*iK*y.pos,m=m)  #position (x.pos , y.pos)
  f <- function(z){1i*log((z-a)*(z-Conj(a))/(z-b)/(z-Conj(b)))}
  
  x <- seq(from=-K,to=K,len=n)
  y <- seq(from=0,to=iK,len=n)
  z <- outer(x,1i*y,"+")
  potential <- f(sn(z,m=m))
  
  get.log.ke <- function(x,y,potential){
    jj <- Re(potential)
    jj.x <- cbind(jj[,-1]-jj[,-ncol(jj)],0)
    jj.y <- rbind(jj[-1,]-jj[-nrow(jj),],0)
    kinetic.energy <- jj.x^2 + jj.y^2
    n.x <- round(n * (x-(-1))/2)
    n.y <- round(n * y)
    return(log(kinetic.energy[n.x , n.y]+0.01))
  }

  return(get.log.ke(x.pos,y.pos,potential))
}

# now fill in code outputs y:
y.elliptic <- rep(NA,nrow(D1.elliptic))
for(i in 1:length(y.elliptic)){
  jj <- D1.elliptic[i,,drop=TRUE]
  y.elliptic[i] <- fluid.speed(jj[1],jj[2],jj[3],jj[4])
}


# Now do the field observations; here the source is known to be at the
# centre of the rectangle:

z.elliptic <- rep(NA,nrow(D2.elliptic))
for(i in 1:length(z.elliptic)){
  jj <- D2.elliptic[i,,drop=TRUE]
  z.elliptic[i] <- fluid.speed(jj[1],jj[2],0,0.5)
}

# Create design matrix plus observations for didactic purposes:
D1 <- round(cbind(D1.elliptic,observation=y.elliptic),2)
D2 <- round(cbind(D2.elliptic,observation=z.elliptic),2)


# create a data vector:
d.elliptic <- c(y.elliptic , z.elliptic)

#now a h1.toy() equivalent:
h1.elliptic <- function(x){
  out <- c(1,x[1])
}

#now a H1.toy() equivalent:
 H1.elliptic <- 
function (D1) 
{
    if (is.vector(D1)) {
        D1 <- t(D1)
    }
    out <- t(apply(D1, 1, h1.elliptic))
    colnames(out)[1] <- "h1.const"
    return(out)
}
                      
h2.elliptic <-
  function(x){
    c(1,x[1]) 
  }

H2.elliptic <-
  function(D2){
    if (is.vector(D2)) {
      D2 <- t(D2)
    }
    out <- t(apply(D2, 1, h2.elliptic))
    colnames(out)[1] <- "h2.const"
    return(out)
  }


#Now an extractor function:
extractor.elliptic <- 
function (D1) 
{
    return(list(x.star = D1[, 1:2, drop = FALSE], t.vec = D1[, 
        3:4, drop = FALSE]))
}

# Now a whole bunch of stuff to define a phi.fun.elliptic()
# and, after that, to call it:
phi.fun.elliptic <- 
function (rho, lambda, psi1, psi1.apriori, psi2, psi2.apriori, 
    theta.apriori, power) 
{
    "pdm.maker.psi1" <- function(psi1) {
        jj.omega_x <- diag(psi1[1:2])
        rownames(jj.omega_x) <- names(psi1[1:2])
        colnames(jj.omega_x) <- names(psi1[1:2])
        jj.omega_t <- diag(psi1[3:4])
        rownames(jj.omega_t) <- names(psi1[3:4])
        colnames(jj.omega_t) <- names(psi1[3:4])
        sigma1squared <- psi1[5]
        return(list(omega_x = jj.omega_x, omega_t = jj.omega_t, 
            sigma1squared = sigma1squared))
    }
    "pdm.maker.psi2" <- function(psi1) {
        jj.omegastar_x <- diag(psi2[1:2])
        sigma2squared <- psi2[3]
        return(list(omegastar_x = jj.omegastar_x, sigma2squared = sigma2squared))
    }
    jj.mean <- theta.apriori$mean
    jj.V_theta <- theta.apriori$sigma
    jj.discard.psi1 <- pdm.maker.psi1(psi1)
    jj.omega_t <- jj.discard.psi1$omega_t
    jj.omega_x <- jj.discard.psi1$omega_x
    jj.sigma1squared <- jj.discard.psi1$sigma1squared
    jj.discard.psi2 <- pdm.maker.psi2(psi2)
    jj.omegastar_x <- jj.discard.psi2$omegastar_x
    jj.sigma2squared <- jj.discard.psi2$sigma2squared
    jj.omega_t.upper <- chol(jj.omega_t)
    jj.omega_t.lower <- t(jj.omega_t.upper)
    jj.omega_x.upper <- chol(jj.omega_x)
    jj.omega_x.lower <- t(jj.omega_x.upper)
    jj.a <- solve(solve(jj.V_theta) + 2 * jj.omega_t, solve(jj.V_theta, 
        jj.mean))
    jj.b <- t(2 * solve(solve(jj.V_theta) + 2 * jj.omega_t) %*% 
        jj.omega_t)
    jj.c <- jj.sigma1squared/sqrt(det(diag(nrow = nrow(jj.V_theta)) + 
        2 * jj.V_theta %*% jj.omega_t))
    jj.A <- solve(jj.V_theta + solve(jj.omega_t)/4)
    jj.A.upper <- chol(jj.A)
    jj.A.lower <- t(jj.A.upper)
    list(rho = rho, lambda = lambda, psi1 = psi1, psi1.apriori = psi1.apriori, 
        psi2 = psi2, psi2.apriori = psi2.apriori, theta.apriori = theta.apriori, 
        power = power, omega_x = jj.omega_x, omega_t = jj.omega_t, 
        omegastar_x = jj.omegastar_x, sigma1squared = jj.sigma1squared, 
        sigma2squared = jj.sigma2squared, omega_x.upper = jj.omega_x.upper, 
        omega_x.lower = jj.omega_x.lower, omega_t.upper = jj.omega_t.upper, 
        omega_t.lower = jj.omega_t.lower, a = jj.a, b = jj.b, 
        c = jj.c, A = jj.A, A.upper = jj.A.upper, A.lower = jj.A.lower)
}

# OK, that's the function defined.  Now to create some jj.* variables
# to call it:

jj.psi1 <- c(rep(1,4),0.3)
names(jj.psi1)[1:4] <- colnames(D1.elliptic)
names(jj.psi1)[5] <- "sigma1squared"

jj.mean.psi1 <- rep(1,5)
names(jj.mean.psi1) <- names(jj.psi1)

jj.sigma.psi1 <- diag(0.1,nrow=5)
rownames(jj.sigma.psi1) <- names(jj.psi1)
colnames(jj.sigma.psi1) <- names(jj.psi1)

jj.psi2 <- c(1,1,0.3)
names(jj.psi2)[1:2] <- colnames(D2.elliptic)
names(jj.psi2)[3] <- "sigma2squared"

jj.mean.psi2 <- rep(1,4)
names(jj.mean.psi2) <- c("x.sink", "y.sink","rho","lambda")

jj.sigma.psi2 <- diag(0.1,4)
rownames(jj.sigma.psi2) <- names(jj.mean.psi2)
colnames(jj.sigma.psi2) <- names(jj.mean.psi2)

jj.mean.th <- c(1,0.5)
names(jj.mean.th) <- colnames(D1.elliptic)[3:4]

jj.sigma.th <- diag(rep(1,2))
rownames(jj.sigma.th) <- colnames(D1.elliptic)[3:4]
colnames(jj.sigma.th) <- colnames(D1.elliptic)[3:4]

# Now call phi.fun.elliptic():
phi.elliptic <-
  phi.fun.elliptic(
                   rho=1,
                   lambda=0.1,
                   psi1=jj.psi1,
                   psi2=jj.psi2,
                   psi1.apriori=list(mean=jj.mean.psi1, sigma=jj.sigma.psi1),
                   psi2.apriori=list(mean=jj.mean.psi2, sigma=jj.sigma.psi2),
                   theta.apriori=list(mean=jj.mean.th, sigma=jj.sigma.th),
                   power=2
                   )

# Now an E.theta.elliptic():
E.theta.elliptic <- 
function (D2 = NULL, H1 = NULL, x1 = NULL, x2 = NULL, phi, give.mean = TRUE) 
{
    if (give.mean) {
        m_theta <- phi$theta.apriori$mean
        return(H1(D1.fun(D2, t.vec = m_theta)))
    }
    else {
        out <- matrix(0, 2,2)
        rownames(out) <- c("h1.const","x")
        colnames(out) <- c("h1.const","x")
        return(out)
    }
}

#Now an Edash.theta.elliptic().  Because the basis vector is not a
#function of theta, this is a bit academic as we can use a function
#that is identical to Edash.theta.toy():

Edash.theta.elliptic <-
function (x, t.vec, k, H1, fast.but.opaque = FALSE, a = NULL, 
    b = NULL, phi = NULL) 
{
    if (fast.but.opaque) {
        edash.mean <- a + crossprod(b, t.vec[k, ])
    }
    else {
        V_theta <- phi$theta.apriori$sigma
        m_theta <- phi$theta.apriori$mean
        omega_t <- phi$omega_t
        edash.mean <- solve(solve(V_theta) + 2 * omega_t, solve(V_theta, 
            m_theta) + 2 * crossprod(omega_t, t.vec[k, ]))
    }
    jj <- as.vector(edash.mean)
    names(jj) <- rownames(edash.mean)
    edash.mean <- jj
    return(H1(D1.fun(x, edash.mean)))
}



# Define a wrapper for equation 8:
# First, calculate the constant to subtract to ensure that
# the support has a maximum of about zero: 

maximum.likelihood.support <- p.eqn8.supp(theta=c(0,1/2), D1=D1.elliptic, D2=D2.elliptic, H1=H1.elliptic, H2=H2.elliptic, d=d.elliptic, include.prior=FALSE, lognormally.distributed=FALSE, return.log=TRUE, phi=phi.elliptic)

support <- function(x){
p.eqn8.supp(theta=x, D1=D1.elliptic, D2=D2.elliptic, H1=H1.elliptic, H2=H2.elliptic, d=d.elliptic, include.prior=FALSE, lognormally.distributed=FALSE, return.log=TRUE, phi=phi.elliptic) - maximum.likelihood.support
}

#define a local function called optim() for aesthetic reasons (ie it 
# improves the appearance of the  call to optim():

optim <- 
  function(par,fn){
    stats::optim(par=par,fn=fn,control=list(fnscale = -1))$par
  }



###################################################
### code chunk number 51: head_D1
###################################################
head(D1)


###################################################
### code chunk number 52: head_D2
###################################################
head(D2)


###################################################
### code chunk number 53: calc_b_sn
###################################################
b <- sn(D1[1,3] + 1i*D1[1,4],m=m) #point corresponding to first line of D1
fsnz2 <- f(sn(z,m=m))


###################################################
### code chunk number 54: code_obs_file
###################################################
png("code_obs.png",width=800,height=800)


###################################################
### code chunk number 55: code_obs_plot
###################################################
view(x,y,fsnz2,nlevels=44,imag.contour=FALSE,real.cont=TRUE,scheme=-1,drawlabels=FALSE,axes=FALSE,xlab="",ylab="")
points(x=K*D1[1,1],y=D1[1,2]*iK,pch=4)
rect(-K,0,K,iK,lwd=3)


###################################################
### code chunk number 56: code_obs_close
###################################################
null <- dev.off()


###################################################
### code chunk number 57: calc_b2
###################################################
b <- sn(0 + 1i*iK/2,m=m) 
fsnzz <- f(sn(z,m=m))


###################################################
### code chunk number 58: true_flow_file
###################################################
png("true_flow.png",width=800,height=800)


###################################################
### code chunk number 59: true_flow_plot
###################################################
view(x,y,fsnzz,nlevels=44,imag.contour=FALSE,real.cont=TRUE,scheme=-1,drawlabels=FALSE,axes=FALSE,xlab="",ylab="")
points(x=K*D2[,1],y=D2[,2]*iK,pch=4)
rect(-K,0,K,iK,lwd=3)


###################################################
### code chunk number 60: true_flow_close
###################################################
null <- dev.off()


###################################################
### code chunk number 61: support
###################################################
support(c(0,1/2)) #centre of the rectangle
support(c(-1,1))   #top left corner


###################################################
### code chunk number 62: mle_calc
###################################################
mle <- optim(c(0,1/2),support)


###################################################
### code chunk number 63: print_mle
###################################################
mle


###################################################
### code chunk number 64: support_of_mle
###################################################
support(mle)


