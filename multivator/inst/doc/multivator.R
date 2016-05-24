### R code from vignette source 'multivator.Rnw'

###################################################
### code chunk number 1: set_seed_chunk
###################################################
set.seed(0)
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: time_saver
###################################################
calc_from_scratch <- TRUE


###################################################
### code chunk number 3: multivator.Rnw:117-120
###################################################
ignore <- require("multivator",quietly=TRUE)
ignore <- require("mvtnorm",quietly=TRUE)
ignore <- require("emulator",quietly=TRUE)


###################################################
### code chunk number 4: show_two_processes
###################################################
show_diff_processes <- function(){
  jj <- seq(from=0, to=10, by=0.4)

  roughness <- c(3,0.5)

  jjm <- matrix(c(jj,jj))
  colnames(jjm) <- "a"
  jjn <- c("foo","bar")

  jjt <- factor(rep(jjn,each=length(jj)) ,levels=jjn)

  mm <- mdm(xold=jjm,types=jjt)

  co <- 1

  hh <- mhp(M=matrix(c(1,co,co,1),2,2), B=array(roughness,c(1,1,2)),names="a",levels=jjn)
  
  
  LoF <- default_LoF(mm)

  set.seed(0)
  o <- obs_maker(mm,hh,LoF,beta=c(1,0,1,0))
  o <- matrix(o,ncol=2)
  matplot(o,type='o',pch=16,ylab="independent variable",xlab="dependent variable")
  abline(0,0)
  
  jjx <- 22   # x coord of LHS of segments
  jjy <- 3.7  # length of segments
  jjd <- 0.1  # separation of segments
  jje <- 0.03 # length of serifs
  jjl <- 5    # length of text
  
  jj1 <- jjx + 2/roughness[1]^0.5
  jj2 <- jjx + 2/roughness[2]^0.5
  
  ## First the segments:
  segments(jjx,jjy-jjd,jj1,jjy-jjd,lty=1,col='black')
  segments(jjx,jjy+jjd,jj2,jjy+jjd,lty=2,col='red'  )

  ## Then the left serif:
  segments(jjx,jjy-jjd-jje,jjx,jjy-jjd+jje,lty=1,col='black')
  segments(jjx,jjy+jjd-jje,jjx,jjy+jjd+jje,lty=2,col='red')
 
  ## Then the right serif:
  segments(jj1,jjy-jjd-jje,jj1,jjy-jjd+jje,lty=1,col='black')
  segments(jj2,jjy+jjd-jje,jj2,jjy+jjd+jje,lty=2,col='red')
  
  text(22-jjl,jjy,"roughness lengths")
}

show_diff_processes()


###################################################
### code chunk number 5: setopts
###################################################
options(digits=3)
set.seed(0) 


###################################################
### code chunk number 6: toy_design_matrix
###################################################
data("mtoys")
head(toy_mm)


###################################################
### code chunk number 7: show_toy_d
###################################################
head(toy_d)


###################################################
### code chunk number 8: showpoint
###################################################
toy_point


###################################################
### code chunk number 9: show_multem
###################################################
(e <- multem(toy_point, toy_expt, toy_mhp, toy_LoF, give = TRUE))


###################################################
### code chunk number 10: sample.from.posterior
###################################################
rmvnorm(n=5,mean=e$mstar,sigma=e$cstar)


###################################################
### code chunk number 11: setup_uni
###################################################
wanted <- types(toy_mm)=='temp'
m_uni <- xold(toy_mm[wanted,])
d_uni <- toy_d[wanted]
s_uni <- diag(B(toy_mhp)[,,"temp"])
x_uni <- xold(toy_point)[1,,drop=FALSE]
f_uni <- toy_LoF[["temp"]]
A_uni <- solve(corr.matrix(m_uni,scales=s_uni))   # NB: A^{-1}


###################################################
### code chunk number 12: use_uni
###################################################
interpolant.quick(
                  x      = x_uni,
                  d      = d_uni,
                  xold   = m_uni,
                  Ainv   = A_uni,
                  scales = s_uni,
                  func   = f_uni,
                  give.Z = TRUE)


###################################################
### code chunk number 13: use_uni_forjj
###################################################
jj_univariate <- 
interpolant.quick(
                  x      = x_uni,
                  d      = d_uni,
                  xold   = m_uni,
                  Ainv   = A_uni,
                  scales = s_uni,
                  func   = f_uni,
                  give.Z = TRUE)


###################################################
### code chunk number 14: mm_maker
###################################################
mm <- toy_mm_maker(81,82,83)
d <- obs_maker(mm, toy_mhp, toy_LoF, toy_beta)
jj_expt <- experiment(mm,d)


###################################################
### code chunk number 15: mhpopt
###################################################
mhp_opt <- optimal_params(jj_expt, toy_LoF, option="b")


###################################################
### code chunk number 16: show_estimated_M
###################################################
M(mhp_opt)


###################################################
### code chunk number 17: show_true_M
###################################################
M(toy_mhp)


###################################################
### code chunk number 18: estimate_toy_mm2
###################################################
est2 <- multem(toy_mm2, toy_expt, toy_mhp, toy_LoF)


###################################################
### code chunk number 19: holdbackhalf
###################################################
par(pty="s")

jj <- multem(toy_mm2, toy_expt, toy_mhp, toy_LoF, give=TRUE)
y <- jj$mstar
sd <- sqrt(diag(jj$cstar))
quan <- 1
cols <- c("red","green","blue")[types(toy_mm2)]

plot(toy_d2,y,xlim=c(-2,16),ylim=c(-2,16),asp=1,col=cols,pch=16, xlab="observed",ylab="emulated")

segments(x0=toy_d2, y0=y-sd*quan, y1=y+sd*quan,col=cols)

abline(0,1)
legend("topleft",c("temp","rain","humidity"),pch=16,col=c("red","green","blue"))


###################################################
### code chunk number 20: simple_function_thing
###################################################
fa <- function(x) sin(5*sum(x)) 
fb <- function(x) 7*sin(5*sum(x)) + sin(20*diff(x))


###################################################
### code chunk number 21: designmatrix
###################################################
# number of observation points:
na <- 33    # observation of 'a'
nb <- 09    # observation of 'b'


###################################################
### code chunk number 22: multivator.Rnw:1096-1097
###################################################
set.seed(0)


###################################################
### code chunk number 23: dowork
###################################################
xa <- latin.hypercube(na,2) # so rows of 'xa' are observation points for 'a'
xb <- xa[seq_len(nb),]
#xb <- latin.hypercube(nb,2)


###################################################
### code chunk number 24: namecols
###################################################
colnames(xa) <- colnames(xb) <- c("x","y")


###################################################
### code chunk number 25: feval
###################################################
a_obs <- apply(xa,1,fa)
b_obs <- apply(xb,1,fb) 


###################################################
### code chunk number 26: thingplot
###################################################
RS_mdm <- mdm(rbind(xa,xb),types=c(rep("a",na),rep("b",nb)))
RS_expt <- experiment(mm=RS_mdm, obs= c(a_obs,b_obs))
RS_opt <- optimal_params(RS_expt, option="b")


###################################################
### code chunk number 27: dfgfdgfdg
###################################################
n <- 20
xnew <- latin.hypercube(n,2,names=c("x","y"))
#xnew <- cbind(x=runif(20),y=runif(20))


###################################################
### code chunk number 28: simple_mdm
###################################################
RS_new_mdm <- mdm(rbind(xnew,xnew),rep(c("a","b"),each=n))
RS_prediction <- multem(x=RS_new_mdm, expt=RS_expt, hp=RS_opt)


###################################################
### code chunk number 29: uni_multi_fig
###################################################
nf <- layout(matrix(1:2,1,2))
par(mai=c(0,0,0,0)+0.8)
par(pty="s")

#first an emulator for the components separately:
sa <- optimal.scales(val=xa, scales.start=c(4,4), d=a_obs)
sb <- optimal.scales(val=xb, scales.start=c(4,4), d=b_obs)

# now use univariate emulation:
ya_emulated <- interpolant.quick(x=xnew, d=a_obs, xold=xa, pos.def.matrix=diag(sa))
yb_emulated <- interpolant.quick(x=xnew, d=b_obs, xold=xb, pos.def.matrix=diag(sb))

# now calculate f() at xnew:
ya_calc <- apply(xnew,1,fa)
yb_calc <- apply(xnew,1,fb)

# extract the multivator predictions:
ya_multivated <- RS_prediction[1:n]
yb_multivated <- RS_prediction[(n+1):(n+n)]

# now the scattergraphs:
#plot(ya_calc,ya_emulated,main="a emulated")
plot(yb_calc,yb_emulated,
     xlab='observed',ylab='predicted',
     xlim=c(-10,10),ylim=c(-10,10),
     main='(a), univariate emulation'
     )
abline(0,1)

plot(yb_calc,yb_multivated,
     xlab='observed',ylab='predicted',
     xlim=c(-10,10),ylim=c(-10,10),
     main='(b), multivariate emulation'
     )
abline(0,1)



###################################################
### code chunk number 30: simple_prediction
###################################################



###################################################
### code chunk number 31: datatemp
###################################################
data("mcneall")
dim(mcneall_temps)


###################################################
### code chunk number 32: showmap
###################################################
showmap(mcneall_temps[,1],FALSE)


###################################################
### code chunk number 33: datamcneall
###################################################
dim(mcneall_pc)
head(mcneall_pc,2)


###################################################
### code chunk number 34: showeigenmaps
###################################################

showmap_works_in_Sweave <- function(z, ...){ #bespoke function needed because showmap() output doesn't play nicely with layout().
  long <- seq(from=2.81,to=357,length.out=64)
  lat  <- c(-79.811531,seq(from=-74.81,to=86,len=30),86.6)
  z <- t(matrix(z,32,64))
  image(x=long,y=lat,z=z,col=terrain.colors(12), ...)
  contour(x=long,y=lat, landmask,level=0.5,drawlabels=FALSE,method='simple',add=TRUE,lwd=1.2,col='black')
}

nf <- layout(matrix(1:4,2,2))
par(mai=c(0,0,0,0)+0.8)
showmap_works_in_Sweave(eigenmaps[,1],main='Principal component 1')
showmap_works_in_Sweave(eigenmaps[,2],main='Principal component 2')
showmap_works_in_Sweave(eigenmaps[,3],main='Principal component 3')
showmap_works_in_Sweave(eigenmaps[,4],main='Principal component 4')


###################################################
### code chunk number 35: optimize_mcneall_hps (eval = FALSE)
###################################################
## jj <- apart(mcneall_pc, 17:20)
## opt_mcneall <- optimal_params(jj, start_hp=opt_mcneall, option='a')


###################################################
### code chunk number 36: show_mcneall_covs
###################################################
(CM <- M(opt_mcneall))


###################################################
### code chunk number 37: showcovs
###################################################
CM/sqrt(tcrossprod(diag(CM)))


