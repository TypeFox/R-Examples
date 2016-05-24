starship <- function(data,optim.method="Nelder-Mead",initgrid=NULL,
	inverse.eps=.Machine$double.eps,param="FMKL",optim.control=NULL,return.data=FALSE) 
{
# call adaptive grid first to find a first minimum
if (is.null(initgrid) ) {
  # default initgrid - depends on parameters
  if (param=="FMKL"| param=="FKML" | param=="fmkl"| param=="fkml" | param=="freimer" | param=="frm") {
    initgrid=list(
      lcvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,0.8,1,1.5), 
      ldvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,0.8,1,1.5))
    }
  if (param=="RS"| param=="rs" | param=="ramberg" | param=="ram") {
    initgrid = list(lcvect=c(0.1,0.2,4,0.8,1,1.5),
      ldvect=c(0.1,0.2,4,0.8,1,1.5))
    }
  if (param=="fm5") {
    initgrid=list(
      lcvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,0.8,1,1.5), 
      ldvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,0.8,1,1.5),
      levect=c(-0.5,-0.25,0,0.25,0.5))
    }
  if (param=="VSK"| param=="vsk" | param=="gpd" | param=="GPD") {
    initgrid=list(
      ldvect=c(-1.5,-.5,0,.2,.4,0.8,1.5,5), # lambda
      lcvect=c(0.3,0.5,0.7)) # delta - restricted to [0,1]
    }
  }
else 	{
  # initgrid taken from arguments - the parameter values in grid are not checked
	}
gridmin <- starship.adaptivegrid(data,initgrid,inverse.eps=inverse.eps,param=param)
# If they haven't sent any control parameters, scale by max(lambda1,1),
# lambda2 (can't be <= 0), don't scale for lambda3, lambda4 or lambda5
if (is.null(optim.control) ) {
	if (param=="fm5") {
		optim.control <- list(parscale=c(max(1,abs(gridmin$lambda[1])), abs(gridmin$lambda[2]),1,1,1))
		}
	else {	
		optim.control <- list(parscale=c(max(1,abs(gridmin$lambda[1])), abs(gridmin$lambda[2]), 1,1))
		}
	}
# else use what they sent - this should allow the user to change other stuff 
# in the control - this will mean that parscale is the default unless they
# change it also
# call optimiser
optimmin <- optim(par=gridmin$lambda,fn=starship.obj,method=optim.method,
	control=optim.control,data=data,param=param,inverse.eps=inverse.eps)
result <- list(lambda=optimmin$par,grid.results=gridmin,optim.results=optimmin,param=param)
if (return.data) {result$data = data}
# Apply starship class to result
class(result) <- "starship"
result$method.code <- "SM"
result$method.name <- "Starship"
# Add names to the lambda element - what names to use will depend on the parameterisation
if (param=="VSK"| param=="vsk" | param=="gpd" | param=="GPD") {
  names(result$lambda) <- c("alpha","beta","delta","lambda")
  } else { # any other parameterisation
  names(result$lambda) <- paste("lambda",1:length(result$lambda),sep="")
  }
result
}

starship.adaptivegrid <- function(data,initgrid,inverse.eps=1e-8,param="FMKL") 
{
data <- sort(data)
quarts <- quantile(data)
nombo <- length(data)
minresponse <- 1000
if (is.null(initgrid)){
  warning("use of starship.adaptivegrid without specifying the argument initgrid is now deprecated.  Continuing with a failsafe initgrid, but please change your code")
  initgrid = list(
    lcvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,0.8,1,1.5), 
    ldvect=c(0,.1,.2,.4,0.8,1),
    levect=c(-0.5,-0.25,0,0.25,0.5))
  }
if (param=="fm5"){
	minlambda <- c(NA,NA,NA,NA,NA)
	# There are other options than copying this code, but they seem to involve doing a param check inside all the for loops, and that seems terribly inefficient to me.  If it doesn't actuallly impact efficiency, let me know - RK
	for (le in initgrid$levect) {
	for (ld in initgrid$ldvect) {
		for (lc in initgrid$lcvect) {
			# calculate expected lambda2 on basis of IQR
			# IQR for 0,1,lc,ld - Changed for fm5
			iqr1 <- qgl(.75,c(0,1,lc,ld,le),param=param) -
qgl(.25,c(0,1,lc,ld,le),param=param)
			# actual IQR
			iqr <- quarts[4] - quarts[2]
			# so estimated lambda 2 from IQR
			lbguess <- iqr1/iqr
			for (lb in (c(0.5,0.7,1,1.5,2)*lbguess) )
				{
				# calculate expected lambda1 on the basis of median
				lavect <- quarts[3] -c(qgl(0.65,c(0,lb,lc,ld,le),param=param),
  qgl(0.55,c(0,lb,lc,ld,le),param=param),qgl(0.5,c(0,lb,lc,ld,le),param=param), 
  qgl(0.45,c(0,lb,lc,ld,le),param=param),qgl(0.35,c(0,lb,lc,ld,le),param=param) )
				for (la in lavect) {
					# calculate uniform g-o-f
					response <- starship.obj(c(la,lb,lc,ld,le),data,inverse.eps,param)
					if (response < minresponse) {
						minresponse <- response
						minlambda <- c(la,lb,lc,ld,le)
						} # new minimum
						# otherwise try the next
					} #lavect
				} #lbvect
			} # lcvct
		} # ldvect
		} # levect
	} # if param fm5 - not indendent in a consistent manner, but there's no room left in 80 columns
else	{
	minlambda <- c(NA,NA,NA,NA)
	for (ld in initgrid$ldvect) {
		for (lc in initgrid$lcvect) {
			# calculate expected lambda2 on basis of IQR
			# IQR for 0,1,lc,ld
			iqr1 <- qgl(.75,0,1,lc,ld,param=param) -
qgl(.25,0,1,lc,ld,param=param)
			# actual IQR
			iqr <- quarts[4] - quarts[2]
			# so estimated lambda 2 from IQR
			lbguess <- iqr1/iqr
			for (lb in (c(0.5,0.7,1,1.5,2)*lbguess) )
				{
				# calculate expected lambda1 on the basis of median
				lavect <- quarts[3] -c(qgl(0.65,0,lb,lc,ld,param=param),
	qgl(0.55,0,lb,lc,ld,param=param),qgl(0.5,0,lb,lc,ld,param=param), 
	qgl(0.45,0,lb,lc,ld,param=param),qgl(0.35,0,lb,lc,ld,param=param) )
				for (la in lavect) {
					# calculate uniform g-o-f
					response <- starship.obj(c(la,lb,lc,ld),data,inverse.eps,param)
					if (response < minresponse) {
						minresponse <- response
						minlambda <- c(la,lb,lc,ld)
						} # new minimum
						# otherwise try the next
					} #lavect
				} #lbvect
			} # lcvct
		} # ldvect
	} # if param not fm5
list(response=minresponse,lambda=minlambda)
}

starship.obj <- function(par,data,inverse.eps,param="fmkl") 
{
l1 <- par[1]; l2 <- par[2];
l3 <- par[3]; l4 <- par[4];
if (param=="fm5") { l5 <- par[5] }
# Check that these are legitimate parameter values.  If not, give a very
# large internal g-o-f measure.  We do this instead of NAs to make the
# optimistations easy to code.  Should investigate using NAs instead
if (!gl.check.lambda(par, param=param,vect=TRUE)) {
	return(54321)
	}
x <- sort(data)
# defining other variables
nombo <- length(x)

xacc <- inverse.eps
# values sent to pgl
# lower & upper limit on u
x1 <- xacc 
x2 <- 1.0 - xacc

u <- pgl(x,par,inverse.eps=inverse.eps,param=param);
# write.table(matrix(c(u),nrow=1),"interim-output/u1.txt",append=T,sep=",",quote=F,col.names=F)
response <- .anddarl(u,nombo);
# write.table(matrix(c(response),nrow=1),"interim-output/response.txt",append=T,sep=",",quote=F,col.names=F)
return(response)
}

# ANDERSON-DARLING TEST
.anddarl <- function(u,nombo)
{
# Anderson-Darling g-o-f assessment for uniform
ad <- c(dim(nombo),dim(1));
for(j in 1:nombo) {ad[j] <- ((2*j-1)/nombo)*(log(u[j])+log(1-u[nombo+1-j]))}
# This is useful in opt to illustrate what its doing
# hist(u)
asum <- -nombo-sum(ad);
asum;
}
# END OF ANDERSON-DARLING TEST
