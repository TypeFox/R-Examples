biviso <- function(y, w=NULL,eps=NULL,eps2=1e-9,ncycle=50000,
                   fatal=TRUE,warn=TRUE) {
#
# Function 'biviso'.  To perform bivariate isotonic regression for simple
# (increasing) linear ordering on both variables.  Uses Applied Statistics
# Algorithm AS 206 (Isotonic regression in two independent variables;
# Gordon Bril, Richard Dykstra, Carolyn Pillers, and Tim Robertson;
# Algorithm AS 206; JRSSC (Applied Statistics), vol. 33, no. 3, pp.
# 352-357, 1984.)

# Check that ncycle makes sense:
if(ncycle!=round(ncycle) | ncycle < 2)
    stop("Argument ncycle must be an integer with value at least 2.\n")

# Check that y is of the right shape:
if(!is.numeric(y) | !is.matrix(y))
	stop("Argument \"y\" must be a numeric matrix.\n")

if(is.null(w)) w <- matrix(1,nrow=nrow(y),ncol=ncol(y)) else {
	if(!isTRUE(all.equal(dim(y),dim(w))))
		stop("Arguments \"y\" and \"w\" must have the same dimension.\n")
}

# Set epsilon:
if(is.null(eps)) eps <- sqrt(.Machine$double.eps)

nr   <- nrow(y)
nc   <- ncol(y)
nd   <- max(nr,nc)
rslt <- .Fortran(
	"smooth",
	nrow=as.integer(nr),
	ncol=as.integer(nc),
	ndim=as.integer(nd),
	x=as.double(y),
	w=as.double(w),
	a=double(4*nr*nc),
	b=double(5*nd),
	ncycle=as.integer(ncycle),
	icycle=integer(1),
	g=double(nr*nc),
	eps1=as.double(eps),
	eps2=as.double(eps2),
	ifault=integer(1),
        fx=double(nd),
        pw=double(nd),
        wi=double(nd),
        wt=double(nd),
        nw=integer(nd),
	PACKAGE="Iso"
)
if(rslt$ifault != 0) {
	if(rslt$ifault == 4 && warn) {
		warning(paste("A near zero weight less than delta=0.00001\n",
                              "was replaced by delta.\n",sep=""))
	} else if(fatal) {
		stop(paste("Failed with ifault = ",rslt$ifault,".\n",sep=""))
        } else if(warn) {
		warning(paste("Algorithm gave ifault = ",rslt$ifault,".\n",sep=""))
	}
}
m <- matrix(rslt$g,nrow=nr,ncol=nc)
attr(m,"icycle") <- rslt$icycle
attr(m,"ifault") <- rslt$ifault
m
}
