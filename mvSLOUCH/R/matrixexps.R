.calc.exptA<-function(t,lAcalcs=NULL,A=NULL){
    exptA<-NA
    if (is.null(lAcalcs)){## precalculate matrices if not done
	if(is.null(A)){print("No A matrix to work with");return(exptA);}
        else{
	    lAcalcs<-.decompEigenA.S(list(A=A,B=NULL),NULL,NA,list(bCalcA=TRUE,bCovCalc=FALSE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)[[1]]
	    exptA<-.calc.exptA(t,lAcalcs)
	}
    }else{    
	if (lAcalcs$decomp){
    	    P<-lAcalcs$eigA$vectors
	    vLambda<-lAcalcs$eigA$values
	    kY<-length(vLambda)
	    exptA<-Re(P%*%diag(exp(t*vLambda),ncol=kY,nrow=kY)%*%lAcalcs$invP) 
	}else{
	    if (lAcalcs$TwoByTwo){exptA<-.calc.exptA.2dim(t,lAcalcs$A)}
	    else{print("Cannot exponentiate this type of matrix right now")}
	}
    }
    exptA[which(abs(exptA)<1e-15)]<-0
    exptA
}

.calc.exptA.2dim<-function(t,A){
## The calculations for a 2x2 matrix comes from 
## Bernstein D., So W.; Some Explicit Formulas for the Matrix Exponential; IEEE Transactions on Automatic Control; Vol 38(8) 1993
    At<-t*A
    a<-At[1,1]
    b<-At[1,2]
    cM<-At[2,1]
    d<-At[2,2]
    D<-(a-d)^2+4*b*cM
    exptA<-NA
    if (D>0){
	delta<-sqrt(D)/2
	chd<-cosh(delta)
	shd<-sinh(delta)
	exptA<-rbind(c(chd+(a-d)*shd/(2*delta),b*shd/delta),c(cM*shd/delta,chd-(a-d)*shd/(2*delta)))
    }
    if (D==0){exptA<-rbind(c(1+(a-d)/2,b),c(cM,1-(a-d)/2))}
    if (D<0){
	delta<-sqrt(abs(D))/2
	sindelta<-sin(delta)
	cosdelta<-cos(delta)
	exptA<-rbind(c(cosdelta+(a-d)*sindelta/(2*delta),b*sindelta/delta),c(cM*sindelta/delta,cosdelta-(a-d)*sindelta/(2*delta)))
    }
    exp((a+d)/2)*exptA
}

.calc.integral.evAStevA<-function(t,S,lAcalcs=NULL,A=NULL){
    Integral<-NA    
    if (is.null(lAcalcs)){## precalculate matrices if not done
    	    if(is.null(A)){print("No A matrix to work with");return(Integral);}
	    else{
		lAcalcs<-.decompEigenA.S(list(A=A,B=NULL),NULL,NA,list(bCalcA=TRUE,bCovCalc=FALSE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)[[1]]
		Integral<-.calc.integral.evAStevA(t,S,lAcalcs,A)
	    }
    }else{
	if (lAcalcs$decomp){ ## our matrix is decomposable
	    vLambda<-lAcalcs$eigA$values
	    P<-lAcalcs$eigA$vectors
	    invP<-lAcalcs$invP
	    kY<-length(vLambda)
	    
	    V<-matrix(0:(kY^2-1),nrow=kY,ncol=kY,byrow=TRUE)
	    V<-apply(V,c(1,2),.CalcVlq2,"vlambda"=vLambda,"t"=t,"k"=kY)
	    
	    invPSP<-invP%*%S%*%t(invP)
	    invPSP[which(Mod(invPSP)<1e-15)]<-0
	    Integral<-Re(P%*%(V*(invPSP))%*%t(P))  ##V* Hadamard product
	    ## the integral has to be real no matter what the decomposition of A is we don't want complex R classes so need to take Real ie 0*0i
	}
	else{
	    if (lAcalcs$TwoByTwo){Integral<-.calc.integral.evAStevA.2dim(t,S,lAcalcs$A)}
	    else{print("Cannot calculate the integral for this class of matrices")}	
	}
    }
    Integral[which(abs(Integral)<1e-15)]<-0
    Integral
}

.calc.integral.evAStevA.2dim<-function(t,S,A){
## The calculations of expM for a 2x2 matrix comes from 
## Bernstein D., So W.; Some Explicit Formulas for the Matrix Exponential; IEEE Transactions on Automatic Control; Vol 38(8) 1993
## then we do a step by step calculation of the integral
## The symbolic calculations to get the derivative were done in Mathematica
    a<-A[1,1]
    b<-A[1,2]
    cM<-A[2,1]
    d<-A[2,2]
    s11<-S[1,1]
    s12<-S[1,2]
    s21<-S[2,1]
    s22<-S[2,2]
    D<-(a-d)^2+4*b*cM ## discriminant ok since sign won't change sgn(D)=sgn(v^2 D)=sgn(disc(vA))
    Int11<-NA
    Int12<-NA
    Int21<-NA
    Int22<-NA

    if (D==0){
	Int11<-(1/(4* (a + d)^3))*(-2*a^2*s11- 4 * a * d *s11-10*d^2*s11+8*b*d*(s12+s21)-8*b^2*s22+exp((a+d)*t)*(2*(a^2*s11+2*a*d*s11+5*d^2*s11-4 * b * d*(s12+s21)+4* b^2 * s22)+2*(a+d)*(a^2*s11+2*a*d*s11-3*d^2*s11+ 4 * b* d*(s12+s21)-4*b^2*s22)*t+(a+d)^2*((a-d)*(a*s11-d*s11+2*b*(s12+s21))+4*b^2*s22)*t^2))
	Int12<-(1/(4* (a + d)^3))*(-2* a^2* s12 - 2* d* (6* a + d)* s12 + 8* cM* (d* s11 - b* s21) + 8* a* b* s22 + exp((a + d)* t)* (2* (a^2 + 6* a *d + d^2)* s12 - 8* a* b* s22 + 2* (a + d)* ((a - d)^2* s12 + 4 *a* b *s22)* t - (a - d)* (a + d)^2* (a* s12 - d* s12 + 2* b* s22)* t^2 + cM* (-8* d* s11 + 8* b* s21 + 8 *(a + d)* (d* s11 - b* s21)* t + 2* (a + d)^2* (a* s11 - d *s11 + 2* b* s21)* t^2)))
	Int21<-(1/(4* (a + d)^3))*(-2* a^2* s21 - 2* d* (6* a + d)* s21 + 8* cM* (d* s11 - b* s12) + 8* a* b* s22 + exp((a + d)* t)* (2* (a^2 + 6 *a* d + d^2)* s21 - 8 *a* b *s22 + 2* (a + d)* ((a - d)^2* s21 + 4 *a *b *s22)* t - (a - d)* (a + d)^2* (a* s21 - d* s21 + 2 *b* s22)* t^2 + cM* (-8* d* s11 + 8* b* s12 + 8 *(a + d)* (d* s11 - b* s12)* t + 2* (a + d)^2* (a* s11 - d* s11 + 2* b* s12)* t^2)))
	Int22<-(1/(4* (a + d)^3))*(-8* cM^2 *s11 + 8* a* cM* (s12 + s21) - 2 *(5* a^2 + 2* a* d + d^2) *s22 +exp((a + d)* t)* (2 *(4*cM^2 * s11 -4* a* cM* (s12 + s21) + (5 *a^2 + 2* a* d + d^2) *s22) - 2 *(a + d) *(4* cM^2* s11 - 4* a* cM* (s12 + s21) + (a - d) *(3* a + d) *s22) *t + (a + d)^2* (4* cM^2 *s11 - 2* cM* (a - d)* (s12 + s21) + (a - d)^2* s22)* t^2))
    }
    if (D>0){
	Int11<- (1/D)* ((2* b* (b* cM - a *d)*exp((a + d)* t) *(-2* cM* s11 + (a - d)* (s12 + s21) + 2* b *s22) - D *(d* (a + d)* s11 - b* (cM* s11 + d* (s12 + s21)) + b^2* s22))/( 2* (a + d)* (-b *cM + a *d)) - (1/(2* b* cM - 2 *a* d))* exp((a + d)* t)* ((a^2* d* s11 + d^3* s11 - b* d* (-3* cM* s11 + d* (s12 + s21)) + a* (-2 *d^2 *s11 + b* (-cM *s11 + d* (s12 + s21)) + b^2* s22) + b^2 *(-2* cM* (s12 + s21) + d* s22)) *cosh(sqrt(D)*t) + sqrt(D)* ((a - d) *d* s11 + b* (-cM *s11 + d *(s12 + s21)) - b^2 *s22)* sinh(sqrt(D) *t)))
	Int12<-	(1/D)* ((D* (cM* d* s11 - b *cM *s21 - 2* a* d* s12 + b *cM* s12 + a* b* s22) + 2 *(-b* cM + a* d)*exp((a + d)*t )* (-a* cM* s11 + cM* d* s11 + a^2* s12 - 2 *b* cM* s21 - 2* a* d *s12 + d^2 *s12 + 2* b* cM* s12 + b* (a - d)* s22))/(2* (a + d)* (-b *cM + a* d)) + (1/(2* b* cM - 2* a* d))*exp((a + d)*t) *((cM* d* (-a + d)* s11 + 2 *b^2* cM *s22 + b* (2* cM^2 *s11 - cM *(a + d)* (s12 + s21) + a* (a - d)* s22))* cosh(sqrt(D)*t) + sqrt(D) *(cM *(-d* s11 + b* (s12 + s21)) - a* b* s22)* sinh(sqrt(D)* t)))
	Int21<- (1/D)* ((D* (cM* d *s11 - b *cM* s12 - 2 *a *d* s21 + b *cM *s21 + a* b* s22) + 2 *(-b* cM + a* d)*exp((a + d)* t)* (-a* cM *s11 + cM* d* s11 + a^2* s21 - 2* b* cM* s12 - 2* a *d *s21 + d^2* s21 + 2 *b* cM* s21 + b *(a - d)* s22))/(2* (a + d)* (-b* cM + a* d)) + (1/(2* b* cM - 2* a* d))*exp((a + d)*t) *((cM* d* (-a + d)* s11 + 2 *b^2* cM *s22 + b* (2* cM^2 *s11 - cM *(a + d)* (s12 + s21) + a* (a - d)* s22))* cosh(sqrt(D)*t) + sqrt(D) *(cM *(-d *s11 + b* (s12 + s21)) - a* b* s22)* sinh(sqrt(D)* t)))
	Int22<-	(1/D)* ((2* cM* (b* cM - a* d)* exp((a + d)* t)* (2* cM* s11 - a *(s12 + s21) +  d* (s12 + s21) - 2 *b* s22) - D*(cM^2 *s11 +  a* (a + d)* s22 - cM *(a* (s12 + s21) + b* s22)))/(  2* (a + d) *(-b *cM + a *d)) -  (1/(2* b* cM - 2 *a* d))* exp((a + d)* t) *((a^3 *s22 - a^2* (cM *(s12 + s21) + 2* d* s22) + cM *(cM* (d* s11 - 2* b* (s12 + s21)) - b* d *s22) + a *(cM^2* s11 + d^2* s22 + cM* (d *(s12 + s21) + 3* b* s22))) *cosh(sqrt(D) *t) +       sqrt(D) *(-cM^2* s11 + a* (-a + d)* s22 +cM *(a* (s12 + s21) - b* s22))* sinh(sqrt(D)* t)))
    }
    if (D<0){
	D<-(-1)*D
	Int11<-	(1/(2* D))* (((-1 + exp((a + d)* t)) *((a - d) *(a* s11 - d *s11 + 2* b *(s12 + s21)) + 4 *b^2* s22))/(a + d) + (1/((a + d)^3 + (a + d)* D))* ((a + d)^2 *((a - d) *(a *s11 - d* s11 + 2 *b* (s12 + s21)) + 4 *b^2* s22) + D* ((a + d) *(-4 *d* s11 + (a + d)* exp((a + d)* t)* s11 + 2* b* (s12 + s21)) + (-1 + exp((a + d) *t))* s11* D ) - (a + d)* exp((a + d)* t)* ((a^2 - d^2)* (a *s11 - d *s11 + 2* b* (s12 + s21)) + 4* b^2* (a + d)* s22 + (a* s11 - 3 *d* s11 + 2* b* (s12 + s21)) *D) *cos(t *sqrt(D)) + (a + d)* exp((a + d) *t) *sqrt(D)* (a^2 *s11 + 2* a* d *s11 - 3* d^2* s11 +4* b* d* (s12 + s21) - 4 *b^2 *s22 + s11* D)* sin(t *sqrt(D))))
	Int12<-	(1/(2* D))* (2* cM* s11 - a *s12 - d *s12 + 2* b* s22 + (4 *(a + d)* (-cM *d* s11 + a* d* s12 + b* cM *s21 - a* b* s22))/((a + d)^2 + D) + (2* cM* d* s11 - 4* b* cM* s21 + a^2* s12 + d^2* s12 - 2* b* d* s22 - 2* a* (cM* s11 + d* s12 - b* s22) - s12* D)/(a + d) + exp((a + d)* t)* ((-2* cM* d* s11 + 4 *b* cM* s21 - a^2* s12 - d^2 *s12 + 2 *b* d *s22 + 2* a* (cM* s11 + d* s12 - b* s22) + s12 *D)/(a + d) + (1/((a + d)^2 + D))* ((a + d)* (2 *cM* d* s11 - 4* b *cM *s21 + a^2 *s12 + d^2 *s12 - 2 *b* d* s22 - 2* a* (cM* s11 + d *s12 - b* s22)) + (-2* cM* s11 + (a + d)* s12 - 2 *b* s22)* D)* cos(t *sqrt(D)) + (sqrt(D)* (4* cM* (d* s11 - b *s21) + (a - d)^2 *s12 + 4* a* b* s22 + s12* D)* sin(t*sqrt(D)))/((a + d)^2 + D)))
	Int21<-	(1/(2* D))* (2* cM* s11 - a *s21 - d *s21 + 2* b* s22 + (4* (a + d)* (-cM *d* s11 + a *d* s21 + b* cM* s12 - a* b* s22))/((a + d)^2 + D) + (2* cM* d* s11 - 4* b* cM* s12 + a^2* s21 + d^2 *s21 - 2* b* d* s22 - 2* a* (cM* s11 + d* s21 - b* s22) - s21* D)/(a + d) + exp((a + d)* t)* ((-2* cM* d* s11 + 4* b* cM* s12 - a^2* s21 - d^2 *s21 + 2* b* d* s22 + 2* a* (cM* s11 + d* s21 - b* s22) + s21 *D)/(a + d) + (1/((a + d)^2 + D))* ((a + d)* (2 *cM* d* s11 - 4* b* cM* s12 + a^2* s21 + d^2* s21 - 2* b* d* s22 - 2* a* (cM* s11 + d *s21 - b* s22)) + (-2* cM* s11 + (a + d)* s21 - 2 *b* s22)* D)* cos(t *sqrt(D)) + (sqrt(D)* (4* cM* (d* s11 - b *s12) + (a - d)^2 *s21 + 4* a* b* s22 + s21* D)* sin(t*sqrt(D)))/((a + d)^2 + D)))
	Int22<-	(1/(2* D))* (-((4* cM^2* s11 - 2* cM* (a - d) *(s12 + s21) + (a - d)^2* s22 + s22* D)/(a + d)) + ((a + d) *(4 *cM^2* s11 - 2* cM* (a - d)* (s12 + s21) + (a - d)^2 *s22) + (2* cM* (s12 + s21) + (-3* a + d) *s22) *D)/((a + d)^2 + D) + exp((a + d) *t)* ((4* cM^2* s11 - 2* cM* (a - d)* (s12 + s21) + (a - d)^2* s22 + s22* D)/(a + d) - (((a + d)* (4 *cM^2 *s11 - 2 *cM* (a - d)* (s12 + s21) + (a - d)^2 *s22) + (2* cM* (s12 + s21) + (-3* a + d) *s22)* D)* cos(t *sqrt(D)))/((a + d)^2 + D) + ( sqrt(D)* (-4* cM^2* s11 + 4 *a* cM * (s12 + s21) + (-3* a^2 + 2* a* d + d^2) *s22 + s22 *D)* sin(t* sqrt(D)))/((a + d)^2 + D)))
    }
    rbind(c(Int11,Int12),c(Int21,Int22))    
}
