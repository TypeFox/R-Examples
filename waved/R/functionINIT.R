#  
#
#   This file contains all the functions to perform
#   the WaveD transform in R.   The top fucntion is WaveD.
#   
#
#
#
# Written for R by Marc Raimondo 
# and Michael Stewart, the University of Sydney, 2006. 
# Copyright (c) 2006. 
#
#  references: ``Wavelet deconvolution in a periodic setting''
#  by  I.M. Johnstone, G. Kerkyacharian, D. Picard and M. Raimondo (2004).[JKPR]
#  The paper is available  by www from Marc Raimondo's web page
#		
#
#              ``Translation Invariant Deconvolution in a periodic setting''
#  by  D.Donoho and M. Raimondo (2005). IJWMIP
#  The paper is available  by www from Marc Raimondo's web page
#		
#
#      Software Manual: "The WaveD Transform in R"
#       by Marc Raimondo and Michael Stewart 
#      available from the Journal of Statistical Software (2007)  
#
#    R-code version by Marc Raimondo and Michael Stewart University of Sydney 2006.
#   
#
#
#
#
#
#########################################################################



WaveD<-function(yobs,g=c(1,rep(0,(length(yobs)-1))),MC=FALSE,SOFT=FALSE,F=find.j1(g,scale(yobs))[2],L=3,deg=3,eta=sqrt(6),thr=maxithresh(yobs,g,eta=eta),label="WaveD")
{
# 1-d Wavelet deconvolution. 
#  Inputs (required)
#    yobs =f*g+Noise   
#    g    Sample of the (known) function g
#   Inputs (optional)
#    L    Lowest resolution level (default=3)
#    F    Finest resolution level (default=Automatic)
#    deg   deg  of the Meyer Wavelet (default=3)
#    eta: smoothing parameter (default=conservative sqrt(6))
#  Outputs 
#    fhat= Estimated function 'f'
#    
#   
#  reference: [JKPR04]
    

################################
 cl <- match.call()

###Preliminary parameter estimation#################

nn<-length(yobs);
J<-log2(nn);
psyJ_fft<-wavelet_YM(J-2,J,deg);
yobs_fft<-fft(yobs);
yobs_w<- FWT_TI(yobs_fft,psyJ_fft);
noise=sqrt(nn)*yobs_w[1:(nn/2)];
s<-mad(noise);


######################################################
##Ordinary deconvolution##############################

g_fft<-fft(g);
x_fft<-yobs_fft/g_fft; 

#################################################Threshold set up

if(length(thr)<F-L && length(thr)>1){thr=c(thr,rep(thr[length(thr)],log2(nn)));
warning("Vector of threshold has not length J-L \n")}
if(length(thr)==1){thr=rep(thr[1],log2(length(yobs)))}


#############################################################################
#############if MC=TRUE WaveD returns only the TI-WaveD estimator############
if(MC){

 # loop to compute up coefficients up to F
 phyJ_fft<-scaling_YM(L,J,deg);
 f_Coarse <- IFWT_TI(x_fft,phyJ_fft,L,thr[1],nn,SOFT=SOFT);
 f_sum<-f_Coarse
#if F=J-1 use fine_YM instead of wavelet_YM
if (F<J-1) (Fmax<-F) else (Fmax<-J-2)
for (j in L:Fmax){
    psyJ_fft<-wavelet_YM(j,J,deg);
	f_detail <- IFWT_TI(x_fft,psyJ_fft,j,thr[j-L+2],nn,SOFT=SOFT);
    f_sum<-f_sum+f_detail;     
 }
if (Fmax==J-2) {
      j=J-1;
      psyJ_fft<-fine_YM(j,J,deg); 
      f_detail = IFWT_TI(x_fft,psyJ_fft,j,thr[j-L+2],nn,SOFT=SOFT);
     f_sum=f_sum+f_detail;  
      y=f_sum} else {y=f_sum}
}


############################################################################
############################################################################
#############if MC=FALSE WaveD returns a full WaveD class output############
else{

#################################################################
#Aux computation for waved outputs
M=find.j1(g,s)[1];
y_thr0_w=FWaveD(yobs,g=g,thr=0,F=F,SOFT=SOFT);
#Note here should be F+1 due to IwaveD which stops at to F-1
y0=IWaveD(y_thr0_w,L,deg,(F+1));
y_TI=WaveD(yobs,g,MC=TRUE,F=F,thr=thr);


 # Compute FWaveD and threshold


y_w=FWaveD(yobs,g=g,thr=thr,F=F,SOFT=SOFT); 

y_res=y_thr0_w-y_w;

#
# Invert the WaveD Transform
#

#Note here should be F+1 due to IwaveD goes to F-1

y=IWaveD(y_w,L,deg,(F+1));

#######################Summary list#########################

if(SOFT){pol='Soft'}else{pol='Hard'}

pv=shapiro.test(yobs_w*sqrt(nn))$p;
if(pv<0.01)
{
print("Warning, WaveD-fit residuals normality-test  P-value=")
print(pv)
print("larger threshold eta-constant may be required")
}


if(length(thr)==1){thr=rep(thr[1],F-L+2)}

output=list(sd=s,noise=noise,j1=find.j1(g,scale(yobs))[2],Coarse=L,degree=deg,Low=L,M=M,thr=thr[1:(F-L+2)],w=y_thr0_w,w.thr=y_w,FWaveD=y_thr0_w,iw=y0,ordinary=y,waved=y_TI,POLICY=pol,residuals=y_res[1:(2^(F+1))],percent=threshsum(y_res,L,F),levels=c(L,L:F),eta=eta,F=F,n=nn,data=yobs,g=g,p=pv,label="WaveD",call=cl)
class(output)="wvd"
output
}



}



#####################################################
#####################################################
#####################################################
##########functions for ordinary WaveD below#########

##################################################
##################################################
PhaseC<-function(l,j)
#
#Phase matrix to compute wavelet coefficients
# in the Fourier Domain
#

{
nb=2^j-1

k=0:nb
k=matrix(k,nrow=1)
l=matrix(l,nrow=1)
mat=t(l)%*%k

phase_mat=exp(2*pi*1i*mat/(2^j));

}

# Written by  Marc Raimondo (the University of Sydney), 2006. 
# Comments? e-mail marcr@maths.usyd.edu.au



########################################
########################################

WaveDjC<-function(y_fft,f2fft,j)
# compute--scaling--coef--at-level-j
# uses Plancherel equality.
#  Usage
#    coef=coefj(f1fft,phij0fft,j,n)
#  Inputs  
#    y_fft=fft(yobs)
#    f2fft=fft(phij0) or fft(psij0) accordingly
#    j=resolution level
#    n=sample size
#  Outputs
#    vector of coefs at resolution level j: length= 2^j

{
n=length(y_fft);
j1=j-1;
j2=j1-1;
x2=2^j1+2^j2+4;

l1=(0:x2);
l2=((n-x2):(n-1));
l_0=c(l1,l2);
l1=l1+1;
l2=l2+1;
l_1=c(l1,l2);

#initial vector without phase
# corresponds to k=0.


y1=y_fft[l_1];
v1=matrix(y1,nrow=1);
v2=f2fft[l_1];
v2=matrix(v2,nrow=1);


#componentwise multiplication
vector1=v1*Conj(v2);

#matrix of phases

phasemat=PhaseC(l_0,j);



#coef is the product
# of vector1 by phase matrix:




y=Re(vector1%*%phasemat)/(n^2);

}


##############################
##############################

WaveDjD<-function(y_fft,f2fft,j)
# compute--wavelet--coef--at-level-j        #changed
# uses Plancherel equality.
#  Usage
#    coef=coefj(f1fft,phij0fft,j,n)
#  Inputs  
#    y_fft=fft(yobs)
#    f2fft=fft(phij0) or fft(psij0) accordingly
#    j=resolution level
#    n=sample size
#  Outputs
#    vector of coefs at resolution level j: length= 2^j

{
n=length(y_fft);
j1=j-1;
j2=j1-1;
j3=j-3;

x1=2^j1-2^j3-1;
x2=2^j+2^j2+8;
l1=(x1:x2);
l2=((n-x2):(n-x1+2));

l_1=c(l1,l2);
l_0=l_1-1;

#initial vector without phase
# corresponds to k=0.


y1=y_fft[l_1];
v1=matrix(y1,nrow=1);
v2=f2fft[l_1];
v2=matrix(v2,nrow=1);


#componentwise multiplication
vector1=v1*Conj(v2);

#matrix of phases

phasemat=PhaseC(l_0,j);



#coef is the product
# of vector1 by phase matrix:




y=Re(vector1%*%phasemat)/(n^2);

}


##############################
##############################

WaveDjF<-function(f1fft,f2fft,j)
# compute--fine wavelet--coef--at-level-j    #changed
# uses Plancherel equality.
#  Usage
#    coef=coefj(f1fft,phij0fft,j,n)
#  Inputs  
#    f1fft=fft(yobs)
#    f2fft=fft(phij0) or fft(psij0) accordingly
#    j=resolution level
#    n=sample size
#  Outputs
#    vector of coefs at resolution level j: length= 2^j

{
n=length(f1fft);
j1=j-1;
j2=j-2;
j3=j-3;



j1=j-1;
j2=j-2;
j3=j-3;
x1=2^j1-2^j3;
x2=2^j+2^j1+2^j3+14;

l_1=(x1:x2);
l_0=l_1-1;



#initial vector without phase
# corresponds to k=0.


y1=f1fft[l_1];
v1=matrix(y1,nrow=1);
v2=f2fft[l_1];
v2=matrix(v2,nrow=1);


#componentwise multiplication
vector1=v1*Conj(v2);

#matrix of phases

phasemat=PhaseC(l_0,j);



#coef is the product
# of vector1 by phase matrix:




y=Re(vector1%*%phasemat)/(n^2);

}


########################################
########################################
FWaveD<-function(y,g=1,L=3,deg=3,F=(log2(length(y))-1),thr=rep(0,log2(length(y))),SOFT=FALSE)
#  FWaveD: Forward WaveD Transform
#  Usage
#    w = FWaveD(y,g,L,deg,F,SOFT)
#  Inputs
#    y    1-d signal; length(y) = 2^J
#    g    sample of known function; length(g)= 2^J
#    L    Coarsest Level of V_0;  L << J
#    deg  degree of polynomial window 2 <= deg <=4
#    F     Finest resolution level
#    SOFT: thresholding policy (default=HARD)
#  Outputs
#    w    1-d estimated wavelet transform of f from x=f*g+noise.
#  Note: if g has no input, returns the Meyer transform
#  Description
#    This algorithm is based on [JKPR] . 

{ x = y;
        nn = length(y);
	J = log2(nn);
     
if(length(thr)<F-L && length(thr)>1 ){thr=c(thr,rep(thr[length(thr)],log2(nn)));
warning("Vector of threshold has not length J-L \n")}
if(length(thr)==1){thr=rep(thr[1],log2(length(y)))}


#############################################################
y_fft = fft(y);

if(length(g)==1){g_fft=rep(1,nn)}else{
g_fft= fft(g);}
##
w=rep(0,nn);
# Perform deconvolution in the Fourier domain   
 
 y_fft=y_fft/g_fft;





#####################################


if(SOFT){
#
#  Compute Coefficients at Coarse Level.
#

waveL0_fft=scaling_YM(L,J,deg);


	w[1:(2^L)] = HardThresh(WaveDjC(y_fft,waveL0_fft,L),thr[1]);
       
  
#
#  Loop to Get Detail Coefficients for levels  j=L,...,F.
#
	for (j in L:F){ 
      waveJ0_fft=wavelet_YM(j,J,deg);
	  w[(2^j+1):(2^(j+1))] = SoftThresh(WaveDjD(y_fft,waveJ0_fft,j),thr[j-L+2]);
  }
	#
#  Calculate Fine Level Detail Coefficients (for j=J-1).
#  if F=J-1
  
  if(F==J-1){
  waveJ_fft=fine_YM(J-1,J,deg);
  co= WaveDjF(y_fft,waveJ_fft,J-1);
	w[(2^(J-1)+1):2^J] = SoftThresh(WaveDjF(y_fft,waveJ_fft,J-1),thr[J-L+1]);
  
  y=w;}else{y=w}
  }else
{
#
#  Compute Coefficients at Coarse Level.
#

waveL0_fft=scaling_YM(L,J,deg);


	w[1:(2^L)] = HardThresh(WaveDjC(y_fft,waveL0_fft,L),thr[1]);
       
  
#
#  Loop to Get Detail Coefficients for levels  j=L,...,F.
#
	for (j in L:F){ 
      waveJ0_fft=wavelet_YM(j,J,deg);
	  w[(2^j+1):(2^(j+1))] = HardThresh(WaveDjD(y_fft,waveJ0_fft,j),thr[j-L+2]);
  }
	#
#  Calculate Fine Level Detail Coefficients (for j=J-1).
# if necessary
if(F==J-1){
  waveJ_fft=fine_YM(J-1,J,deg);
  co= WaveDjF(y_fft,waveJ_fft,J-1); 
	w[(2^(J-1)+1):2^J] = HardThresh(WaveDjF(y_fft,waveJ_fft,J-1),thr[J-L+1]);
  
  y=w}else{y=w};
  }
	
    

}


##################
###################

projVj<-function(beta,n,deg)
# projection onto Vj
#  Usage
#    y = projVj(beta,n)
# 
#  Inputs
#    beta=wavelet coef at resolution level j
#    n=sample size
#    deg=degree of Meyer wavelet
	
{
nj = length(beta);
j = log2(nj);
	 
      
#computing fine sampling of Psy(j,0)_fft 

     psyJ_fft=scaling_YM(j,log2(n),deg);
     
     #dyadic grid at resolution level j
     dyad_grid=(1:(2^j))/2^j;
     
     #index of non-zero wavelet coef
     dyad_grid=floor(n*dyad_grid);
     dyad_grid=dyad_grid-dyad_grid[1]+1;
     #zero padding on fine sampling
    x_w=rep(0,n);
    #imbedding non-zero wavelet coef on fine grid 
    x_w[dyad_grid]=beta;
     
     
   #computing projection at level j
   #in the Fourier domain (convolution)
   
   x_w_fft=fft(x_w);
   
   y_fft=x_w_fft*psyJ_fft;
   
   y=Re(fft(y_fft,inverse=TRUE));
   
   
     
 }   
 


####################################
####################################

projWj<-function(beta,n,deg)
# projection onto Wj
#  Usage
#    y = projWj(beta,n)
# 
#  Inputs
#    beta=wavelet coef at resolution level j
#    n=sample size
#    deg=degree of Meyer wavelet
	
{
nj = length(beta);
j = log2(nj);
	 
      
#computing fine sampling of Psy(j,0)_fft 

     psyJ_fft=wavelet_YM(j,log2(n),deg);
     
     #dyadic grid at resolution level j
     dyad_grid=(1:(2^j))/2^j;
     
     #index of non-zero wavelet coef
     dyad_grid=floor(n*dyad_grid);
     dyad_grid=dyad_grid-dyad_grid[1]+1;
     #zero padding on fine sampling
    x_w=rep(0,n);
    #imbedding non-zero wavelet coef on fine grid 
    x_w[dyad_grid]=beta;
     
     
   #computing projection at level j
   #in the Fourier domain (convolution)
   
   x_w_fft=fft(x_w);
   
   y_fft=x_w_fft*psyJ_fft;
   
   y=Re(fft(y_fft,inverse=TRUE));
   
   
     
 }   
 

#############################
#############################

projFj<-function(beta,n,deg)
# projection onto Fj
#  Usage
#    y = projFj(beta,n)
# 
#  Inputs
#    beta=wavelet coef at resolution level j
#    n=sample size
#    deg=degree of Meyer wavelet
	
{
nj = length(beta);
j = log2(nj);
	 
      
#computing fine sampling of Psy(j,0)_fft 

     psyJ_fft=fine_YM(j,log2(n),deg);
     
     #dyadic grid at resolution level j
     dyad_grid=(1:(2^j))/2^j;
     
     #index of non-zero wavelet coef
     dyad_grid=floor(n*dyad_grid);
     dyad_grid=dyad_grid-dyad_grid[1]+1;
     #zero padding on fine sampling
    x_w=rep(0,n);
    #imbedding non-zero wavelet coef on fine grid 
    x_w[dyad_grid]=beta;
     
     
   #computing projection at level j
   #in the Fourier domain (convolution)
   
   x_w_fft=fft(x_w);
   
   y_fft=x_w_fft*psyJ_fft;
   
   y=Re(fft(y_fft,inverse=TRUE));
   
   
     
 }   
 


################################################
#################################################
IWaveD<-function(w,C=3,deg=3,F=log2(length(w)))
# IWaveD -- Inverse WaveD Transform (periodized Meyer Wavelet)
#  
#  Inputs
#    w   1-d wavelet transform, length(wc) = 2^J.
#    L    Coarsest Level of V_0;  L << J
#    deg  degree of polynomial window 2 <= deg <=4
#  Outputs
#    x    1-d reconstructed signal; length(x) = 2^J
#
#  Description
#    The IWaveD transform is obtained by the command
#        f = IFWaveD_PL(w,C,deg,F)
#    to reconstruct x, use the IWaveD_YM.
#

{ nn = length(w);
        wc = w;
    
if (F==log2(nn)){ J=F} else{J=F+1};
 
#
#  Reconstruct Projection at Coarse Level.
#
f_sum=0;
	beta = w[1:(2^C)];
  
	f_sum =projVj(beta,nn,deg);
	
 
#
#  Loop to Get Projections at detail levels j=C,...,J-2.
#
	for (j in C:(J-2)){
		alpha = w[dyad(j)];
    
		detail_proj =projWj(alpha,nn,deg);
		f_sum=f_sum+detail_proj;
	}
 
#
#  Calculate Projection for fine detail level, j=J-1.
#

if (F==log2(nn)) {
	alpha = w[dyad(J-1)];
   
	fine_detail = projFj(alpha,nn,deg);
	f_sum=f_sum+fine_detail}



	y = f_sum/nn

}

#########################################################
#########################################################
#########################################################
##########functions for TI-WaveD below###################

#
#
# This file contains all the fucntions to perform TIWaveD
#  (Tranlation Invariant Wavelet Deconvolution)
#
#

#######################################################################

findONE<-function(x)
#
#find location of >=1
#

{
 n1<-sum(x>=1)
(sort.list(-(x>=1)*1))[1:n1]
}

########################################################################

findZERO<-function(x)
#
#find location of <=0
#

{
 n1<-sum(x<=0)
(sort.list(-(x<=0)*1))[1:n1]
}
#########################################################################


MeyerWindow<-function(xi,deg)
{
#
# WindowMeyer -- auxiliary window function for Meyer wavelets.
#  
#   nu = WindowMeyer(xi,deg)
#  Inputs
#    xi     abscissa values for window evaluation
#    deg    degree of the polynomial defining Nu on [0,1]
#          1 <= deg <= 3
#  Outputs
#    nu     polynomial of degree 'deg' if x in [0,1]
#           1 if x > 1 and 0 if x < 0.
#
#
#




if (deg==3){
				 nu = xi^4 * ( 35 - 84 * xi + 70 * xi^2 - 20 * xi^3)
		} 
else if (deg==2){
	nu = xi^3*(10 -15* xi + 6*xi^2)} 

else if (deg == 1){
		nu = xi^2* (3 - 2*xi)}

else if (deg == 0){
		nu = xi}

ix0<-findZERO(xi)
if (length(ix0)>0) {nu[ix0]<-rep(0,length(ix0))}

ix1<-findONE(xi)
if (length(ix1)>0) {nu[ix1]<-rep(1,length(ix1))}

y=nu
}




########################################################################################

phyHAT<-function(x,deg)
{
# auxillary function for computing Meyer scaling function
# in the Fourier domain

aux1<-MeyerWindow((3*abs(x)-1),deg)


cos_part<-cos(pi*aux1/2);

y1<-cos_part*(abs(x)>1/3 & abs(x)<=2/3);
y2<-rep(1,length(x))*(abs(x)<=1/3);
y<-y1+y2;


  

}


################################################################################


psyHAT<-function(x,deg)
{
# auxillary function for computing Meyer wavelet function
# in the Fourier domain

aux1<-MeyerWindow((3*abs(x)-1),deg)
aux2<-MeyerWindow((3*abs(x)/2-1),deg)



sin_part<-(exp(-1i*pi*x))*sin(pi*aux1/2);
cos_part<-(exp(-1i*pi*x))*cos(pi*aux2/2);

y1<-sin_part*(abs(x)>1/3 & abs(x)<=2/3);
y2<-cos_part*(abs(x)>2/3 & abs(x)<=4/3);

y=y1+y2;

}

###########################################################################3
scaling_YM<-function(j,j_max,deg)
{
#generate phyHAT(j,0)--the fft of the scaling function--in the periodic setting.
# deg=deg of the Meyer wavelet (1,2,3).

 
# Written for R by Marc Raimondo 
# and Michael Stewart, the University of Sydney, 2006. 
# Copyright (c) 2006. 
#
#  reference: ``Translation Invariant Deconvolution in a periodic setting''
#  by  D.Donoho and M. Raimondo (2005). IJWMIP
#  The paper is available  by www from 
#		http://www.usyd.edu.au/u/marcr/
#


omega_pos<-seq(0,10/(2*pi),by=1/(2^j))
n_pos<-length(omega_pos);

nn<-2^j_max;
aux1=phyHAT(omega_pos,deg);

#rotate and take zero out

psyHAT_LEFT=rep(0,nn);

psyHAT_LEFT[1:n_pos]<-aux1;
aux2<-psyHAT_LEFT;



aux2<-aux2[-1]

aux2<-c(aux2, 0);
psyHAT_RIGHT<-Conj(rot90(aux2));



y<-(psyHAT_RIGHT+psyHAT_LEFT)*2^(-j/2)*nn;



 }   

###########################################################################

wavelet_YM<-function(j,j_max,deg)
{
#generate psyHAT(j,0)--the fft of the wavelet function--in the periodic setting.
# deg=deg of the Meyer wavelet (1,2,3).




omega_pos<-seq(0,10/(2*pi),by=1/(2^j))
n_pos<-length(omega_pos);

nn<-2^j_max;

omega_neg<-seq(-10/(2*pi),-1/2^j,,by=1/(2^j))

n_omega<-2*length(omega_pos)-1;

n<-2^j_max-n_omega;

aux1<-psyHAT(omega_pos,deg);

#rotate and take zero out




psyHAT_LEFT=rep(0,nn);

psyHAT_LEFT[1:n_pos]<-aux1;

aux2<-psyHAT_LEFT;



aux2<-aux2[-1]

aux2<-c(aux2, 0);
psyHAT_RIGHT<-Conj(rot90(aux2));



y<-(psyHAT_RIGHT+psyHAT_LEFT)*2^(-j/2)*nn;



 }   

##############################################################################

fine_YM<-function(j,j_max,deg)
{
#generate psyHAT--the fft of the wavelet function--in the periodic setting.
# deg=deg of the Meyer wavelet (1,2,3).
# at the finest resolution level
 



omega_pos<-seq(0,10/(2*pi),by=1/(2^j));

n_pos<-2^j-2^(j-3);

nn<-2^j_max;


aux1<-psyHAT(omega_pos,deg);

#rotate and take zero out



n2<-n_pos+2^(j-3)+1;
psyHAT_LEFT=rep(0,n2);

psyHAT_LEFT[1:n_pos]<-aux1[1:n_pos];
psyHAT_LEFT[(n_pos+1):(n_pos+2^(j-3)+1)]<--1;
aux2<-psyHAT_LEFT;
aux2<-aux2[-1]



psyHAT_RIGHT<-Conj(rot90(psyHAT_LEFT));

psyHAT_RIGHT<-psyHAT_RIGHT[-1];

psyHAT_RIGHT<-psyHAT_RIGHT[-length(psyHAT_RIGHT)]

y<-c(psyHAT_LEFT, psyHAT_RIGHT)*2^(-j/2)*nn;



 }  

########################################################################

rot90<-function(x)
{
#
#is the 90 degree counterclockwise rotation of matrix x.
#


if (is.null(dim(x))){
n1=length(x)
x[n1:1]
}
else {
n1=dim(x)[1]
t(x)[n1:1,]
}

}


#################################################################

HardThresh<-function(y,t)
{
# HardThresh -- Apply Hard Threshold 
#  Usage 
#    x = HardThresh(y,t)
#  Inputs 
#    y     Noisy Data 
#    t     Threshold
#  Outputs 
#    x     y 1_{|y|>t}
#
	x   = y * (abs(y) > t);



  
# Written for R by Marc Raimondo 
# and Michael Stewart, the University of Sydney, 2006. 
#  
#
# Reference: David Donoho and Iain Johnstone JASA 1994.
# 
#
}

#################################################################

SoftThresh<-function(y,t)
{
# SoftThresh -- Apply Hard Threshold 
#  Usage 
#    x = HardThresh(y,t)
#  Inputs 
#    y     Noisy Data 
#    t     Threshold
#  Outputs 
#    x     y  sign(y)(|y|-t)_+
#

        aux=(abs(y)-t);
aux=(aux+abs(aux))/2;
 
	x   = sign(y) * aux;



  
# Written for R by Marc Raimondo 
# and Michael Stewart, the University of Sydney, 2006. 
#  
#
# Reference: David Donoho and Iain Johnstone JASA 1994.
# 
#
}















##################################################################################
MultiThresh1<-function(s,g,L,eta)
{
# MultiThresh --Find optimal threshold by level-by-level. 
#
#  Inputs 
#    s    Noise-sd (or estimate). 
#    g_fft  noisy  Sample of g_fft 
#    L    Lowest resolution level=C (Coarsest)
#    eta  smoothing parameter (default=1)
#  Outputs 
#    y= vector of thresholds from level=L to J-1.
#
#  reference: ``Wavelet deconvolution in a periodic setting''
#  by  I.M. Johnstone, G. Kerkyacharian, D. Picard and M. Raimondo (2004).[JKPR]
#  The paper is available (as a prepint) by www from 
#		http://rome:www.usyd.edu.au/u/marcr/
   	  

#compute optimal threshold for deconvolution
# based on the Fourier transform of g.
# For the last resolution level we use spline extrapolation.

n=length(g);
g_fft<-fft(g);
F<-log(n)/log(2);

thr<-rep(0,F-L-1);
  for (j in L:(F-2)){
thr[j-L+1] <- (sum(abs(g_fft[dyad(j)])^(-2))+sum(abs(g_fft[n-dyad(j)+1])^(-2)))/2^j ;
  }

aux1<-s*sqrt(thr)*sqrt(log(n)/n)*(eta/sqrt(4*pi));


#linear extrapolation for Jmax
n1<-length(aux1);

d1<-(aux1[n1]-aux1[1])/(n1-1);
aux2<-(n1+1)*d1;
y<-c(aux1, aux2);



}    


#################################################################
################################################################
#
#
#
#
#

FWT_TI<-function(f_fft,psyJ_fft)
{
#
#Inputs:  f_fft=fft(f)
#         psyJ_fft=fft(Ondelette au niveau 'lev')
#         lev=resolution level 
#         thr=threshold 
#          nn=length(f)
#
# Output: wavelet coefficients at resolution level lev (non-ordered in time)
#
 
nn<-length(f_fft);
  Aj_fft<-f_fft*psyJ_fft;
  
#need to scale *sqrt(pi) as in Matlab

  Aj<-Re(fft(Aj_fft,inverse=TRUE))/(nn*nn);
 
  
    
   
#  Written  by Marc Raimondo, the University of Sydney, 
#  Comments? e-mail marcr@maths.usyd.edu.au
 
 }   
     
########################################################33
###########################################################
scale<-function(yobs,L=3,deg=3)
{
# scale-estimation-for-Meyer-Wavelet-coefficients
# Can be used for -direct-and--indirect-noisy--observations
#  Inputs  
#    yobs    1-d signal; length(x) = 2^J
#    L    Coarsest Level of V_0;  L << J
#    deg  degree of polynomial window 2 <= deg <=4
# 
#  Outputs
#    y=estimated noise standard deviation. 


n<-length(yobs);
J<-round(log(n)/log(2),1);



j<-J-1;
 psyJ_fft<-fine_YM(j,J,deg);
# psyJ_fft<-wavelet_YM(j-1,J,deg);
x_fft<-fft(yobs);

yobs_w<- FWT_TI(x_fft,psyJ_fft);




y<-mad(yobs_w*sqrt(n));
 

#y<-sqrt(var(yobs_w*sqrt(n)));

#  Written by Marc Raimondo, the University of Sydney, 2003.     
#  Comments? e-mail marcr@maths.usyd.edu.au
 
}


##############################################################
#

IFWT_TI<-function(f_fft,psyJ_fft,lev,thr,nn,SOFT=FALSE)
{
#
#Inputs:  f_fft=fft(f)
#         psyJ_fft=fft(Ondelette au niveau 'lev')
#         lev=resolution level 
#         thr=threshold 
#          nn=length(f)

 if (SOFT){
  Aj_fft<-f_fft*psyJ_fft;
  
  Aj<-Re(fft(Aj_fft,inverse<-TRUE))/(nn^2);
  Bj<-SoftThresh(Aj,thr);
  Bj_fft<-fft(Bj);
  fj_fft<-Bj_fft*Conj(psyJ_fft);
  J <- log2(nn);
  y<-Re(fft(fj_fft,inverse=TRUE))/(nn*2^(J-lev));
  }else{
  Aj_fft<-f_fft*psyJ_fft;
  
  Aj<-Re(fft(Aj_fft,inverse<-TRUE))/(nn^2);
  Bj<-HardThresh(Aj,thr);
  Bj_fft<-fft(Bj);
  fj_fft<-Bj_fft*Conj(psyJ_fft);
  J <- log2(nn);
  y<-Re(fft(fj_fft,inverse=TRUE))/(nn*2^(J-lev));
  }

    
   

 
 }  

##################################################
#
# 
   dyad<-function(j)
{
# dyad -- Index entire j-th dyad of 1-d wavelet xform
#  Usage
#    ix = dyad(j);
#  Inputs
#    j     integer
#  Outputs
#    ix    list of all indices of wavelet coeffts at j-th level

    i <- (2^(j)+1):(2^(j+1)) ;

}
# Original version by David L. Donoho
# Copyright (c) 1993. 
# Written for R by Marc raimondo (University of Sydney)    
    

#############################################################################################
#
#
#

fftshift<-function(x)
{
#
# for visualizing the Fourier transform with  the zero-frequency component in the middle of the spectrum.
#

m=length(x);
 p = ceiling(m/2);
    new_index = c(((p+1):m),(1:p));
y=x[new_index];

}

###################################################################


##########################################################################

  speczoom<-function(y_test,fenetre)
{
#Plot of spectrum with zero frequency in the middle
# of the window
#input is y_test (kconvolution kernel)
# fenetre is (window size)

n=length(y_test);
y_fft=fft(y_test);
a=n/2-fenetre/2+1;
b=n/2+fenetre/2+1;

spectrum=fftshift(y_fft);
spectrum=spectrum[a:b];
x=(a:b)-n/2-1; 
#subplot(211)
#plot((1:n)/n,y_test,'-.b');
#subplot(212);
#plot(x,log(abs(dum)),'-.r');

plot(x,log(abs(spectrum)),type='l',xlim=c(-500,500),ylim=c(-6,0));

y=log(abs(spectrum));

#min(abs(dum))



}
#
###########################################################
###########################################################

fftshift<-function(y)
{
#
# Shift zero-frequency component to center of spectrum.
#

n=length(y);
n2=floor(n/2);
ind1=(1:n2);
ind2=((n2+1):n);
aux1=(y[ind1]);
aux2=y[ind2];
y=c(aux2,aux1);






}
#
# Written for R by Marc raimondo (University of Sydney)    
    
########################################################
########################################################
stoptime<-function(g,sigma)
{
# compute the stoping time
# in the Fourier domain using noisy egein values
#Input: 
# g=noisy fft of convolution kernel
# sigma=estimated noise's sd
#Output
# [M,J_1]=stime(g,sigma)
# where M is frequency cut-off
# and   J_1 is corresponding resolution level cut-off
n=length(g);
#bandwidth
h=n/2;

aux=log(abs(g[1:h]));
threshold_fft=log((sigma/sqrt(n)))+ log((log(sqrt(n)/sigma))^(1/2));


aux1=1*(aux[1:h]- threshold_fft)>0;

vs1=sort(aux1);
ind1=sort.list(aux1);
freq_cut=ind1[1];
resolution_cut=floor(log2(freq_cut))-1;

if (min(aux1)>0){
    y=c(n/2,(log2(n)-1))}
else {
y=c(freq_cut,resolution_cut);
}
}
    
#
############################################################
#############################################################

find.j1<-function(g,sigma)
# compute the stoping time
# in the Fourier domain using noisy eigen values
# --scaled to band-limited wavelet with 2^{j/2}=sqrt{l} 
# for more details see [CR05]
#Input: 
# g =convolution kernel
# sigma=estimated noise's sd
#Output
# [M,J_1]=stime(g,sigma)
# where M is frequency cut-off
# and   J_1 is corresponding resolution level cut-off
######################################################
#
#
{
g_fft=fft(g);
n=length(g_fft);
aux2=g_fft[(1:(n/2))];
aux3=sqrt((1:(n/2)));
aux4=aux2/aux3;
aux5=rot90(aux4);
aux6=c(aux4,aux5);
y=stoptime(aux6,sigma);
}
 


######################################################
######################################################

multires <- function(wcUntrimmed,lowest=3,coarse=3,highestplot=NULL,descending=FALSE,sc=1)
#this function plots MRA
{
  L <- log2(length(wcUntrimmed))
  if((L%%1)>0) 
    warning("Vector of wavelet coefficients has length not a power of 2\n")
  wc <- wcUntrimmed[-(1:(2^coarse))] 
  #print(rbind(length(wcUntrimmed),length(wc)))
  highest <- ceiling(L)-1
  if(is.null(highestplot))  highestplot <- highest
  reslev <- rep(coarse:highest,2^(coarse:highest)) 
  #print(cbind(wc,reslev))
  num <- numeric()
  for(i in (min(reslev):max(reslev))){num <- c(num,2*(1:(2^i))-1)}
  #cbind(num,den)
  M <- 2*max(abs(wc))/sc;
 den <- 2^(reslev+1)
  ind <- (reslev>=lowest)&(reslev<=highestplot)
  x <- (num/den)[ind]
  y1 <- reslev[ind]
  if(descending) y1 <- (highest+lowest-y1)
  leny1=length(y1)
  lenwcind <- length(wc[ind])
  #print(cbind(leny1,lenwcind))
  y2 <- y1 + wc[ind]/M
  
  cbind(x,y1,y2) 
  plot(c(x[1],x[1]),c(y1[1],y2[1]),type="l",xlim=range(x),ylim=range(c(y1,y2)),lab=c(5,highestplot-lowest+1,7),xlab="",ylab="Resolution Level")
  for(i in 2:length(x)){
	lines(c(x[i],x[i]),c(y1[i],y2[i]))
  }
}


####################################################
####################################################

maxithresh<-function(data,g,L=2,F=(log2(length(data))-1),eta=sqrt(6))
{
# maxithresh --Find optimal threshold level-by-level. 
#
#  Inputs 
#    data: y=f*g+noise. 
#    g: convolution kernel (sample of) 
#    L:  here should be Lowest set to L-1, to match C in FOurier domain
#    F: Finest resolution level
#    eta  smoothing parameter (default=sqrt(6))
#  Outputs 
#    y= vector of thresholds from level=L to F.
#
#  reference: ``Wavelet deconvolution in a periodic setting''
#  by  I.M. Johnstone, G. Kerkyacharian, D. Picard and M. Raimondo (2004).[JKPR]
#  The paper is available (as a prepint) by www from 
#		http://www.maths.usyd.edu.au/u/marcr/
   	  
#compute optimal threshold for deconvolution
# based on the Fourier transform of g.


s=scale(data);
n=length(g);
g_fft<-fft(g);


thr<-rep(0,F-L);

if (F==(log2(n)-1))
{ for (j in L:(F-1)){
thr[j-L+1] <- (sum(abs(g_fft[dyad(j)])^(-2))+sum(abs(g_fft[n-dyad(j)+1])^(-2)))/2^j ;
  }
aux1<-s*sqrt(thr)*sqrt(log(n)/n)*(eta/sqrt(4*pi));


#linear extrapolation for Jmax
n1<-length(aux1);

d1<-(aux1[n1]-aux1[1])/(n1-1);
aux2<-(n1+1)*d1;
y<-c(aux1, aux2);
} else { for (j in L:(F)){
thr[j-L+1] <- (sum(abs(g_fft[dyad(j)])^(-2))+sum(abs(g_fft[n-dyad(j)+1])^(-2)))/2^j ;
  }
aux1<-s*sqrt(thr)*sqrt(log(n)/n)*(eta/sqrt(4*pi));
y<-aux1;
}


    
}

#####################################
#####################################

threshsum<-function(w.res,L=3,F=(log2(length(w.res))-1))
#  
#  Usage
#    out = threshsum(w.res,L=3,F=log2(length(w))-1)
#  Inputs
#    w.res:residual after thresholding
#     
#    L    Coarsest Level of V_0;  L << J
#    F     Finest resolution level
#  Outputs
#   percentage  of thresholding per level 
#   j=C,L,L+1,...,F (incl.scaling)

{
  

res=rep(0,F-L+1);


#scaling coef



               sum1=sum(abs(w.res[1:(2^L)])<=0)
                res[1]=sum1/2^L;

       
  
#
#  Loop to Get Detail Coefficients for levels  j=L,...,F.
#
	for (j in L:F){ 
      sum1=sum(abs(w.res[dyad(j)])<=0)
                res[j-L+2]=sum1/2^j;
	  
  }
 


  
  y=1-res;
  
	
    

}
#########################################
#########################################
plotspec<-function(g,s)
{
#This function plots 
#the log spectrum and the optimal noise threshold
#which is used to find_j1

n=length(g);
fen=n/2-2;
g_fft=fft(g);
sigma=s;

speczoom(g,n)
aux2=g_fft[(1:(n/2))];
aux3=sqrt((1:(n/2)));
aux4=aux2/aux3;
aux5=rot90(aux4);
aux6=c(aux4,aux5);


#####

threshold_fft=log((sigma/sqrt(n)))+ log((log(sqrt(n)/sigma))^(1/2));
fen=n/2-2;
aux3=log(sqrt((1:fen/2)));
aux5=rot90(aux3);

aux6=c(aux5, aux3);
aux7=aux6+threshold_fft;

lines((-fen:(fen-1)),aux7,lty=2)
}

####################################################################
####################################################################

##########################################################
##########################################################
summary.wvd<-function(object,...)
{
#this function gives a  summary  
#of objects of class 'WaveD'
#as produced by the WaveD function
# y.wvd=WaveD(y,g)

y.wvd=object





lab=y.wvd$label
n=y.wvd$n;
t = (1:n)/n;


 cat("\nCall:\n", deparse(y.wvd$call), "\n\n", sep = "")



cat('Degree of Meyer wavelet =',y.wvd$deg,', Coarse resolution level=',y.wvd$L)
cat('\n')
cat('Sample size =',y.wvd$n,', Maximum resolution level=',log2(y.wvd$n)-1,'.')

cat('\n')
cat('WaveD optimal Fourier freq=',y.wvd$M,   '; WaveD optimal fine resolution level j1=',y.wvd$j1)
cat('\n')
if(y.wvd$j1!=y.wvd$F){
cat('Warning: WaveD finest resolution level has been set to F=',y.wvd$F)
cat('\n')}

opt.thr=round(maxithresh(y.wvd$data,y.wvd$g,eta=y.wvd$eta),3);
opt.thr=opt.thr[1:length(y.wvd$lev)];

aux=prod(opt.thr==round(y.wvd$thr,3));


if(aux){threshold='Maxiset threshold'}else{threshold='Manual threshold'}
cat('The choice of the threshold is:',threshold)
cat('\n')
cat('Thresholding policy=',y.wvd$POL, '.   Threshold constant gamma=',round(y.wvd$eta,3)) 
cat('\n')

cat('\n')
mat.thr=matrix(rep(0,3*length(y.wvd$thr)),ncol=3)
mat.thr[,2]=round(y.wvd$thr,3);
mat.thr[,3]=round(y.wvd$per,3);

nb=y.wvd$F-y.wvd$L+2;
aux=rep(0,nb);
aux[1]=round(max(abs(y.wvd$w[1:2^y.wvd$C])),3)
for (j in 2:nb){ 
aux[j]=round(max(abs(y.wvd$w[dyad(y.wvd$L+j-2)])),3)
}

mat.thr[,1]=aux;
dimnames(mat.thr) <- list(paste("level ",y.wvd$lev,"   ",sep=""),c('Max|w|','Threshold',' % of thresholding'))
print(mat.thr)

cat('\n')
cat('Noise-proxy statistics:\n')
cat(c('Estimated standard deviation= ',round(y.wvd$s,3)))
cat('\n')
cat(c('Shapiro test for normality, P=',round(y.wvd$p,5)))
cat('\n')
}
##############################################################
#############################################################
plot.wvd<-function(x,...)
{
#this fucntion gives (generic) plot  
#of objects of class 'wvd'
#as produced by the WaveD function
y.wvd=x

n=y.wvd$n;
t = (1:n)/n;
par(mfrow=c(2,2))
plot(t,y.wvd$data,type='l',main='Observations')
multires(y.wvd$w.thr,highestplot=y.wvd$F)
title(main='Thresholded FWaveD transform ')
plotspec(y.wvd$g,y.wvd$s)

title(main=cbind('Fourier domain, optimal resolution level F=',y.wvd$j1))
plot(t,y.wvd$waved,type='l',main='TI-WaveD estimate')
#readline()


#par(mfrow=c(2,2))
#multires(y.wvd$w,hi=y.wvd$F)
#title(main='FWaveD transform (without thresholding)')
#multires(y.wvd$w.thr,hi=y.wvd$F)
#title(main='Thresholded FWaveD transform ')
#plot(t,y.wvd$iw,type='l',main='Inverse WaveD transform (no thresholding)')
#plot(t,y.wvd$ordinary,type='l',main='Ordinary WaveD estimate')

cat('Hit enter to see the next plot\n')

readline()
par(mfrow=c(2,2))
par(pty='s')
noise.scaled=y.wvd$noise/y.wvd$s;
plot(noise.scaled,type='l',main=cbind('Noise proxy. Shapiro test for normality, P=',round(y.wvd$p,3)))
plot(seq(-4,4,le=n),dnorm(seq(-4,4,le=n)),type='l',lwd=2,lty=2,xlab=' ',ylab=' ',main='N(0,1) pdf:dashed curve, estimated density: plain curve')
lines(density(noise.scaled),lwd=2)

qqnorm(noise.scaled,cex=0.4)
qqline(noise.scaled,lwd=2)
boxplot(noise.scaled,main='Boxplot',horizontal=TRUE)
par(pty='m')
}


#   Written by Marc Raimondo and Michael Stewart,
# the University of Sydney, 2006. 
#   Copyright (c) 2005. 
#   Comments? e-mail marcr@maths.usyd.edu.au
#
# Setting up noisy-blur-signal model parameters. 
#          n=sample size
#          a=box-car width (preferrably  a quadratic irrational e.g. a=1/sqrt(5)
#          al, be=shape and scale of Gamma pdf. 
#           sigma_xxx=noise level.






##############################################################
# below are aux function to perform example.waved
#
#
#
#   This file contains all the functions to perform
#   example.wave.R, the script that implement the experiment 
#   for wavelet deconvolution. 
#
#
#
# Written for R by Marc Raimondo 
# and Michael Stewart, the University of Sydney, 2006. 
# Copyright (c) 2006. 
#
#  reference: ``Wavelet deconvolution in a periodic setting''
#  by  I.M. Johnstone, G. Kerkyacharian, D. Picard and M. Raimondo (2004).[JKPR]
#
#  The paper is available  by www from 
#		http://www.usyd.edu.au/u/marcr/
#


#################################################################
make.lidar<-function(n){
#
#Make artificial lidar signal
#
t<-(1:n)/n
1*(t>0.15)*(t<0.65)+1*(t>0.28)*(t<0.48)+(133.33*t-106.66)*(t>0.8)*(t<0.815)+(-133.33*t+110.6639)*(t>0.815)*(t<0.83)+(133.33*t-119.997)*(t>0.9)*(t<0.915)+(-133.33*t+123.9969)*(t>0.915)*(t<0.93)



}

########################################################3



##########################################################

BlurSignal<-function(f,g)
{
#
#compute the convolution of f and g
#

nf=length(f)
f.fft<-fft(f)
g.fft<-fft(g)
fg.fft<-f.fft*g.fft
Re(fft(fg.fft,inverse=TRUE))/nf



}
########################################################
###################################################################
make.doppler<-function(n)
#
#Make doppler signal
#
{
t<-(1:n)/n
 sig = (sqrt(t*(1-t)))*sin((2*pi*1.05)/(t+0.05));
}

#
# As in [DJ94]
#



################################
################################


#######################################
########################################
dyadjk<-function(j,k)
{
# dyadjk -- return Index of wavelet coefficient (j,k)-- 
#  Usage
#    ix = dyad(j,k);
#  Inputs
#    j     integers
#  Outputs
#     Index of wavelet coefficient (j,k) in the vector of wavelet coef.

ix=dyad(j)[1]+k;
   


}
# 


##############################################
#
#
# Below is the function that generates
# data examples of [RS]
#
################################################




waved.example<-function(pr=TRUE,gr=TRUE)
{

if(pr){
n=2^11 
al = 0.5
be=0.25
seed=11;
sigma.med=0.05
t = (1:n)/n

} 

if(!pr){
cat('Please enter the sample size (must be a power of 2, default n=2048, should be>256 to run all figures) \n')
n <- readline()
if(n==''){n=2048} 
n=as.numeric(n)
cat('\n')

cat('Please enter the noise sd ( default sd=0.05) \n')
sigma.med <- readline()
if(sigma.med==''){sigma.med=0.05}
sigma.med=as.numeric(sigma.med)

cat('\n')
cat('Please enter the Degree of Ill-Posedness (default=0.5) \n')
al <- readline()
if(al==''){al=0.5}
al=as.numeric(al)
cat('\n')

cat('Please enter the scale parameter of the blurring kernel (default=0.25) \n')
be <- readline()
if(be==''){be=0.25}
be=as.numeric(be)

cat('Please enter the seed number  (must be >0, default=11) \n')
seed <- readline()
if(seed==''){seed=11}
seed=as.numeric(seed)



cat('\n')

}

  


############################################################
t=(1:n)/n
temp<-dgamma(t,shape=al,scale=be);
aux1<-max(temp)
temp2<-temp/aux1
aux2<-sum(temp2)
GAMMA<-temp2/aux2
GAMMA=GAMMA;
LIDAR<-make.lidar(n)/1.7
DOPPLER=make.doppler(n)*1.9;
lidar.blur=BlurSignal(LIDAR,GAMMA)
doppler.blur=BlurSignal(DOPPLER,GAMMA)
set.seed(seed)
noise=rnorm(n);
noiseT=rt(n,2);
sdnoise=sqrt(var(noise));
noiseT=noiseT/sdnoise;

lidar.noisy=lidar.blur+sigma.med*noise
doppler.noisy=doppler.blur+sigma.med*noise

##t noise data
lidar.noisyT=lidar.blur+sigma.med*noiseT
doppler.noisyT=doppler.blur+sigma.med*noiseT


GAMMA_fft=fft(GAMMA);
GAMMA_fft_noisy=GAMMA_fft+sigma.med*fft(noise)/sqrt(n);
GAMMA_noisy=Re(fft(GAMMA_fft_noisy,inverse=TRUE))/n;
g=GAMMA;
g.noisy=GAMMA_noisy;


##t-noise
GAMMA_fft_noisyT=GAMMA_fft+sigma.med*fft(noiseT)/sqrt(n);
GAMMA_noisyT=Re(fft(GAMMA_fft_noisyT,inverse=TRUE))/n;
g.noisyT=GAMMA_noisyT;


############################
##
## Graphics display below if gr=TRUE



if(gr){
cat('-------------------------------------------------------\n')
	cat('Initializing noisy-blurred signals model:\n')
cat('sample size n =',n)
cat('\n')
cat('noise  sd = ',sigma.med)
cat('\n')
#
cat('Convolution kernel g:\n')
cat('gamma-distribution with shape paremeter=',al)
cat('\n')
cat('and scale parameter= ',be)
cat('\n')
cat('(effective) Degree of Ill-Posedness (DIP)= ',al)
cat('\n')
cat('The seed number has been set to ',seed)
cat('\n')



cat('\n')
cat('Blurred Signals to Noise Ratios:\n')


bsnrl<-round(10*log(var(lidar.blur)/sigma.med^2)/log(10),1)

cat('Lidar   BSNR(dB) =',bsnrl)
cat('\n')



bsnrl<-round(10*log(var(doppler.blur)/sigma.med^2)/log(10),1)

cat('Doppler   BSNR(dB) =',bsnrl)
cat('\n')




cat('-------------------------------------------------------\n')
###run figures below



if(pr){
cat('This is  the same set up as in: \n')
cat('The WaveD Transform in R by Raimondo and Stewart [RS].\n')
cat('To change  model parameters use data.demo=waved.example(F) \n')
}
 cat('\n')


if(!pr){
cat('This is not the same set up as in: \n')
cat('The WaveD Transform in R by Raimondo and Stewart [RS].\n')
cat('To get the same  model parameters as in [RS]  use data.demo=waved.example(T) \n')
cat('Warning: the figures/captions of [RS] may not be relevant when using new model parameters. \n')

}
 cat('\n')

########################################################

cat('Would you like to see the output? (Y/N) \n')
answer <- readline()
pr2=switch(answer,y=,Y=TRUE,FALSE)




cat('-------------------------------------------------------\n')

if(pr2){

######################################################
cat('Figure 1 shows the original LIDAR and doppler signals\n')
par(mfrow=c(1,2))


plot(t,LIDAR,type='l')
plot(t,DOPPLER,type='l')
cat('Hit enter to see the next plot\n')
readline()


cat('-------------------------------------------------------\n')

######################################################
cat('Figure 2 shows the blurred LIDAR and doppler signals\n')
par(mfrow=c(1,2))
plot(t,lidar.blur,type='l')
plot(t,doppler.blur,type='l')
cat('Hit enter to see the next plot\n')
readline()


cat('-------------------------------------------------------\n')
####################################################
cat('Figure 3 shows the noisy-blurred LIDAR and doppler signals\n')


par(mfrow=c(1,2))


plot(t,lidar.noisy,type='l')
plot(t,doppler.noisy,type='l')



##########################3
#cat('Hit enter to execute the WaveD command\n')
#readline()



Fmax=floor(log(length(lidar.noisy))/log(2))-1



lidar.wvd=WaveD(lidar.blur,g,F=6,thr=0)
lidar.noisy.wvd=WaveD(lidar.noisy,g,F=6,thr=0)


##############################



cat('Hit enter to see the next plot\n')
readline()




par(mfrow=c(2,2))



multires(lidar.wvd$w,lowest=3,highestplot=6)
multires(lidar.noisy.wvd$w,lowest=3,highestplot=6)


plot(t,lidar.wvd$iw,type='l')
 plot(t,lidar.noisy.wvd$iw,type='l')




cat('-------------------------------------------------------\n')
################################################

cat('Figure 4 of [RS] shows a wavelet decomposition of the LIDAR signal\n')
cat('Left: wavelet transform from noise free data, right...from noisy data\n')
cat('Top plots are wavelet coefficients according to time and resolution level\n')
cat('Bottom plots are corresponding inverse wavelet transforms\n')
cat('Left plots illustrate: large wavelet coefficeint nearby LIDAR discontinuities \n')
cat('Right plots illustrate: noise effect on wavelet coefficient increases with resolution level \n')



cat('Hit enter to see the next plot\n')
readline()


##########################################
cat('-------------------------------------------------------\n')
cat('Figure 5 of [RS] shows the effect of a large threshold (right) and  small threshold (left)  \n')




par(mfrow=c(1,2))
plot(t,WaveD(lidar.noisy,g,F=6,thr=0.2)$ord,type='l')
plot(t,WaveD(lidar.noisy,g,F=6,thr=0.02)$ord,type='l')

cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')


###########################################################
readline()

cat('Figure 6 of [RS] shows the effect of the Maxiset threshold, left plots: no thresholding, right plots: Maxiset Thresholding \n')


par(mfrow=c(2,2))



lidar.maxi.wvd=WaveD(lidar.noisy,g);

multires(lidar.maxi.wvd$w,lowest=3,highestplot=6)
multires(lidar.maxi.wvd$w.thr,lowest=3,highestplot=6)


 plot(t,WaveD(lidar.noisy,g,F=6,thr=0)$ord,type='l')
 plot(t,WaveD(lidar.noisy,g,F=6)$ord,type='l')





cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')
readline()
##############################################

cat('Figure 7 of [RS]: the Maxiset threshold do not prevent noise in high resolution level...  \n')
cat('Computing the full WaveD transform takes  approx.15sec...\n')
lidar.Fmax.wvd=WaveD(lidar.noisy,g,F=Fmax)
par(mfrow=c(1,2))
multires(lidar.Fmax.wvd$w.thr)
plot(t,lidar.Fmax.wvd$ord,type='l')

cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')
readline()




#######################################

par(mfrow=c(1,2))
cat('Figure 8 of [RS]: To prevent noise in high resolution levels \n')
cat(' WaveD is fitted with an  adaptive level selection in the Fourier domain \n')

cat('(left) noise free eigen values: maximum Fourier freq=', lidar.maxi.wvd$M,'\n')
cat('(left) noise free eigen values: maximum Resolution level=', lidar.maxi.wvd$F,'\n')

lidar.NEV.wvd=WaveD(lidar.noisy,g.noisy)

cat('(right) noisy  eigen values: maximum Fourier freq=', lidar.NEV.wvd$M,'\n')
cat('(right) noisy  eigen values: maximum Resolution level=', lidar.NEV.wvd$F,'\n')

plotspec(lidar.maxi.wvd$g,lidar.maxi.wvd$s)
plotspec(g.noisy,lidar.maxi.wvd$s)


cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')
readline()

#######################################

cat('Figure 9 of [RS]: cycle-spining improves the WaveD fit (right) \n')
cat('Ordinary WaveD (left), Translation Invariant WaveD (right)\n')

par(mfrow=c(1,2))
plot(t,lidar.maxi.wvd$ord,type='l')
plot(t,lidar.maxi.wvd$waved,type='l')

cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')


##########################################
readline()

cat('Figure [10] of [RS]: WaveD with Hard thresholding (left), soft thresholding (right)\n')

lidar.soft.wvd=WaveD(lidar.noisy,g,SOFT=TRUE)
plot(t,lidar.maxi.wvd$ord,type='l')
plot(t,lidar.soft.wvd$ord,type='l')


cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')
readline()



#############################################

cat('Figure [11],[12] of [RS]: WaveD analysis of the blurred Doppler in Gaussian noise\n')
cat('Illustration of the plot function for wvd objects\n')
cat('The first 4 plots illustrate the WaveD fit \n')
cat('The last 4 plots show a residual/noise analysis after a WaveD fit\n')




 doppler.wvd=WaveD(doppler.noisy,g);
plot(doppler.wvd)




cat('Hit enter to read the summary of the doppler.wvd object \n')
readline()
cat('-------------------------------------------------------\n')


summary(doppler.wvd)









############################################


cat('Hit enter to see the next plot\n')
cat('-------------------------------------------------------\n')
readline()





print('Figure 14,13: WaveD fit to lidar data with non-Gaussian noise (default setting)')
cat('Hit enter to plot  the lidarT.wvd object \n')
readline()

lidarT.wvd=WaveD(lidar.noisyT,g)



plot(lidarT.wvd)
print('The TI-WaveD estimate exhibits a large noise residual  even after thresholding')




cat('Hit enter to read the summary of the lidarT.wvd object \n')
readline()
cat('-------------------------------------------------------\n')

summary(lidarT.wvd)
print('The P-value of the Shapiro test=0 suggesting a non-Gaussian noise')

cat('Hit enter to see the next plot \n')
readline()


#############################################################


print('Figure 15:  Improved WaveD fit to lidar data with non-Gaussian noise')


lidarTS.wvd=WaveD(lidar.noisyT,g,SOFT=TRUE)
lidarT8.wvd=WaveD(lidar.noisyT,g,eta=sqrt(8))
lidarT12.wvd=WaveD(lidar.noisyT,g,eta=sqrt(12))

par(mfrow=c(2,2))

plot(t,lidarTS.wvd$ord,type='l',lwd=2,main='(a)')
lines(t,LIDAR,lwd=2,lty=2)
plot(t,lidarT8.wvd$ord,type='l',lwd=2,main='(b)')
lines(t,LIDAR,lwd=2,lty=2)
plot(t,lidarT8.wvd$waved,type='l',lwd=2,main='(c)')
lines(t,LIDAR,lwd=2,lty=2)
plot(t,lidarT12.wvd$waved,type='l',lwd=2,main='(d)')
lines(t,LIDAR,lwd=2,lty=2)




cat('Hit enter to see the next figure\n')
readline()
cat('-------------------------------------------------------\n')





#############################################################################
cat('Figure 16: LIDAR WaveD estimation with noisy eigen values\n')


plot(lidar.NEV.wvd)



############################################


cat('This is the last Figure.\n')

}
}

output=list(lidar.noisy=lidar.noisy, lidar.noisyT=lidar.noisyT,doppler.noisyT=doppler.noisyT,lidar.blur=lidar.blur, doppler.noisy=doppler.noisy, 
doppler.blur=doppler.blur,t=t,n=n,g=GAMMA,lidar=LIDAR,doppler=DOPPLER,seed=seed,sigma=sigma.med,
g.noisy=GAMMA_noisy,g.noisyT=g.noisyT,dip=be,k.scale=al)






}

