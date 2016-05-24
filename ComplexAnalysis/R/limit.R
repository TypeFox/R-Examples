limit <-
function(f,z0,z=NULL,track=FALSE){
#input check
if(is.finite(z0) & length(z)<1){stop("Missing arugment z")}
#edit specific directions
if(is.character(z)){
if(z=="left"){z<-z0-1};if(z=="right"){z<-z0+1}
if(z=="up"){z<-z0+1i};if(z=="down"){z<-z0-1i}
if(z=="upleft"){z<-z0-1+1i};if(z=="upright"){z<-z0+1+1i}
if(z=="downleft"){z<-z0-1-1i};if(z=="downright"){z<-z0+1-1i}
if(is.character(z)){stop("Invalid argument z")}}
# initialisation
n<-1                 #counting
record<-NULL         #initialising evaluation record
eps<-1e-3            #epsilon (differential)
divisor<-2           #multiplicative increment to eps
conv<-FALSE          #notifying convergence
#amendment to the divisor if limit to infinity
if(is.infinite(z0)){divisor<-0.2;eps<-1;
if(is.nan(Re(z0))){z<-sign(Im(z0))*1i}else{z<-sign(z0)};z0<-0}
# evaluation
repeat{
  L0<-z0+eps*exp(1i*Arg(z-z0));L<-f(L0)
  if(n>1){delta1<-Re(L-record[(n-1),4]);delta2<-Im(L-record[(n-1),4])}else{delta1<-Inf;delta2<-Inf}
  if(is.nan(L)==TRUE){L<-Inf};record<-rbind(record,c(n,eps,L0,L,delta1,delta2))
  if(is.nan(delta1) | is.nan(delta2) | is.infinite(eps)){break}
  if(n>1){if(is.infinite(delta1) | is.infinite(delta2)) {break}}

  S1<-Re(record[n,5]);S2<-Re(record[n,6])
  S3<-Re(record[n-1,5]);S4<-Re(record[n-1,6])
  Im.zero<-round(abs(Im(L)),3)
  Re.zero<-round(abs(Re(L)),3)
  if(n>3){
  if(Re.zero>0 & Im.zero==0 & abs(Re(delta1))<1e-6){if(sign(S1)!=sign(S3) | abs(S1)>abs(S3)) {conv<-TRUE;break}}
  if(Im.zero>0 & Re.zero==0 & abs(Re(delta2))<1e-6){if(sign(S2)!=sign(S4) | abs(S2)>abs(S4)) {conv<-TRUE;break}}
  if(Im.zero>0 & Re.zero>0){if(abs(Re(delta1))<1e-6 | abs(Re(delta2))<1e-6){if(sign(S1)!=sign(S3) | sign(S2)!=sign(S4) | abs(S1)>abs(S3) | abs(S2)>abs(S4)){conv<-TRUE;break}}}
  }
  eps<-eps/divisor
  n<-n+1
}
if(conv==TRUE){value<-record[n,4]}else{record<-record[1:(n-1),];value<-record[nrow(record),4]}
value<-unique(value)
colnames(record)<-c("count","epsilon","input","output","delta1","delta2")
if(track==TRUE){return(list(value=value,record=record))}else{return(value=value)}
}
