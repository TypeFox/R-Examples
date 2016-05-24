PMLE.Clayton.Weibull <-
function(l.trunc,x.trunc,GOF=TRUE){
l=l.trunc
x=x.trunc
  
B_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  RR = (S == rep(1,length(S)))
  S[RR] = 1-10^-8
  (1-S)^-alpha+u^-alpha-1
}

B_a_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  RR = (S == rep(1,length(S)))
  S[RR] = 1-10^-8
  -log(1-S)/(1-S)^alpha-log(u)/u^alpha
}

B_L_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  -alpha*S*V/(1-S)^(alpha+1)
}

B_X_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX-1)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  -alpha*lamdaL*S*(nuL/nuX)*V*log(1-u)/(lamdaX^2*(1-S)^(alpha+1))
}

B_nL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  -alpha*lamdaL*S*V*(1/nuX)*log(-1/lamdaX*log(1-u))/(1-S)^(alpha+1)
}

B_nX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  -alpha*lamdaL*S*V*(-nuL/nuX^2)*log(-1/lamdaX*log(1-u))/(1-S)^(alpha+1)
}

B_aa_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  RR = (S == rep(1,length(S)))
  S[RR] = 1-10^-8
  (log(1-S))^2/(1-S)^alpha+(log(u))^2/u^alpha
}

B_LL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_L = B_L_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-log(1-u)/lamdaX)^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  ((alpha+1)*S/(1-S)+1)*B_L*(-V)
}

B_XX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_X = B_X_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-log(1-u)/lamdaX)^(nuL/nuX-1)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (lamdaL*nuL*(-alpha-1)*S*V*log(1-u)/(lamdaX*nuX*(1-S))-
     lamdaL*nuL/(lamdaX*nuX)*V*log(1-u)-(nuL/nuX-1)-2)*B_X/lamdaX
}

B_nLnL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_nL = B_nL_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-log(1-u)/lamdaX)^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((-alpha-1)*S/(1-S)-1)*lamdaL*V+1)*B_nL*log(-log(1-u)/lamdaX)/nuX
}

B_nXnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_nX = B_nX_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-log(1-u)/lamdaX)^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*lamdaL*nuL/nuX^2*V*log(-log(1-u)/lamdaX)-
     nuL*log(-log(1-u)/lamdaX)/nuX^2-2/nuX)*B_nX
}

B_aL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-log(1-u)/lamdaX)^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (alpha*log(1-S)-1)/(1-S)^(alpha+1)*V*S
}

B_aX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX-1)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (alpha*log(1-S)-1)/(1-S)^(alpha+1)*
    lamdaL*nuL/(lamdaX^2*nuX)*V*log(1-u)*S
}

B_anL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (alpha*log(1-S)-1)/(1-S)^(alpha+1)*lamdaL*V*log(-log(1-u)/lamdaX)/nuX*S
}

B_anX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  -(alpha*log(1-S)-1)/(1-S)^(alpha+1)*lamdaL*V*log(-log(1-u)/lamdaX)*nuL/nuX^2*S
}

B_LX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_L = B_L_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*lamdaL*V-1)*nuL/(nuX*lamdaX)*B_L
}

B_LnL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_L = B_L_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*(-lamdaL*V)+1)/nuX*log(-1/lamdaX*log(1-u))*B_L
}

B_LnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_L = B_L_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*(-lamdaL*V)+1)*(-nuL/nuX^2)*log(-1/lamdaX*log(1-u))*B_L
}

B_XnL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_X = B_X_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*(1/nuX*log(-1/lamdaX*log(1-u)))*
     (-lamdaL*V)+1/nuX*log(-1/lamdaX*log(1-u))+1/nuL)*B_X
}

B_XnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_X = B_X_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*(-nuL/nuX^2)*log(-1/lamdaX*log(1-u))*
     (-lamdaL*V)-nuL/nuX^2*log(-1/lamdaX*log(1-u))-1/nuX)*B_X
}

B_nLnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B_nL = B_nL_func(u)
  S = exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX))
  V = (-1/lamdaX*log(1-u))^(nuL/nuX)
  RR = (S == rep(1,length(S)))
  RV = (V == rep(Inf,length(V)))
  S[RR] = 1-10^-8
  V[RV] = 10^8
  (((alpha+1)*S/(1-S)+1)*((-nuL/nuX^2)*log(-1/lamdaX*log(1-u)))*
     (-lamdaL*V)+((-nuL/nuX^2)*log(-1/lamdaX*log(1-u))-1/nuX))*B_nL
}

######
H_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  u^(-alpha-1)*B^(-1/alpha-1)
}
########

H_a_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  H = H_func(u)
  B = B_func(u)
  B_a = B_a_func(u)
  (H*(-log(u)+log(B)/alpha^2+(-1/alpha-1)*B_a/B))*alpha
}

alpha_tuda = log(8)
lamdaL_tuda = log(2)
lamdaX_tuda = log(1)
nuL_tuda = log(2)
nuX_tuda = log(1)
integrate(H_a_func,0,1)
a=integrate(H_func,0,1)$value
alpha_tuda = alpha_tuda+10^-8
b=integrate(H_func,0,1)$value
(b-a)/10^-8


########
H_L_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_L = B_L_func(u)
  ((-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_L)*lamdaL
}

########
H_X_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_X = B_X_func(u)
  ((-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_X)*lamdaX
}

#########
H_nL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nL = B_nL_func(u)
  ((-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_nL)*nuL
}

#########
H_nX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nX = B_nX_func(u)
  ((-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_nX)*nuX
}

#########
H_aa_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  H = H_func(u)
  H_a = H_a_func(u)
  B = B_func(u)
  B_a = B_a_func(u)
  B_aa = B_aa_func(u)
  H_a+((-log(u)+log(B)/alpha^2+(-1/alpha-1)*B_a/B)*H_a/alpha+
         H*(2*B_a/(alpha^2*B)-2*log(B)/alpha^3+
              (-1/alpha-1)*(B_aa*B-B_a^2)/B^2))*alpha^2
}

###########
H_LL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  H_L = H_L_func(u)
  B = B_func(u)
  B_L = B_L_func(u)
  B_LL = B_LL_func(u)
  H_L+((-1/alpha-2)*B_L/B*H_L/lamdaL+(-1/alpha-1)*B_LL/(u^(alpha+1)*B^(1/alpha+2)))*lamdaL^2
}

############
H_XX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_X = B_X_func(u)
  B_XX = B_XX_func(u)
  H_X = H_X_func(u)
  H_X+((-1/alpha-2)*B_X/B*H_X/lamdaX+(-1/alpha-1)*B_XX/(u^(alpha+1)*B^(1/alpha+2)))*lamdaX^2
}

###########
H_nLnL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nL = B_nL_func(u)
  B_nLnL = B_nLnL_func(u)
  H_nL = H_nL_func(u)
  H_nL+((-1/alpha-2)*B_nL/B*H_nL/nuL+u^(-alpha-1)*(-1/alpha-1)*B^(-1/alpha-2)*B_nLnL)*nuL^2
}

H_nXnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nX = B_nX_func(u)
  B_nXnX = B_nXnX_func(u)
  H_nX = H_nX_func(u)
  H_nX+((-1/alpha-2)*B_nX/B*H_nX/nuX+u^(-alpha-1)*(-1/alpha-1)*B^(-1/alpha-2)*B_nXnX)*nuX^2
}

###########
H_aL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_a = B_a_func(u)
  B_L = B_L_func(u)
  B_aL = B_aL_func(u)
  H_a = H_a_func(u)
  ((-1/alpha-1)*B_L/B*H_a/alpha+u^(-alpha-1)*B^(-1/alpha-1)*
     (B_L/(B*alpha^2)+(-1/alpha-1)*(B_aL*B-B_L*B_a)/B^2))*lamdaL*alpha
}

###########
H_aX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_a = B_a_func(u)
  B_X = B_X_func(u)
  B_aX = B_aX_func(u)
  H_a = H_a_func(u)
  ((-1/alpha-1)*B_X/B*H_a/alpha+u^(-alpha-1)*B^(-1/alpha-1)*
     (B_X/(B*alpha^2)+(-1/alpha-1)*(B_aX*B-B_X*B_a)/B^2))*lamdaX*alpha
}

###########
H_anL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_a = B_a_func(u)
  B_nL = B_nL_func(u)
  B_anL = B_anL_func(u)
  H_a = H_a_func(u)
  ((-1/alpha-1)*B_nL/B*H_a/alpha+u^(-alpha-1)*B^(-1/alpha-1)*
     (B_nL/(B*alpha^2)+(-1/alpha-1)*(B_anL*B-B_nL*B_a)/B^2))*nuL*alpha
}

###########
H_anX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_a = B_a_func(u)
  B_nX = B_nX_func(u)
  B_anX = B_anX_func(u)
  H_a = H_a_func(u)
  ((-1/alpha-1)*B_nX/B*H_a/alpha+u^(-alpha-1)*B^(-1/alpha-1)*
     (B_nX/(B*alpha^2)+(-1/alpha-1)*(B_anX*B-B_nX*B_a)/B^2))*nuX*alpha
}

##########
H_LX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_X = B_X_func(u)
  B_LX = B_LX_func(u)
  H_L = H_L_func(u)
  ((-1/alpha-2)*B_X/B*H_L/lamdaL+(-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_LX)*lamdaL*lamdaX
}

##########
H_LnL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nL = B_nL_func(u)
  B_LnL = B_LnL_func(u)
  H_L = H_L_func(u)
  ((-1/alpha-2)*B_nL/B*H_L/lamdaL+(-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_LnL)*lamdaL*nuL
}

###########
H_LnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nX = B_nX_func(u)
  B_LnX = B_LnX_func(u)
  H_L = H_L_func(u)
  ((-1/alpha-2)*B_nX/B*H_L/lamdaL+(-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_LnX)*lamdaL*nuX
}

##############
H_XnL_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nL = B_nL_func(u)
  B_XnL = B_XnL_func(u)
  H_X = H_X_func(u)
  ((-1/alpha-2)*B_nL/B*H_X/lamdaX+(-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_XnL)*lamdaX*nuL
}

##############
H_XnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nX = B_nX_func(u)
  B_XnX = B_XnX_func(u)
  H_X = H_X_func(u)
  ((-1/alpha-2)*B_nX/B*H_X/lamdaX+(-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_XnX)*lamdaX*nuX
}

##############
H_nLnX_func = function(u){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  B = B_func(u)
  B_nX = B_nX_func(u)
  B_nLnX = B_nLnX_func(u)
  H_nL = H_nL_func(u)
  ((-1/alpha-2)*B_nX/B*H_nL/nuL+(-1/alpha-1)*u^(-alpha-1)*B^(-1/alpha-2)*B_nLnX)*nuL*nuX
}

##################
A_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  FL=1-exp(-lamdaL*l^nuL)
  FX^(-alpha)+FL^(-alpha)-1
}

A_a_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  FL=1-exp(-lamdaL*l^nuL)
  -log(FL)/(FL^alpha)-log(FX)/(FX^alpha)
}

A_L_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FL=1-exp(-lamdaL*l^nuL)
  -alpha*l^nuL*(1-FL)/ (FL^(alpha+1))
}

A_X_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  -alpha*x^nuX*(1-FX)/ (FX^(alpha+1))
}

A_nL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FL=1-exp(-lamdaL*l^nuL)
  -alpha*lamdaL*l^nuL*log(l)*(1-FL)/(FL^(alpha+1))
}

A_nX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  -alpha*lamdaX*x^nuX*log(x)*(1-FX)/(FX^(alpha+1))
}

A_aa_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  FL=1-exp(-lamdaL*l^nuL)
  (log(FL))^2/(FL^alpha)+(log(FX))^2/(FX^alpha)
}

##################
pdf_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  fx= nuX*lamdaX*x^(nuX-1)*exp(-lamdaX*x^nuX)
  fy= nuL*lamdaL*l^(nuL-1)* exp(-lamdaL*l^nuL)
  FX=1-exp(-lamdaX*x^nuX)
  FL=1-exp(-lamdaL*l^nuL)
  (1+alpha)*fy*fx*((FL)*(FX))^(-alpha-1)/A^(1/alpha+2)
}

##################
pdf_a_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  FL=1-exp(-lamdaL*l^nuL)
  (1/(1+alpha)-log(FL)-log(FX)+1/alpha^2*log(A)+(-1/alpha-2)*A_a/A)*alpha
}

##################
pdf_L_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  ((-alpha-1)*(l^nuL*(1-FL))/(FL)+(-1/alpha-2)*A_L/A+1/lamdaL-l^nuL)*lamdaL
}

##################
pdf_X_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  ((-alpha-1)*(x^nuX*(1-FX))/(FX)+(-1/alpha-2)*A_X/A+1/lamdaX-x^nuX)*lamdaX
}

#################
pdf_nL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  ((-alpha-1)*(lamdaL*l^nuL*log(l)*(1-FL))/(FL)+
     (-1/alpha-2)*A_nL/A+(1-lamdaL*l^nuL)*log(l)+1/nuL)*nuL
}

##################
pdf_nX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  ((-alpha-1)*(lamdaX*x^nuX*log(x)*(1-FX))/(FX)+
     (-1/alpha-2)*A_nX/A+(1-lamdaX*x^nuX)*log(x)+1/nuX)*nuX
}

###############
pdf_aa_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_aa = A_aa_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_a = pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_a+(-1/(1+alpha)^2-2/alpha^3*log(A)+2/alpha^2*A_a/A+
           (-1/alpha-2)*(A_aa*A-A_a^2)/A^2)*alpha^2
}

################
A_LL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  ((-alpha-1)*l^nuL*(1-FL)/(FL)-l^nuL)*A_L
}

pdf_LL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_LL = A_LL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x) 
  pdf_L = pdf_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  pdf_L+((alpha+1)*(l^(2*nuL)*(1-FL))/(FL)^2+
           (-1/alpha-2)*(A_LL*A-A_L^2)/A^2-1/lamdaL^2)*lamdaL^2
}

A_XX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  ((-alpha-1)*x^nuX*(1-FX)/(FX)-x^nuX)*A_X
}

pdf_XX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_XX = A_XX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x) 
  pdf_X = pdf_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  pdf_X+((alpha+1)*(x^(2*nuX)*(1-FX))/(FX)^2+
           (-1/alpha-2)*(A_XX*A-A_X^2)/A^2-1/lamdaX^2)*lamdaX^2
}

#############
A_nLnL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  (((alpha+1)*(1-FL)/(FL)+1)*
     (-lamdaL*l^nuL*log(l))+log(l))*A_nL
}

pdf_nLnL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nLnL = A_nLnL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nL = pdf_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  
  pdf_nL+((-alpha-1)*lamdaL*l^nuL*log(l)*(1-FL)/(FL)*
            (((1-FL)/(FL)+1)*(-lamdaL*l^nuL*log(l))+log(l))+
            (-1/alpha-2)*(A_nLnL*A-A_nL^2)/A^2-lamdaL*l^nuL*(log(l))^2-1/nuL^2)*nuL^2
}

##############
A_nXnX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  
  (((alpha+1)*(1-FX)/(FX)+1)*
     (-lamdaX*x^nuX*log(x))+log(x))*A_nX
}

pdf_nXnX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nXnX = A_nXnX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nX = pdf_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  
  pdf_nX+((-alpha-1)*lamdaX*x^nuX*log(x)*(1-FX)/(FX)*
            (((1-FX)/(FX)+1)*(-lamdaX*x^nuX*log(x))+log(x))+
            (-1/alpha-2)*(A_nXnX*A-A_nX^2)/A^2-lamdaX*x^nuX*(log(x))^2-1/nuX^2)*nuX^2
}

A_aL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FL=1-exp(-lamdaL*l^nuL)
  (1-alpha*log(FL))*-l^nuL*(1-FL)/(FL)^(alpha+1)
}

pdf_aL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_aL = A_aL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  (-l^nuL*(1-FL)/(FL)+A_L/(alpha^2*A)+
     (-1/alpha-2)*(A_aL*A-A_a*A_L)/A^2)*alpha*lamdaL
}

##################
A_aX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  (1-alpha*log(FX))*-x^nuX*(1-FX)/(FX)^(alpha+1)
}

pdf_aX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_aX = A_aX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  (-x^nuX*(1-FX)/(FX)+A_X/(alpha^2*A)+
     (-1/alpha-2)*(A_aX*A-A_a*A_X)/A^2)*alpha*lamdaX
}

###############
A_anL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FL=1-exp(-lamdaL*l^nuL)
  (1-alpha*log(FL))*-lamdaL*l^nuL*log(l)*(1-FL) /
    (FL)^(alpha+1)
}

pdf_anL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_anL = A_anL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  (-lamdaL*l^nuL*log(l)*(1-FL)/(FL)+A_nL/(alpha^2*A)+
     (-1/alpha-2)*(A_anL*A-A_a*A_nL)/A^2)*alpha*nuL
}

###############
A_anX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  FX=1-exp(-lamdaX*x^nuX)
  FL=1-exp(-lamdaL*l^nuL)
  (1-alpha*log(FX))*-lamdaX*x^nuX*log(x)*(1-FX)/(FX)^(alpha+1)
}

pdf_anX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_a = A_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_anX = A_anX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  (-lamdaX*x^nuX*log(x)*(1-FX)/(FX)+A_nX/(alpha^2*A)+
     (-1/alpha-2)*(A_anX*A-A_a*A_nX)/A^2)*alpha*nuX
}

##############
A_LnL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  (((alpha+1)*(1-FL)/(FL)+1)*(-lamdaL*l^nuL)+1)*log(l)*A_L
}

pdf_LnL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_LnL = A_LnL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FL=1-exp(-lamdaL*l^nuL)
  ((-alpha-1)*l^nuL*log(l)*(1-FL)/(FL)*(((1-FL)/(FL)+1)*(-lamdaL*l^nuL)+1)+
     (-1/alpha-2)*(A_LnL*A-A_L*A_nL)/A^2-l^nuL*log(l))*lamdaL*nuL
}

##############
pdf_LnX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  (-1/alpha-2)*(-A_L*A_nX)/A^2*lamdaL*nuX
}

############
A_XnX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  (((alpha+1)*(1-FX)/(FX)+1)*(-lamdaX*x^nuX)+1)*log(x)*A_X
}

pdf_XnX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_XnX = A_XnX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  FX=1-exp(-lamdaX*x^nuX)
  ((-alpha-1)*x^nuX*log(x)*(1-FX)/(FX)*(((1-FX)/(FX)+1)*(-lamdaX*x^nuX)+1)+
     (-1/alpha-2)*(A_XnX*A-A_X*A_nX)/A^2-x^nuX*log(x))*lamdaX*nuX
}

##############
pdf_XnL_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  (-1/alpha-2)*(-A_X*A_nL)/A^2*lamdaX*nuL
}

################
pdf_nLnX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nL = A_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_nX = A_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  (-1/alpha-2)*(-A_nX*A_nL)/A^2*nuL*nuX
}

#############
pdf_LX_func = function(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x){
  alpha = exp(alpha_tuda)
  lamdaL = exp(lamdaL_tuda)
  lamdaX = exp(lamdaX_tuda)
  nuL = exp(nuL_tuda)
  nuX = exp(nuX_tuda)
  A = A_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_L = A_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  A_X = A_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  (-1/alpha-2)*(-A_X*A_L)/A^2*lamdaL*lamdaX
}

score_func = function(er){
  prob = integrate(H_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_a = integrate(H_a_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_L = integrate(H_L_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_X = integrate(H_X_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_nL = integrate(H_nL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_nX = -prob_nL
  
  pdf_a = pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_L = pdf_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_X = pdf_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nL = pdf_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nX = pdf_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  
  score = c(-n*prob_a/prob+sum(pdf_a),-n*prob_L/prob+sum(pdf_L),-n*prob_X/prob+sum(pdf_X),
            -n*prob_nL/prob+sum(pdf_nL),-n*prob_nX/prob+sum(pdf_nX))
  score
}


hessian_func = function(er){
  prob = integrate(H_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_a = integrate(H_a_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_L = integrate(H_L_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_X = integrate(H_X_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_nL = integrate(H_nL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_nX = -prob_nL
  prob_aa = integrate(H_aa_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_LL = integrate(H_LL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_XX = integrate(H_XX_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_nLnL = integrate(H_nLnL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_nXnX = prob_nLnL
  prob_aL = integrate(H_aL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_aX = integrate(H_aX_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_anL = integrate(H_anL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_anX = -prob_anL
  prob_LX = integrate(H_LX_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_LnL = integrate(H_LnL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_LnX = -prob_LnL
  prob_XnL = integrate(H_XnL_func, lower = 0, upper = 1,rel.tol = er)$value
  prob_XnX = -prob_XnL
  prob_nLnX = -prob_nLnL
  
  pdf_a = pdf_a_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_L = pdf_L_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_X = pdf_X_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nL = pdf_nL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nX = pdf_nX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_aa = pdf_aa_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_LL = pdf_LL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_XX = pdf_XX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nLnL = pdf_nLnL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nXnX = pdf_nXnX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_aL = pdf_aL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_aX = pdf_aX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_anL = pdf_anL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_anX = pdf_anX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_LX = pdf_LX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_LnL = pdf_LnL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_LnX = pdf_LnX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_XnL = pdf_XnL_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_XnX = pdf_XnX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  pdf_nLnX = pdf_nLnX_func(alpha_tuda,lamdaL_tuda,lamdaX_tuda,nuL_tuda,nuX_tuda,l,x)
  
  hessian = matrix(c(
    -n/prob^2*(prob_aa*prob-prob_a^2)+n*mean(pdf_aa,na.rm = T),
    -n/prob^2*(prob_aL*prob-prob_a*prob_L)+n*mean(pdf_aL,na.rm = T),
    -n/prob^2*(prob_aX*prob-prob_a*prob_X)+n*mean(pdf_aX,na.rm = T),
    -n/prob^2*(prob_anL*prob-prob_a*prob_nL)+n*mean(pdf_anL,na.rm = T),
    -n/prob^2*(prob_anX*prob-prob_a*prob_nX)+n*mean(pdf_anX,na.rm = T),
    -n/prob^2*(prob_aL*prob-prob_a*prob_L)+n*mean(pdf_aL,na.rm = T),
    -n/prob^2*(prob_LL*prob-prob_L^2)+n*mean(pdf_LL,na.rm = T),
    -n/prob^2*(prob_LX*prob-prob_X*prob_L)+n*mean(pdf_LX,na.rm = T),
    -n/prob^2*(prob_LnL*prob-prob_L*prob_nL)+n*mean(pdf_LnL,na.rm = T),
    -n/prob^2*(prob_LnX*prob-prob_L*prob_nX)+n*mean(pdf_LnX,na.rm = T),
    -n/prob^2*(prob_aX*prob-prob_a*prob_X)+n*mean(pdf_aX,na.rm = T),
    -n/prob^2*(prob_LX*prob-prob_X*prob_L)+n*mean(pdf_LX,na.rm = T),
    -n/prob^2*(prob_XX*prob-prob_X^2)+n*mean(pdf_XX,na.rm = T),
    -n/prob^2*(prob_XnL*prob-prob_X*prob_nL)+n*mean(pdf_XnL,na.rm = T),
    -n/prob^2*(prob_XnX*prob-prob_X*prob_nX)+n*mean(pdf_XnX,na.rm = T),
    -n/prob^2*(prob_anL*prob-prob_a*prob_nL)+n*mean(pdf_anL,na.rm = T),
    -n/prob^2*(prob_LnL*prob-prob_L*prob_nL)+n*mean(pdf_LnL,na.rm = T),
    -n/prob^2*(prob_XnL*prob-prob_X*prob_nL)+n*mean(pdf_XnL,na.rm = T),
    -n/prob^2*(prob_nLnL*prob-prob_nX^2)+n*mean(pdf_nLnL,na.rm = T),
    -n/prob^2*(prob_nLnX*prob-prob_nL*prob_nX)+n*mean(pdf_nLnX,na.rm = T),
    -n/prob^2*(prob_anX*prob-prob_a*prob_nX)+n*mean(pdf_anX,na.rm = T),
    -n/prob^2*(prob_LnX*prob-prob_L*prob_nX)+n*mean(pdf_LnX,na.rm = T),
    -n/prob^2*(prob_XnX*prob-prob_X*prob_nX)+n*mean(pdf_XnX,na.rm = T),
    -n/prob^2*(prob_nLnX*prob-prob_nL*prob_nX)+n*mean(pdf_nLnX,na.rm = T),
    -n/prob^2*(prob_nXnX*prob-prob_nX^2)+n*mean(pdf_nXnX,na.rm = T)
  ),5,5,byrow = TRUE)
  hessian
}


rel_tol_func = function(CI){
  
  e =.Machine$double.eps^0.25*10^sample(-4:2,1) 
  repeat{
    EE = ( 
      try(integrate(H_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_a_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_L_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_X_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_nL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &        
        try(integrate(H_aa_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_LL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_XX_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_nLnL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_aL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_aX_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_anL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_LX_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_LnL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0 &
        try(integrate(H_XnL_func, lower = 0, upper = 1,rel.tol = e)$value*0,silent = T) == 0)
    if( EE == T | CI==3 ){break}
    if( EE == F ){ e =.Machine$double.eps^0.25*10^CI ; CI=CI+1} 
    
  }
  e
}

g_func = function(lamdaX,nuX){
  gamma(1+1/exp(nuX))/(exp(lamdaX)^(1/exp(nuX)))
}
g_lamdaX_func = function(lamdaX,nuX){
  (-1/exp(nuX))*gamma(1+1/exp(nuX))/(exp(lamdaX)^(1/exp(nuX)))
}
g_nuX_func = function(lamdaX,nuX){
  (digamma(1+1/exp(nuX))*gamma(1+1/exp(nuX))*(-exp(-nuX))+
     gamma(1+1/exp(nuX))*lamdaX*exp(-nuX))/(exp(lamdaX)^(1/exp(nuX)))
  
}


O_B_func = function(u){
  (1-exp(-lamdaL*(-1/lamdaX*log(1-u))^(nuL/nuX)))^-alpha+u^-alpha-1
}
O_h_func = function(u){
  O_B = O_B_func(u)
  u^(-alpha-1)*O_B^(-1/alpha-1)
}
O_A_func = function(alpha,lamdaL,lamdaX,nuL,nuX,l,x){
  (1-exp(-lamdaL*l^nuL))^-alpha+(1-exp(-lamdaX*x^nuX))^-alpha-1
}
O_pdf_func = function(alpha,lamdaL,lamdaX,nuL,nuX,l,x){
  O_A = O_A_func(alpha,lamdaL,lamdaX,nuL,nuX,l,x)
  
  (1+alpha)*nuL*lamdaL*l^(nuL-1)*exp(-lamdaL*l^nuL)*
    nuX*lamdaX*x^(nuX-1)*exp(-lamdaX*x^nuX)*
    ((1-exp(-lamdaL*l^nuL))*(1-exp(-lamdaX*x^nuX)))^(-alpha-1)/O_A^(1/alpha+2)
}

n=length(l)
par_new = matrix(c(NA,NA,NA,NA,NA),5,1)
par_old = matrix(c(NA,NA,NA,NA,NA),5,1)
par_A = c()
par_L = c()
par_X = c()
par_nL = c()
par_nX = c()


tau=cor(l,x,method = "kendall")
par_old[1,1] = log(2*tau/(1-tau))
par_old[2,1] = log(1/mean(l))
par_old[3,1] = log(1/mean(x))
par_old[4,1] = log(1)
par_old[5,1] = log(1)

AI = 1
BI = 1
alpha_tuda = par_old[1,1];lamdaL_tuda = par_old[2,1];lamdaX_tuda = par_old[3,1];
nuL_tuda = par_old[4,1];nuX_tuda = par_old[5,1]
par_A[AI] = par_old[1,1];par_L[AI] = par_old[2,1];par_X[AI] = par_old[3,1];
par_nL[AI] = par_old[4,1];par_nX[AI] = par_old[5,1]

e = rel_tol_func(-6)
repeat{    
  par_new = par_old-solve(hessian_func(e))%*%score_func(e)      
  error1 = abs(exp(par_new[1,1])-exp(par_old[1,1]))
  error2 = abs(exp(par_new[2,1])-exp(par_old[2,1]))
  error3 = abs(exp(par_new[3,1])-exp(par_old[3,1]))
  error4 = abs(exp(par_new[4,1])-exp(par_old[4,1]))
  error5 = abs(exp(par_new[5,1])-exp(par_old[5,1]))
  Error1 = max(error1,error2,error3,error4,error5)
  alpha_tuda = par_new[1,1]
  lamdaL_tuda = par_new[2,1]
  lamdaX_tuda = par_new[3,1]
  nuL_tuda = par_new[4,1]
  nuX_tuda = par_new[5,1]
  e =rel_tol_func(-6)
  Error2 = try(max(eigen(hessian_func(e))$value),silent = T)
  Error3 = ( exp(par_new[1,1]) <= 10^-6 | exp(par_new[1,1]) > 20 | Error1 > 2 
             |min(exp(par_new[2:5,1])) <= 10^-8 |AI/100 == floor(AI/100)  
             | try(Error2*0,silent=T)!=0 )     
  if( Error1 <10^-4 & Error2 < 0 ){break}  
  if( AI> 100000 ){break} 
  if( Error3 == T ){   
    BI=BI+1
    repeat{
      par_new[1,1] = log(2*tau/(1-tau)*exp(runif(1,-1,1)))
      par_new[2,1] = log(1/mean(l)*exp(runif(1,-0.5,0.5)))
      par_new[3,1] = log(1/mean(x)*exp(runif(1,-0.5,0.5)))
      par_new[4,1] = runif(1,-0.5,0.5)
      par_new[5,1] = runif(1,-0.5,0.5)
      alpha_tuda = par_new[1,1]
      lamdaL_tuda = par_new[2,1]
      lamdaX_tuda = par_new[3,1]
      nuL_tuda = par_new[4,1]
      nuX_tuda = par_new[5,1]
      e =rel_tol_func(-6)
      Error2 = try(max(eigen(hessian_func(e))$value),silent = T)  
      if( exp(par_new[1,1]) > 10^-6 & exp(par_new[1,1]) < 20 & 
            try(Error2*0,silent=T)==0 ){break}    
    }
  }
  par_old = par_new
  AI = AI+1
  par_A[AI] = par_old[1,1];par_L[AI] = par_old[2,1];par_X[AI] = par_old[3,1];
  par_nL[AI] = par_old[4,1];par_nX[AI] = par_old[5,1] 
  e = rel_tol_func(-6) 
}
alpha_tuda = par_new[1,1]
lamdaL_tuda = par_new[2,1]
lamdaX_tuda = par_new[3,1]
nuL_tuda = par_new[4,1]
nuX_tuda = par_new[5,1]



##Table 4 (weibull)
alpha = exp( par_new[1,1])
lamdaL = exp( par_new[2,1])
lamdaX = exp( par_new[3,1])
nuL = exp( par_new[4,1])
nuX = exp( par_new[5,1])
mean_L = gamma(1+1/nuL)/(lamdaL^(1/nuL))
mean_X = gamma(1+1/nuX)/(lamdaX^(1/nuX))

TT = solve(-hessian_func(e))
se_alpha=exp(alpha_tuda)*sqrt(TT[1,1])
se_lamdaL=exp(lamdaL_tuda)*sqrt(TT[2,2])
se_lamdaX=exp(lamdaX_tuda)*sqrt(TT[3,3])
se_nuL=exp(nuL_tuda)*sqrt(TT[4,4])
se_nuX=exp(nuX_tuda)*sqrt(TT[5,5])
se_mu=sqrt(t(matrix(c(0,0,g_lamdaX_func(lamdaX_tuda,nuX_tuda),0,
                      g_nuX_func(lamdaX_tuda,nuX_tuda)),5,1))%*%TT%*%
             matrix(c(0,0,g_lamdaX_func(lamdaX_tuda,nuX_tuda),0,
                      g_nuX_func(lamdaX_tuda,nuX_tuda)),5,1))
alpha_res = c(EST = round(alpha,4),SE = round(se_alpha,4))
lambda_L_res = c(EST = round(lamdaL,6),SE = round(se_lamdaL,6))
lambda_X_res = c(EST = round(lamdaX,7),SE = round(se_lamdaX,7))
nu_L_res = c(EST = round(nuL,4),SE = round(se_nuL,4))
nu_X_res = c(EST = round(nuX,4),SE = round(se_nuX,4))
mu_res = c(EST = round(mean_X,4),SE = round(se_mu,4))
LL=-n*log(integrate(O_h_func, lower = 0, upper = 1)$value)+
  sum(log(O_pdf_func(alpha,lamdaL,lamdaX,nuL,nuX,l,x)))
AIC_res = round(-2*LL+2*5,2)
BIC_res = round(-2*LL+5*log(n),2)

C.test=K.test=NULL
F_par=F_emp=prop=NULL

if(GOF==TRUE){

F.func=function(ll,xx){
  Fll=1-exp(-lamdaL*ll^nuL)
  Fxl=1-exp(-lamdaX*ll^nuX)
  Fxx=1-exp(-lamdaX*xx^nuX)
  F1 = integrate(H_func, lower = 0, upper = Fxl)$value
  HH_func = function(u){
    B=Fll^(-alpha)+u^(-alpha)-1
    u^(-alpha-1)*B^(-1/alpha-1)
  }
  F2 = integrate(HH_func, lower = Fxl, upper = Fxx)$value
  F1+F2
}

prop=F.func(Inf,Inf)
F_par=F_emp=numeric(n)
for(i in 1:n){
  F_par[i]=F.func(l[i],x[i])/prop
  F_emp[i]=mean( (l<=l[i])&(x<=x[i]) )
}
C.test=sum( (F_emp-F_par)^2 )
K.test=max( abs( F_emp-F_par ) )

plot(F_emp,F_par,xlab="F_empirical",ylab="F_parametric",xlim=c(0,1),ylim=c(0,1))
lines(x = c(0,1), y = c(0,1))

}

list(n=n, alpha = alpha_res,lambda_L = lambda_L_res,lambda_X = lambda_X_res,
     nu_L = nu_L_res,nu_X = nu_X_res,mu = mu_res,
     ML=LL,AIC=AIC_res,BIC=BIC_res,Iteration=AI,
     c=prop,C=C.test,K=K.test,F_empirical=F_emp,F_parametric=F_par)
}


