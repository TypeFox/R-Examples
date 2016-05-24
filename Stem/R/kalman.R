`kalman` <-
function (z, coordinates, p, n, d, r, phi_j, max.iter, precision, covariates, Gdiag, Sigmaetadiag, cov.spat) {


zz   = ts(z)
dist = as.matrix(dist(coordinates,diag=TRUE)) #distance matrix

####################
###Model definition 
###if you want to check the model use phi_j=phi_start
####################
SSmodel  = list(z	= zz,
		Fmat 	= phi_j$K,
		Gmat 	= phi_j$G,
		Vmat 	= phi_j$sigma2omega * cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist),
		Wmat 	= phi_j$Sigmaeta,
		m0   		= t(phi_j$m0),
		C0   		= phi_j$C0,
		phi  		= phi_j,
		XXX  		= covariates,
		beta 	= phi_j$beta,
		flag.cov 	= TRUE,
		n		= n,
		p		= p,
		d		= d,
		m 		= NA,
		C		= NA,
		loglik		= NA
		)

####################
###kalman filtering and smoothing
####################
mod1.filter   	= filtering(SSmodel) 
mod1.smoother 	= smoothing(mod1.filter) 


####calculating  B0 e P_0_n=C_0* using kalman filter, smoother and initial values
R1    = SSmodel$Gmat %*% SSmodel$C0 %*% t(SSmodel$Gmat) + SSmodel$Wmat  #P_1^0
B0    = SSmodel$C0 %*% SSmodel$Gmat %*% 	solve(R1)
P_0_n = SSmodel$C0 + B0 %*% (mod1.smoother$C[[1]] - R1 ) %*% t(B0)

########################
###Define the elements for B_function (output of mod1$filter)
###Use B_function (see functions.R) ---> list of B_t (t=1,...,n)
########################
nobs 	= n
m    	= mod1.filter$m
C    	= mod1.filter$C
B    	= list()

for (tt in (nobs-1):1) {            
      nextstep = B_function(
			m     = matrix(m[tt,],nrow=1),
			C     = C[[tt]],
			Gmatx = SSmodel$Gmat,
			Wtx   = SSmodel$Wmat,
			mx    = matrix(m[tt+1,],nrow=1),
			Cx    = C[[tt+1]])

	m[tt,]  = nextstep$ms #equal to mod1.smoother$m
	C[[tt]] = nextstep$Cs #equal to mod1.smoother$C
	B[[tt]] = nextstep$B  #what I need
}

#B_n=C_n * t(Gmat) * solve(R_n+1)
#where  R_n+1=Gmat * C_n * t(Gmat) +Wmat

B[[nobs]] = mod1.filter$C[[nobs]] %*% t(SSmodel$Gmat) %*%
			solve(SSmodel$Gmat %*% mod1.filter$C[[nobs]] %*%
			t(SSmodel$Gmat)+SSmodel$Wmat)

################################
###Ricursion for LAG ONE COVARIANCE SMOOTHER
################################
CCC = list() #empty list for lag one covariance values

###FIRST STEP: Calculate C*_{n,n-1} (start of the iterative procedure --> CCC[[n]])
###Pn_n-1=R_n=G_n%*%C_n-1%*%t(G_n)+W_n
###A_n=R_n%*%F_n%*%solve(t(F_n)%*%R_n%*%F_n+V_n)
R_n 			= SSmodel$Gmat %*% mod1.filter$C[[nobs-1]] %*%t(SSmodel$Gmat)+ SSmodel$Wmat
A_n 			= R_n%*%SSmodel$Fmat %*%
			solve(t(SSmodel$Fmat) %*% R_n %*% SSmodel$Fmat + SSmodel$Vmat)
CCC[[nobs]] 	= (diag(p)-A_n %*% t(SSmodel$Fmat)) %*% SSmodel$Gmat %*% mod1.filter$C[[nobs-1]]

################################
###loop for calculating the other elements of CCC using cov_lagone fucntion (see function.R)
################################
for (tt in (nobs):2) {
	passiCCC =cov_lagone(
		C_t_minus_1 	= mod1.filter$C[[tt-1]],
		B_t_minus_1 	= B[[tt-1]],
		Gmat        		= SSmodel$Gmat,
		if (tt>=3)  
			B_t_minus_2= B[[tt-2]] else B_t_minus_2 = B0,
		CCCx       		 = CCC[[tt]])

    	CCC[[tt-1]] = passiCCC$cov
}

################################################################
#Parameters estimates and calculate of the elements of Q function
################################################################
###PARAMETER 1: m0=y_0_n
y_0_n = t(SSmodel$m0) + B0 %*% (t(matrix(mod1.smoother$m[1,],nrow=1))- SSmodel$Gmat %*% t(SSmodel$m0))
m0_j  = y_0_n

####Q element n.2 (NB: y0_n=m0)
Q_addendo2 = Q_function_addendo2(C0=phi_j$C0, m_0=m0_j, P_0_n=P_0_n, y_0_n=y_0_n)

###S11, S10 and S00
S11list = list()
for (tt in 1:nobs) {
    S11list[[tt]] = t(matrix(mod1.smoother$m[tt,],nrow=1)) %*% matrix(mod1.smoother$m[tt,],nrow=1) +  mod1.smoother$C[[tt]]
}
S11 = sumMatrices(S11list)

S00list = list()
S00list[[1]] = m0_j %*% t(m0_j) + P_0_n
for (tt in 2:nobs) {
    S00list[[tt]] = t(matrix(mod1.smoother$m[tt-1,],nrow=1)) %*% matrix(mod1.smoother$m[tt-1,],nrow=1) +  mod1.smoother$C[[tt-1]]
}
S00 = sumMatrices(S00list)

S10list = list()
S10list[[1]] = t(matrix(mod1.smoother$m[1,],nrow=1)) %*% t(m0_j) + CCC[[1]]
for (tt in 2:nobs) {
    S10list[[tt]] = t(matrix(mod1.smoother$m[tt,],nrow=1)) %*% matrix(mod1.smoother$m[tt-1,],nrow=1) +  CCC[[tt]]
}
S10 = sumMatrices(S10list)

################################
###PARAMETER 2: Sigmaeta
#Sigmaeta_j = (S11 - S10 %*% t(phi_j$G) - phi_j$G %*% t(S10) + phi_j$G %*% S00 %*% t(phi_j$G))/nobs
Sigmaeta_j = (S11 - S10 %*% solve(S00) %*% t(S10))/nobs
if(Gdiag & Sigmaetadiag & p>1) Sigmaeta_j = diag(diag(Sigmaeta_j))

if (Sigmaetadiag) {
   Sigmaeta_j = matrix(0,p,p)
   num = (S11 - S10 %*% t(phi_j$G) - phi_j$G %*% t(S10) + phi_j$G %*% S00 %*% t(phi_j$G))
   for (i in 1:p) { Sigmaeta_j[i,i] = num[i,i] / nobs }
}

################################
###PARAMETER 3: G
G_j = S10 %*% solve(S00)
if(Gdiag & Sigmaetadiag & p>1)  G_j = diag(diag(G_j))

if (Gdiag) {
   G_j = matrix(0,p,p)
   num = solve(Sigmaeta_j) %*% S10
   den = solve(Sigmaeta_j) %*% S00
   for (i in 1:p) { G_j[i,i] = num[i,i] / den[i,i]  }
}


###Q element n.3
Q_addendo3 = Q_function_addendo3(Sigmaeta=Sigmaeta_j, G=G_j, n=n, S11=S11, S00=S00, S10=S10)

#############################
###PARAMETER 4: sigma2omega
###sigma2omega=sigma2epsilon if logb=0 (exp(logb)=1)
BB_list = list()
for (tt in 1:nobs) {
	BB_list[[tt]] = 	t(SSmodel$Fmat) %*% mod1.smoother$C[[tt]] %*% SSmodel$Fmat +
				((t(matrix(zz[tt,],nrow=1)) - covariates[,,tt]%*%phi_j$beta - t(SSmodel$Fmat) %*% t(matrix(mod1.smoother$m[tt,],nrow=1)))
				%*% t(t(matrix(zz[tt,],nrow=1)) -covariates[,,tt]%*%phi_j$beta- t(SSmodel$Fmat) %*% t(matrix(mod1.smoother$m[tt,],nrow=1)))   )
}
BB = sumMatrices(BB_list)

D = solve(cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist)) %*% BB
sigma2omega_j = sum(diag(D))/(n*d)

#############################
###PARAMETER 5: beta coefficients
############################
#\sum_t X_t^\prime \Sigma_e^-1 v_t
#v_t=z_t-K_t y_t
Sigmae_inversa = solve(sigma2omega_j * cov.spat(d=d , logb=phi_j$logb , logtheta=phi_j$logtheta , dist=dist))

vvt_list = list ()
for (tt in 1:nobs) {
    vt 		=  matrix(unlist(z[tt,]),ncol=1) - t(SSmodel$Fmat) %*% mod1.smoother$m[tt,]
    vvt_list[[tt]] 	= t(covariates[,,tt]) %*%  Sigmae_inversa  %*% vt
}
v = sumMatrices(vvt_list)


MM_list = list()
for (tt in 1:nobs) {
    MM_list[[tt]] 	= t(covariates[,,tt]) %*% Sigmae_inversa %*% covariates[,,tt]
}
MM = sumMatrices(MM_list)

if(det(MM) != 0) {beta_j = solve(MM) %*% v}
if(det(MM)  < 10^(-7)) {cat("Error in beta estimation! The matrix can not be inverted!!!!")}

#############################
###PARAMETER 6: theta e logb
##exponential spatial covariance function
##Newton raphson algorithm
#############################
n_iter_NR = 1
n_iter_Hess.list = c()
convergence_NR = FALSE

logb     = phi_j$logb
logtheta = phi_j$logtheta

while(!convergence_NR && n_iter_NR < 50) {

        Q_addendo1_old = Q_function_addendo1(sigma2omega=sigma2omega_j ,n=n, Sigmastar=do.call(cov.spat, list(d=d, logb=logb, logtheta=logtheta, dist=dist)),B=BB)
        Q_prev = Q_addendo1_old+Q_addendo2+Q_addendo3

        cond.hessiana = FALSE
        n_iter_Hess   = 1
        logtheta.iniz = logtheta

        while(!cond.hessiana && n_iter_Hess < 30) {
        cov.spat.mat = do.call(cov.spat, list(logb=logb,d=d,logtheta=logtheta,dist=dist))
	
		derivata_prima_logtheta 	= d1_Q(
								n	= n,
								X	= cov.spat.mat,
								d1_X  = d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
								sigma2omega = sigma2omega_j,
								B       = BB)
		derivata_seconda_logtheta 	= d2_Q(
								n	= n,
								X      = cov.spat.mat,
								d1_X	= d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
								d2_X = d2_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
								sigma2omega = sigma2omega_j,
								B      = BB)
		
		derivata_prima_logb    		= d1_Q(
								n	 = n,
								X       = cov.spat.mat,
								d1_X  = d1_Sigmastar_logb.exp(logb=logb,d=d),
								sigma2omega = sigma2omega_j,
								B       = BB)

		derivata_seconda_logb    	= d2_Q(
								n	= n,
								X 	= cov.spat.mat,
								d1_X	= d1_Sigmastar_logb.exp(logb=logb,d=d),
								d2_X	= d2_Sigmastar_logb.exp(logb=logb,d=d),
								sigma2omega = sigma2omega_j,
								B	= BB)
	
		derivata_mista        		=  d12_Q(
								n	= n,
								X	= cov.spat.mat,
								d1_X_theta  = d1_Sigmastar_logtheta.exp(logtheta=logtheta,dist=dist),
								d1_X_logb   = d1_Sigmastar_logb.exp(logb=logb,d=d),
								sigma2omega = sigma2omega_j,
								B	= BB)


		hessiana  = matrix( c(derivata_seconda_logtheta, derivata_mista, derivata_mista, derivata_seconda_logb),2,2)
		cond.hessiana = det(hessiana) > 10^(-3)

		if(!cond.hessiana) {
			kk1=10
			kk2=10
			logtheta_vec = log(seq((0.01*exp(logtheta)),(10*exp(logtheta)),length=kk1))
			logb_vec  = log(seq((0.01*exp(logb)),(10*exp(logb)),length=kk2))
			QQ=matrix(NA,kk1,kk2)

			for(i in 1:kk1){
				for(j in 1:kk2){
					QQ[i,j] = Q_function_addendo1(sigma2omega=sigma2omega_j , n=n ,
					Sigmastar=cov.spat(d=d , logb=logb_vec[j] , logtheta=logtheta_vec[i], dist=dist),B=BB)
				}
			}
               
               val.col.min = apply(QQ,2,min)
               colonna = which.min(val.col.min)
               righe = apply(QQ,2,which.min)
               riga = righe[colonna]

               logtheta = logtheta_vec[riga]
               logb     = logb_vec[colonna]
		}

	n_iter_Hess = n_iter_Hess + 1
        }

        Q_addendo1_prev =	 Q_function_addendo1(sigma2omega=sigma2omega_j, n=n,
					Sigmastar=cov.spat(d=d, logb=logb, logtheta=logtheta, dist=dist),B=BB)

	gradiente = matrix( c(derivata_prima_logtheta, derivata_prima_logb),2,1)

        logtheta.logb_old = matrix(c(logtheta , logb),2,1)
        delta = solve(hessiana) %*% gradiente
        logtheta.logb_new = logtheta.logb_old - delta

        dist_rel_num = sqrt(t(logtheta.logb_new - logtheta.logb_old) %*% (logtheta.logb_new - logtheta.logb_old))
        dist_rel_den = sqrt(t(logtheta.logb_old) %*% logtheta.logb_old)
        dist_rel = dist_rel_num / dist_rel_den

	convergence_NR = dist_rel < precision
        logb  = logtheta.logb_new[2,]
	logtheta = logtheta.logb_new[1,]

        n_iter_Hess.list[[n_iter_NR]] = n_iter_Hess - 1
        n_iter_NR = n_iter_NR + 1
        cat(paste("***NR Algorithm - iteration n.",n_iter_NR-1),"\n")


	} #end while loop


logtheta_j 	= logtheta
logb_j     	= logb


Q_addendo1 = Q_function_addendo1(sigma2omega=sigma2omega_j, n=n, Sigmastar=cov.spat(d=d, logb=logb_j, logtheta=logtheta_j, dist=dist),B=BB)
Q_new      	 = Q_addendo1 + Q_addendo2 + Q_addendo3


###########################
###New parameter vector (k e C0 don't change)
############################

phi_jj = list(
                 loglik	        	= mod1.filter$loglik,
                 K      	    	= phi_j$K,
                 sigma2omega	= sigma2omega_j,
                 logtheta      	= logtheta_j,
                 logb       		= logb_j,
                 beta       		= beta_j,
                 G          		= G_j,
                 Sigmaeta   	= Sigmaeta_j,
                 m0         		= m0_j,
                 C0         		= phi_j$C0)

return(list(phi    = phi_jj,
	        Q_prev = Q_prev,
                Q_new  = Q_new,
                 n_iter_NR =  n_iter_NR - 1,
                 m.smoother = mod1.smoother$m))
                 #m.filter   = mod1.filter$m,
                 #c.smoother = mod1.smoother$C,
                 #c.filter   = mod1.filter$C))              
}

