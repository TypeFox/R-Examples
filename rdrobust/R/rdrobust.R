### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014
### version 0.6  17Jun2014
### version 0.61 03Sep2014
### version 0.7  14Oct2014
### version 0.8  04Feb2015
### version 0.9  28Mar2016

rdrobust = function(y, x, covs = NULL, fuzzy=NULL, cluster=NULL, c=0, p=1, q=2, deriv=0, h=NULL, b=NULL, rho=NULL, 
                    scalepar=1, kernel="tri", bwselect="mserd", scaleregul=1, sharpbw=FALSE, vce="nn",  nnmatch=3, level=95, all=FALSE, subset = NULL) {
  
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
  
  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  if (!is.null(covs)){
    if (!is.null(subset))  covs <- covs[subset]
    na.ok <- na.ok & complete.cases(covs)
  } 
  
  if (!is.null(fuzzy)){
    if (!is.null(subset)) fuzzy <- fuzzy[subset]
    na.ok <- na.ok & complete.cases(fuzzy)
  } 

  x = matrix(x[na.ok])
  y = matrix(y[na.ok])
  
  if (!is.null(covs))    covs    = matrix(covs)[na.ok, , drop = FALSE]
  if (!is.null(fuzzy))   fuzzy   = matrix(  fuzzy[na.ok])
  if (!is.null(cluster)) cluster = matrix(cluster[na.ok])
  
  if (vce=="nn") {
    order_x = order(x)
    x = x[order_x,,drop=FALSE]
    y = y[order_x,,drop=FALSE]
    if (!is.null(covs))    covs    =  matrix(covs)[order_x,,drop=FALSE]
    if (!is.null(fuzzy))   fuzzy   =   fuzzy[order_x,,drop=FALSE]
    if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
  }

  kernel = tolower(kernel)
  #bwselect = toupper(bwselect)
  vce = tolower(vce)
  
  X_l=x[x<c,,drop=FALSE];  X_r=x[x>=c,,drop=FALSE]
  Y_l=y[x<c,,drop=FALSE];  Y_r=y[x>=c,,drop=FALSE]
  x_min = min(x);  x_max = max(x)
  N_l = length(X_l);   N_r = length(X_r)
  range_l = abs(max(X_l)-min(X_l));   range_r = abs(max(X_r)-min(X_r))
  N = N_r + N_l
  quant = -qnorm(abs((1-(level/100))/2))
  
  if (deriv>0 & p<=deriv) {
   p = deriv + 1
   q = p+1
  }


  #####################################################   CHECK ERRORS
  exit=0
    if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
      print("kernel incorrectly specified")
      exit = 1
    }
    
  if  (bwselect!="mserd" & bwselect!="msetwo" & bwselect!="msesum" & bwselect!="msecomb1" & bwselect!="msecomb2" & bwselect!="cerrd" & bwselect!="certwo" & bwselect!="cersum" & bwselect!="cercomb1" & bwselect!="cercomb2" & bwselect!=""){
    print("bwselect incorrectly specified")  
    exit = 1
  }
  
  if (bwselect=="CCT" | bwselect=="IK" | bwselect=="CV" | bwselect=="cct" | bwselect=="ik" | bwselect=="cv"){
    print("bwselect options IK, CCT and CV have been depricated. Please see help for new options")  
    exit = 1
  }
  
  if (vce!="nn" & vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){ 
    print("vce incorrectly specified")
    exit = 1
  }
    
    if (c<=x_min | c>=x_max){
      print("c should be set within the range of x")
      exit = 1
    }
    
    if (p<=0 | q<=0 | deriv<0 | nnmatch<=0 ){
      print("p,q,deriv and matches should be positive integers")
      exit = 1
    }
    
    if (p>=q & q>0){
      print("p should be set higher than q")
      exit = 1
    }
    
    if (deriv>p & deriv>0 ){
      print("deriv should be set higher than p")
      exit = 1
    }
    
    p_round = round(p)/p;  q_round = round(q)/q;  d_round = round(deriv+1)/(deriv+1);  m_round = round(nnmatch)/nnmatch
    
    if (p_round!=1 | q_round!=1 | d_round!=1 | m_round!=1 ){
      print("p,q,deriv and matches should be integer numbers")
      exit = 1
    }
    

    if (level>100 | level<=0){
      print("level should be set between 0 and 100")
      exit = 1
    }
    
    if (!is.null(rho)){  
       if (rho<0){
          print("rho should be greater than 0")
          exit = 1
        }
    }
  
    if (exit>0) stop()
    if (!is.null(h)) bwselect = "Manual"
    if (!is.null(h) & is.null(rho) & is.null(b)) {
      rho = 1
      b = h
    }
    if (!is.null(h) & !is.null(rho) ) b = h/rho
    
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
  }   else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
  }   else  {
    kernel_type = "Triangular"
  }

    ############################################################################################
    #print("Preparing data.") 
    
    if (is.null(h)) {
      rdbws=rdbwselect(y=y, x=x,  covs=covs, cluster=cluster, fuzzy=fuzzy, c=c, deriv=deriv, p=p, q=q, vce=vce, bwselect=bwselect, kernel=kernel, scaleregul=scaleregul, sharpbw=sharpbw)
      h_l = rdbws$bws[1]
      h_r = rdbws$bws[2]
      b_l = rdbws$bws[3]
      b_r = rdbws$bws[4]
      
      if (!is.null(rho)) {
        b_l = h_l/rho
    		b_r = h_r/rho
      }
      
    } else{
      if (length(h)==1) h_l = h_r = h
      if (length(h)==2) {
        h_l = h[1]
        h_r = h[2]
      }
      if (is.null(b)) {
        b_l = h_l
        b_r = h_r
      } else {
        if (length(b)==1) b_l = b_r = b
        if (length(b)==2) {
          b_l = b[1]
          b_r = b[2]
        }  
      }  
    }

  w_h_l = rdrobust_kweight(X_l,c,h_l,kernel);	w_h_r = rdrobust_kweight(X_r,c,h_r,kernel)
  w_b_l = rdrobust_kweight(X_l,c,b_l,kernel);	w_b_r = rdrobust_kweight(X_r,c,b_r,kernel)
  ind_h_l = w_h_l> 0;		ind_h_r = w_h_r> 0
  ind_b_l = w_b_l> 0;		ind_b_r = w_b_r> 0
  N_h_l = sum(ind_h_l);  N_b_l = sum(ind_b_l)
  N_h_r = sum(ind_h_r);	N_b_r = sum(ind_b_r)
  
  #if (N_h_l<5 | N_h_r<5 | N_b_l<5 | N_b_r<5){
  #  stop("Not enough observations to perform calculations")
  #  exit(1)
  #}
  
  ind_l = ind_b_l; ind_r = ind_b_r
  if (h_l>b_l) ind_l = ind_h_l   
  if (h_r>b_r) ind_r = ind_h_r   
  
  eN_l = sum(ind_l); eN_r = sum(ind_r)
  eY_l  = Y_l[ind_l,,drop=FALSE];	eY_r  = Y_r[ind_r,,drop=FALSE]
  eX_l  = X_l[ind_l,,drop=FALSE];	eX_r  = X_r[ind_r,,drop=FALSE]
  W_h_l = w_h_l[ind_l];	W_h_r = w_h_r[ind_r]
  W_b_l = w_b_l[ind_l];	W_b_r = w_b_r[ind_r]
  

  edups_l = edups_r = edupsid_l= edupsid_r = 0	
  if (vce=="nn") {
    for (i in 1:eN_l) {
      edups_l[i]=sum(eX_l==eX_l[i])
    }
    for (i in 1:eN_r) {
      edups_r[i]=sum(eX_r==eX_r[i])
    }
    for (i in 1:eN_l) {
      edupsid_l[i:(i+edups_l[i]-1)]=1:edups_l[i]
      i=i+edups_l[i]-1
    }
    for (i in 1:eN_r) {
      edupsid_r[i:(i+edups_r[i]-1)]=1:edups_r[i]
      i=i+edups_r[i]-1
    }
  }          
          
  u_l = (eX_l-c)/h_l;	u_r = (eX_r-c)/h_r
  R_q_l = matrix(NA,eN_l,(q+1)); R_q_r = matrix(NA,eN_r,(q+1))
  for (j in 1:(q+1))  {
    R_q_l[,j] = (eX_l-c)^(j-1);  R_q_r[,j] = (eX_r-c)^(j-1)
  }
  R_p_l = R_q_l[,1:(p+1)]; R_p_r = R_q_r[,1:(p+1)]

  #display("Computing RD estimates.")
  L_l = crossprod(R_p_l*W_h_l,u_l^(p+1)); L_r = crossprod(R_p_r*W_h_r,u_r^(p+1)) 
  invG_q_l  = qrXXinv((sqrt(W_b_l)*R_q_l));	invG_q_r  = qrXXinv((sqrt(W_b_r)*R_q_r))
  invG_p_l  = qrXXinv((sqrt(W_h_l)*R_p_l));	invG_p_r  = qrXXinv((sqrt(W_h_r)*R_p_r))
  e_p1 = matrix(0,(q+1),1); e_p1[p+2]=1
  e_v  = matrix(0,(p+1),1); e_v[deriv+1]=1
  Q_q_l = t(t(R_p_l*W_h_l) - h_l^(p+1)*(L_l%*%t(e_p1))%*%t(t(invG_q_l%*%t(R_q_l))*W_b_l))
  Q_q_r = t(t(R_p_r*W_h_r) - h_r^(p+1)*(L_r%*%t(e_p1))%*%t(t(invG_q_r%*%t(R_q_r))*W_b_r))
  D_l = eY_l; D_r = eY_r
  eC_l=eC_r=eT_l=eT_r=eZ_l=eZ_r=NULL
  dZ = dT =g_l=g_r= 0
  
  if (!is.null(fuzzy)) {
    dT=1
    T_l  = fuzzy[x<c,,drop=FALSE];  eT_l  = T_l[ind_l,,drop=FALSE]
    T_r  = fuzzy[x>=c,,drop=FALSE]; eT_r  = T_r[ind_r,,drop=FALSE]
    D_l  = cbind(D_l,eT_l); D_r = cbind(D_r,eT_r)
  }
  
  if (!is.null(covs)) {
    dZ = ncol(covs)
    Z_l  = covs[x<c,,drop=FALSE];	  eZ_l = Z_l[ind_l,,drop=FALSE]
    Z_r  = covs[x>=c,,drop=FALSE];	eZ_r = Z_r[ind_r,,drop=FALSE]
    D_l  = cbind(D_l,eZ_l); D_r = cbind(D_r,eZ_r)
    U_p_l = crossprod(R_p_l*W_h_l,D_l); U_p_r = crossprod(R_p_r*W_h_r,D_r)
  }
              
  if (!is.null(cluster)) {
    C_l  = cluster[x<c,,drop=FALSE]; C_r= cluster[x>=c,,drop=FALSE]
    eC_l  = C_l[ind_l];	     eC_r  = C_r[ind_r]
    g_l = length(unique(eC_l));	g_r = length(unique(eC_r))
  }
                                                     
  beta_p_l = invG_p_l%*%crossprod(R_p_l*W_h_l,D_l); beta_q_l = invG_q_l%*%crossprod(R_q_l*W_b_l,D_l); beta_bc_l = invG_p_l%*%crossprod(Q_q_l,D_l) 
  beta_p_r = invG_p_r%*%crossprod(R_p_r*W_h_r,D_r); beta_q_r = invG_q_r%*%crossprod(R_q_r*W_b_r,D_r); beta_bc_r = invG_p_r%*%crossprod(Q_q_r,D_r)
  beta_p  = beta_p_r  - beta_p_l
  beta_q  = beta_q_r  - beta_q_l
  beta_bc = beta_bc_r - beta_bc_l

  if (is.null(covs)) {	
  tau_cl = tau_Y_cl = scalepar*factorial(deriv)*beta_p[(deriv+1),1]
  tau_bc = tau_Y_bc = scalepar*factorial(deriv)*beta_bc[(deriv+1),1]
  s_Y = 1
  if (!is.null(fuzzy)) {
     tau_T_cl =  factorial(deriv)*beta_p[(deriv+1),2]
     tau_T_bc = 	factorial(deriv)*beta_bc[(deriv+1),2]
     tau_cl = tau_Y_cl/tau_T_cl
     s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2))
     B_F = c(tau_Y_cl-tau_Y_bc , tau_T_cl-tau_T_bc)
     tau_bc = tau_cl - t(s_Y)%*%B_F
     sV_T = c(0 , 1)
  }	
  } else {	
    ZWD_p_l  = crossprod(eZ_l*W_h_l,D_l)
    ZWD_p_r  = crossprod(eZ_r*W_h_r,D_r)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU_p_l =  crossprod(U_p_l[,colsZ],invG_p_l%*%U_p_l) 
    UiGU_p_r =  crossprod(U_p_r[,colsZ],invG_p_r%*%U_p_r) 
    ZWZ_p_l = ZWD_p_l[,colsZ] - UiGU_p_l[,colsZ] 
    ZWZ_p_r = ZWD_p_r[,colsZ] - UiGU_p_r[,colsZ]     
    ZWY_p_l = ZWD_p_l[,1:(1+dT)] - UiGU_p_l[,1:(1+dT)] 
    ZWY_p_r = ZWD_p_r[,1:(1+dT)] - UiGU_p_r[,1:(1+dT)]     
    ZWZ_p = ZWZ_p_r + ZWZ_p_l
    ZWY_p = ZWY_p_r + ZWY_p_l
    gamma_p = chol2inv(chol(ZWZ_p))%*%ZWY_p
    s_Y = c(1 ,  -gamma_p[,1])
    if (is.null(fuzzy)) {
        tau_cl = t(s_Y)%*%beta_p[(deriv+1),]
        tau_bc = t(s_Y)%*%beta_bc[(deriv+1),]
    } else {
      s_T  = c(1,    -gamma_p[,2])
      sV_T = c(0, 1, -gamma_p[,2])
      tau_Y_cl = factorial(deriv)*t(s_Y)%*%c(beta_p[(deriv+1),1], beta_p[(deriv+1),colsZ])
      tau_T_cl = factorial(deriv)*t(s_T)%*%c(beta_p[(deriv+1),2], beta_p[(deriv+1),colsZ])
      tau_Y_bc = factorial(deriv)*t(s_Y)%*%c(beta_bc[(deriv+1),1],beta_bc[(deriv+1),colsZ])
      tau_T_bc = factorial(deriv)*t(s_T)%*%c(beta_bc[(deriv+1),2],beta_bc[(deriv+1),colsZ])
      tau_cl = tau_Y_cl/tau_T_cl
      B_F = c(tau_Y_cl-tau_Y_bc , tau_T_cl-tau_T_bc)
      s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2))
      tau_bc = tau_cl - t(s_Y)%*%B_F
      s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2) , -(1/tau_T_cl)*gamma_p[,1] + (tau_Y_cl/tau_T_cl^2)*gamma_p[,2])
    }
  }

#*display("Computing variance-covariance matrix.")
  
  hii_l=hii_r=predicts_p_l=predicts_p_r=predicts_q_l=predicts_q_r=0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_p_l=R_p_l%*%beta_p_l
    predicts_p_r=R_p_r%*%beta_p_r
    predicts_q_l=R_q_l%*%beta_q_l
    predicts_q_r=R_q_r%*%beta_q_r
    if (vce=="hc2" | vce=="hc3") {
      hii_l=matrix(NA,eN_l,1)	
      for (i in 1:eN_l) {
        hii_l[i] = R_p_l[i,]%*%invG_p_l%*%(R_p_l*W_h_l)[i,]
      }
      hii_r=matrix(NA,eN_r,1)	
      for (i in 1:eN_r) {  
        hii_r[i] = R_p_r[i,]%*%invG_p_r%*%(R_p_r*W_h_r)[i,]
      }
    }
  }
  						
	res_h_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_p_l, hii_l, vce, nnmatch, edups_l, edupsid_l, p+1)
	res_h_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_p_r, hii_r, vce, nnmatch, edups_r, edupsid_r, p+1)
	if (vce=="nn") {
			res_b_l = res_h_l;	res_b_r = res_h_r
	} 	else {
			res_b_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_q_l, hii_l, vce, nnmatch, edups_l, edupsid_l, q+1)
			res_b_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_q_r, hii_r, vce, nnmatch, edups_r, edupsid_r, q+1)
  }
			                       
	V_Y_cl_l = invG_p_l%*%rdrobust_vce(dT+dZ, s_Y, R_p_l*W_h_l, res_h_l, eC_l)%*%invG_p_l
	V_Y_cl_r = invG_p_r%*%rdrobust_vce(dT+dZ, s_Y, R_p_r*W_h_r, res_h_r, eC_r)%*%invG_p_r
	V_Y_bc_l = invG_p_l%*%rdrobust_vce(dT+dZ, s_Y, Q_q_l,       res_b_l, eC_l)%*%invG_p_l
	V_Y_bc_r = invG_p_r%*%rdrobust_vce(dT+dZ, s_Y, Q_q_r,       res_b_r, eC_r)%*%invG_p_r
	V_tau_cl = factorial(deriv)^2*(V_Y_cl_l+V_Y_cl_r)[deriv+1,deriv+1]
	V_tau_rb = factorial(deriv)^2*(V_Y_bc_l+V_Y_bc_r)[deriv+1,deriv+1]
	se_tau_cl = scalepar*sqrt(V_tau_cl);	se_tau_rb = scalepar*sqrt(V_tau_rb)

	if (!is.null(fuzzy)) {
		V_T_cl_l = invG_p_l%*%rdrobust_vce(dT+dZ, sV_T, R_p_l*W_h_l, res_h_l, eC_l)%*%invG_p_l
		V_T_cl_r = invG_p_r%*%rdrobust_vce(dT+dZ, sV_T, R_p_r*W_h_r, res_h_r, eC_r)%*%invG_p_r
		V_T_bc_l = invG_p_l%*%rdrobust_vce(dT+dZ, sV_T, Q_q_l, res_b_l, eC_l)%*%invG_p_l
		V_T_bc_r = invG_p_r%*%rdrobust_vce(dT+dZ, sV_T, Q_q_r, res_b_r, eC_r)%*%invG_p_r
		V_T_cl = factorial(deriv)^2*(V_T_cl_l+V_T_cl_r)[deriv+1,deriv+1]
		V_T_rb = factorial(deriv)^2*(V_T_bc_l+V_T_bc_r)[deriv+1,deriv+1]
		se_tau_T_cl = sqrt(V_T_cl);	se_tau_T_rb = sqrt(V_T_rb)
	}
  
  tau = c(tau_cl, tau_bc, tau_bc)
  se  = c(se_tau_cl,se_tau_cl,se_tau_rb)
  t  =  tau/se
  pv = 2*pnorm(-abs(t))
  ci = matrix(NA,nrow=3,ncol=2)
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("Lower","Upper")
  ci[1,] = c(tau[1] - quant*se[1], tau[1] + quant*se[1])
  ci[2,] = c(tau[2] - quant*se[2], tau[2] + quant*se[2])
  ci[3,] = c(tau[3] - quant*se[3], tau[3] + quant*se[3])
    
  if (!is.null(fuzzy)) {  
      tau_T = c(tau_T_cl, tau_T_bc, tau_T_bc)
      se_T = c(se_tau_T_cl, se_tau_T_cl, se_tau_T_rb)
      t_T = tau_T/se_T
      pv_T = 2*pnorm(-abs(t_T))
      ci_T=matrix(NA,nrow=3,ncol=2)
      ci_T[1,] = c(tau_T[1] - quant*se_T[1], tau_T[1] + quant*se_T[1])
      ci_T[2,] = c(tau_T[2] - quant*se_T[2], tau_T[2] + quant*se_T[2])
      ci_T[3,] = c(tau_T[3] - quant*se_T[3], tau_T[3] + quant*se_T[3])
  }

  #print("Estimation Completed.") 

    coef=matrix(tau,3,1)
    se  =matrix(se,3,1)
    z   =matrix(t,3,1)
    pv  =matrix(pv,3,1)
    ci=ci
  


  bws=matrix(c(h_l,b_l,h_r,b_r),2,2)
  rownames(coef)=rownames(se)=rownames(se)=rownames(z)=rownames(pv)=c("Conventional","Bias-Corrected","Robust")
  colnames(coef)="Coeff"
  colnames(se)="Std. Err."
  colnames(z)="z"
  colnames(pv)="P>|z|"
  colnames(bws)=c("left","right")
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("CI Lower","CI Upper")
    
  tabl1.str=matrix(NA,4,1)
  rownames(tabl1.str)=c("Number of Obs", "NN Matches", "BW Type", "Kernel Type")
  dimnames(tabl1.str) <-list(c("Number of Obs", "NN Matches", "BW Type", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=N
  tabl1.str[2,1]=nnmatch
  tabl1.str[3,1]=bwselect
  tabl1.str[4,1]=kernel_type
  
  tabl2.str=matrix(NA,6,2)
  colnames(tabl2.str)=c("Left","Right")
  rownames(tabl2.str)=c("Number of Obs","Order Loc Poly (p)","Order Bias (q)","BW Loc Poly (h)","BW Bias (b)","rho (h/b)")
  tabl2.str[1,]=formatC(c(N_l,N_r),digits=0, format="f")
  tabl2.str[2,]=formatC(c(p,p),digits=0, format="f")
  tabl2.str[3,]=formatC(c(q,q),digits=0, format="f")
  tabl2.str[4,]=formatC(c(h_l,h_r),digits=4, format="f")
  tabl2.str[5,]=formatC(c(b_l,b_r),digits=4, format="f")
  tabl2.str[6,]=formatC(c(h_l/b_l,h_r/b_r),digits=4, format="f")
  
  tabl3.str=matrix("",2,6)
  colnames(tabl3.str)=c("Coef","Std. Err.","z","P>|z|","CI Lower","CI Upper")
  rownames(tabl3.str)=c("Conventional", "Robust")
  tabl3.str[1,1]  =formatC(coef[1],digits=4, format="f")
  tabl3.str[1,2]  =formatC(se[1],  digits=4, format="f")
  tabl3.str[1,3]  =formatC(z[1],   digits=4, format="f")
  tabl3.str[1,4]  =formatC(pv[1],  digits=4, format="f")
  tabl3.str[1,5:6]=formatC(ci[1,], digits=4, format="f")
  tabl3.str[2,4]  =formatC(pv[3],  digits=4, format="f")
  tabl3.str[2,5:6]=formatC(ci[3,] ,digits=4, format="f")
  
  if (all==TRUE){
    tabl3.str=formatC(cbind(coef,se,z,pv,ci),digits=4, format="f")                   
    colnames(tabl3.str)=c("Coef","Std. Err.","z","P>|z|","CI Lower","CI Upper")
  }

  out=list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,tabl3.str=tabl3.str,coef=coef,bws=bws,se=se,z=z,pv=pv,ci=ci,p=p,q=q,h=h,b=b,rho=rho,N=N,N_l=eN_l,N_r=eN_r)
  out$call <- match.call()
  class(out) <- "rdrobust"
  return(out)
}

#rdrobust <- function(y,x, ...) UseMethod("rdrobust")

#rdrobust.default <- function(y,x,  ...){
#  est <- rdrobustEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "rdrobust"
#  est
#}

print.rdrobust <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nSummary:\n")
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F)
  cat("\nEstimates:\n")
  print(x$tabl3.str,quote=F)
}

summary.rdrobust <- function(object,...) {
  TAB <- cbind(Estimate    =object$coef,
               "Std. Error"=object$se,
               "z"         =object$z,
               "Pr(>|z|)"  =object$pv,
               "95% CI"    =object$ci)
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.rdrobust"
  res
}

#print.summary.rdrobust <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coef)
#  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
#}