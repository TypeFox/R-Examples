### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014
### version 0.6  17Jun2014
### version 0.61 03Sep2014
### version 0.7  14Oct2014
### version 0.8  04Feb2015
### version 0.9  28Mar2016

rdbwselect = function(y, x, covs = NULL, fuzzy = NULL, cluster = NULL, c=0, p=1, q=2, deriv=0, 
                      kernel="tri", bwselect="mserd", scaleregul=1, sharpbw=FALSE, vce="nn",  nnmatch=3, all=FALSE, subset = NULL){
  
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
  if (!is.null(fuzzy))   fuzzy   = matrix(fuzzy[na.ok])
  if (!is.null(cluster)) cluster = matrix(cluster[na.ok])
  
  if (vce=="nn") {
    order_x = order(x)
    x = x[order_x,,drop=FALSE]
    y = y[order_x,,drop=FALSE]
    if (!is.null(covs))    covs    =  matrix(covs)[order_x,,drop=FALSE]
    if (!is.null(fuzzy))   fuzzy   =   fuzzy[order_x,,drop=FALSE]
    if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
  }
  
  ### reescaling
  y_sd = sd(y)
	y = y/y_sd
	x_sd = sd(x)
  x = x/x_sd
	c_orig = c
  c = c/x_sd
	
  #x_sd = sd(x)
  x_iq = quantile(x,.75) - quantile(x,.25)
  
  X_l = x[x<c];   
  x_l_min = min(X_l)
  x_l_max = max(X_l)
  range_l = abs(c-x_l_min)
  X_r = x[x>=c]
  x_r_min = min(X_r)
  x_r_max = max(X_r)
  range_r = abs(c-x_r_max)

  Y_l = y[x<c];    Y_r = y[x>=c]
  N_l = length(X_l);   N_r = length(X_r)
  x_min=min(x);  x_max=max(x)
  N = N_r + N_l
    
  if (deriv>0 & p==deriv) {
    p = deriv + 1
    q = p+1
  }
  
    exit=0
    #################  ERRORS
    if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
      print("kernel incorrectly specified")
      exit = 1
    }
    
    if  (bwselect!="mserd" & bwselect!="msetwo" & bwselect!="msesum" & bwselect!="msecomb1" & bwselect!="msecomb2"  & bwselect!="cerrd" & bwselect!="certwo" & bwselect!="cersum" & bwselect!="cercomb1" & bwselect!="cercomb2" & bwselect!=""){
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
      print("p, q, deriv and matches should be positive integers")
      exit = 1
    }
    
    if (p>=q){
      print("q should be set higher than p")
      exit = 1
    }
    
    if (deriv>=p){
      print("deriv should be set higher than p")
      exit = 1
    }
    
    p_round = round(p)/p;    q_round = round(q)/q;    d_round = round(deriv+1)/(deriv+1);    m_round = round(nnmatch)/nnmatch
        
    if (p_round!=1 | q_round!=1 | d_round!=1 | m_round!=1 ){
      print("p,q,deriv and matches should be integer numbers")
      exit = 1
    }
    
    if (exit>0) stop()
    
  
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
    C_c=2.34
  }  else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
    C_c=1.843
  }   else  {
    kernel_type = "Triangular"
    C_c=2.576
  }
  


  
  #***********************************************************************

  dZ=Z_l=Z_r=T_l=T_r=C_l=C_r=Cind_l=Cind_r=g_l=g_r=dups_l=dups_r=dupsid_l=dupsid_r=NULL
    
  if (vce=="nn") {
      for (i in 1:N_l) {
        dups_l[i]=sum(X_l==X_l[i])
      }
      for (i in 1:N_r) {
        dups_r[i]=sum(X_r==X_r[i])
      }
      
      for (i in 1:N_l) {
        dupsid_l[i:(i+dups_l[i]-1)]=1:dups_l[i]
        i=i+dups_l[i]-1
      }
      for (i in 1:N_r) {
        dupsid_r[i:(i+dups_r[i]-1)]=1:dups_r[i]
        i=i+dups_r[i]-1
      }
  }
  
  
  
  if (!is.null(covs)) {
    dZ = ncol(covs)
    Z_l  = covs[x<c,,drop=FALSE];  Z_r  = covs[x>=c,,drop=FALSE]

  }
  perf_comp=FALSE
  if (!is.null(fuzzy)) {
    T_l  = fuzzy[x<c,,drop=FALSE];  T_r  = fuzzy[x>=c,,drop=FALSE]; 
    if (var(T_l)==0 | var(T_r)==0) perf_comp=TRUE
    if (perf_comp==TRUE | sharpbw==TRUE) {
      T_l = T_r = NULL
      }
    }
   
  if (!is.null(cluster)) {
    C_l  = cluster[x<c,,drop=FALSE]; C_r= cluster[x>=c,,drop=FALSE]
    g_l = length(unique(C_l));	g_r = length(unique(C_r))
  }
                                                                              
    #***********************************************************************
    c_bw = C_c*min(c(1,x_iq/1.349))*N^(-1/5)
  
    #*** Step 1: d_bw
    C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range_l, 0, vce, nnmatch, kernel, dups_l, dupsid_l)
    C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range_r, 0, vce, nnmatch, kernel, dups_r, dupsid_r)
  
    #if (C_d_l$V=="NaN" | C_d_l$B=="NaN" | C_d_l$R=="NaN") C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=c_bw, 0, vce, nnmatch, kernel, dups_l, dupsid_l)
    #if (C_d_r$V=="NaN" | C_d_r$B=="NaN" | C_d_r$R=="NaN") C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=c_bw, 0, vce, nnmatch, kernel, dups_r, dupsid_r)
    #if (C_d_l$V==. | C_d_l$B==. | C_d_l$R==.) display("Invertibility problem in the computation of preliminary bandwidth below the threshold")  
    #if (C_d_r$V==. | C_d_r$B==. | C_d_r$R==.) display("Invertibility problem in the computation of preliminary bandwidth above the threshold")  
    #if (C_d_l$V==0 | C_d_l$B==0) display("Not enough variability to compute the preliminary bandwidth below the threshold. Range defined by bandwidth: ")  
    #if (C_d_r$V==0 | C_d_r$B==0) display("Not enough variability to compute the preliminary bandwidth above the threshold. Range defined by bandwidth: ")  
    
    
    #*** TWO
    if  (bwselect=="msetwo" |  bwselect=="certwo" | bwselect=="msecomb2" | bwselect=="cercomb2"  | all=="TRUE")  {		
      d_bw_l = (  C_d_l$V              /   C_d_l$B^2             )^C_d_l$rate
      d_bw_r = (  C_d_r$V              /   C_d_r$B^2             )^C_d_l$rate
      C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
      b_bw_l = (  C_b_l$V              /   (C_b_l$B^2 + scaleregul*C_b_l$R)        )^C_b_l$rate
      C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
      b_bw_r = (  C_b_r$V              /   (C_b_r$B^2 + scaleregul*C_b_r$R)        )^C_b_l$rate
      C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
      h_bw_l = (  C_h_l$V              /   (C_h_l$B^2 + scaleregul*C_h_l$R)         )^C_h_l$rate
      C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
      h_bw_r = (  C_h_r$V              /   (C_h_r$B^2 + scaleregul*C_h_r$R)         )^C_h_l$rate
                  
      #if (C_b_l$V==0 | C_b_l$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) below the threshold. Range defined by bandwidth = %f\n", d_bw_l)  
      #if (C_b_r$V==0 | C_b_r$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) above the threshold. Range defined by bandwidth = %f\n", d_bw_r)  
      #if (C_h_l$V==0 | C_h_l$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) below the threshold. Range defined by bandwidth = %f\n", b_bw_l) 
      #if (C_h_r$V==0 | C_h_r$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) above the threshold. Range defined by bandwidth = %f\n", b_bw_r) 
    }
  
#  *** SUM
  if  (bwselect=="msesum" | bwselect=="cersum" |  bwselect=="msecomb1" | bwselect=="msecomb2" |  bwselect=="cercomb1" | bwselect=="cercomb2"  |  all=="TRUE")  {
    d_bw_s = ( (C_d_l$V + C_d_r$V)  /  (C_d_r$B + C_d_l$B)^2 )^C_d_l$rate
    C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
    C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
    b_bw_s = ( (C_b_l$V + C_b_r$V)  /  ((C_b_r$B + C_b_l$B)^2 + scaleregul*(C_b_r$R+C_b_l$R)) )^C_b_l$rate
    C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
    C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
    h_bw_s = ( (C_h_l$V + C_h_r$V)  /  ((C_h_r$B + C_h_l$B)^2 + scaleregul*(C_h_r$R + C_h_l$R)) )^C_h_l$rate

    #if (C_b_l$V==0 | C_b_l$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) below the threshold. Range defined by bandwidth = %f\n", d_bw_s)  
    #if (C_b_r$V==0 | C_b_r$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) above the threshold. Range defined by bandwidth = %f\n", d_bw_s)  
    #if (C_h_l$V==0 | C_h_l$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) below the threshold. Range defined by bandwidth = %f\n", b_bw_s) 
    #if (C_h_r$V==0 | C_h_r$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) above the threshold. Range defined by bandwidth = %f\n", b_bw_s) 
}


#                     *** RD
if  (bwselect=="mserd" | bwselect=="cerrd" | bwselect=="msecomb1" | bwselect=="msecomb2" | bwselect=="cercomb1" | bwselect=="cercomb2" | bwselect=="" | all=="TRUE" ) {
  d_bw_d = ( (C_d_l$V + C_d_r$V)  /  (C_d_r$B - C_d_l$B)^2 )^C_d_l$rate
  C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
  C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
  b_bw_d = ( (C_b_l$V + C_b_r$V)  /  ((C_b_r$B - C_b_l$B)^2 + scaleregul*(C_b_r$R + C_b_l$R)) )^C_b_l$rate
  C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l)
  C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r)
  h_bw_d = ( (C_h_l$V + C_h_r$V)  /  ((C_h_r$B - C_h_l$B)^2 + scaleregul*(C_h_r$R + C_h_l$R)) )^C_h_l$rate
    
  #if (C_b_l$V==0 | C_b_l$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) below the threshold. Range defined by bandwidth = %f\n", d_bw_d)  
  #if (C_b_r$V==0 | C_b_r$B==0) printf("{err}Not enough variability to compute the bias bandwidth (b) above the threshold. Range defined by bandwidth = %f\n", d_bw_d)  
  #if (C_h_l$V==0 | C_h_l$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) below the threshold. Range defined by bandwidth = %f\n", b_bw_d) 
  #if (C_h_r$V==0 | C_h_r$B==0) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h) above the threshold. Range defined by bandwidth = %f\n", b_bw_d) 
}	
                     
#if (C_b_l$V==. | C_b_l$B==. | C_b_l$R==.) printf("{err}Invertibility problem in the computation of bias bandwidth (b) below the threshold") 			
#if (C_b_r$V==. | C_b_r$B==. | C_b_r$R==.) printf("{err}Invertibility problem in the computation of bias bandwidth (b) above the threshold")  
#if (C_h_l$V==. | C_h_l$B==. | C_h_l$R==.) printf("{err}Invertibility problem in the computation of loc. poly. bandwidth (h) below the threshold") 
#if (C_h_r$V==. | C_h_r$B==. | C_h_r$R==.) printf("{err}Invertibility problem in the computation of loc. poly. bandwidth (h) above the threshold") 
         
                     
if  (bwselect=="mserd" | bwselect=="cerrd" | bwselect=="msecomb1" | bwselect=="msecomb2" | bwselect=="cercomb1" | bwselect=="cercomb2" | bwselect=="" | all=="TRUE" ) {
  h_mserd = x_sd*h_bw_d
  b_mserd = x_sd*b_bw_d
}	
if  (bwselect=="msesum" | bwselect=="cersum" |  bwselect=="msecomb1" | bwselect=="msecomb2" |  bwselect=="cercomb1" | bwselect=="cercomb2"  |  all=="TRUE")  {
  h_msesum = x_sd*h_bw_s
  b_msesum = x_sd*b_bw_s
}
if  (bwselect=="msetwo" |  bwselect=="certwo" | bwselect=="msecomb2" | bwselect=="cercomb2"  | all=="TRUE")  {		
  h_msetwo_l = x_sd*h_bw_l
  h_msetwo_r = x_sd*h_bw_r
  b_msetwo_l = x_sd*b_bw_l
  b_msetwo_r = x_sd*b_bw_r
}
if  (bwselect=="msecomb1" | bwselect=="cercomb1" | all=="TRUE" ) {
  h_msecomb1 = min(c(h_mserd,h_msesum))
  b_msecomb1 = min(c(b_mserd,b_msesum))
}
if  (bwselect=="msecomb2" | bwselect=="cercomb2" |  all=="TRUE" ) {
  h_msecomb2_l = median(c(h_mserd,h_msesum,h_msetwo_l))
  h_msecomb2_r = median(c(h_mserd,h_msesum,h_msetwo_r))
  b_msecomb2_l = median(c(b_mserd,b_msesum,b_msetwo_l))
  b_msecomb2_r = median(c(b_mserd,b_msesum,b_msetwo_r))
}
cer_h = N^(-(p/((3+p)*(3+2*p))))
cer_b = N^(-(q/((3+q)*(3+2*q))))
	if  (bwselect=="cerrd" | all=="TRUE" ){
		h_cerrd = h_mserd*cer_h
		b_cerrd = b_mserd*cer_b
	}
	if  (bwselect=="cersum" | all=="TRUE" ){
		h_cersum = h_msesum*cer_h
		b_cersum=  b_msesum*cer_b
		}
	if  (bwselect=="certwo" | all=="TRUE" ){
		h_certwo_l   = h_msetwo_l*cer_h
		h_certwo_r   = h_msetwo_r*cer_h
		b_certwo_l   = b_msetwo_l*cer_b
		b_certwo_r   = b_msetwo_r*cer_b
		}
	if  (bwselect=="cercomb1" | all=="TRUE" ){
		h_cercomb1 = h_msecomb1*cer_h
		b_cercomb1 = b_msecomb1*cer_b
		}
	if  (bwselect=="cercomb2" | all=="TRUE" ){
		h_cercomb2_l = h_msecomb2_l*cer_h
		h_cercomb2_r = h_msecomb2_r*cer_h
		b_cercomb2_l = b_msecomb2_l*cer_b
		b_cercomb2_r = b_msecomb2_r*cer_b
	}


if (all=="FALSE"){
  results = matrix(NA,1,4)
  colnames(results)=c("h_l","h_r","b_l","b_r")
  #rownames(results)=c("CCT","IK","CV")
  rownames(results)=bwselect
  if  (bwselect=="mserd" | bwselect=="") results[1,] = c(h_mserd,      h_mserd,      b_mserd,      b_mserd)
  if  (bwselect=="msetwo")   results[1,] = c(h_msetwo_l,   h_msetwo_r,   b_msetwo_l,   b_msetwo_r)
  if  (bwselect=="msesum")   results[1,] = c(h_msesum,     h_msesum,     b_msesum,     b_msesum)
  if  (bwselect=="msecomb1") results[1,] = c(h_msecomb1,   h_msecomb1,   b_msecomb1,   b_msecomb1)
  if  (bwselect=="msecomb2") results[1,] = c(h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r) 
  if  (bwselect=="cerrd")    results[1,] = c(h_cerrd,      h_cerrd,      b_cerrd,      b_cerrd)
  if  (bwselect=="certwo")   results[1,] = c(h_certwo_l,   h_certwo_r,   b_certwo_l,   b_certwo_r)
  if  (bwselect=="cersum")   results[1,] = c(h_cersum,     h_cersum,     b_cersum,     b_cersum)
  if  (bwselect=="cercomb1") results[1,] = c(h_cercomb1,   h_cercomb1,   b_cercomb1,   b_cercomb1)
  if  (bwselect=="cercomb2") results[1,] = c(h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r)
}

  if (all=="TRUE"){
    bwselect="All"
    results = matrix(NA,10,4)
    colnames(results)=c("h_l","h_r","b_l","b_r")
    #rownames(results)=c("CCT","IK","CV")
    rownames(results)=c("mserd","msetwo","msesum","msecomb1","msecomb2","cerrd","certwo","cersum","cercomb1","cercomb2") 
    results[1,] =c(h_mserd,      h_mserd,      b_mserd,      b_mserd)
    results[2,] =c(h_msetwo_l,   h_msetwo_r,   b_msetwo_l,   b_msetwo_r)
    results[3,] =c(h_msesum,     h_msesum,     b_msesum,     b_msesum)
    results[4,] =c(h_msecomb1,   h_msecomb1,   b_msecomb1,   b_msecomb1)
    results[5,] =c(h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r) 
    results[6,] =c(h_cerrd,      h_cerrd,      b_cerrd,      b_cerrd)
    results[7,] =c(h_certwo_l,   h_certwo_r,   b_certwo_l,   b_certwo_r)
    results[8,] =c(h_cersum,     h_cersum,     b_cersum,     b_cersum)
    results[9,] =c(h_cercomb1,   h_cercomb1,   b_cercomb1,   b_cercomb1)
    results[10,]=c(h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r)
  }
  
  tabl1.str=matrix(NA,4,1)
  dimnames(tabl1.str) <-list(c("BW Selector", "Number of Obs", "NN Matches", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=bwselect
  tabl1.str[2,1]=N
  tabl1.str[3,1]=nnmatch
  tabl1.str[4,1]=kernel_type
  
  tabl2.str=matrix(NA,3,2)
  colnames(tabl2.str)=c("Left","Right")
  rownames(tabl2.str)=c("Number of Obs","Order Loc Poly (p)","Order Bias (q)")
  tabl2.str[1,]=formatC(c(N_l,N_r),digits=0, format="f")
  tabl2.str[2,]=formatC(c(p,p),digits=0, format="f")
  tabl2.str[3,]=formatC(c(q,q),digits=0, format="f")

  bws=results
  out = list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,bws=bws,bws,bwselect=bwselect,kernel=kernel_type,p=p,q=q)
  out$call <- match.call()
  class(out) <- "rdbwselect"
  return(out)
}


#rdbwselect <- function(y,x, ...) UseMethod("rdbwselect")

#rdbwselect.default <- function(y,x, ...){
#  est <- rdbwselectEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "rdbwselect"
#  est
#}

print.rdbwselect <- function(x,...){
  cat("Call:\n")
  print(x$call)
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F) 
  cat("\n")
  print(x$bws)  
}

summary.rdbwselect <- function(object,...) {
  TAB <- object$bws
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.rdbwselect"
  res
}

#print.summary.rdbwselect <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
#}