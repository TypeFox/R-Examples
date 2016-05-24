# mutistagecor.R author: xuefei mi, first version 17-04-2013, for selectiongain package v2.0.6

# mutistagecor.R author: xuefei mi, second version 11-08-2013, for selectiongain package since v2.0.18

# the package are totally rewritten

# mutistagecor.R author: xuefei mi, 3th version 03-01-2014, for selectiongain package since v2.0.30

#M for type of testers are added

# mutistagecor.R author: xuefei mi, 4th version 11-01-2014, for selectiongain package since v2.0.30

# parental effect was considered

# v45, for longin 2015, genomic selection in wheat paper. 2015-03-01

# v46, add 3 method for parenet-bs1-gca etc. 2015-03-06

#  mutistagecor.R author: xuefei mi, 2015-10-08, for selectiongain package since v2.0.35, 40, 47

# defualt values of Vline in covwithoutmas was added


# v48, integrated bug fixing of gsonly scheme in TAG 2015 paper. marked by ####PleaseCheck



`multistagecor` <-function(maseff=NA,VGCAandE=c(0,0,0,0,0),VSCA=c(0,0,0,0),VLine=c(0,0,0,0,0), ecoweight=c(0,rep(1,length(T-1))/length(T-1)), rhop=0,  T, L,M=rep(1,length(T)),Rep,index=FALSE,covtype=c("LonginII"),detail=FALSE)
{

# we had two optition of VGCA,  either VGCA or VGCAandE  
V=VGCAandE
VarianceType=V[1]
  
  if (VarianceType=="VC2" )
  {
    Vg=1
    Vgl=0.5
    Vgy=0.5
    Vgly=1
    Ve=2
  }
  else if (VarianceType=="VC1")
  {
    Vg=1
    Vgl=0.25
    Vgy=0.25
    Vgly=0.5
    Ve=1
  }
  else if (VarianceType=="VC3")
  {
    Vg=1
    Vgl=1
    Vgy=1
    Vgly=2
    Ve=4
  }
 
  else  if (length(VGCAandE)==5& all(is.numeric(VGCAandE)))
  { 
    Vg=V[1]
    Vgl=V[2]
    Vgy=V[3]
    Vgly=V[4]
    Ve=V[5]
  }else
  {
   stop("the value of VGCAandE has to be numeric with length 5, or VC1,VC2,VC3 as defined in Longin's paper")
  }
 		
 
 
  dim = length(L)+1
  
  if (length(L) != length(Rep))
  {
    warning("Location and Replicates do not have same length", call. = FALSE)
  }
  
   if (length(L) != length(T))
  {
    warning("Location and testers do not have same length", call. = FALSE)
  }
   
  if (length(L) != length(M) & covtype=="LonginII-M")
  {
    warning("Location and type of testers do not have same length", call. = FALSE)
  }
   

  if (length(L) ==1 && index==TRUE)
  {
    warning("selection index need n >1", call. = FALSE)
  }
  

  # preparation for the VSCA
  
  
  VSCAType=VSCA[1]
  
  if (VSCAType=="VC2.2" )
  {
    Vs=0.5
    Vsl=0.25
    Vsy=0.25
    Vsly=0.5

  }
  else if (VSCAType=="VC2.1")
  {
    Vg=0.5
    Vgl=0.125
    Vgy=0.125
    Vgly=0.25

  }else if(all(is.numeric(VSCA))& length (VSCA)==4)
  {
    Vs=VSCA[1]
    Vsl=VSCA[2]
    Vsy=VSCA[3]
    Vsly=VSCA[4]

   }else
   {
     stop("the value of VSCA has to be numeric with length 4, or VC2.1,VC2.2, as defined in Longin's paper")
   }
  

  VLineType=VLine[1]
  
  if (VLineType=="VC4" )
  {
    VL=1
    VLl=0.15
    VLy=0.15
    VLly=0.50
    VLe=0.50

  }
  else if (VLineType=="VC5")
  {
     VL=1
    VLl=0.30
    VLy=0.30
    VLly=1
    VLe=1

  } else if (VLineType=="VC6")
  {
     VL=1
    VLl=0.60
    VLy=0.60
    VLly=2
    VLe=2

  }else if(all(is.numeric(VLine))& length (VLine)==5)
  {
    VL=VLine[1]
    VLl=VLine[2]
    VLy=VLine[3]
    VLly=VLine[4]
    VLe=VLine[5]

   }else
   {
     stop("the value of VLine has to be numeric with length 5, or VC4,VC5,VC6, as defined in Longin's paperII")
   }


# function for calculating cov without markers


covwithoutmas  <-function(V=c(0,0,0,0,0),VSCA=c(0,0,0,0),VLine=c(0,0,0,0,0), ecoweight,rhop, T, L,M=length(T),Rep,index=FALSE,covtype=c("LonginII")) 
  ####PleaseCheck : I think in this line "VLineVLine" actually corresponds to "VLine"
  # modifyed by mi at 2015-11-05, yes this is a typo
  {
 
# we had two optition of VGCA,  either VGCA or VGCAandE  

VarianceType=V[1]
  
  if (VarianceType=="VC2" )
  {
    Vg=1
    Vgl=0.5
    Vgy=0.5
    Vgly=1
    Ve=2
  }
  else if (VarianceType=="VC1")
  {
    Vg=1
    Vgl=0.25
    Vgy=0.25
    Vgly=0.5
    Ve=1
  }
  else if (VarianceType=="VC3")
  {
    Vg=1
    Vgl=1
    Vgy=1
    Vgly=2
    Ve=4
  }
 
  else  if (length(VGCAandE)==5& all(is.numeric(VGCAandE)))
  { 
    Vg=V[1]
    Vgl=V[2]
    Vgy=V[3]
    Vgly=V[4]
    Ve=V[5]
  }else
  {
   stop("the value of VGCAandE has to be numeric with length 5, or VC1,VC2,VC3 as defined in Longin's paper")
  }
 		
 
 
  dim = length(L)+1
  
  if (length(L) != length(Rep))
  {
    warning("Location and Replicates do not have same length", call. = FALSE)
  }
  
   if (length(L) != length(T))
  {
    warning("Location and testers do not have same length", call. = FALSE)
  }
   
  if (length(L) != length(M) & covtype=="LonginII-M")
  {
    warning("Location and type of testers do not have same length", call. = FALSE)
  }
   

  if (length(L) ==1 && index==TRUE)
  {
    warning("selection index need n >1", call. = FALSE)
  }
  

  # preparation for the VSCA
  
  
  VSCAType=VSCA[1]
  
  if (VSCAType=="VC2.2" )
  {
    Vs=0.5
    Vsl=0.25
    Vsy=0.25
    Vsly=0.5

  }
  else if (VSCAType=="VC2.1")
  {
    Vg=0.5
    Vgl=0.125
    Vgy=0.125
    Vgly=0.25

  }else if(all(is.numeric(VSCA))& length (VSCA)==4)
  {
    Vs=VSCA[1]
    Vsl=VSCA[2]
    Vsy=VSCA[3]
    Vsly=VSCA[4]

   }else
   {
     stop("the value of VSCA has to be numeric with length 4, or VC2.1,VC2.2, as defined in Longin's paper")
   }
  

  VLineType=VLine[1]
  
  if (VLineType=="VC4" )
  {
    VL=1
    VLl=0.15
    VLy=0.15
    VLly=0.50
    VLe=0.50

  }
  else if (VLineType=="VC5")
  {
     VL=1
    VLl=0.30
    VLy=0.30
    VLly=1
    VLe=1

  } else if (VLineType=="VC6")
  {
     VL=1
    VLl=0.60
    VLy=0.60
    VLly=2
    VLe=2

  }else if(all(is.numeric(VLine))& length (VLine)==5)
  {
    VL=VLine[1]
    VLl=VLine[2]
    VLy=VLine[3]
    VLly=VLine[4]
    VLe=VLine[5]

   }else
   {
     stop("the value of VLine has to be numeric with length 5, or VC4,VC5,VC6, as defined in Longin's paperII")
   }

   
  
# calculate the covariance  

if (covtype=="LonginII")
 { 
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1]    + Ve/L[i-1]/Rep[i-1]/T[i-1]
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1],T[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1],T[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-M")
 { 
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental")
 { 
  
	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5
	
  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)
  cov[1,2]<-covlg*a2+a1*VL
  cov[2,1]=cov[1,2]
	
	cov[2,2]=VL+VLy+ (VLl+VLly)/L[1]+VLe/L[1]/Rep[1]/T[1]
	 for (i in 3:dim)
  {
	cov[2,i]=covlg
	cov[i,2]=cov[2,i]
	
	cov[1,i]=a1*covlg+a2*Vg
	cov[i,1]=cov[1,i]
	
	# covVLl,Vgl are assumed to be 0, if not 0, add the code here
	
	}
	
  for (i in 3:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]
	
		
  }

  for (i in 3:dim)
  {
    for (j in 3:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-gca")
 { 
  
	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5
	
  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)
  
  cov[1,1]<-Vg
  
  cov[1,2]=covlg
  cov[2,1]=cov[1,2]
	
	cov[2,2]=VL+VLy+ (VLl+VLly)/L[1]+VLe/L[1]/Rep[1]/T[1]
	 for (i in 3:dim)
  {
	cov[2,i]=covlg
	cov[i,2]=cov[2,i]
	
	cov[1,i]=Vg
	cov[i,1]=cov[1,i]
	
	# covVLl,Vgl are assumed to be 0, if not 0, add the code here
	
	}
	
  for (i in 3:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]
	
		
  }

  for (i in 3:dim)
  {
    for (j in 3:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-perse")
 { 
  
	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5
	
  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)
  
  cov[1,1]=VL
  
  cov[1,2]=VL
  cov[2,1]=cov[1,2]
	
	cov[2,2]=VL+VLy+ (VLl+VLly)/L[1]+VLe/L[1]/Rep[1]/T[1]
	 for (i in 3:dim)
  {
	cov[2,i]=covlg
	cov[i,2]=cov[2,i]
	
	cov[1,i]=covlg
	cov[i,1]=cov[1,i]
	
	# covVLl,Vgl are assumed to be 0, if not 0, add the code here
	
	}
	
  for (i in 3:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]
	
		
  }

  for (i in 3:dim)
  {
    for (j in 3:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-BS1")
 { 
  
	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5
	
  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)
  cov[1,2:dim]=a1*covlg+a2*Vg
  cov[2:dim,1]=cov[1,2:dim]
	
	# if might need to be modifyed if dim>3
	
	
	
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]
	
		
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-BS1-perse")
 { 
  
	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5
	
  cov= diag(dim)*(VL)
  cov[1,2:dim]=a1*covlg
  cov[2:dim,1]=cov[1,2:dim]
	
	# if might need to be modifyed if dim>3
	
	
	
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]
	
		
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="Heffner")
 {
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
	LT=L*T
	LTR=L*T*Rep
	
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy/(i-1)+(Vgl+Vgly)/sum(L[1:(i-1)])+  (Vs + Vsy)/sum(T[1:(i-1)]) + (Vsl+Vsly)/sum(LT[1:(i-1)])   + Ve/sum(LTR[1:(i-1)])
  
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {    
        cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+Vs/max(T[i-1],T[j-1])+
        Vsl/max(T[i-1],T[j-1])/max(L[i-1],L[j-1])
        cov[j,i]=cov[i,j]
      }
     }
  } 
 }else 
 {
  stop("covtype is not specified, see Longin's paperII")
 }
  
# calculate the optimal selection index = G^-1 /P


cov

}

# main function begins

        dim = length(L)+1 
        Vg=1

	if (all(is.numeric(VGCAandE)))
	{ 
	   if (length(VGCAandE)==5)
           { 
              Vg=VGCAandE[1] 
	   }else
           {
             stop("the value of VGCAandE has to be numeric with length 5, or VC1,VC2,VC3 as defined in Longin's paper")
           }
 		
        }
        
        if (length(L)!=length(Rep) | length(T)!=length(L))
        {
            stop("T, L and Rep must have the same length")   
        }
        
        
        if (is.na(maseff))
	 {
              tempb="empty"
	      cov=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep, covtype=covtype)
				           a1=ecoweight[1]
	                 a2=ecoweight[2]
	                 covlg=rhop*(Vg*VL)^0.5

				
              if (index==TRUE)
             {                
                P=cov[-1,-1]
                G=matrix(rep(1,(dim-1)^2), nrow = dim-1, ncol=dim-1, byrow=TRUE) * Vg
                if (covtype=="LonginII-Parental" | covtype== "LonginII-Parental-gca" | covtype== "LonginII-Parental-perse")
									 {
									  G[1,1]=VL
										# genotypic covariance of stage 1
										
									  covlg=rhop*(Vg*VL)^0.5
									   for (i in 2:c(dim-1))
                     {
	                     G[1,i]=covlg
	                     G[i,1]=G[1,i]
	
	                     # covVLl,Vgl are assumed to be 0, if not 0, add the code here
	
	                    }
									  
									 }
									 
									  
									
									
                
                for (i in 2:c(dim-1))
                {
                   tempb=matrix(c(a1,rep((1-a1)/(dim-2),i-1)),nrow=1,ncol=i,byrow=TRUE)
									 
                   tempp=P[1:i,1:i]
                   tempg=G[1:i,1:i]
                   tempb= t((solve(tempp)%*% tempg )  %*% t( tempb))
                   
                   tempb=tempb / sum(tempb)
                   	if (covtype=="LonginII-Parental-BS1")
									 {
									  	if (tempb[1]!=0)
									  	{
									  	tempb=tempb / tempb[1]
									  	}else
									  	{
									  	tempb=tempb / sum(tempb)
									  	}
									  }
                   
                   P[i,i]=  tempb %*% tempp %*% t(tempb)
                   P[1:(i-1),i]= tempb %*% tempp[,1:(i-1)]
                   P[i,1:(i-1)]=t(P[1:(i-1),i])
                }

# calculate the covariance and give output
                   cov[2:dim,2:dim]=P
									
									if (covtype=="LonginII-Parental")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                   # 	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]
									
									# formulation have to be improved when dim >3
									
									 cov[1,3:dim]=a1*b1*VL+a2*b2*Vg+(a1*b2+a2*b1)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }else if (covtype=="LonginII-Parental-gca")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                   # 	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]
									
									# formulation have to be improved when dim >3
									
									 cov[1,3:dim]=b2*Vg+(b1)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }else if (covtype=="LonginII-Parental-perse")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                   # 	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]
									
									# formulation have to be improved when dim >3
									
									 cov[1,3:dim]=b1*VL+(b2)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }else 	if (covtype=="LonginII-Parental-BS1")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                    #	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]
									
									# formulation have to be improved when dim >3
									   
										cov[1,2]=a2*Vg+a1*covlg
									  cov[2,1]=cov[1,2]
									
									 cov[1,3:dim]=(a2*b1+a2*b2)*Vg+(a1*b2+a1*b1)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }
									
									
              }else if (index!=TRUE)
							{
							 
							}
	
                   output= list(cov2cor(cov),tempb,cov)
	   } else if(!is.na(maseff))
       {
             if (length(maseff)!=1)
             {
                 stop("if maseff is not NA the then the length of it has to be 1.")   
             }
              if (covtype=="LonginII-Parental")
             {
                cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
								covtype=covtype)
               	dim=dim-1
                corp=cov	
								covlg=rhop*(Vg*VL)^0.5
	                cormas= diag(dim+1)
		            cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
		            cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
		            cormas[3:c(dim+1),1]=corp[2:c(dim),1]
	            	cormas[1,2]=maseff^2*corp[1,3]
		            cormas[2,1]=maseff^2*corp[3,1]
	            	cormas[2,2]=maseff^2*Vg
	              cormas[1,1]=corp[1,1]
		          
                    cormas[2,3]=maseff^2*corp[2,3] 
	 	           cormas[3,2]=cormas[2,3]  

	          	for (i in 4:c(dim+1))
	      	   {
		           cormas[2,i]=maseff^2*Vg
	 	           cormas[i,2]=cormas[2,i]
	           } 
                
                
            
             
             tempb="empty"
             if (index==TRUE)
               {
                warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
               }
             }else if(covtype=="LonginII-Parental-gca")
             {
               
                 cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
                                    covtype=covtype)
                 dim=dim-1
                 corp=cov							
                 cormas= diag(dim+1)
                 cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
                 cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
                 cormas[3:c(dim+1),1]=corp[2:c(dim),1]
                 cormas[1,2]=maseff^2*Vg
                 cormas[2,1]=maseff^2*Vg
                 cormas[2,2]=maseff^2*Vg
                 cormas[1,1]=corp[1,1]
                 
                 cormas[2,3]=maseff^2*corp[1,2] 
                 cormas[3,2]=cormas[2,3]  
                 
                 for (i in 4:c(dim+1))
                 {
                   cormas[2,i]=maseff^2*Vg
                   cormas[i,2]=cormas[2,i]
                 } 
                 
                 
              
               tempb="empty"
               if (index==TRUE)
               {
                 warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
               }
               
             }else if(covtype=="LonginII-Parental-perse")
             {
               
               cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
                                  covtype=covtype)
               dim=dim-1
               corp=cov  						
               cormas= diag(dim+1)
               cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
               cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
               cormas[3:c(dim+1),1]=corp[2:c(dim),1]
               cormas[1,2]=maseff^2*corp[1,3] 
               cormas[2,1]=maseff^2*corp[1,3] 
               cormas[2,2]=maseff^2*Vg
               cormas[1,1]=corp[1,1]
               
               cormas[2,3]=maseff^2*corp[1,3] 
               cormas[3,2]=cormas[2,3]  
               
               for (i in 4:c(dim+1))
               {
                 cormas[2,i]=maseff^2*Vg
                 cormas[i,2]=cormas[2,i]
               } 
               
               
            

             tempb="empty"
              if (index==TRUE)
              {
                warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
              }

               
             }else if(covtype=="LonginII-Parental-BS1")
             {
               
               cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
                                  covtype=covtype)
               dim=dim-1
               corp=cov    					
               cormas= diag(dim+1)
               cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
               cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
               cormas[3:c(dim+1),1]=corp[2:c(dim),1]
               cormas[1,2]=maseff^2*corp[1,3] 
               cormas[2,1]=maseff^2*corp[1,3] 
               cormas[2,2]=maseff^2*Vg
               cormas[1,1]=corp[1,1]
               
               cormas[2,3]=maseff^2
               cormas[3,2]=cormas[2,3]  
               
               for (i in 4:c(dim+1))
               {
                 cormas[2,i]=maseff^2*Vg
                 cormas[i,2]=cormas[2,i]
               } 
               
               
               
               
               tempb="empty"
               if (index==TRUE)
               {
                 warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
               }
               
               
             }else if (index!=TRUE & covtype!="LonginII-Parental"& covtype!="LonginII-Parental-gca"& covtype!="LonginII-Parental-perse"& covtype!="LonginII-Parental-BS1")
             {
                cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
								covtype=covtype)
               	dim=dim-1
                corp=cov
							
	              cormas= diag(dim+1)
		            cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
		            cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
		            cormas[3:c(dim+1),1]=corp[2:c(dim),1]
	            	cormas[1,2]=maseff^2*corp[2,1]  ####PleaseCheck  : if the Length of L is 2 and we are using mas, corp will be of dimentions 2x2. Then corp[3,1] does not exist. The construction of the 
					                                                   # cov matrix (cov=corp) indicates that for covtype=LonginII corp[3,1]=corp[2,1]=corp[1,1]=Vg So we can replace corp[3,1] by 
																	   # corp[2,1] or corp[1,1] or Vg . I replace it by corp[2,1] and the packege worked. Is this change correct?
	            	                     # modifyed by mi 2015.11.05, yes you can use corp[2,1]   
	            	
		            cormas[2,1]=cormas[1,2]
	            	cormas[2,2]=maseff^2*Vg
	              cormas[1,1]=corp[1,1]
		
	          	for (i in 3:c(dim+1))
	      	   {
		           cormas[2,i]=maseff^2*Vg
	 	           cormas[i,2]=cormas[2,i]
	           } 
              


             }
             
	           tempb="empty"
             output= list(cov2cor(cormas),tempb,cormas)
         }
if (detail==TRUE)
{
output
}else
{
output[[1]]
}

}












