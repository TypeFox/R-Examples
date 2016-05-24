KSADCHRY=function(e,tip,B,alpha){ 
 
  M=length(e)
  Ff=0
	observed=0	
	expected=0
    if (tip==1){ #Excursion test
		if (B==16){	
			Ff[1]=.196
			Ff[2]=.392
			Ff[3]=.575
			Ff[4]=.732
			Ff[5]=1
			observed[1]=sum(e[]==0)/M
			observed[2]=observed[1]+sum(e[]==1)/M
			observed[3]=observed[2]+sum(e[]==2)/M
			observed[4]=observed[3]+sum(e[]==3)/M
			observed[5]=observed[4]+sum(e[]>=4)/M		
		}else if (B==32){
			Ff[1]=.279
			Ff[2]=.54
			Ff[3]=.75
			Ff[4]=1
			observed[1]=sum(e[]<=1)/M
			observed[2]=observed[1]+sum((e[]>=2) & (e[]<=3))/M
			observed[3]=observed[2]+sum((e[]>=4) & (e[]<=5))/M
			observed[4]=observed[3]+sum((e[]>=6))/M
		} else if (B==64){
			Ff[1]=.198
			Ff[2]=.390
			Ff[3]=.563
			Ff[4]=.769
			Ff[5]=1
			observed[1]=sum(e[]<=1)/M
			observed[2]=observed[1]+sum((e[]>=2) & (e[]<=3))/M
			observed[3]=observed[2]+sum((e[]>=4) & (e[]<=5))/M
			observed[4]=observed[3]+sum((e[]>=6) & (e[]<=8))/M
			observed[5]=observed[4]+sum((e[]>=9))/M
		}else if (B==128){
			Ff[1]=.21
			Ff[2]=.410
			Ff[3]=.587
			Ff[4]=.771
			Ff[5]=1
			observed[1]=sum(e[]<=2)/M	
			observed[2]=observed[1]+sum((e[]>=3) & (e[]<=5))/M
			observed[3]=observed[2]+sum((e[]>=6) & (e[]<=8))/M
			observed[4]=observed[3]+sum((e[]>=9) & (e[]<=12))/M
			observed[5]=observed[4]+sum((e[]>=13))/M
		}else if (B==256){
			Ff[1]=.198
			Ff[2]=.387
			Ff[3]=.593
			Ff[4]=.804
			Ff[5]=1
			observed[1]=sum(e[]<=3)/M	
			observed[2]=observed[1]+sum((e[]>=4) & (e[]<=7))/M
			observed[3]=observed[2]+sum((e[]>=8) & (e[]<=12))/M
			observed[4]=observed[3]+sum((e[]>=13) & (e[]<=19))/M
			observed[5]=observed[4]+sum((e[]>=20))/M
		}	
	} else if (tip==2){
        if (B==64){ #Height test
			Ff[1]=0.252
			Ff[2]=0.479
			Ff[3]=0.661
			Ff[4]=0.834
			Ff[5]=1
			observed[1]=sum(e[]<=6)/M	
			observed[2]=observed[1]+sum((e[]>=7) & (e[]<=8))/M
			observed[3]=observed[2]+sum((e[]>=9) & (e[]<=10))/M
			observed[4]=observed[3]+sum((e[]>=11) & (e[]<=13))/M
			observed[5]=observed[4]+sum((e[]>=14))/M
		}else if (B==128){
			Ff[1]=0.18
			Ff[2]=0.421
			Ff[3]=0.63
			Ff[4]=0.813
			Ff[5]=1
			observed[1]=sum(e[]<=8)/M	
			observed[2]=observed[1]+sum((e[]>=9) & (e[]<=11))/M
			observed[3]=observed[2]+sum((e[]>=12) & (e[]<=14))/M
			observed[4]=observed[3]+sum((e[]>=15) & (e[]<=18))/M
			observed[5]=observed[4]+sum((e[]>=19))/M
		}else if (B==256){
			Ff[1]=0.196
			Ff[2]=0.426
			Ff[3]=0.62
			Ff[4]=0.815
			Ff[5]=1
			observed[1]=sum(e[]<=12)/M	
			observed[2]=observed[1]+sum((e[]>=13) & (e[]<=16))/M
			observed[3]=observed[2]+sum((e[]>=17) & (e[]<=20))/M
			observed[4]=observed[3]+sum((e[]>=21) & (e[]<=26))/M
			observed[5]=observed[4]+sum((e[]>=27))/M
		}else if (B==512){
			Ff[1]=0.18
			Ff[2]=0.385
			Ff[3]=0.599
			Ff[4]=0.794
			Ff[5]=1
			observed[1]=sum(e[]<=17)/M
			observed[2]=observed[1]+sum((e[]>=18) & (e[]<=22))/M
			observed[3]=observed[2]+sum((e[]>=23) & (e[]<=28))/M
			observed[4]=observed[3]+sum((e[]>=29) & (e[]<=36))/M
			observed[5]=observed[4]+sum((e[]>=37))/M
		}else if (B==1024){
			Ff[1]=0.196
			Ff[2]=0.399
			Ff[3]=0.599
			Ff[4]=0.803
			Ff[5]=1
			observed[1]=sum(e[]<=25)/M	
			observed[2]=observed[1]+sum((e[]>=26) & (e[]<=32))/M
			observed[3]=observed[2]+sum((e[]>=33) & (e[]<=40))/M
			observed[4]=observed[3]+sum((e[]>=41) & (e[]<=52))/M
			observed[5]=observed[4]+sum((e[]>=53))/M
         	}
    } else if (tip==3){ #Expansion
		if (B==32){
			Ff[1]=0.307
			Ff[2]=0.473
			Ff[3]=0.739
			Ff[4]=1	
			observed[1]=sum(e[]<=6)/M	
			observed[2]=observed[1]+sum((e[]==7))/M
			observed[3]=observed[2]+sum((e[]>=8) & (e[]<=9))/M
			observed[4]=observed[3]+sum((e[]>=10))/M
		}else if (B==64){
			Ff[1]=0.188
			Ff[2]=0.422
			Ff[3]=0.634
			Ff[4]=0.84
			Ff[5]=1
			observed[1]=sum(e[]<=8)/M	
			observed[2]=observed[1]+sum((e[]>=9) & (e[]<=10))/M
			observed[3]=observed[2]+sum((e[]>=11) & (e[]<=12))/M
			observed[4]=observed[3]+sum((e[]>=13) & (e[]<=15))/M
			observed[5]=observed[4]+sum((e[]>=16))/M
		}else if (B==128){
			Ff[1]=0.196
			Ff[2]=0.361
			Ff[3]=0.596
			Ff[4]=0.81
			Ff[5]=1
			observed[1]=sum(e[]<=12)/M	
			observed[2]=observed[1]+sum((e[]>=13) & (e[]<=14))/M
			observed[3]=observed[2]+sum((e[]>=15) & (e[]<=17))/M
			observed[4]=observed[3]+sum((e[]>=18) & (e[]<=21))/M
			observed[5]=observed[4]+sum((e[]>=22))/M
		}
	}

	N=length(Ff)
	
	expected2=0
	expected2=round(Ff*2000)
    obs=0
	obs2=0
	obs2=round(observed*2000)
	obs[1]=obs2[1]
	expected[1]=expected2[1]
	for (i in 2:length(observed)){
		obs[i]=obs2[i]-obs2[i-1]
		expected[i]=expected2[i]-expected2[i-1]
	}
	

 	z=rep(0:(length(expected)-1),expected)
 	e=rep(0:(length(obs)-1),obs) 
  
	test=kSamples::ad.test(e,z=z,method="simulated",dist=FALSE,Nsim=1000)
	ADtest=test$ad
	if (ADtest[1,3]<alpha){ 
		sonucAD=0
	}
	else
	{ 
		sonucAD=1
	}
	KStest=matrix(0,1,2)
	test2=ks.test(e,z)

	KStest[1]=test2$statistic
	KStest[2]=test2$p.value

	if (KStest[2]<alpha){ 
		sonucKS=0
	}
	else
	{ 
		sonucKS=1
	}
  
  KKtest=0
	KKtest[1]=sum(((obs-expected)^2)/expected)
  KKtest[2]=1-pchisq(KKtest[1],N-1)
	if (KKtest[2]<alpha){ 
	  sonucKK=0 
	}
	else
	{ 
	  sonucKK=1
	}
   
 	result=list(sonucAD=sonucAD,ADtest=ADtest,sonucKS=sonucKS,KStest=KStest,sonucKK=sonucKK,KKtest=KKtest)
	return(result)
}