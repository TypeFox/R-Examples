uncond=function(a,b,c,d,n,alpha){
		if (a+b+c+d!=n){return(paste("Caution: a+b+c+d is not equal to n."))}	
		if (a+b+c+d==n){
		z=qnorm(1-alpha/2)
		psidach=(b+c)/n
		tetadach=(b-c)/n
		dec=c()
		dec2=c()
		teta_v=seq(-0.9995,0.9995,0.0005)
		i=0
		if (b>0 & c>0 & (a+d)>0){
			for (teta in teta_v ){
			i=i+1
				if (teta==0){psiteta=psidach
				
				if ((a+d)*(log(1-psiteta)-log(1-psidach))+b*(log(psiteta)-log(psidach+tetadach))+c*(log(psiteta)-log(psidach-tetadach))>= -z^2/2){
				dec[i]=1}
				if ((a+d)*(log(1-psiteta)-log(1-psidach))+b*(log(psiteta)-log(psidach+tetadach))+c*(log(psiteta)-log(psidach-tetadach))< -z^2/2){
				dec[i]=0}
				}
				if(teta!=0){
				B=1/2*(psidach)+1/2*teta*(tetadach)
				C=tetadach*teta-((a+d)/n)*teta^2
				psiteta=B+sqrt(B^2-C)
				
				if ((a+d)*(log(1-psiteta)-log(1-psidach))+b*(log(psiteta+teta)-log(psidach+tetadach))+c*(log(psiteta-teta)-log(psidach-tetadach))>= -z^2/2){
				dec[i]=1
				}
				if ((a+d)*(log(1-psiteta)-log(1-psidach))+b*(log(psiteta+teta)-log(psidach+tetadach))+c*(log(psiteta-teta)-log(psidach-tetadach))< -z^2/2){
				dec[i]=0}
				}
			}
			ma_all=cbind(teta_v,dec)
			like_l=min(ma_all[ma_all[,2]==1,1])
			like_u=max(ma_all[ma_all[,2]==1,1])
			cint=c(like_l,like_u)
		}
		
		# second scenario
		if(c==0 & b>0 & (a+d)>0){ 
			for (teta in teta_v ){
			i=i+1
			psiteta=max(teta,b/n-(1-b/n)*teta)
			if (psiteta!=1){	
				if ((a+d)*(log(1-psiteta)-log(1-psidach))+b*(log(psiteta+teta)-log(psidach+tetadach))>= -z^2/2){
				dec2[i]=1
				}
				if ((a+d)*(log(1-psiteta)-log(1-psidach))+b*(log(psiteta+teta)-log(psidach+tetadach))< -z^2/2){
				dec2[i]=0}
				
				}
				}
				
				ma_all=cbind(teta_v,dec2)
				like_l=min(ma_all[ma_all[,2]==1,1])
				like_u=max(ma_all[ma_all[,2]==1,1])
				cint=c(like_l,like_u)
		}
			
		# third scenario
		if (a==0 & d==0){
			psiteta=1
			for (teta in teta_v ){
			i=i+1
			if (b*(log(psiteta+teta)-log(psidach+tetadach))+c*(log(psiteta-teta)-log(psidach-tetadach))>= -z^2/2){
				dec2[i]=1
				}
				if (b*(log(psiteta+teta)-log(psidach+tetadach))+c*(log(psiteta-teta)-log(psidach-tetadach))< -z^2/2){
				dec2[i]=0}
				}
				ma_all=cbind(teta_v,dec2)
				like_l=min(ma_all[ma_all[,2]==1,1])
				like_u=max(ma_all[ma_all[,2]==1,1])
			    cint=c(like_l,like_u)
				
		}
		# fourth scenario
		if (c==0 & b==0 & (a+d)>0){
			for (teta in teta_v ){
			i=i+1
			psiteta=abs(teta)
			if (psiteta!=1){
				if ((a+d)*(log(1-psiteta)-log(1-psidach))>= -z^2/2){
				dec2[i]=1
				}
				if ((a+d)*(log(1-psiteta)-log(1-psidach))< -z^2/2){
				dec2[i]=0}
				}
				}
				ma_all=cbind(teta_v,dec2)
				like_l=min(ma_all[ma_all[,2]==1,1])
				like_u=max(ma_all[ma_all[,2]==1,1])
			    cint=c(like_l,like_u)
				
		}
		c(like_l,like_u)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = tetadach)
		class(rval) <- "htest"
		return(rval)
		}}