BICcalc <-
function(conv,G,p,mod,q,logl,n,gauss){
	if(conv==0){
		BIC <- -Inf
	}
	if(conv==1){
		freepar <- G-1 + G*p

		if(substring(mod,1,1)=="U"){
			freepar <- freepar + G*(p*q - q*(q-1)/2)
		}
		if(substring(mod,1,1)=="C"){
			freepar <- freepar + (p*q - q*(q-1)/2)
		}
		if(substring(mod,2,2)=="U"){
			if(substring(mod,3,3)=="U"){
				freepar <- freepar + G*p
			}
			if(substring(mod,3,3)=="C"){
				freepar <- freepar + G
			}
		}
		if(substring(mod,2,2)=="C"){
			if(substring(mod,3,3)=="U"){
				freepar <- freepar + p
			}
			if(substring(mod,3,3)=="C"){
				freepar <- freepar + 1
			}
		}
		if(substring(mod,3,3)=="1"){
			freepar <- freepar + (p*q-q*(q-1)/2) + G + (p-1)
		}
		if(substring(mod,3,3)=="2"){
			freepar <- freepar + G*(p*q-q*(q-1)/2) + G + (p-1)
		}
		if(substring(mod,3,3)=="3"){
			freepar <- freepar + (p*q-q*(q-1)/2) + 1 + G*(p-1)
		}
		if(substring(mod,3,3)=="4"){
			freepar <- freepar + G*(p*q-q*(q-1)/2) + 1 + G*(p-1)
		}
		if(substring(mod,4,4)=="U"){
			freepar <- freepar + G
		}
		if(substring(mod,4,4)=="C"){
			freepar <- freepar + 1
		}
		if(gauss){freepar <- freepar - 1}
		BIC <- 2 * max(logl) - freepar*log(n)
	}
	BIC
}
