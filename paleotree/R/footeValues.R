#' Calculates Values for Foote's Inverse Survivorship Analyses
#'
#' This function calculates the intermediary values needed for
#' fitting Foote's inverse survivorship analyses, as listed in the
#' table of equations in Foote (2003), with the analyses themselves
#' described further in Foote (2001) and Foote (2005).

#' @details Although most calculations in this function agree
#' with the errata for Foote's 2003 table (see references), there were some additional
#' corrections for Prob(D|FL) made as part of a personal communication in 2013
#' between the package author and Michael Foote.

#' @inheritParams inverseSurv

#' @param p Instantaneous origination/branching rate of taxa. Under 
#' a continuous model, assumed to be \emph{per interval}, or equal
#' to the product of interval lengths and the rates per lineage time
#' units for each interval. Under a pulsed mode (p_cont=FALSE), p is a 
#' per-interval 'rate' which can exceed 1 (because diversity can
#' more than double; Foote, 2003a). Given as a vector with length
#' equal to the number of intervals, so a different value may be
#' given for each separate interval. Must be the same length as
#' q and r. 

#' @param q Instantaneous extinction rate of taxa. Under 
#' a continuous model, assumed to be \emph{per interval}, or
#' equal to the product of interval lengths and the rates per lineage
#' time units for each interval. Under a pulsed mode (q_cont=FALSE), q is a  
#' per-interval 'rate' but which cannot be observed to exceed 1
#' (because you can't have more taxa go extinct than exist). Given as
#' a vector with length equal to the number of intervals, so a 
#' different value may be given for each separate interval. 
#' Must be the same length as p and r. 

#' @param r Instantaneous sampling rate of taxa, assumed to be
#' \emph{per interval}, or equal to the product of interval lengths
#' and the rates per lineage time units for each interval. Given as
#' a vector with length equal to the number of intervals, so a 
#' different value may be given for each separate interval. 
#' Must be the same length as p and q. 

#' @param PA_n The probability of sampling a taxon after the last interval 
#' included in a survivorship study. Usually zero for extinct groups, 
#' although more logically has the value of 1 when there are still extant
#' taxa (i.e., if the last interval is the Holocene and the group is
#' still alive, the probability of sampling them later is probably 1...).
#' Should be a value of 0 to 1.

#' @param PB_1 The probability of sampling a taxon before the first interval 
#' included in a survivorship study. Should be a value of 0 to 1.

#' @return Returns a matrix with number of rows equal to the number of intervals 
#' (i.e. the length of p, q and r) and named columns representing the different
#' values calculated by the function: "Nb", "Nbt", "NbL", "NFt", "NFL", "PD_bt", 
#' "PD_bL", "PD_Ft", "PD_FL", "PA", "PB", "Xbt", "XbL", "XFt" and "XFL".

#' @author David W. Bapst, with advice from Michael Foote.

#' @references
#' Foote, M. 2001. Inferring temporal patterns of preservation, origination, and 
#' extinction from taxonomic survivorship analysis. \emph{Paleobiology} 27(4):602-630.
#'
#' Foote, M. 2003a. Origination and Extinction through the Phanerozoic: A New
#' Approach. \emph{The Journal of Geology} 111(2):125-148.
#'
#' Foote, M. 2003b. Erratum: Origination and Extinction through the Phanerozoic:
#' a New Approach. \emph{The Journal of Geology} 111(6):752-753.
#'
#' Foote, M. 2005. Pulsed origination and extinction in the marine realm.
#' \emph{Paleobiology} 31(1):6-20.

#' @examples
#' #very simple example with three intervals, same value for all parameters
#' 
#' #example rates (for the most part)
#' rate<-rep(0.1,3)                  
#' #all continuous
#' footeValues(rate,rate,rate)	
#' #origination pulsed
#' footeValues(rate,rate,rate,p_cont=FALSE)		 
#' #extinction pulsed
#' footeValues(rate,rate,rate,q_cont=FALSE) 	 
#' #all pulsed
#' footeValues(rate,rate,rate,p_cont=FALSE,q_cont=FALSE) 

#' @export
footeValues<-function(p,q,r,PA_n=0,PB_1=0,p_cont=TRUE,q_cont=TRUE,Nb=1){
	#set '03 Table Values
	# TEST p, q, r
	if (length(p)!=length(q)){stop("p is not same length as q!")}
	if (length(p)!=length(r)){stop("p is not same length as r!")}
	if (length(r)!=length(q)){stop("q is not same length as r!")}
	#following relies on separate p, q, r of all intervals
	#assumes interval length is 1; i.e. the rates have been rescaled accordingly
	n<-length(p)
	#Nbt
	if(q_cont){
		Nbt<-Nb*exp(-q)
	}else{
		Nbt<-Nb*(1-q)
		}
	#NbL
	if(q_cont){
		NbL<-Nb*(1-exp(-q))
	}else{
		NbL<-Nb*q
		}
	#NFt
	if(p_cont){
		if(q_cont){
			NFt<-Nb*exp(p-q)*(1-exp(-p))
		}else{
			NFt<-Nb*exp(p)*(1-q)*(1-exp(-p))
			}
	}else{
		if(q_cont){
			NFt<-Nb*p*exp(-q)
		}else{
			NFt<-Nb*p*(1-q)
			}
		}
	#NFL
	NFL<-numeric()
	if(p_cont){
		if(q_cont){
			for(i in 1:n){
				if(p[i]==q[i]){
					NFL[i]<-Nb*(exp(-q[i])+p[i]-1)
				}else{
					NFL[i]<-Nb*(((q[i]*exp(p[i]-q[i]))+((p[i]-q[i])*exp(-q[i]))-p[i])/(p[i]-q[i]))
					}
				}
		}else{
			NFL<-Nb*q*(exp(p)-1)
			}
	}else{
		if(q_cont){
			NFL<-Nb*p*(1-exp(-q))
		}else{
			NFL<-Nb*p*q
			}
		}
	#PD_bt
	PD_bt<-1-exp(-r)
	#PD_bL
	if(q_cont){
		PD_bL<-(((r+q*exp(-(q+r)))/(q+r))-exp(-q))/(1-exp(-q))
	}else{
		PD_bL<-1-exp(-r)
			}
	#PD_Ft
	if(p_cont){
		PD_Ft<-(((r+p*exp(-(p+r)))/(p+r))-exp(-p))/(1-exp(-p))
	}else{
		PD_Ft<-1-exp(-r)
		}
	#PD_FL
	#these are corrected equations which do not agree with the erratum from Foote (2003b)
		#mostly, but not entirely, from personal communications with Mike Foote
	if(p_cont){
		if(q_cont){
			PD_FL<-numeric()
			for(i in 1:n){
				if(p[i]==q[i]){
					PD_FL[i]<-(Nb*p[i]/NFL[i])*((r[i]/(p[i]+r[i]))-((1-exp(-p[i]))/p[i])
						+(p[i]*(1-exp(-(p[i]+r[i])))/((p[i]+r[i])^2)))
				}else{
					PD_FL[i]<-(Nb/NFL[i])*(((p[i]*r[i]*(exp(p[i]-q[i])-1))/((q[i]+r[i])*(p[i]-q[i])))
						+((p[i]*q[i]*exp(-(q[i]+r[i]))*(exp(p[i]+r[i])-1))/((p[i]+r[i])*(q[i]+r[i])))
						-(exp(-q[i])*(exp(p[i])-1)))
					}
				}
		}else{
			PD_FL<-(-((p*(-exp(-r))-(exp(p)*r)+p+r)/((exp(p)-1)*(p+r))))
			}
	}else{
		if(q_cont){
			PD_FL<-(-((q*(-exp(-r))-(exp(q)*r)+q+r)/((exp(q)-1)*(q+r))))
		}else{
			PD_FL<-1-exp(-r)
			}
		}
	#PA for a given i
	PA<-numeric()
	PA[n]<-PA_n
	for(i in 1:(n-1)){
		if(q_cont){
			PA_i<-0
			for(k in (i+1):n){
				if((i+1)<=(k-1)){m<-(i+1):(k-1)}else{m<-NA}
				PA_i<-PA_i+ifelse(is.na(m[1]),1,exp(-sum(q[m])))*(1-exp(-q[k]))*(1-(
					ifelse(is.na(m[1]),1,exp(-sum(r[m])))*(1-PD_bL[k])))
				}
			PA[i]<-PA_i+(exp(-sum(q[(i+1):n])))*(1-exp(-sum(r[(i+1):n]))*(1-PA[n]))
		}else{
			PA_i<-0
			for(k in (i+1):n){
				if((i+1)<=(k-1)){m<-(i+1):(k-1)}else{m<-NA}
				PA_i<-PA_i+ifelse(is.na(m[1]),1,prod(1-q[m]))*q[k]*(1-(
					ifelse(is.na(m[1]),1,exp(-sum(r[m])))*(1-PD_bL[k])))
				}
			PA[i]<-PA_i+(prod(1-q[(i+1):n])*(1-exp(-sum(r[(i+1):n]))*(1-PA[n])))
			}
		}
	#PB
	PB<-numeric()
	PB[1]<-PB_1
	for(i in 2:n){
		if(p_cont){
			PB_i<-0
			for(k in 1:(i-1)){
				if((k+1)<=(i-1)){m<-(k+1):(i-1)}else{m<-NA}
				PB_i<-PB_i+((ifelse(is.na(m[1]),1,exp(-sum(p[m])))*(1-exp(-p[k])))*(1-
					ifelse(is.na(m[1]),1,exp(-sum(r[m])))*(1-PD_Ft[k])))
				}
			PB[i]<-PB_i+(exp(-sum(p[1:(i-1)])))*(1-exp(-sum(r[1:(i-1)]))*(1-PB[1]))
		}else{
			PB_i<-0
			for(k in 1:(i-1)){
				if((k+1)<=(i-1)){m<-(k+1):(i-1)}else{m<-NA}
				PB_i<-PB_i+((ifelse(is.na(m[1]),1,prod(1/(1+p[m])))*(p[k]/(1+p[k])))*(1-
					ifelse(is.na(m[1]),1,exp(-sum(r[m])))*(1-PD_Ft[k])))
				}
			PB[i]<-PB_i+(prod(1/(1+p[1:(i-1)]))*(1-exp(-sum(r[1:(i-1)]))*(1-PB[1])))
			}
		}
	#Xbt
	Xbt<-Nbt*PB*PA
	#XbL
	XbL<-NbL*PB*PD_bL+Nbt*PB*PD_bt*(1-PA)
	#XFt
	XFt<-NFt*PA*PD_Ft+Nbt*PA*PD_bt*(1-PB)	
	#XFL
	XFL<-(NFL*PD_FL)+(NbL*(1-PB)*PD_bL)+(NFt*(1-PA)*PD_Ft)+(Nbt*(1-PB)*PD_bt*(1-PA))
	res<-data.frame(cbind(Nb=rep(Nb,n),Nbt,NbL,NFt,NFL,PD_bt,PD_bL,PD_Ft,PD_FL,PA,PB,Xbt,XbL,XFt,XFL))
	return(res)
	}