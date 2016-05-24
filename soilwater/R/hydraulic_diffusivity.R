NULL
#'
#' Soil Water Retention Curve 'swc', Hydraulic Conductivity 'khy' , Soil Water Capacity 'cap' , Soil Water (Hydraulic) Diffusivity 'diffusivity'
#' 
#' @title Soil water Retantion Curve and Unsaturated Hydraulic Conductivity
#' 
#' 
#' @param psi soil wwater pressure head 
#' @param alpha inverse of a length - scale parameters in Van Genuchten Formula
#' @param n shape parameter in Van Genuchten Formula
#' @param m shape parameter in Van Genuchten Formula. Default is \code{1-1/n}  
#' @param theta_sat saturated water content
#' @param theta_res residual water content
#' @param ksat saturated hydraulic conductivity
#' @param lambda,b lambda and b exponents in Brook and Corey formula. It is used in case \code{type_swc} and/or \code{type_khy} are equal to \code{BrooksAndCorey}.
#' @param psi_s psi_s value (capillary fringe) in Brook and Corey formula. It is used in case \code{type_swc} and/or \code{type_khy} are equal to \code{BrooksAndCorey}.
#' @param v exponent in Mualem Formula for Hydraulic Conductivity
#' @param saturation_index logical index, If \code{TRUE} (Default) the function \code{swc()} returns soil water content, otherwise a saturation index between 0 and 1.
#' @param type_swc type of Soil Water Retention Curve. Default is \code{"VanGenuchten"} and actually the only implemented type 
#' @param type_khy type of Soil Hydraulic Conductivity Curve. Default is \code{"Mualem"} and actually the only implemented type
#' @param ... further arguments which are passed to \code{swc()} and \code{khy()}
#' 
#' 
#' @rdname hydraulic_diffusivity
#' @export
#' @examples
#' 
#' library(soilwater)
#' soiltype <- c("sand","silty-sand","loam","clay")
#' theta_sat <- c(0.44,0.39,0.51,0.48)
#' theta_res <- c(0.02,0.155,0.04,0.10)
#' alpha <- c(13.8,6.88,9.0,2.7) # 1/meters
#' n <- c(2.09,1.881,1.42,1.29) 
#' m <- 1-1/n
#' v <- array(0.5,length(soiltype))
#' ks <- c(1.5e-1,1e-4*3600,3.3e-2,4.1e-4)/3600   # meters/seconds
#' 
#' psi <- -(1:2000)/1000
#' 
#' D <- as.data.frame(array(0.1,c(length(psi),length(soiltype))))
#' names(D) <- soiltype

#' for (it in names(D)) {
#'   
#'   i=which(names(D)==it)
#'   D[,i] <- diffusivity(psi=psi,
#'             v=v[i],ksat=ks[i],alpha=alpha[i],
#'             n=n[i],m=m[i],theta_sat=theta_sat[i],
#'             theta_res=theta_res[i])
#'  
#' }

#'# plot diffusivity on log scale 
#' lty <- 1:length(names(D) )
#' 
#' plot(psi,D[,1],lty=lty[1],main="Diffusvity vs psi",xlab="psi [m]",
#' ylab="D [m^2/s]",type="l",ylim=range(D),ylog=TRUE)
#' for (i in 2:ncol(D)) {
#'   lines(psi,D[,i],lty=lty[i]) 
#' }
#' legend("topleft",lty=lty,legend=names(D))
#'Dinv <- 1/D 
#' 
#' # pot diffusivity on log scale 
#' lty <- 1:length(names(D) )
#' 
#' plot(psi,Dinv[,1],lty=lty[1],main="1/Diffusvity vs psi",
#' xlab="psi [m]",ylab="1/D [s/m^2]",type="l",ylim=range(Dinv),ylog=TRUE)
#' for (i in 2:ncol(Dinv)) {
#'   lines(psi,Dinv[,i],lty=lty[i]) 
#' }
#' legend("topright",lty=lty,legend=names(D))




# script to genereta plots of soil water hydralic diffusivity 

#  soil water retention curves 

swc <- function (psi=0.5,alpha=1.0,n=1.5,m=1-1/n,theta_sat=0.4,theta_res=0.05,psi_s=-1/alpha,lambda=m*n,saturation_index=FALSE,type_swc=c("VanGenuchten","BrooksAndCorey"),...) {
  
	if (length(type_swc)>0) type_swc <- type_swc[1]
	
	if (type_swc=="VanGenuchten"){
		
		psiloc <- psi
		psi[psi>0] <- 0.0
		
		sat_index <- (1+(-alpha*psi)^n)^(-m)
		
	}else if (type_swc=="BrooksAndCorey"){
		
		psiloc <- psi 
		psi[psi>psi_s] <- psi_s
		psirel <- psi/psi_s
		sat <- psirel^(-1/lambda)
		
	} else {
		return(NA)
	}
	theta <- sat_index*(theta_sat-theta_res)+theta_res
	out <- theta
	if (saturation_index) out<- sat_index
	
	return(out)
  
}

# hydraulic conductivity 
NULL
#'
#'
#' @rdname hydraulic_diffusivity
#' @export

#'  

khy <- function (psi=0.5,v=0.5,ksat=0.01,alpha=1.0,n=1.5,m=1-1/n,theta_sat=0.4,theta_res=0.05,psi_s=-1/alpha,lambda=m*n,b=NA,type_swc="VanGenuchten",type_khy=c("Mualem","BrooksAndCorey"),...) {
  
  if (length(type_khy)>1) type_khy <- type_khy[1]
	
  if (type_khy=="Mualem") {
  s <- soilwater::swc(psi=psi,alpha=alpha,m=m,n=n,theta_sat=theta_sat,theta_res=theta_res,psi_s=psi_s,lambda=lambda,type_swc=type_swc,saturation_index=TRUE,...)
  
  k <- ksat*s^v*(1-(1-s^(1/m))^m)^2
}else if (type_khy=="BrooksAndCorey") {
	
	s <- soilwater::swc(psi=psi,alpha=alpha,m=m,n=n,theta_sat=theta_sat,theta_res=theta_res,psi_s=psi_s,lambda=lambda,type_swc=type_swc,saturation_index=TRUE,...)
	k <- ksat*s^b
	
} else {	
	k<- NA
}
  return(k)
}
NULL
#'
#' @rdname hydraulic_diffusivity
#' @export
#'  

# hydraulic capacity 
# it is calculated by dariving swc 

cap <- function (psi=0.5,alpha=1.0,n=1.5,m=1-1/n,theta_sat=0.4,theta_res=0.05,type_swc="VanGenuchten",...) {
 
	if (type_swc!="VanGenuchten") return(NA)
	
	psiloc <- psi
    psi[psi>0] <- 0.0
    
    out <- alpha*(theta_sat-theta_res)*(1+(-alpha*psi)^n)^(-m-1)*m*n*(-alpha*psi)^(n-1)
    
    return(out)
  
}


NULL
#'
#' @rdname hydraulic_diffusivity
#' @export
#'  
# hydraulic diffusivity 
diffusivity <- function (psi=0.5,v=0.5,ksat=0.01,alpha=1.0,n=1.5,m=1-1/n,theta_sat=0.4,theta_res=0.05,...) {
  
  k <- soilwater::khy(psi=psi,v=v,ksat=ksat,alpha=alpha,n=n,m=m,theta_sat=theta_sat,theta_res=theta_res,...)
  C <- soilwater::cap(psi=psi,alpha=alpha,n=n,theta_sat=theta_sat,theta_res=theta_res,...)
  
  return(k/C) 
  
}

# soil used 
