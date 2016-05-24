#' Wetland model.
#' 
#' The wetland model function.
#'
#' @param t time
#' @param init init
#' @param parameters model parameters
#' @param nr number of raster map rows
#' @param nc number of raster map columns
#'
#' @return Model solver
#'
#' @keywords solver
#'
#' @export
#' 
#' @examples
#' ## Not run spdynmod()

##if (require(deSolve) == F) {install.packages('deSolve',repos='http://cran.r-project.org');if (require(deSolve) == F) print ('Error: deSolve package is not installed on your machine')}

############# MODEL
spdynmod<-function(t,init,parameters,nr,nc) { 

rpath = paste(find.package('spdynmod'),'/extdata',sep='')

nr<-get('nr')
nc<-get('nc')
NN<-get('NN')

#nc<<-NULL
#nr<<-NULL
#NN<-NULL

#r<- raster(paste(rpath,'/mc84_1_reclass.asc',sep=''))
#nr<<-dim(r)[1]
#nc<<-dim(r)[2]

fak<-NULL
fak<- raster::raster(paste(rpath,'/log_cr10_acum_rm_t1_aver.asc',sep=''))

avd<-NULL
avd<- raster::raster(paste(rpath,'/ramblas_cr10_dist_t1_ave.asc',sep=''))

fa_avd<-NULL
fa_avd<-fak+(1-avd)

dr1<-NULL
dr2<-NULL
dr1<- raster::raster(paste(rpath,'/rambla11_cr10_dist_t1.asc',sep='')) 
dr2<- raster::raster(paste(rpath,'/rambla22_cr10_dist_t1.asc',sep='')) 


	#NN<<-nr*nc
	sm <- matrix(nrow=nr,ncol=nc,init[1:NN])
	es <- matrix(nrow=nr,ncol=nc,init[(NN+1):(2*NN)])
	rb <- matrix(nrow=nr,ncol=nc,init[((2*NN)+1):(3*NN)])
	baresoil <- matrix(nrow=nr,ncol=nc,init[((3*NN)+1):(4*NN)])
	drambla1<-raster::as.matrix(dr1) 
	drambla2<-raster::as.matrix(dr2) 
	fa_avd2<-raster::as.matrix(fa_avd)

	Time <- t
	tprb <- parameters['tprb']
	tpsm <- parameters['tpsm']

	IPRH <- inputData(t, 'IPRH')
	eframbla1 <- 1-drambla1
	eframbla2 <- 1-drambla2

	daguarb <- (eframbla1)^5 + (eframbla2)^5
	daguasm <- fa_avd2^4 

	edarb <- (1+IPRH)*daguarb 
	edasm <- (1+IPRH)*daguasm 

	tpss2sm <- IPRH/2
	tpbs2sm <- IPRH/0.2
	tpss2rb <- IPRH/0.4540899989
	tpbs2rb <- IPRH/0.1
	tpsm2rb <- IPRH/3

	tactsm2rb <- edarb*tpsm2rb*tprb
	tactss2rb <- edarb*tpss2rb*tprb
	tactbs2sm <- edasm*tpbs2sm*tpsm
	tactss2sm <- edasm*tpss2sm*tpsm
	tactbs2rb <- edarb*tpbs2rb*tprb

	ss2rb <-  rb * tactss2rb * es * ( 1 - ( rb / 25 ) )
	ss2sm <-  sm * tactss2sm * es * ( 1 - ( sm / 25 ) )
	sm2rb <-  rb * tactsm2rb * sm * ( 1 - ( rb / 25 ) )
	bs2sm <-  sm * tactbs2sm * baresoil * ( 1 - ( sm / 25 ) )
	bs2rb <-  rb * tactbs2rb * baresoil * ( 1 - ( rb / 25 ) )

	 dsm = ss2sm  + bs2sm  - sm2rb 
	 des = - ss2sm  - ss2rb 
	 drb = sm2rb  + ss2rb  + bs2rb 
	 dbaresoil =  - bs2sm  - bs2rb 

######### BEGIN DISPERSAL
#sm_disp<-0
		w22<-which(sm > 1)
#print(length(w22))
		try(if (length(w22)>0) { # ERROR to be checked with random initial state variables!
			wg22<-sapply(w22,neigh_cell)
#print(wg22)
				disp222bs<-sapply(wg22,function(x) x[which(baresoil[x]>=1)]) 
				which(disp222bs>0)->ind
				as.numeric(disp222bs[ind])->xs

				if (length(ind)>0) {
						dbaresoil[xs]<-dbaresoil[xs]-0.25
						dsm[xs]<-dsm[xs]+0.25
						#sm_disp<-3


				} 

				disp222ss<-sapply(wg22,function(x) x[which(es[x]>=1)]) 
				which(disp222ss>0)->ind
				as.numeric(disp222ss[ind])->xs

				if (length(ind)>0) {
						des[xs]<-des[xs]-0.25
						dsm[xs]<-dsm[xs]+0.25
						#sm_disp<-4

	} 
})
######### END DISPERSAL

	 list(c(dsm,des,drb,dbaresoil))
}

#system('aplay -t wav hangout_dingtone.wav')
