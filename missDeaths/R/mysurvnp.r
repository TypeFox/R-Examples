################
# observable survival
########################
my.survnp <- function(tt,st,R,ratetable,konc,conf.int=0.95){
	#estimates event-free (EFS) and net survival
	#tt ... cas
	#st ... censoring indicator (1=dogodek, nic sicer)
	#R ... spremenljivke za populacijsko tabelo
	#ratetable ... uporabljena populacijska tabela
	#konc ... maximum potential observation time for this individual
	if(missing(konc)) konc <- rep(max(tt),length(tt))				#ce ni podan koncni datum krnjenja - das  vsem na zadnjega opazenega
	
	if(any((tt>konc)&(st==1)))stop("Events are observed after administrative censoring \n")		#ce imamo za nekoga podan konc, ne more imeti dogodka po tem
	if(any((tt!=konc)&(st==0)))stop("Only administrative censoring is allowed - max times should be reported for all individuals \n")			#za nekoga imamo cas krnjenja, ki ni enak koncu
	
	data <- data.frame(time=round(tt,8),cens=st,konc=konc)
	data <- data[order(tt),]						#uredim po casih
	R <- R[order(tt),]							#uredim po casih		
	#ttun <- sort(unique(c(data$time,data$konc)))				#unique casi dogodkov ali censoringa
	ttun <- sort(unique(c(data$time,data$konc,1:max(data$time,data$konc))))	#unique casi dogodkov ali censoringa - za net je vazno, da jih je cimvec, za ostale je vseeno

	fk = (attr(ratetable, 'type') != 1)
	#fk <- (attributes(ratetable)$factor != 1)				#spremenljivke iz popul. tabel, ki se vecajo v casu
	nfk <- length(fk)							#stevilo spr. v popul. tabelah
	
	difft <- c(ttun[1],diff(ttun))						#sprememba casa na posameznem intervalu

	le.net <- vari.net <- rep(NA,length(ttun))				#hazard, stevilo ljudi ob vsakem casu in prirastek variance
	le.efs <- vari.efs <- rep(NA,length(ttun))				#hazard, stevilo ljudi ob vsakem casu in prirastek variance
	Yosebni <- rep(1,nrow(data))						#Yosebni je verjetnost, da je se ziv v populaciji
	Yosebnin <- Spicum  <- rep(1,nrow(data))				#Spicum potrebujemo za net survival, na zacetku 1 za vse, Yosebnin je za net

	for(it in 1:length(ttun)){
		#population survival probabilities
		Spi <- relsurv:::exp.prep(R,difft[it],ratetable)			#populacijska verjetnost, da zivijo se naprej
		Spicum <- Spicum*Spi						#unconditional (od zacetka do konca), potrebujem za net survival
		dLpi <- -log(Spi)
	
		#at risk indicator
		Yosebni <- Yosebni*Spi						#Yosebni se zmanjsa - nekateri umrejo
		
		#number of events
		Nei <- data$time==ttun[it]&data$cens==1				#individual event indicator
		who.event <- which(Nei==1)					#kdo ima dogodek ob tem casu
		Ne <- length(who.event)						#stevilo dogodkov ob tem casu
		Npi <- Yosebni*dLpi
		Np <- sum(Npi)

		#net number of events
		Neisp <- Nei/Spicum
		Nesp <- sum(Neisp)						#stevilo dogodkov ob tem casu v hyp world
		
		#hazard - net
		Yoni <- sum(Yosebni/Spicum)					#number at risk in hyp world
		le.net[it] <- Nesp/Yoni						#trenutno tveganje
		vari.net[it] <- sum(Nei/Spicum^2)/Yoni^2			#prirastek k varianci

		#hazard - event free
		le.efs[it] <- (Ne+Np)/sum(Yosebni)				#trenutno tveganje
		vari.efs[it] <- sum((Nei + Npi)^2)/sum(Yosebni)^2		#prirastek k varianci

		#preparing for the next time interval
		who.cens <- which(data$konc==ttun[it])				#kdo je krnjen oz. konca opazovanje ob tem casu
		Yosebni[who.event] <- Yosebni[who.event]-1			#odstejem 1 pri tistemu, ki ima dogodek
		Yosebni[who.cens] <- 0						#dokoncno vrzem ven vse, ki imajo konc (tudi ce je imel prej event)			
		R <- R+  matrix(difft[it] * t(fk),ncol=3,nrow=nrow(R),byrow=T)	#povecam starost in koledarski cas za vse
		
	}
	se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval

	#NET
	vari.net <- c(0,cumsum(vari.net))
	surv.net <- exp(-c(0,cumsum(le.net)))
	
	#log scale conf. interval
	lower.net <- exp(as.vector(log(surv.net) - sqrt(vari.net)*se.fac))
	upper.net <- exp(as.vector(log(surv.net) + sqrt(vari.net)*se.fac))
	
	#EFS
	vari.efs <- c(0,cumsum(vari.efs))
	surv.efs <- exp(-c(0,cumsum(le.efs)))

	#log scale conf. interval
	lower.efs <- exp(as.vector(log(surv.efs) - sqrt(vari.efs)*se.fac))
	upper.efs <- exp(as.vector(log(surv.efs) + sqrt(vari.efs)*se.fac))
		
	
	out <- list(time=c(0,ttun),
	surv.net=surv.net,std.err.net=sqrt(vari.net),lower.net=lower.net,upper.net=upper.net,
	surv.efs=surv.efs,std.err.efs=sqrt(vari.efs),lower.efs=lower.efs,upper.efs=upper.efs
	)
	
	out
}