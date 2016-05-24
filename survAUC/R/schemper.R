###########################################################
###             Schemper-Henderson estimator
###########################################################

### Programm fuer Cox-Modell (aus Lusa et al. 2007)
### liefert Objekt mit Zeitpunkten wie survfit ("timep")
### und prediction error curve ("Mhat")
### -> mit seval, etc. weiterverwertbar

### train.fit = fit von cph()



schemper <- function(train.fit, traindata, newdata)
{
	if(!inherits(train.fit,"rms"))
		stop("\nThe Cox model has to be estimated via the cph function of the rms package.\n")
    f.Mt <- function(tempo, tutti.tempi, stima.surv, tempi.evento,
					 Stj, ind.censura, num.sogg)
		{
			Stj1 <- unique(Stj[tempi.evento == tempo])
			primo <- rep(1 - Stj1, num.sogg)
			primo[tutti.tempi <= tempo] <- 0
			secondo <- Stj1 * (1 - ind.censura)
			secondo[tutti.tempi > tempo] <- 0
			terzo <- ind.censura * (((1 - Stj1) * Stj1)/stima.surv + Stj1 * (1 - Stj1/stima.surv))
			terzo[tutti.tempi > tempo] <- 0
			terzo[is.na(terzo)] <- 0
			ris <- primo + secondo + terzo
			return(sum(ris)/num.sogg)
		}
    f.Mt.cox <- function(tempo, tutti.tempi, stima.surv, tempi.evento,
						 Stj0, ind.censura, num.sogg, lin.pred)
		{
			Stj00 <- unique(Stj0[tempi.evento == tempo])
			Stj1 <- Stj00^exp(lin.pred)
			primo <- 1 - Stj1
			primo[tutti.tempi <= tempo] <- 0
			secondo <- Stj1 * (1 - ind.censura)
			secondo[tutti.tempi > tempo] <- 0
			terzo <- ind.censura * (((1 - Stj1) * Stj1)/stima.surv + Stj1 * (1 - Stj1/stima.surv))
			terzo[tutti.tempi > tempo] <- 0
			terzo[is.na(terzo)] <- 0
			ris <- primo + secondo + terzo
			return(sum(ris)/num.sogg)
		}
    f.assegna.surv <- function(tempo, tempi.eventi)
		{
			if (any(tempo == tempi.eventi)) {
				pos <- (c(1:length(tempi.eventi)) * as.numeric(tempo == tempi.eventi))
				pos <- pos[pos != 0]
			}
			else {
				tmp <- (tempo - tempi.eventi)
				if (all(tmp < 0))
					pos <- NA
				else {
					tmp[tmp < 0] <- Inf
					pos <- order(tmp)[1]
				}
			}
			return(pos)
		}
    tsurv <- as.numeric(newdata$time)
    surv  <- as.numeric(newdata$status)
    lin.pred <- predict(train.fit, newdata,"lp")
    num.sogg <- length(tsurv)
    km <- survfit(Surv(tsurv, surv) ~ 1)
    km.fit <- survfit(Surv(time, status) ~ 1, data=traindata)
    tempi.eventi <- km$time[km$n.event != 0]
    pos.surv <- apply(as.matrix(tsurv), 1, f.assegna.surv, tempi.eventi)
    surv.tj <- approx(km.fit$time,km.fit$surv,xout =tempi.eventi , method = "constant", f = 0, yleft=1,yright=min(km.fit$surv, na.rm=T))$y
    surv.tot.km <- (surv.tj)[pos.surv]
    ind.censura <- as.numeric(!as.logical(surv))
    Mt <- apply(as.matrix(tempi.eventi), 1, f.Mt, tsurv, surv.tot.km,
				tempi.eventi, surv.tj, ind.censura, num.sogg)
    numero.eventi <- km$n.event[km$n.event != 0]
    surv0.tj.cox <- approx(train.fit$time,train.fit$surv,xout =tempi.eventi , method = "constant", f = 0, yleft=1,yright=min(km.fit$surv, na.rm=T))$y
    surv0.tot.cox <- (surv0.tj.cox)[pos.surv]
    surv.tot.cox <- surv0.tot.cox^exp(lin.pred)
    Mtx <- apply(as.matrix(tempi.eventi), 1, f.Mt.cox, tsurv,
				 surv.tot.cox, tempi.eventi, surv0.tj.cox, ind.censura,
				 num.sogg, lin.pred)
    Gkm <- survfit(Surv(tsurv, ind.censura) ~ 1)
    tempi.censure <- Gkm$time[Gkm$n.event != 0]
    if (!length(tempi.censure))
	cens.tot.km <- rep(1, length(tempi.eventi))
    else {
        pos.surv.censure <- apply(as.matrix(tempi.eventi), 1,
								  f.assegna.surv, tempi.censure)
        cens.tot.km <- (Gkm$surv[Gkm$n.event != 0])[pos.surv.censure]
        cens.tot.km[tempi.eventi < min(Gkm$time[Gkm$n.event !=
									   0])] <- 1
    }
    pesi <- numero.eventi/cens.tot.km
    peso.tot <- sum(pesi)
    D <- sum(Mt * pesi)/peso.tot
    Dx <- sum(Mtx * pesi)/peso.tot
    V <- (D - Dx)/D
    return(list(Model = train.fit$call, D = D, Dx = Dx, V = V,  Mhat=Mtx, Mhat.0=Mt, timep=tempi.eventi))
}
