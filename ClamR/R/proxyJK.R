proxyJK <- function(x, y, dx)
  {
###############   dx is one whole period, the
    ########   actual increment will be half of dx
      PLOT = FALSE
    JOUT = list()
    pmids = vector()
    omids = vector()
    delw = vector()
    xout = seq(from = min(x), to = max(x), length = 4 * length(x))
    hxy = smooth.spline(x, y)
    predy = predict(hxy, xout)
    for (i in 1:length(x)) {
        t1 = x[i]
        t2 = t1 + dx/2
        if (t2 > max(x)) 
		break
        flag = x >= t1 & x <= t2
        ax = x[flag]
        ay = y[flag]
        zx = predy$x[predy$x >= t1 & predy$x <= t2]
        zy = predy$y[predy$x >= t1 & predy$x <= t2]
        fzed = fft(zy)
        deltax = zx[2] - zx[1]
        N = length(zy)
        m = floor(N/2) + 1
        f = seq(from = 0, to = 0.5, length = m) * (1/(deltax))
        Phs = Arg(fzed[2])
        Bphs = sqrt(Mod(fzed[2]))
        ef = f[2]
        dC2 = ay
        myEx = ax
        mid = mean(ax)
        Pos = mean(ay)
        Amp = diff(range(ay))/2
        Prd = 1/ef
        xin = c(Phs, Pos, Amp, Prd)
        FOUT = proxyA(ax, ay, xin)
        PhsALL = FOUT$par[1]
        PosALL = FOUT$par[2]
        AmpALL = FOUT$par[3]
        PrdALL = FOUT$par[4]
        wimidALL = (AmpALL/2) * sin((mid - PhsALL) * 2 * pi/PrdALL) + 
		PosALL
        ALLPAR = c(PhsALL, PosALL, AmpALL, PrdALL, wimidALL)
        wi = (AmpALL/2) * sin((myEx - PhsALL) * 2 * pi/PrdALL) + 
		PosALL
        if (PLOT == TRUE) {
            plot(ax, ay, main = i)
            lines(ax, wi, col = rgb(0.5, 0, 0), lwd = 2)
        }
        Kjack = length(ax)
        PhsPseudo = vector()
        PosPseudo = vector()
        AmpPseudo = vector()
        PrdPseudo = vector()
        wimidPseudo = vector()
        PhsTILDE = vector()
        PosTILDE = vector()
        AmpTILDE = vector()
        PrdTILDE = vector()
        wimidTILDE = vector()
        for (MM in 1:Kjack) {
            bx = ax[-MM]
            by = ay[-MM]
            KOUT = proxyA(bx, by, xin)
            PhsTEMP = KOUT$par[1]
            PosTEMP = KOUT$par[2]
            AmpTEMP = KOUT$par[3]
            PrdTEMP = KOUT$par[4]
            temmwimid = (AmpTEMP/2) * sin((mid - PhsTEMP) * 2 * 
										  pi/PrdTEMP) + PosTEMP
            PhsTILDE[MM] = PhsTEMP
            PosTILDE[MM] = PosTEMP
            AmpTILDE[MM] = AmpTEMP
            PrdTILDE[MM] = PrdTEMP
            wimidTILDE[MM] = temmwimid
            PhsPseudo[MM] = Kjack * PhsALL - (Kjack - 1) * PhsTEMP
            PosPseudo[MM] = Kjack * PosALL - (Kjack - 1) * PosTEMP
            AmpPseudo[MM] = Kjack * AmpALL - (Kjack - 1) * AmpTEMP
            PrdPseudo[MM] = Kjack * PrdALL - (Kjack - 1) * PrdTEMP
            wimidPseudo[MM] = Kjack * wimidALL - (Kjack - 1) * 
			temmwimid
        }
        PSTILDE = cbind(PhsTILDE, PosTILDE, AmpTILDE, PrdTILDE, 
						wimidTILDE, PhsPseudo, PosPseudo, AmpPseudo, PrdPseudo, 
						wimidPseudo)
        if (FALSE) {
            for (MM in 1:Kjack) {
                PhsTEMP = PhsTILDE[MM]
                PosTEMP = PosTILDE[MM]
                AmpTEMP = AmpTILDE[MM]
                PrdTEMP = PrdTILDE[MM]
                wimidTEMP = (AmpTEMP/2) * sin((mid - PhsTEMP) * 
											  2 * pi/PrdTEMP) + PosTEMP
                wi = (AmpTEMP/2) * sin((ax - PhsTEMP) * 2 * pi/PrdTEMP) + 
				PosTEMP
                lines(ax, wi, col = rgb(1, 0.7, 0.7))
                px = ax[MM]
                py = wi[MM]
                text(px, py, labels = MM, xpd = TRUE)
            }
        }
        if (FALSE) {
            plot(ax, ay)
            for (MM in 1:Kjack) {
                PhsTEMP = PhsPseudo[MM]
                PosTEMP = PosPseudo[MM]
                AmpTEMP = AmpPseudo[MM]
                PrdTEMP = PrdPseudo[MM]
                wimidTEMP = (AmpTEMP/2) * sin((mid - PhsTEMP) * 
											  2 * pi/PrdTEMP) + PosTEMP
                wi = (AmpTEMP/2) * sin((ax - PhsTEMP) * 2 * pi/PrdTEMP) + 
				PosTEMP
                lines(ax, wi, col = "red")
            }
            for (MM in 1:Kjack) {
                plot(ax, ay)
                bx = ax[-MM]
                by = ay[-MM]
                points(bx, by, col = "red")
                TOUT = proxyA(bx, by, xin)
                PhsTEMP = TOUT$par[1]
                PosTEMP = TOUT$par[2]
                AmpTEMP = TOUT$par[3]
                PrdTEMP = TOUT$par[4]
                wimidTEMP = (AmpTEMP/2) * sin((mid - PhsTEMP) * 
											  2 * pi/PrdTEMP) + PosTEMP
                wi = (AmpTEMP/2) * sin((ax - PhsTEMP) * 2 * pi/PrdTEMP) + 
				PosTEMP
                lines(ax, wi, col = "red")
                locator(1)
            }
        }
        Phs = mean(PhsPseudo)
        Pos = mean(PosPseudo)
        Amp = mean(AmpPseudo)
        Prd = mean(PrdPseudo)
        wimid = mean(wimidPseudo)
        nnm1 = Kjack * (Kjack - 1)
        PhsV = sum((PhsPseudo - Phs)^2)/nnm1
        PosV = sum((PosPseudo - Pos)^2)/nnm1
        AmpV = sum((AmpPseudo - Amp)^2)/nnm1
        PrdV = sum((PrdPseudo - Prd)^2)/nnm1
        wimidV = sum((wimidPseudo - wimid)^2)/nnm1
        wi = (Amp/2) * sin((myEx - Phs) * 2 * pi/Prd) + Pos
        if (PLOT == TRUE) {
            lines(myEx, wi, col = "purple")
            lines(zx, zy, lty = 2, col = "green")
            locator(1)
        }
        wimid = (Amp/2) * sin((mid - Phs) * 2 * pi/Prd) + Pos
        JKest = c(Phs, Pos, Amp, Prd, wimid)
        JKvar = c(PhsV, PosV, AmpV, PrdV, wimidV)
        omids[i] = mid
        pmids[i] = wimidALL
        tval = qt(0.95, (Kjack - 1))
        delw[i] = tval * sqrt(wimidV)
        JOUT[[i]] = list(par = ALLPAR, mid = mid, ax = ax, predmid = wimid, 
						 JKest = JKest, JKvar = JKvar, PSTILDE = PSTILDE)
    }
    return(list(JOUT = JOUT, omids = omids, pmids = pmids, delw = delw, x=x, y=y, predy=predy))
}
