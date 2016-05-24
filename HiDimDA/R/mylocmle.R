# Modified version of the locmle function from the locfde package, that makes an extra test on lines 57 and 59 to ensure that V, mu and 
# bett[2] are composed by valid finite numbers. 
# Other than that this file is identical to the original locfns.R file as written by Bradley Efron, Brit B. Turnbull, and Balasubramanian Narasimhan.

##
## We need the splines library, so we load it.
##

#.onLoad <- function(lib, pkg) require(splines)
require(splines) # Code modification  !!

mylocmle <-
function(z, xlim , Jmle = 35, d = 0, s = 1, ep = 1/100000, sw = 0, Cov.in)
{
        ## uses z-values in [-xlim,xlim] to find mles for p0,del0,sig0
        ## Jmle number of iterations, beginning at (del0,sig0)=(d,s)
        ## sw=1 returns correlation matrix
        N = length(z)
        if (missing(xlim)) {
          if (N>500000) b = 1
          else b=4.3 * exp(-0.26*log(N,10)) 
          xlim=c(median(z),b*diff(quantile(z)[c(2,4)])/(2*qnorm(.75)))
        }
        aorig=xlim[1]-xlim[2]
        borig=xlim[1]+xlim[2]
        z0=z[which(z>=aorig & z<=borig)]
        N0 = length(z0)
        Y = c(mean(z0), mean(z0^2))
        that = N0/N
        ######## find mle estimates ###########################
        for(j in 1:Jmle) {
                bet = c(d/s^2, -1/(2 * s^2))
                aa = (aorig - d)/s
                bb = (borig - d)/s
                H0 = pnorm(bb) - pnorm(aa)
                fa = dnorm(aa)
                fb = dnorm(bb)
                H1 = fa - fb
                H2 = H0 + aa * fa - bb * fb
                H3 = (2 + aa^2) * fa - (2 + bb^2) * fb
                H4 = 3 * H0 + (3 * aa + aa^3) * fa - (3 * bb + bb^3) * fb
                H = c(H0, H1, H2, H3, H4)
                r = d/s
                I = matrix(rep(0, 25), 5)
                for(i in 0:4)
                        I[i + 1, 0:(i + 1)] = choose(i, 0:i)
                u1 = s^(0:4)
                II = pmax(row(I) - col(I), 0)
                II = r^II
                I = u1 * (I * II)
                E = as.vector(I %*% H)/H0
                E1 = E[2]
                E2 = E[3]
                E3 = E[4]
                E4 = E[5]
                mu = c(E1, E2)
                V = matrix(c(E2 - E1^2, E3 - E1 * E2, E3 - E1 * E2, E4 -
                        E2^2), 2)
                if ( !all(is.finite(mu)) || !all(is.finite(V)) ) break	# Code modification  !!
		qrV <- qr(V)						# Code modification  !!
		if (qrV$rank!=nrow(V)) break				# Code modification  !!
#                bett = bet + solve(V, Y - mu)/(1 + 1/j^2)
                bett = bet + solve(qrV, Y - mu)/(1 + 1/j^2)		# Code modification  !!
                if (!is.finite(bett[2])) break				# Code modification  !!
#               if(bett[2]>0)bett=bet+.1*solve(V, Y - mu)/(1 + 1/j^2)
                if(bett[2]>0)bett=bet+.1*solve(qrV, Y - mu)/(1 + 1/j^2) # Code modification  !!
                if (!is.finite(bett[2])) break
                else if (bett[2]>=0) break
                d =  - bett[1]/(2 * bett[2])
                s = 1/sqrt(-2 * bett[2])
                if(sum((bett - bet)^2)^0.5 < ep)
                        break
        }
#        if (is.na(bett[2])) {
        if (!all(is.finite(mu)) || !all(is.finite(V)) || qrV$rank!=nrow(V) || !is.finite(bett[2])) {   # Code modification  !!
          mle = rep(NA, 6)
          Cov.lfdr = NA
          Cor = matrix(NA, 3, 3)
        }
        else if (bett[2] >=0)  {
          mle = rep(NA, 6)
          Cov.lfdr = Cov.out = NA
          Cor = matrix(NA, 3, 3)
        }
        else {
          aa = (aorig - d)/s
          bb = (borig - d)/s
          H0 = pnorm(bb) - pnorm(aa)
          p0 = that/H0
          ############  sd calcs ###########################
          J = s^2 * matrix(c(1, 0, 2 * d, s), 2)
#          JV = J %*% solve(V)
          JV = J %*% solve(qrV)		# Code modification  !!
          JVJ = JV %*% t(J)
          mat2 = cbind(0, JVJ/N0)
          mat1 = c((p0 * H0 * (1 - p0 * H0))/N, 0, 0)
          mat = rbind(mat1, mat2)
          h = c(H1/H0, (H2 - H0)/H0)
          matt = c(1/H0,  - (p0/s) * t(h))
          matt = rbind(matt, cbind(0, diag(2)))
          C = matt %*% (mat %*% t(matt))
          mle = c(p0, d, s, diag(C)^0.5)
          if(sw == 1) {
                sd = mle[4:6]
                Co = C/outer(sd, sd)
                dimnames(Co) = list(c("p0", "d", "s"), c("p0", "d", "s"))
                Cor = Co[c(2, 3, 1), c(2, 3, 1)]
          }
          if (!missing(Cov.in)) {
            i0 = which(Cov.in$x > aa & Cov.in$x < bb)
            Cov.out = loccov(N, N0, p0, d, s, Cov.in$x, Cov.in$X, Cov.in$f,
                             JV, Y, i0, H, h, Cov.in$sw)
          }
        }
        names(mle) = c("p0", "del0", "sig0", "sd.p0", "sd.del0", "sd.sig0")
        mle = mle[c(2, 3, 1, 5, 6, 4)]
        out = list(mle=mle)
        if (sw==1) {
          Cor = list(Cor = Cor)
          out = c(out, Cor)
        }
        if (!missing(Cov.in)) {
          if (Cov.in$sw == 2) {
            pds. = list(pds.=Cov.out)
            out = c(out, pds.)
          }
          else if (Cov.in$sw == 3) {
            Ilfdr = list(Ilfdr=Cov.out)
            out = c(out, Ilfdr)
          }
          else {
            Cov.lfdr = list(Cov.lfdr=Cov.out)
            out = c(out, Cov.lfdr)
          }
        }
        if ((sw==1) | !missing(Cov.in)) return(out)
        else return(mle)
}

loccov = function(N, N0, p0, d, s, x, X, f, JV, Y, i0, H, h, sw) {
  M = rbind(1, x - Y[1], x^2 - Y[2])
  if (sw==2) {
    K = length(x)
    K0 = length(i0)    
    toprow = c(1 - N0/N, -t(h) %*% JV / s)
    botrow = cbind(0, JV / p0)
    mat = rbind(toprow, botrow)
    M0 = M[,i0]
    dpds.dy0 = mat %*% M0 / N / H[1]
    dy0.dy = matrix(0, K0, K)
    dy0.dy[,i0] = diag(1, K0)
    dpds.dy = dpds.dy0 %*% dy0.dy
    rownames(dpds.dy) = c("p", "d", "s")
    return(dpds.dy)
  }
  else {
    xstd = (x - d)/s
    U = cbind(xstd - H[2]/H[1], xstd^2 - H[3]/H[1])
    M[,-i0] = 0
    dl0plus.dy = cbind(1 - N0/N, U %*% JV / s) %*% M /N/H[1]/p0  
    G <- t(X) %*% (f * X)
#    dl.dy = X %*% solve(G) %*% t(X)
    if (!all(is.finite(G))) return(NA)		# Code modification  !!
    qrG <- qr(G)				# Code modification  !!
    if (qrG$rank!=nrow(G)) return(NA)	 	# Code modification  !!		
    dl.dy = X %*% solve(qrG) %*% t(X)		# Code modification  !!
    dlfdr.dy = dl0plus.dy - dl.dy
    if (sw==3) return(dlfdr.dy)
    else {
      Cov.lfdr = dlfdr.dy %*% (f * t(dlfdr.dy))
      return(Cov.lfdr)
    }
  }
}

loccov2 = function(X, X0, i0, f, ests, N) {
        d = ests[1]
        s = ests[2]
        p0 = ests[3]
        theo = I(ncol(X0)==1)
        Xtil <- X[i0,]
        X0til <- X0[i0,]
        G <- t(X) %*% (f * X)
    	if (!all(is.finite(G))) return(NA)			# Code modification  !!
    	qrG <- qr(G)						# Code modification  !!
    	if (qrG$rank!=nrow(G)) return(NA)			# Code modification  !!
        G0 <- t(X0til) %*% X0til
    	if (!all(is.finite(G0))) return(NA)			# Code modification  !!
    	qrG0 <- qr(G0)						 # Code modification  !!
    	if (qrG0$rank!=nrow(G0)) return(NA)			 # Code modification  !!		
#        B0 <- X0 %*% (solve(G0) %*% t(X0til)) %*% Xtil		 
        B0 <- X0 %*% (solve(qrG0) %*% t(X0til)) %*% Xtil	 # Code modification  !!
        C <- B0 - X
#        Ilfdr = C %*% solve(G, t(X))
#        Cov <- C %*% solve(G) %*% t(C)
        Ilfdr = C %*% solve(G, t(X))				# Code modification  !!
        Cov <- C %*% solve(G) %*% t(C)				# Code modification  !!
        if (theo)
          D = matrix(1,1,1)
        else
          D = matrix(c(1, 0, 0, d, s^2, 0, s^2 + d^2, 2 * d * s^2, s^3), 3)
#	gam. = solve(G0, t(X0til)) %*% (Xtil %*% solve(G, t(X)))
	gam. = solve(qrG0, t(X0til)) %*% (Xtil %*% solve(qrG, t(X)))	 # Code modification  !!
	pds. = D %*% gam.
        if (theo) pds. = rbind(pds., matrix(0, 2, nrow(X)))
	pds.[1,] = pds.[1,] - 1/N
	m1 = pds. %*% f
	m2 = pds.^2 %*% f
	stdev = as.vector(sqrt(m2 - m1^2/N))
	stdev[1] = p0 * stdev[1]
        pds.[1,] = p0 * pds.[1,]
        rownames(pds.) = c("p", "d", "s")
        list(Ilfdr = Ilfdr, pds.=pds., stdev=stdev, Cov=Cov)
}

