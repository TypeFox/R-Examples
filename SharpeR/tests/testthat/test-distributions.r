# Copyright 2013 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of SharpeR.
#
# SharpeR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SharpeR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SharpeR.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2013.01.10
# Copyright: Steven E. Pav, 2013-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
is.sorted <- function(xs,pragma=c("ascending","descending")) {
	pragma <- match.arg(pragma)
	retv <- switch(pragma,
								 ascending = !is.unsorted(xs),
								 descending = !is.unsorted(rev(xs)))
	return(retv)
}
# are the data approximately U(0,1)?
is.appx.uniform <- function(xs,slop=1.0) {
	q.check <- seq(0,1,length.out=21)
	retv <- max(abs(quantile(ecdf(xs),q.check) - q.check)) < slop * max(diff(q.check))
	return(retv)
}

set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("do they run?")#FOLDUP

test_that("sr functions",{#FOLDUP
	zeta <- 1.1
	ope <- 252
	df <- 5*ope - 1

	# sample:
	set.char.seed("dee9af9b-cb59-474f-ac3b-acd60faa8ba2")
	rvs <- rsr(500, df=df, zeta=zeta, ope=ope)

	dens <- dsr(rvs, df=df, zeta=zeta, ope=ope)
	expect_true(all(dens >= 0))
	dens <- dsr(rvs, df=df, zeta=zeta, ope=ope, log=TRUE)

	pvs <- psr(rvs, df=df, zeta=zeta, ope=ope)
	expect_true(all(pvs >= 0) && all(pvs <= 1))
	pvs <- psr(rvs, df=df, zeta=zeta, ope=ope, lower.tail=FALSE)
	pvs <- psr(rvs, df=df, zeta=zeta, ope=ope, lower.tail=FALSE, log.p=TRUE)

	qvs <- qsr(ppoints(100), df=df, zeta=zeta, ope=ope)
	qvs <- qsr(ppoints(100), df=df, zeta=zeta, ope=ope, lower.tail=FALSE)
	qvs <- qsr(log(ppoints(100)), df=df, zeta=zeta, ope=ope, lower.tail=FALSE, log.p=TRUE)

	# once again, without optional zeta
	set.char.seed("7166310f-7462-4890-bdd0-0df9ca5d97bd")
	rvs <- rsr(500, df=df, ope=ope)

	dens <- dsr(rvs, df=df, ope=ope)
	expect_true(all(dens >= 0))
	dens <- dsr(rvs, df=df, ope=ope, log=TRUE)

	pvs <- psr(rvs, df=df, ope=ope)
	expect_true(all(pvs >= 0) && all(pvs <= 1))
	pvs <- psr(rvs, df=df, ope=ope, lower.tail=FALSE)
	pvs <- psr(rvs, df=df, ope=ope, lower.tail=FALSE, log.p=TRUE)

	qvs <- qsr(ppoints(100), df=df, ope=ope)
	qvs <- qsr(ppoints(100), df=df, ope=ope, lower.tail=FALSE)
	qvs <- qsr(log(ppoints(100)), df=df, ope=ope, lower.tail=FALSE, log.p=TRUE)
	
	# sentinel:
	expect_true(TRUE)
})#UNFOLD

test_that("sropt functions",{#FOLDUP
	ngen <- 128
	ope <- 252
	df1 <- 8
	df2 <- ope * 10
	drag <- 0
	zeta.s <- 1.1

	# sample:
	set.char.seed("fc90523d-abb4-4543-8476-d48e3a6f28a3")
	rvs <- rsropt(ngen, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag)

	dens <- dsropt(rvs, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag)
	expect_true(all(dens >= 0))
	dens <- dsropt(rvs, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, log=TRUE)
	expect_false(all(dens >= 0))

	pvs <- psropt(rvs, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag)
	expect_true(all(pvs >= 0) && all(pvs <= 1))
	pvs <- psropt(rvs, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, lower.tail=FALSE)
	pvs <- psropt(rvs, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, log.p=TRUE)
	pvs <- psropt(rvs, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, lower.tail=FALSE, log.p=TRUE)

	qvs <- qsropt(ppoints(100), df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag)
	qvs <- qsropt(ppoints(100), df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, lower.tail=FALSE)
	qvs <- qsropt(log(ppoints(100)), df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, lower.tail=FALSE, log.p=TRUE)
	qvs <- qsropt(log(ppoints(100)), df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, drag=drag, lower.tail=FALSE, log.p=TRUE)

	# once again, without optional zeta.s
	set.char.seed("7166310f-7462-4890-bdd0-0df9ca5d97bd")
	rvs <- rsropt(ngen, df1=df1, df2=df2, ope=ope, drag=drag)

	dens <- dsropt(rvs, df1=df1, df2=df2, ope=ope, drag=drag)
	expect_true(all(dens >= 0))
	dens <- dsropt(rvs, df1=df1, df2=df2, ope=ope, drag=drag, log=TRUE)
	expect_false(all(dens >= 0))

	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, drag=drag)
	expect_true(all(pvs >= 0) && all(pvs <= 1))
	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, drag=drag, lower.tail=FALSE)
	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, drag=drag, log.p=TRUE)
	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, drag=drag, lower.tail=FALSE, log.p=TRUE)

	qvs <- qsropt(ppoints(100), df1=df1, df2=df2, ope=ope, drag=drag)
	qvs <- qsropt(ppoints(100), df1=df1, df2=df2, ope=ope, drag=drag, lower.tail=FALSE)
	qvs <- qsropt(log(ppoints(100)), df1=df1, df2=df2, ope=ope, drag=drag, lower.tail=FALSE, log.p=TRUE)
	qvs <- qsropt(log(ppoints(100)), df1=df1, df2=df2, ope=ope, drag=drag, lower.tail=FALSE, log.p=TRUE)

	# once again, with optional drag
	set.char.seed("7166310f-7462-4890-bdd0-0df9ca5d97bd")
	rvs <- rsropt(ngen, df1=df1, df2=df2, ope=ope)

	dens <- dsropt(rvs, df1=df1, df2=df2, ope=ope)
	expect_true(all(dens >= 0))
	dens <- dsropt(rvs, df1=df1, df2=df2, ope=ope, log=TRUE)
	expect_false(all(dens >= 0))

	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope)
	expect_true(all(pvs >= 0) && all(pvs <= 1))
	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, lower.tail=FALSE)
	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, log.p=TRUE)
	pvs <- psropt(rvs, df1=df1, df2=df2, ope=ope, lower.tail=FALSE, log.p=TRUE)

	qvs <- qsropt(ppoints(100), df1=df1, df2=df2, ope=ope)
	qvs <- qsropt(ppoints(100), df1=df1, df2=df2, ope=ope, lower.tail=FALSE)
	qvs <- qsropt(log(ppoints(100)), df1=df1, df2=df2, ope=ope, lower.tail=FALSE, log.p=TRUE)
	qvs <- qsropt(log(ppoints(100)), df1=df1, df2=df2, ope=ope, lower.tail=FALSE, log.p=TRUE)

	# sentinel:
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

context("distribution functions: basic monotonicity")#FOLDUP

test_that("psr/qsr monotonicity",{#FOLDUP
	set.char.seed("1ccb4a05-fd09-43f7-a692-80deebfd67f4")
	
	# psr
	ps <- seq(0.1,0.9,length.out=9)
	for (df in c(256,1024)) {
		for (snr in c(0,1)) {
			for (ope in c(1,52)) {
				for (lp in c(TRUE,FALSE)) {
					if (lp) { checkps <- log(ps) } else { checkps <- ps }
					for (lt in c(TRUE,FALSE)) {
						qs <- qsr(checkps, df, snr, ope, lower.tail=lt, log.p=lp)
						if (lt) { 
							expect_true(is.sorted(qs,pragma="ascending"))
						} else {
							expect_true(is.sorted(qs,pragma="descending"))
						}
						pret <- psr(qs, df, snr, ope, lower.tail=lt, log.p=lp)
						expect_equal(checkps,pret)
					}
				}
			}
		}
	}
})#UNFOLD

test_that("pT2/qT2 monotonicity",{#FOLDUP
	set.char.seed("f28b5cd4-dcdb-4e8e-9398-c79fb9038351")
	
	# pT2
	ps <- seq(0.1,0.9,length.out=9)
	for (df1 in c(2,4,8)) {
		for (df2 in c(256,1024)) {
			for (delta2 in c(0,0.02)) {
				for (lp in c(TRUE,FALSE)) {
					if (lp) { checkps <- log(ps) } else { checkps <- ps }
					for (lt in c(TRUE,FALSE)) {
						qs <- qT2(checkps, df1, df2, delta2, lower.tail=lt, log.p=lp)
						if (lt) { 
							expect_true(is.sorted(qs,pragma="ascending"))
						} else {
							expect_true(is.sorted(qs,pragma="descending"))
						}
						pret <- pT2(qs, df1, df2, delta2, lower.tail=lt, log.p=lp)
						expect_equal(checkps,pret)
					}
				}
			}
		}
	}
})#UNFOLD

test_that("psropt/qsropt monotonicity",{#FOLDUP
	set.char.seed("22ad9afb-49c4-4f37-b32b-eab413b32750")
	
	# psropt
	ps <- seq(0.1,0.9,length.out=9)
	for (df1 in c(2,4,8)) {
		for (df2 in c(256,1024)) {
			for (snrstar in c(0,0.05)) {
				for (ope in c(1,2)) {
					for (drag in c(0,0.1)) {
						for (lp in c(TRUE,FALSE)) {
							if (lp) { checkps <- log(ps) } else { checkps <- ps }
							for (lt in c(TRUE,FALSE)) {
								qs <- qsropt(checkps, df1, df2, snrstar, ope, drag, lower.tail=lt, log.p=lp)
								if (lt) { 
									expect_true(is.sorted(qs,pragma="ascending"))
								} else {
									expect_true(is.sorted(qs,pragma="descending"))
								}
								pret <- psropt(qs, df1, df2, snrstar, ope, drag, lower.tail=lt, log.p=lp)
								expect_equal(checkps,pret)
							}
						}
					}
				}
			}
		}
	}
})#UNFOLD

test_that("plambdap/qlambdap monotonicity",{#FOLDUP
	set.char.seed("a5ae65f3-257a-47a3-af8e-46b34dfcebb0")
	
	# plambdap
	ps <- seq(0.1,0.9,length.out=9)
	for (df in c(4,8,16)) {
		for (tstat in c(-1,0,1)) {
			for (lp in c(TRUE,FALSE)) {
				if (lp) { checkps <- log(ps) } else { checkps <- ps }
				for (lt in c(TRUE,FALSE)) {
					qs <- qlambdap(checkps, df, tstat, lower.tail=lt, log.p=lp)
					if (lt) { 
						expect_true(is.sorted(qs,pragma="ascending"))
					} else {
						expect_true(is.sorted(qs,pragma="descending"))
					}
					pret <- plambdap(qs, df, tstat, lower.tail=lt, log.p=lp)
					expect_equal(checkps,pret,tolerance=1e-4)
				}
			}
		}
	}
})#UNFOLD

test_that("pco_sropt/qco_sropt monotonicity",{#FOLDUP
	set.char.seed("161f0496-0229-4013-a65e-ff7c8b236f4a")
	
	# co_sropt
	ps <- seq(0.05,0.95,length.out=9)
	for (df1 in c(2,4,8)) {
		for (df2 in c(256,1024)) {
			for (delta2 in c(0.72,1.2)) {  # this is the observed SR stat
				for (lp in c(TRUE,FALSE)) {
					if (lp) { checkps <- log(ps) } else { checkps <- ps }
					for (lt in c(TRUE,FALSE)) {
						qs <- qco_sropt(checkps, df1, df2, z.s=delta2, ope=1,
														lower.tail=lt, log.p=lp)
						if (lt) { 
							expect_true(is.sorted(qs,pragma="ascending"))
						} else {
							expect_true(is.sorted(qs,pragma="descending"))
						}
						pret <- pco_sropt(qs, df1, df2, z.s=delta2, ope=1,
															lower.tail=lt, log.p=lp)
						expect_true(ifelse(lp,all(pret <= 0),
															 all(0 <= pret) && all(pret <= 1)))
						expect_equal(checkps,pret,tolerance=0.001)
					}
				}
			}
		}
	}
})#UNFOLD
#UNFOLD

context("distribution functions: parameter monotonicity")#FOLDUP

test_that("psr/qsr parameter monotonicity",{#FOLDUP
	set.char.seed("adc578c1-381c-428a-baae-8d5607732176")
	
	# psr:
	# snrs
	snrs <- seq(-1,1,length.out=11)
	ps <- 0.5
	for (df in c(256,1024)) {
		for (ope in c(1,12)) {
			for (lp in c(TRUE,FALSE)) {
				if (lp) { checkps <- log(ps) } else { checkps <- ps }
				for (lt in c(TRUE,FALSE)) {
					qs <- qsr(checkps, df, snrs, ope, lower.tail=lt, log.p=lp)
					expect_true(is.sorted(qs,pragma="ascending"))
				}
			}
		}
	}
	# ope
	# in this case the nct ncp is decreasing so the bias is as well...
	ope <- c(1,2,4,12,52,253)
	ps <- 0.5
	for (df in c(256,1024)) {
		for (snr in c(0,1)) {
			for (lp in c(TRUE,FALSE)) {
				if (lp) { checkps <- log(ps) } else { checkps <- ps }
				for (lt in c(TRUE,FALSE)) {
					qs <- qsr(checkps, df, snr, ope, lower.tail=lt, log.p=lp)
					expect_true(is.sorted(qs,pragma="descending"))
				}
			}
		}
	}
	# df
	# in this case the bias should decrease in the df...
	# see also Ghosh, B. K. "Some monotonicity theorems for χ 2, 
	# F and t distributions with applications." Journal of the Royal 
	# Statistical Society. Series B (Methodological) (1973): 480-492.
	df <- c(24,52,256,512,1024)
	ps <- 0.5
	for (ope in c(52,256)) {
		for (snr in c(0,1)) {
			for (lp in c(TRUE,FALSE)) {
				if (lp) { checkps <- log(ps) } else { checkps <- ps }
				for (lt in c(TRUE,FALSE)) {
					qs <- qsr(checkps, df, snr, ope, lower.tail=lt, log.p=lp)
					expect_true(is.sorted(qs,pragma="descending"))
				}
			}
		}
	}
	
})#UNFOLD

# 2FIX: other distributions monotonicity wrt parameters...
	## pT2
	#ps <- seq(0.1,0.9,length.out=9)
	#for (df1 in c(2,4,8)) {
		#for (df2 in c(256,1024)) {
			#for (delta2 in c(0,0.02)) {
				#for (lp in c(TRUE,FALSE)) {
					#if (lp) { checkps <- log(ps) } else { checkps <- ps }
					#for (lt in c(TRUE,FALSE)) {
						#qs <- qT2(checkps, df1, df2, delta2, lower.tail=lt, log.p=lp)
						#if (lt) { 
							#expect_true(!is.unsorted(qs))
						#} else {
							#expect_true(!is.unsorted(rev(qs)))
						#}
						#pret <- pT2(qs, df1, df2, delta2, lower.tail=lt, log.p=lp)
						#expect_equal(checkps,pret)
					#}
				#}
			#}
		#}
	#}
	## psropt
	#ps <- seq(0.1,0.9,length.out=9)
	#for (df1 in c(2,4,8)) {
		#for (df2 in c(256,1024)) {
			#for (snrstar in c(0,0.05)) {
				#for (ope in c(1,2)) {
					#for (drag in c(0,0.1)) {
						#for (lp in c(TRUE,FALSE)) {
							#if (lp) { checkps <- log(ps) } else { checkps <- ps }
							#for (lt in c(TRUE,FALSE)) {
								#qs <- qsropt(checkps, df1, df2, snrstar, ope, drag, lower.tail=lt, log.p=lp)
								#if (lt) { 
									#expect_true(!is.unsorted(qs))
								#} else {
									#expect_true(!is.unsorted(rev(qs)))
								#}
								#pret <- psropt(qs, df1, df2, snrstar, ope, drag, lower.tail=lt, log.p=lp)
								#expect_equal(checkps,pret)
							#}
						#}
					#}
				#}
			#}
		#}
	#}
	## plambdap
	#ps <- seq(0.1,0.9,length.out=9)
	#for (df in c(4,8,16)) {
		#for (tstat in c(-1,0,1)) {
			#for (lp in c(TRUE,FALSE)) {
				#if (lp) { checkps <- log(ps) } else { checkps <- ps }
				#for (lt in c(TRUE,FALSE)) {
					#qs <- qlambdap(checkps, df, tstat, lower.tail=lt, log.p=lp)
					#if (lt) { 
						#expect_true(!is.unsorted(qs))
					#} else {
						#expect_true(!is.unsorted(rev(qs)))
					#}
					#pret <- plambdap(qs, df, tstat, lower.tail=lt, log.p=lp)
					#expect_equal(checkps,pret,tolerance=1e-4)
				#}
			#}
		#}
	#}
#UNFOLD

context("distribution functions: sanity checks")#FOLDUP

test_that("qlambdap sensible",{#FOLDUP
	set.char.seed("f54698f1-ec37-49a4-8463-d4209f25afbc")
	
	df <- 128
	true.ncp <- 3
	ngen <- ceiling(THOROUGHNESS * 512)
	alpha.tol = 0.05 + 0.05 / THOROUGHNESS

	tvals <- rt(ngen,df,true.ncp)

	for (p in c(0.05,0.25,0.5,0.75,0.95)) {
		tstat <- sapply(tvals,function(t) { return(qlambdap(p,df,t)) })
		expect_equal(mean(tstat >= true.ncp),p,tolerance=alpha.tol)
	}

	# edge cases
	expect_true(Inf == qlambdap(1,df,1,lower.tail=TRUE))
	expect_true(-Inf == qlambdap(1,df,1,lower.tail=FALSE))
	expect_true(-Inf == qlambdap(0,df,1,lower.tail=TRUE))
	expect_true(Inf == qlambdap(0,df,1,lower.tail=FALSE))

	expect_true(1 == plambdap(Inf,df,1,lower.tail=TRUE))
	expect_true(1 == plambdap(-Inf,df,1,lower.tail=FALSE))
	expect_true(0 == plambdap(-Inf,df,1,lower.tail=TRUE))
	expect_true(0 == plambdap(Inf,df,1,lower.tail=FALSE))
})#UNFOLD
#UNFOLD

context("distribution functions: random generation")#FOLDUP

ngen <- ceiling(THOROUGHNESS * 1024)
alpha.slop = 1 + 0.5 / THOROUGHNESS

test_that("psr/rsr uniform generation",{#FOLDUP
	set.char.seed("6728bbac-7b37-4b26-bd59-f7f074dc3bc3")
	
	# psr
	for (df in c(256,1024)) {
		for (snr in c(0,1)) {
			for (ope in c(1,52)) {
				rvs <- rsr(ngen, df, snr, ope)
				aps <- psr(rvs,  df, snr, ope)
				expect_true(is.appx.uniform(aps,slop=alpha.slop))
			}
		}
	}
})#UNFOLD

test_that("pT2/rT2 uniform generation",{#FOLDUP
	set.char.seed("66ff5482-e5e5-4443-881d-6a9e6fd6fc8d")
	
	# pT2
	for (df1 in c(2,4,8)) {
		for (df2 in c(256,1024)) {
			for (delta2 in c(0,0.02)) {
				rvs <- rT2(ngen, df1, df2, delta2)
				aps <- pT2(rvs,  df1, df2, delta2)
				expect_true(is.appx.uniform(aps,slop=alpha.slop))
			}
		}
	}
})#UNFOLD

test_that("psropt/rsropt uniform generation",{#FOLDUP
	set.char.seed("b7cd77c2-f197-411c-8d14-a37ab4694944")
	
	# psropt
	for (df1 in c(2,4,8)) {
		for (df2 in c(256,1024)) {
			for (snrstar in c(0,0.05)) {
				for (ope in c(1,2)) {
					for (drag in c(0,0.1)) {
						rvs <- rsropt(ngen, df1, df2, snrstar, ope, drag)
						aps <- psropt(rvs,  df1, df2, snrstar, ope, drag)
						expect_true(is.appx.uniform(aps,slop=alpha.slop))
					}
				}
			}
		}
	}
})#UNFOLD

test_that("plambdap/rlambdap uniform generation",{#FOLDUP
	set.char.seed("05a8caf6-35cf-40da-9c7f-57045c4809d4")
	
	for (df in c(256,1024)) {
		for (tstat in c(0,1,2)) {
			rvs <- rlambdap(ngen, df, tstat)
			aps <- plambdap(rvs,  df, tstat)
			expect_true(is.appx.uniform(aps,slop=alpha.slop))
		}
	}
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
