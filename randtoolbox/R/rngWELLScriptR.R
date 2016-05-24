int.to.bits <- function(x)
{
	x[x < 0] <- x[x < 0] + 2^32
	y <- matrix(nrow=32, ncol=length(x))
	for (j in 1:32) {
		y[j,] <- floor(x/2^(32-j)) %% 2
	}
	y
}

hex.to.bits <- function(x)
{
	y <- strsplit(x, "")[[1]]
	y <- match(y, c(0:9, letters[1:6])) - 1
	z <- matrix(nrow=4, ncol=length(y))
	for (j in 1:4)
		z[j,] <- floor(y/2^(4-j)) %% 2
	c(z)
}

bits.to.real <- function(x)
{
	sum(x * 2^-(1:32))
}

getWELLState <- function()
{
	if (RNGkind()[1] == "user-supplied") {
		out <- get.description()
		if (out$name == "WELL") {
			return(int.to.bits(out$state))
		}
	}
	stop("None of the WELL generators is set")
}

rngWELLScriptR <- function(n, s, generator, includeState=FALSE)
{
	params <- getParameters(generator)
	k <- params$k
	w <- params$w
	r <- params$r
	p <- params$p
	mask0 <- params$mask0
	mask1 <- params$mask1
	m1 <- params$m1
	m2 <- params$m2
	m3 <- params$m3
	T0 <- params$T0
	T1 <- params$T1
	T2 <- params$T2
	T3 <- params$T3
	T4 <- params$T4
	T5 <- params$T5
	T6 <- params$T6
	T7 <- params$T7
	Temper  <- params$Temper
	x <- double(n)
	for (i in seq.int(length.out=n)) {
		z0 <- (mask1 * s[, r] + mask0 * s[, r-1]) %% 2
		z1 <- (T0 %*% s[, 1] + T1 %*% s[, m1+1]) %% 2
		z2 <- (T2 %*% s[, m2+1] + T3 %*% s[, m3+1]) %% 2
		z3 <- (z1 + z2) %% 2
		z4 <- (T4 %*% z0 + T5 %*% z1 + T6 %*% z2 + T7 %*% z3) %% 2
		s <- cbind(z4, z3, s[, 2:(r-1)])
		if (is.null(Temper)) {
			x[i] <- bits.to.real(s[, 1])
		} else {
			x[i] <- bits.to.real((Temper %*% s[, 1]) %% 2)
		}
	}
	if (includeState) {
		list(x=x, state=s)
	} else {
		x
	}
}

getParameters <- function(generator)
{
	mask <- function(s, w, p)
	{
		c(rep(s, times=w-p), rep(1-s, times=p))
	}

	M0 <- function(w)
	{
		matrix(0, nrow=w, ncol=w)
	}

	M1 <- function(w)
	{
		diag(w)
	}

	M2 <- function(tt, w)
	{
		A <- rbind(M0(w), M1(w), M0(w))
		A[w + 1:w - tt,]
	}

	M3 <- function(tt, w)
	{
		(M1(w) + M2(tt, w)) %% 2
	}

	M4 <- function(a, w)
	{
		A <- M2(1, w)
		A[,w] <- a
		A
	}

	M5 <- function(tt, b, w) # the paper uses reversed sign of t, see errata
	{
		(M1(w) + diag(b) %*% M2(tt, w)) %% 2
	}

	M6 <- function(qq, s, tt, a, w)
	{
		A <- M2( - qq, w) + M2(w - qq, w)
		A[s+1, ] <- 0
		A[, tt+1] <- A[, tt+1] + a
		A %% 2
	}

	MTemp <- function(b, cc, w)
	{
		(M5(-15, cc, w) %*% M5(-7, b, w)) %% 2
	}

	Temper <- NULL

	switch(generator,
	"512a" =
	{
		k <- 512
		w <- 32
		r <- 16
		p <- 0
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 13
		m2 <- 9
		m3 <- 5
		a1 <- hex.to.bits("da442d24")
		T0 <- M3(-16, w)
		T4 <- M3(-2, w)
		T1 <- M3(-15, w)
		T5 <- M3(-18, w)
		T2 <- M3(11, w)
		T6 <- M2(-28, w) # in the paper is M3(-28)
		T3 <- M0(w)
		T7 <- M5(-5, a1, w)
	},
	"521a" =
	{
		k <- 521
		w <- 32
		r <- 17
		p <- 23
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 13
		m2 <- 11
		m3 <- 10
		T0 <- M3(-13, w)
		T4 <- M3(-13, w)
		T1 <- M3(-15, w)
		T5 <- M2(1, w)
		T2 <- M1(w)
		T6 <- M0(w)
		T3 <- M2(-21, w)
		T7 <- M3(11, w)
	},
	"521b" =
	{
		k <- 521
		w <- 32
		r <- 17
		p <- 23
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 11
		m2 <- 10
		m3 <- 7
		T0 <- M3(-21, w)
		T4 <- M3(13, w)
		T1 <- M3(6, w)
		T5 <- M2(-10, w)
		T2 <- M0(w)
		T6 <- M2(-5, w)
		T3 <- M3(-13, w)
		T7 <- M3(13, w)
	},
	"607a" =
	{
		k <- 607
		w <- 32
		r <- 19
		p <- 1
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 16
		m2 <- 15
		m3 <- 14
		T0 <- M3(19, w)
		T4 <- M3(18, w)
		T1 <- M3(11, w)
		T5 <- M1(w)
		T2 <- M3(-14, w)
		T6 <- M0(w)
		T3 <- M1(w)
		T7 <- M3(-5, w)
	},
	"607b" =
	{
		k <- 607
		w <- 32
		r <- 19
		p <- 1
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 16
		m2 <- 8
		m3 <- 13
		T0 <- M3(-18, w)
		T4 <- M3(-24, w)
		T1 <- M3(-14, w)
		T5 <- M3(5, w)
		T2 <- M0(w)
		T6 <- M3(-1, w)
		T3 <- M3(18, w)
		T7 <- M0(w)
	},
	"800a" =
	{
		k <- 800
		w <- 32
		r <- 25
		p <- 0
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 14
		m2 <- 18
		m3 <- 17
		T0 <- M1(w)
		T4 <- M3(16, w)
		T1 <- M3(-15, w)
		T5 <- M2(20, w)
		T2 <- M3(10, w)
		T6 <- M1(w)
		T3 <- M3(-11, w)
		T7 <- M3(-28, w)
	},
	"800b" =
	{
		k <- 800
		w <- 32
		r <- 25
		p <- 0
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 9
		m2 <- 4
		m3 <- 22
		a2 <- hex.to.bits("d3e43ffd")
		T0 <- M3(-29, w)
		T4 <- M1(w)
		T1 <- M2(-14, w)
		T5 <- M3(10, w)
		T2 <- M1(w)
		T6 <- M4(a2, w)
		T3 <- M2(19, w)
		T7 <- M3(-25, w)
	},
	"1024a" =
	{
		k <- 1024
		w <- 32
		r <- 32
		p <- 0
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 3
		m2 <- 24
		m3 <- 10
		T0 <- M1(w)
		T4 <- M3(-11, w)
		T1 <- M3(8, w)
		T5 <- M3(-7, w)
		T2 <- M3(-19, w)
		T6 <- M3(-13, w)
		T3 <- M3(-14, w)
		T7 <- M0(w)
	},
	"1024b" =
	{
		k <- 1024
		w <- 32
		r <- 32
		p <- 0
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 22
		m2 <- 25
		m3 <- 26
		a3 <- hex.to.bits("8bdcb91e")
		T0 <- M3(-21, w)
		T4 <- M3(-14, w)
		T1 <- M3(17, w)
		T5 <- M3(-21, w)
		T2 <- M4(a3, w)
		T6 <- M1(w)
		T3 <- M3(15, w)
		T7 <- M0(w)
	},
	"19937a" =
	{
		k <- 19937
		w <- 32
		r <- 624
		p <- 31
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 70
		m2 <- 179
		m3 <- 449
		T0 <- M3(-25, w)
		T4 <- M1(w)
		T1 <- M3(27, w)
		T5 <- M3(-9, w)
		T2 <- M2(9, w)
		T6 <- M3(-21, w)
		T3 <- M3(1, w)
		T7 <- M3(21, w)
	},
	"19937c" =
	{
		k <- 19937
		w <- 32
		r <- 624
		p <- 31
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 70
		m2 <- 179
		m3 <- 449
		T0 <- M3(-25, w)
		T4 <- M1(w)
		T1 <- M3(27, w)
		T5 <- M3(-9, w)
		T2 <- M2(9, w)
		T6 <- M3(-21, w)
		T3 <- M3(1, w)
		T7 <- M3(21, w)
		b <- hex.to.bits("e46e1700")
		cc <- hex.to.bits("9b868000")
		Temper <- MTemp(b, cc, w)
	},
	"19937b" =
	{
		k <- 19937
		w <- 32
		r <- 624
		p <- 31
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 203
		m2 <- 613
		m3 <- 123
		T0 <- M3(7, w)
		T4 <- M3(-19, w)
		T1 <- M1(w)
		T5 <- M2(-11, w)
		T2 <- M3(12, w)
		T6 <- M3(4, w)
		T3 <- M3(-10, w)
		T7 <- M3(-10, w)
	},
	"21701a" =
	{
		k <- 21701
		w <- 32
		r <- 679
		p <- 27
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 151
		m2 <- 327
		m3 <- 84
		a4 <- hex.to.bits("86a9d87e")
		T0 <- M1(w)
		T4 <- M3(27, w)
		T1 <- M3(-26, w)
		T5 <- M3(-11, w)
		T2 <- M3(19, w)
		T6 <- M6(15, 27, 10, a4, w) # in the paper is M6(15, 10, 27, a4), see errata
		T3 <- M0(w)
		T7 <- M3(-16, w)
	},
	"23209a" =
	{
		k <- 23209
		w <- 32
		r <- 726
		p <- 23
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 667
		m2 <- 43
		m3 <- 462
		T0 <- M3(28, w)
		T4 <- M3(21, w)
		T1 <- M1(w)
		T5 <- M3(-17, w)
		T2 <- M3(18, w)
		T6 <- M3(-28, w)
		T3 <- M3(3, w)
		T7 <- M3(-1, w)
	},
	"23209b" =
	{
		k <- 23209
		w <- 32
		r <- 726
		p <- 23
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 610
		m2 <- 175
		m3 <- 662
		a5 <- hex.to.bits("a8c296d1")
		a6 <- hex.to.bits("5d6b45cc")
		T0 <- M4(a5, w)
		T4 <- M3(-26, w)
		T1 <- M1(w)
		T5 <- M1(w)
		T2 <- M6(15, 15, 30, a6, w) # in the paper is M6(15, 30, 15, a6), see errata
		T6 <- M0(w)
		T3 <- M3(-24, w)
		T7 <- M3(16, w)
	},
	"44497a" =
	{
		k <- 44497
		w <- 32
		r <- 1391
		p <- 15
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 23
		m2 <- 481
		m3 <- 229
		a7 <- hex.to.bits("b729fcec")
		T0 <- M3(-24, w)
		T4 <- M1(w)
		T1 <- M3(30, w)
		T5 <- M3(20, w)
		T2 <- M3(-10, w)
		T6 <- M6(9, 5, 14, a7, w) # in the paper is M6(9, 14, 5, a7), see errata
		T3 <- M2(-26, w)
		T7 <- M1(w)
	},
	"44497b" =
	{
		k <- 44497
		w <- 32
		r <- 1391
		p <- 15
		mask0 <- mask(0, w, p)
		mask1 <- mask(1, w, p)
		m1 <- 23
		m2 <- 481
		m3 <- 229
		a7 <- hex.to.bits("b729fcec")
		T0 <- M3(-24, w)
		T4 <- M1(w)
		T1 <- M3(30, w)
		T5 <- M3(20, w)
		T2 <- M3(-10, w)
		T6 <- M6(9, 5, 14, a7, w) # in the paper is M6(9, 14, 5, a7), see errata
		T3 <- M2(-26, w)
		T7 <- M1(w)
		b <- hex.to.bits("93dd1400")
		cc <- hex.to.bits("fa118000")
		Temper <- MTemp(b, cc, w)
	},
	"default" =
	{
		stop("unsupported parameters")
	})
	stopifnot(r*w - p == k)
	stopifnot(r == ceiling(k/32))
	list(k=k, w=w, r=r, p=p, mask0=mask0, mask1=mask1, m1=m1, m2=m2, m3=m3,
		T0=T0, T1=T1, T2=T2, T3=T3, T4=T4, T5=T5, T6=T6, T7=T7, Temper=Temper)
}

