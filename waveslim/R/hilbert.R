########################################################################

dwt.hilbert <- function(x, wf, n.levels=4, boundary="periodic", ...) {
  switch(boundary,
         "reflection" =  x <- c(x, rev(x)),
         "periodic" = invisible(),
         stop("Invalid boundary rule in dwt.hilbert"))
  N <- length(x)
  J <- n.levels
  if(N/2^J != trunc(N/2^J)) 
    stop("Sample size is not divisible by 2^J")
  if(2^J > N) 
    stop("Wavelet transform exceeds sample size in dwt")

  dict <- hilbert.filter(wf)
  L <- dict$length; storage.mode(L) <- "integer"
  h0 <- dict$lpf[[1]]; storage.mode(h0) <- "double"
  g0 <- dict$lpf[[2]]; storage.mode(g0) <- "double"
  h1 <- dict$hpf[[1]]; storage.mode(h1) <- "double"
  g1 <- dict$hpf[[2]]; storage.mode(g1) <- "double"

  y <- vector("list", J+1)
  names(y) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  x.h <- x.g <- x
  for(j in 1:J) {
    W <- V <- numeric(N/2^j)
    out.h <- .C("dwt", as.double(x.h), as.integer(N/2^(j-1)), L, h1, h0, 
                W = W, V = V, PACKAGE = "waveslim")[6:7]
    out.g <- .C("dwt", as.double(x.g), as.integer(N/2^(j-1)), L, g1, g0, 
                W = W, V = V, PACKAGE = "waveslim")[6:7]
    y[[j]] <- complex(real = out.h$W, imaginary = out.g$W)
    x.h <- out.h$V
    x.g <- out.g$V
  }
  y[[J+1]] <- complex(real = x.h, imaginary = x.g)
  attr(y, "wavelet") <- wf
  attr(y, "levels") <- n.levels
  attr(y, "boundary") <- boundary
  return(y)
}

########################################################################

dwt.hilbert.nondyadic <- function(x, ...) {
  M <- length(x)
  N <- 2^(ceiling(log(M, 2)))
  xx <- c(x, rep(0, N - M))
  y <- dwt.hilbert(xx, ...)
  
  J <- length(y) - 1
  for(j in 1:J) {
    y[[j]] <- y[[j]][1:trunc(M/2^j)]
  }

  return(y)
}

########################################################################

idwt.hilbert <- function(y) {
  switch(attributes(y)$boundary,
    "reflection" =  x <- c(x, rev(x)),
    "periodic" = invisible(),
    stop("Invalid boundary rule in dwt.dbp"))
  J <- attributes(y)$levels
 
  dict <- hilbert.filter(attributes(y)$wavelet)
  L <- dict$length; storage.mode(L) <- "integer"
  h <- dict$hpf; storage.mode(h) <- "double"
  g <- dict$lpf; storage.mode(g) <- "double"

  jj <- paste("s", J, sep="")
  X <- y[[jj]]
  for(j in J:1) {
    jj <- paste("d", j, sep="")
    XX <- numeric(2 * length(y[[jj]]))
    X <- .C("idwt", y[[jj]], as.double(X), as.integer(length(X)),
            L, h, g, XX=XX, PACKAGE="waveslim")$XX
  }
  return(X)
}

########################################################################

modwt.hilbert <- function(x, wf, n.levels=4, boundary="periodic", ...) {
  switch(boundary,
         "reflection" =  x <- c(x, rev(x)),
         "periodic" = invisible(),
         stop("Invalid boundary rule in modwt"))
  N <- length(x)
  storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N) stop("wavelet transform exceeds sample size in modwt")

  dict <- hilbert.filter(wf)
  L <- dict$length; storage.mode(L) <- "integer"
  h0 <- dict$lpf[[1]] / sqrt(2); storage.mode(h0) <- "double"
  g0 <- dict$lpf[[2]] / sqrt(2); storage.mode(g0) <- "double"
  h1 <- dict$hpf[[1]] / sqrt(2); storage.mode(h1) <- "double"
  g1 <- dict$hpf[[2]] / sqrt(2); storage.mode(g1) <- "double"

  y <- vector("list", J+1)
  names(y) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  W <- V <- numeric(N)

  x.h <- x.g <- x  
  for(j in 1:J) {
    out.h <- .C("modwt", as.double(x.h), N, as.integer(j), L, h1, h0,
                W = W, V = V, PACKAGE="waveslim")[7:8]
    out.g <- .C("modwt", as.double(x.g), N, as.integer(j), L, g1, g0,
                W = W, V = V, PACKAGE="waveslim")[7:8]
    y[[j]] <- complex(real = out.h$W, imaginary = out.g$W)
    x.h <- out.h$V
    x.g <- out.g$V
  }
  y[[J+1]] <- complex(real = x.h, imaginary = x.g)
  attr(y, "wavelet") <- wf
  attr(y, "boundary") <- boundary
  attr(y, "levels") <- n.levels
  return(y)
}

########################################################################

imodwt.hilbert <- function(y) {
  if(attributes(y)$boundary != "periodic")
    stop("Invalid boundary rule in imodwt")
  J <- length(y) - 1

  dict <- hilbert.filter(attributes(y)$wavelet)
  L <- dict$length
  ht <- dict$hpf / sqrt(2)
  gt <- dict$lpf / sqrt(2)

  jj <- paste("s", J, sep="")
  X <- y[[jj]]; N <- length(X)
  XX <- numeric(N)
  for(j in J:1) {
    jj <- paste("d", j, sep="")
    X <- .C("imodwt", y[[jj]], X, as.integer(N), as.integer(j), 
      as.integer(L), ht, gt, XX, PACKAGE="waveslim")[[8]]
  }
  return(X)
}

########################################################################

hilbert.filter <- function(name) {
  select.K3L3 <- function() {
    L <- 12
    h0 <- c(1.1594353e-04, -2.2229002e-03, -2.2046914e-03, 4.3427642e-02,
            -3.3189896e-02, -1.5642755e-01, 2.8678636e-01, 7.9972652e-01,
            4.9827824e-01, 2.4829160e-02, -4.2679177e-02, -2.2260892e-03)
    h1 <- qmf(h0)
    g0 <- c(1.6563361e-05, -5.2543406e-05, -6.1909121e-03, 1.9701141e-02,
            3.2369691e-02, -1.2705043e-01, -1.5506397e-02, 6.1333712e-01,
            7.4585008e-01, 2.1675412e-01, -4.9432248e-02, -1.5582624e-02)
    g1 <- qmf(g0)
    return(list(length = L, hpf = list(h1, g1), lpf = list(h0, g0)))
  }
    select.K3L5 <- function() {
    L <- 12
    h0 <- c(5.4258791e-06, -2.1310518e-04, -2.6140914e-03,  1.0212881e-02,
            3.5747880e-02, -4.5576766e-02, 3.9810341e-03, 5.3402475e-01,
            7.8757164e-01, 2.6537457e-01, -1.3008915e-01, -5.9573795e-02,
            1.2733976e-02, 2.8641011e-03, -2.2992683e-04, -5.8541759e-06)
    h1 <- qmf(h0)
    g0 <- c(4.9326174e-07, 3.5727140e-07, -1.1664703e-03, -8.4003116e-04,
            2.8601474e-02, 9.2509748e-03, -7.4562251e-02, 2.2929480e-01,
            7.6509138e-01, 5.8328559e-01, -4.6218010e-03, -1.2336841e-01,
            -6.2826896e-03, 9.5478911e-03, 4.6642226e-05, -6.4395935e-05)
    g1 <- qmf(g0)
    return(list(length = L, hpf = list(h1, g1), lpf = list(h0, g0)))
  }
  select.K4L2 <- function() {
    L <- 12
    h0 <- c(-1.7853301e-03, 1.3358873e-02, 3.6090743e-02, -3.4722190e-02,
            4.1525062e-02, 5.6035837e-01, 7.7458617e-01, 2.2752075e-01,
            -1.6040927e-01, -6.1694251e-02, 1.7099408e-02, 2.2852293e-03)
    h1 <- qmf(h0)
    g0 <- c(-3.5706603e-04, -1.8475351e-04, 3.2591486e-02, 1.3449902e-02,
            -5.8466725e-02, 2.7464308e-01, 7.7956622e-01, 5.4097379e-01,
            -4.0315008e-02, -1.3320138e-01, -5.9121296e-03, 1.1426146e-02)
    g1 <- qmf(g0)
    return(list(length = L, hpf = list(h1, g1), lpf = list(h0, g0)))
  }
  select.K4L4 <- function() {
    L <- 16
    h0 <- c(2.5734665593981519e-05, -6.6909066441298817e-04, 
            -5.5482443985275260e-03, 1.3203474646343588e-02,
            3.8605327384848696e-02, -5.0687259299773510e-02,
            8.1364447220208733e-03, 5.3021727476690994e-01,
            7.8330912249663232e-01, 2.7909546754271131e-01,
            -1.3372674246928601e-01, -6.9759509629953295e-02,
            1.6979390952358446e-02, 5.7323570134311854e-03,
            -6.7425216644469892e-04, -2.5933188060087743e-05)
    h1 <- qmf(h0)
    g0 <- c(2.8594072882201687e-06, 1.9074538622058143e-06,
            -2.9903835439216066e-03, -1.9808995184875909e-03,
            3.3554663884350758e-02, 7.7023844121478988e-03,
            -7.7084571412435535e-02, 2.3298110528093252e-01,
            7.5749376288995063e-01, 5.8834703992067783e-01,
            5.1708789323078770e-03, -1.3520099946241465e-01,
            -9.1961246067629732e-03, 1.5489641793018745e-02,
            1.5569563641876791e-04, -2.3339869254078969e-04)
    g1 <- qmf(g0)
    return(list(length = L, hpf = list(h1, g1), lpf = list(h0, g0)))
  }
  select.K5L7 <- function() {
    L <- 24
    h0 <- c(-2.5841959496364648e-10, 6.0231243121018760e-10,
            2.1451486802217960e-06,  -4.9989222844980982e-06,
            -2.2613489535132104e-04, 5.1967501391358343e-04,
            3.4011963595840899e-03,  -7.1996997688061597e-03,
            -1.7721433874932836e-02, 3.5491112173858148e-02,
            3.0580617312936355e-02,  -1.3452365188777773e-01,
            2.1741748603083836e-03, 5.8046856094922639e-01,
            7.4964083145768690e-01, 2.6775497264154541e-01,
            -7.9593287728224230e-02,  -4.3942149960221458e-02,
            1.9574969406037097e-02, 8.8554643330725387e-03,
            -7.2770446614145033e-04,  -3.1310992841759443e-04,
            1.4045333283124608e-06, 6.0260907100656169e-07)
    h1 <- qmf(h0)
    g0 <- c(-3.8762939244546978e-09, 2.9846463282743695e-07,
            5.6276030758515370e-06, -7.7697066311187957e-05,
            -2.1442686434841905e-04, 2.1948612668324223e-03,
            9.5408758453423542e-04, -1.7149735951945008e-02,
            1.5212479104581677e-03, 5.6600564413983846e-02,
            -4.8900162376504831e-02, -1.3993440493611778e-01,
            2.7793346796113222e-01, 7.6735603850281364e-01,
            5.4681951651005178e-01, 3.6275855872448776e-02,
            -8.8224410289407154e-02, 3.2821708368951431e-05,
            1.7994969189524142e-02, 1.8662128501760204e-03,
            -7.8622878632753014e-04, -5.8077443328549205e-05,
            3.0932895975646042e-06, 4.0173938067104100e-08)
    g1 <- qmf(g0)
    return(list(length = L, hpf = list(h1, g1), lpf = list(h0, g0)))
  }
  select.K6L6 <- function() {
    L <- 24
    h0 <- c(1.4491207137947255e-09  -3.4673992369566253e-09,
            -6.7544152844875963e-06, 1.6157040144070828e-05,
            4.0416340595645441e-04, -9.4536696039781878e-04,
            -4.2924086033924620e-03, 9.0688042722858742e-03,
            1.8690864167884680e-02, -3.7883945370993717e-02,
            -2.7337592282061701e-02, 1.3185812419468312e-01,
            -2.1034481553730465e-02, -5.9035515013747486e-01,
            -7.4361804647499452e-01, -2.5752016951708306e-01,
            9.2725410672739983e-02, 4.9100676534870831e-02,
            -2.4411085480175867e-02, -1.1190458223944993e-02,
            1.7793885751382626e-03, 7.4715940333597059e-04,
            -6.2392430013359510e-06, -2.6075498267775052e-06)
    h1 <- qmf(h0)
    g0 <- c(1.8838569279331431e-08, -1.1000360697229965e-06,
            -1.4600820117782769e-05, 1.6936567299204319e-04,
            2.6967189953984829e-04, -3.1633669438102655e-03,
            -7.2081460313487946e-04, 1.9638595542490079e-02,
            -3.0968325940269846e-03, -5.6722348677476261e-02,
            5.2260784738219289e-02, 1.2763836788794369e-01,
            -2.9566169882112192e-01, -7.6771793937333599e-01,
            -5.3818432160802543e-01, -2.4023872575927138e-02,
            9.9019132161496132e-02, -1.2059411664071501e-03,
            -2.2693488886969308e-02, -1.8724943382560243e-03,
            1.7270823778712107e-03, 1.5415480681200776e-04,
            -1.1712464100067407e-05, -2.0058075590596196e-07)
    g1 <- qmf(g0)
    return(list(length = L, hpf = list(-h1, -g1), lpf = list(-h0, -g0)))
  }

  switch(name,
         "k3l3" = select.K3L3(),
         "k3l5" = select.K3L5(),
         "k4l2" = select.K4L2(),
         "k4l4" = select.K4L4(),
         "k5l7" = select.K5L7(),
         "k6l6" = select.K6L6(),
         stop("Invalid selection for hilbert.filter"))
}

########################################################################

phase.shift.hilbert <- function(x, wf) {
  coe <- function(g)
    sum(0:(length(g)-1) * g^2) / sum(g^2)

  J <- length(x) - 1

  h0 <- hilbert.filter(wf)$lpf[[1]]
  h1 <- hilbert.filter(wf)$hpf[[1]]
  
  for(j in 1:J) {
    ph <- round(2^(j-1) * (coe(h0) + coe(h1)) - coe(h0), 0)
    Nj <- length(x[[j]])
    x[[j]] <- c(x[[j]][(ph+1):Nj], x[[j]][1:ph])
  }

  ph <- round((2^J-1) * coe(h0), 0)
  J <- J + 1
  x[[J]] <- c(x[[J]][(ph+1):Nj], x[[J]][1:ph])

  return(x)
}

########################################################################

modwpt.hilbert <- function(x, wf, n.levels=4, boundary="periodic") {
  N <- length(x)
  storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N) stop("wavelet transform exceeds sample size in modwpt")

  dict <- hilbert.filter(wf)
  L <- dict$length; storage.mode(L) <- "integer"
  h0 <- dict$lpf[[1]] / sqrt(2); storage.mode(h0) <- "double"
  g0 <- dict$lpf[[2]] / sqrt(2); storage.mode(g0) <- "double"
  h1 <- dict$hpf[[1]] / sqrt(2); storage.mode(h1) <- "double"
  g1 <- dict$hpf[[2]] / sqrt(2); storage.mode(g1) <- "double"

  y <- vector("list", sum(2^(1:J)))
  yn <- length(y)
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  names(y) <- paste("w", crystals1, ".", crystals2, sep="")

  W <- V <- numeric(N)
  storage.mode(W) <- storage.mode(V) <- "double"
  for(j in 1:J) {
    ## cat(paste("j =", j, fill=T))
    index <- 0
    jj <- min((1:yn)[crystals1 == j])
    for(n in 0:(2^j / 2 - 1)) {
      index <- index + 1
      if(j > 1)
        x <- y[[(1:yn)[crystals1 == j-1][index]]]
      else
        x <- complex(real=x, imaginary=x)

      if(n %% 2 == 0) {
        zr <- .C("modwt", as.double(Re(x)), N, as.integer(j), L, h1, h0, 
                 W = W, V = V, PACKAGE="waveslim")[7:8]
        zc <- .C("modwt", as.double(Im(x)), N, as.integer(j), L, g1, g0, 
                 W = W, V = V, PACKAGE="waveslim")[7:8]
        y[[jj + 2*n + 1]] <- complex(real=zr$W, imaginary=zc$W)
        y[[jj + 2*n]] <- complex(real=zr$V, imaginary=zc$V)
      }
      else {
        zr <- .C("modwt", as.double(Re(x)), N, as.integer(j), L, h1, h0, 
                 W = W, V = V, PACKAGE="waveslim")[7:8]
        zc <- .C("modwt", as.double(Im(x)), N, as.integer(j), L, g1, g0, 
                 W = W, V = V, PACKAGE="waveslim")[7:8]
        y[[jj + 2*n]] <- complex(real=zr$W, imaginary=zc$W)
        y[[jj + 2*n + 1 ]] <- complex(real=zr$V, imaginary=zc$V)
      }
    }
  }
  attr(y, "wavelet") <- wf
  return(y)
}

########################################################################

phase.shift.hilbert.packet <- function(x, wf) {
  coe <- function(g)
    sum(0:(length(g)-1) * g^2) / sum(g^2)

  dict <- hilbert.filter(wf)
  h0 <- dict$lpf[[1]]; h1 <- dict$hpf[[1]]
  g0 <- dict$lpf[[2]]; g1 <- dict$hpf[[2]]

  xn <- length(x)
  N <- length(x[[1]])
  J <- trunc(log(xn,2))
  jbit <- vector("list", xn)
  jbit[[1]] <- FALSE; jbit[[2]] <- TRUE
  crystals1 <- rep(1:J, 2^(1:J))

  for(j in 1:J) {
    jj <- min((1:xn)[crystals1 == j])
    for(n in 0:(2^j - 1)) {
      if(j > 1) {
        jp <- min((1:xn)[crystals1 == j-1])
        if(n %% 4 == 0 | n %% 4 == 3)
          jbit[[jj + n]] <- c(jbit[[jp + floor(n/2)]], FALSE)
        else
          jbit[[jj + n]] <- c(jbit[[jp + floor(n/2)]], TRUE)
      }
      Sjn0 <- sum((1 - jbit[[jj + n]]) * 2^(0:(j-1)))
      Sjn1 <- sum(jbit[[jj + n]] * 2^(0:(j-1)))
      ph <- round(Sjn0 * coe(h0) + Sjn1 * coe(h1), 0)
      x[[jj + n]] <- c(x[[jj + n]][(ph+1):N], x[[jj + n]][1:ph])
    }
  }
  return(x)
}

########################################################################

modhwt.coh <- function(x, y, f.length = 0) {
  filt <- rep(1, f.length + 1)
  filt <- filt / length(filt)

  J <- length(x) - 1
  coh <- vector("list", J)
  for(j in 1:J) {
    co.spec <- filter(Re(x[[j]] * Conj(y[[j]])), filt) 
    quad.spec <- filter(-Im(x[[j]] * Conj(y[[j]])), filt)
    x.spec <- filter(Mod(x[[j]])^2, filt)
    y.spec <- filter(Mod(y[[j]])^2, filt)
    coh[[j]] <- (co.spec^2 + quad.spec^2) / x.spec / y.spec
  }
  coh
}

########################################################################

modhwt.phase <- function(x, y, f.length = 0) {
  filt <- rep(1, f.length + 1)
  filt <- filt / length(filt)

  J <- length(x) - 1
  phase <- vector("list", J)
  for(j in 1:J) {
    co.spec <- filter(Re(x[[j]] * Conj(y[[j]])), filt) 
    quad.spec <- filter(-Im(x[[j]] * Conj(y[[j]])), filt)
    phase[[j]] <- Arg(co.spec - 1i * quad.spec)
  }
  phase
}

########################################################################

modhwt.coh.seasonal <- function(x, y, S=10, season=365) {
  J <- length(x) - 1
  coh <- shat <- vector("list", J)
  for(j in 1:J) {
    xj <- x[[j]]
    yj <- y[[j]]

    ## Cospectrum
    co <- matrix(Re(xj * Conj(yj)), ncol=season, byrow=TRUE)
    co.spec <- c(apply(co, 2, mean, na.rm=TRUE))
    gamma.c <- my.acf(as.vector(co))
    omega.c <- sum(gamma.c[c(1, rep(seq(season+1, S*season, by=season),
                                    each=2))])
    
    ## Quadrature spectrum
    quad <- matrix(-Im(xj * Conj(yj)), ncol=season, byrow=TRUE)
    quad.spec <- c(apply(quad, 2, mean, na.rm=TRUE))
    gamma.q <- my.acf(as.vector(quad))
    omega.q <- sum(gamma.q[c(1, rep(seq(season+1, S*season, by=season),
                                    each=2))])

    gamma.cq <- my.ccf(as.vector(co), as.vector(quad))
    omega.cq <- sum(gamma.cq[S*season + seq(-S*season+1, S*season, by=season)])

    ## Autospectrum(X)
    autoX <- matrix(Mod(xj)^2, ncol=season, byrow=TRUE)
    x.spec <- c(apply(autoX, 2, mean, na.rm=TRUE))

    ## Autospectrum(Y)
    autoY <- matrix(Mod(yj)^2, ncol=season, byrow=TRUE)
    y.spec <- c(apply(autoY, 2, mean, na.rm=TRUE))

    shat[[j]] <- 4 * (co.spec*omega.c + quad.spec * omega.q +
                      2*co.spec*quad.spec*omega.cq) / x.spec^2 / y.spec^2
    coh[[j]] <- (co.spec^2 + quad.spec^2) / x.spec / y.spec
  }
  list(coh = coh, var = shat)
}

########################################################################

modhwt.phase.seasonal <- function(x, y, season=365) {
  J <- length(x) - 1
  phase <- vector("list", J)
  for(j in 1:J) {
    co.spec <- Re(x[[j]] * Conj(y[[j]]))
    co.spec <- c(apply(matrix(co.spec, ncol=season, byrow=TRUE), 2,
                       mean, na.rm=TRUE))
    quad.spec <- -Im(x[[j]] * Conj(y[[j]]))
    quad.spec <- c(apply(matrix(quad.spec, ncol=season, byrow=TRUE), 2,
                         mean, na.rm=TRUE))
    phase[[j]] <- Arg(co.spec - 1i * quad.spec)
  }
  phase
}

