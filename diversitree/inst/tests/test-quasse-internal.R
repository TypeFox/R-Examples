source("helper-diversitree.R")

context("QuaSSE (internal)")

## Imports that are generally hidden.
quasse.extent <- diversitree:::quasse.extent
expand.pars.quasse <- diversitree:::expand.pars.quasse
make.pde.quasse.fftC <- diversitree:::make.pde.quasse.fftC
make.pde.quasse.fftR <- diversitree:::make.pde.quasse.fftR
make.pde.quasse.mol <- diversitree:::make.pde.quasse.mol

make.branches.quasse.fftC <- diversitree:::make.branches.quasse.fftC
make.branches.quasse.fftR <- diversitree:::make.branches.quasse.fftR
make.branches.quasse.mol <- diversitree:::make.branches.quasse.mol

## Basic control list.
control.fft <- list(tc=1.3,
                    dt.max=1/20,
                    nx=1024,
                    dx=0.01,
                    r=4L,
                    xmid=0,
                    w=5,
                    flags=0L,
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

lambda <- sigmoid.x
mu <- constant.x
diffusion <- 0.01
sd <- 1/20
len <- 2 # Integrate down a branch length of 2

for ( drift in c(0, .01) ) {
  args <- list(lambda=1:4, mu=5, drift=6, diffusion=7)
  pars <- c(.1, .2, 0, 2.5, .03, drift, diffusion)

  ext.fft <- quasse.extent(control.fft, drift, diffusion)

  ndat <- ext.fft$ndat[2]
  control.mol <- modifyList(control.fft, list(nx=ndat, method="mol"))
  ext.mol <- quasse.extent(control.mol, drift, diffusion)

  ## Initial conditions:
  vars.fft <- matrix(0, control.fft$nx, 2)
  vars.fft[seq_len(ndat),2] <- dnorm(ext.fft$x[[2]], 0, sd)

  vars.mol <- matrix(0, control.mol$nx, 2)
  vars.mol[seq_len(ndat),2] <- dnorm(ext.mol$x[[2]], 0, sd)

  ## TEST
  expect_that(vars.mol,          equals(vars.fft[seq_len(ndat),]))
  expect_that(vars.mol, is_identical_to(vars.fft[seq_len(ndat),]))

  pars.fft <-
    expand.pars.quasse(lambda, mu, args, ext.fft, pars)
  pars.mol <-
    expand.pars.quasse(lambda, mu, args, ext.mol, pars)

  ## TEST
  ## after 5 is padding, which can be dropped (not relevant for MOL)
  expect_that(pars.fft$lo[1:5],          equals(pars.mol$lo[1:5]))
  expect_that(pars.fft$lo[1:5], is_identical_to(pars.mol$lo[1:5]))

  ## Bail here if no FFTW support, even though we could do most of this.
  if (!check.fftC(FALSE)) {
    next
  }
  pde.fftC <- with(control.fft, make.pde.quasse.fftC(nx, dx, dt.max, 2L, flags))
  pde.fftR <- with(control.fft, make.pde.quasse.fftR(nx, dx, dt.max, 2L))
  pde.mol <- with(control.mol, make.pde.quasse.mol(ndat, dx, 2L, atol, rtol))

  ans.fftC <- pde.fftC(vars.fft, len, pars.fft$lo, 0)
  ans.fftR <- pde.fftR(vars.fft, len, pars.fft$lo, 0)

  ## TEST
  expect_that(ans.fftC, equals(ans.fftR))

  if ( drift == 0 ) {
    ans.mol <- pde.mol(vars.mol, len, pars.mol$lo, 0)
    expect_that(ans.mol[[1]], equals(ans.fftC[[1]], tolerance=1e-5))
    expect_that(ans.mol[[2]], equals(ans.fftC[[2]][seq_len(ndat),],
                                     tolerance=0.0004))
  }

  branches.fftC <- make.branches.quasse.fftC(control.fft)
  branches.fftR <- make.branches.quasse.fftR(control.fft)
  branches.mol  <- make.branches.quasse.mol(control.mol)

  vars.hi.fft <- matrix(0, control.fft$nx*control.fft$r, 2)
  vars.hi.fft[seq_len(ext.fft$ndat[1]),2] <-
    dnorm(ext.fft$x[[1]], 0, sd)

  vars.hi.mol <- matrix(0, control.mol$nx*control.mol$r, 2)
  vars.hi.mol[seq_len(ext.mol$ndat[1]),2] <-
    dnorm(ext.mol$x[[1]], 0, 1/20)

  ans.b.fftC <- branches.fftC(vars.hi.fft, len, pars.fft, 0)
  ans.b.fftR <- branches.fftR(vars.hi.fft, len, pars.fft, 0)

  ## TEST:
  expect_that(ans.b.fftC, equals(ans.b.fftR))

  if ( drift == 0 ) {
    ans.b.mol <- branches.mol(vars.hi.mol, len, pars.mol, 0)
    expect_that(c(ans.b.fftC[1],
                  matrix(ans.b.fftC[-1], ncol=2)[seq_len(ndat),]),
                equals(ans.b.mol, tolerance=0.00015))
  }
}
