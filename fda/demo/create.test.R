#  run some tests on basis functions

rangeval = c(0,1)
nbasis   = 5

#  -------------------------  Bspline basis  ---------------

library("splines")

if(exists('basisobj'))remove(basisobj)
basisobj = create.bspline.basis()
basisobj

remove(basisobj)
basisobj = create.bspline.basis(rangeval)
basisobj

remove(basisobj)
basisobj = create.bspline.basis(rangeval, nbasis)
basisobj

norder = 3
remove(basisobj)
basisobj = create.bspline.basis(rangeval, nbasis, norder)
basisobj

remove(basisobj)
nbasis = 8
basisobj = create.bspline.basis(rangeval, nbasis, 3,
                       c(0, 0.2, 0.45, 0.5, 0.55, 0.8, 1))
basisobj

remove(basisobj)
basisobj = create.bspline.basis(rangeval, nbasis)

plot(basisobj)

eval.penalty(basisobj, 2)

remove(basisobj)
basisobj = create.bspline.basis(rangeval, nbasis,
                                dropind=c(1,nbasis))

eval.penalty(basisobj, 2)

remove(basisobj)
basisobj = create.bspline.basis(rangeval, nbasis)
basisobj

evalarg = seq(0,1,0.1)
basismat = eval.basis(evalarg, basisobj)
list1 = list(evalarg=evalarg, basismat=basismat)

basisobj$basisvalues[[1]] = list1
class(basisobj$basisvalues[[1]])
basisobj$basisvalues[[1]]

#  -----------------   constant basis  ------------------

remove(basisobj)
basisobj = create.constant.basis()
basisobj

remove(basisobj)
basisobj = create.constant.basis(rangeval)
basisobj

plot(basisobj)

#  ---------------------  exponential basis  ---------------

remove(basisobj)
basisobj = create.exponential.basis()
basisobj

remove(basisobj)
basisobj = create.exponential.basis(rangeval)
basisobj

remove(basisobj)
nbasis = 5
basisobj = create.exponential.basis(rangeval, nbasis)
basisobj

remove(basisobj)
basisobj = create.exponential.basis(rangeval, nbasis,
                                    (-2:2))
basisobj

plot(basisobj)

eval.penalty(basisobj, 2)

remove(basisobj)
basisobj = create.exponential.basis(rangeval, nbasis,
                                    (-2:2), 3)

eval.penalty(basisobj, 2)

#  --------------------  fourier basis  ----------------

rangeval = c(-pi,pi)

remove(basisobj)
basisobj = create.fourier.basis()
basisobj

remove(basisobj)
basisobj = create.fourier.basis(rangeval)
basisobj

remove(basisobj)
basisobj = create.fourier.basis(rangeval, nbasis)
basisobj

remove(basisobj)
basisobj = create.fourier.basis(rangeval, nbasis, pi)
basisobj

plot(basisobj)

eval.penalty(basisobj, 2)

remove(basisobj)
basisobj = create.fourier.basis(rangeval, nbasis, pi, 1)

eval.penalty(basisobj, 2)

#  ----------------------  monomial basis  ----------------

rangeval = c(-1,1)

remove(basisobj)
basisobj = create.monomial.basis()
basisobj

remove(basisobj)
basisobj = create.monomial.basis(rangeval)
basisobj

remove(basisobj)
basisobj = create.monomial.basis(rangeval, nbasis)
basisobj

eval.penalty(basisobj, 2)

remove(basisobj)
basisobj = create.monomial.basis(rangeval, nbasis, dropind=1)

eval.penalty(basisobj, 2)

#  ---------------------  polygonal basis  ------------------

rangeval = c(0,1)

remove(basisobj)
basisobj = create.polygonal.basis()
basisobj

remove(basisobj)
basisobj = create.polygonal.basis(rangeval)
basisobj

remove(basisobj)
basisobj = create.polygonal.basis(rangeval, c(0, 0.5, 1))
basisobj

plot(basisobj)

eval.penalty(basisobj, 0)
eval.penalty(basisobj, 1)

remove(basisobj)
basisobj = create.polygonal.basis(rangeval, c(0,0.25,0.5,0.75,1))
basisobj

remove(basisobj)
basisobj = create.polygonal.basis(rangeval, c(0,0.25,0.5,0.75,1), c(1,5))
basisobj

eval.penalty(basisobj, 0)
eval.penalty(basisobj, 1)

#  -------------------------  power basis  --------------------

rangeval = c(1e-1,1)

remove(basisobj)
basisobj = create.power.basis()
basisobj

remove(basisobj)
basisobj = create.power.basis(rangeval)
basisobj

remove(basisobj)
basisobj = create.power.basis(rangeval, nbasis, seq(0,2,0.5))
basisobj

plot(basisobj)

eval.penalty(basisobj, 2)

remove(basisobj)
basisobj = create.power.basis(rangeval, nbasis, seq(0,2,0.5), c(1,5))
basisobj

plot(basisobj)

eval.penalty(basisobj, 2)


