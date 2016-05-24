# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the neldermead component
# originally written by Michael Baudin for Scilab :
# "Nelder-Mead User's Manual", 2010, Consortium Scilab - Digiteo,
# Michael Baudin, http://wiki.scilab.org/The_Nelder-Mead_Component

neldermead <- function(optbase, method, simplex0, simplex0method, 
  simplex0length, simplexsize0, simplexopt, historysimplex, coords0, rho, chi, 
  gamma, sigma, tolfstdeviation, tolfstdeviationmethod, tolsimplexizeabsolute, 
  tolsimplexizerelative, tolsimplexizemethod, toldeltafv, tolssizedeltafvmethod,
  simplex0deltausual, simplex0deltazero, restartsimplexmethod, restartmax, 
  restarteps, restartstep, restartnb, restartflag, restartdetection, 
  kelleystagnationflag, kelleynormalizationflag, kelleystagnationalpha0, 
  kelleyalpha, startupflag, boxnbpoints, boxnbpointseff, boxineqscaling, 
  checkcostfunction, scalingsimplex0, guinalphamin, boxboundsalpha, 
  boxtermination, boxtolf, boxnbmatch, boxkount, boxreflect, tolvarianceflag, 
  tolabsolutevariance, tolrelativevariance, variancesimplex0, mymethod, 
  myterminate, myterminateflag, greedy, output, exitflag){
  
  newobj<-list(optbase=optimbase(),method='variable',simplex0=simplex(),
    simplex0method='axes',simplex0length=1.0,simplexsize0=0,simplexopt=NULL,
    historysimplex=list(),coords0=NULL,rho=1.0,chi=2.0,gamma=.5,sigma=.5,
    tolfstdeviation=0,tolfstdeviationmethod=FALSE,tolsimplexizeabsolute=0,
    tolsimplexizerelative=.Machine$double.eps,tolsimplexizemethod=TRUE,
    toldeltafv=.Machine$double.eps,tolssizedeltafvmethod=FALSE,
    simplex0deltausual=.05,simplex0deltazero=.0075,restartsimplexmethod='oriented',
    restartmax=3,restarteps=.Machine$double.eps,restartstep=1,restartnb=0,
    restartflag=FALSE,restartdetection='oneill',kelleystagnationflag=FALSE,
    kelleynormalizationflag=TRUE,kelleystagnationalpha0=1.e-4,kelleyalpha=1.e-4,
    startupflag=FALSE,boxnbpoints='2n',boxnbpointseff=0,boxineqscaling=.5,
    checkcostfunction=TRUE,scalingsimplex0='tox0',guinalphamin=1.e-5,
    boxboundsalpha=1.e-6,boxtermination=FALSE,boxtolf=1.e-5,boxnbmatch=5,
    boxkount=0,boxreflect=1.3,tolvarianceflag=FALSE,tolabsolutevariance=0,
    tolrelativevariance=.Machine$double.eps,variancesimplex0=0,mymethod=NULL,
    myterminate=NULL,myterminateflag=FALSE,greedy=FALSE,output=list(),
    exitflag=FALSE)
  
  # Initial optimbase object
  if (!missing(optbase))
    newobj$optbase <- neldermead.set(this=newobj,key='optbase',value=optbase)
  
  # Possible values "variable", "fixed"
  if (!missing(method))
    newobj <- neldermead.set(this=newobj,key='method',value=method)
  
  # Initial simplex
  if (!missing(simplex0))
    newobj <- neldermead.set(this=newobj,key='simplex0',value=simplex0)
  
  # Possible values: "axes", "spendley", "pfeffer"
  if (!missing(simplex0method))
    newobj <- neldermead.set(this=newobj,key='simplex0method',value=simplex0method)
  if (!missing(simplex0length))
    newobj <- neldermead.set(this=newobj,key='simplex0length',value=simplex0length)
  
  
  # Initial size of the simplex, for the tolerance on the simplex size
  if (!missing(simplexsize0))
    newobj <- neldermead.set(this=newobj,key='simplexsize0',value=simplexsize0)
  
  # The optimum simplex, after one optimization process
  if (!missing(simplexopt))
    newobj <- neldermead.set(this=newobj,key='simplexopt',value=simplexopt)
  
  # The history of the simplex optimization
  if (!missing(historysimplex))
    newobj <- neldermead.set(this=newobj,key='historysimplex',value=historysimplex)
  
  # The coordinates of the initial simplex, given by the user
  if (!missing(coords0))
    newobj <- neldermead.set(this=newobj,key='coords0',value=coords0)
  
  # Reflection factor: rho
  if (!missing(rho))
    newobj <- neldermead.set(this=newobj,key='rho',value=rho)
  
  # Expansion factor: chi
  if (!missing(chi))
    newobj <- neldermead.set(this=newobj,key='chi',value=chi)
  
  # Contraction factor: gamma
  if (!missing(gamma))
    newobj <- neldermead.set(this=newobj,key='gamma',value=gamma)
  
  # Shrinkage factor: sigma
  if (!missing(sigma))
    newobj <- neldermead.set(this=newobj,key='sigma',value=sigma)
  
  # The tolerance for the standard deviation
  if (!missing(tolfstdeviation))
    newobj <- neldermead.set(this=newobj,key='tolfstdeviation',value=tolfstdeviation)
  
  # Possible values: TRUE, FALSE
  if (!missing(tolfstdeviationmethod))
    newobj <- neldermead.set(this=newobj,key='tolfstdeviationmethod',value=tolfstdeviationmethod)
  
  # The absolute tolerance for the simplex size
  if (!missing(tolsimplexizeabsolute))
    newobj <- neldermead.set(this=newobj,key='tolsimplexizeabsolute',value=tolsimplexizeabsolute)
  
  # The relative tolerance for the simplex size
  if (!missing(tolsimplexizerelative))
    newobj <- neldermead.set(this=newobj,key='tolsimplexizerelative',value=tolsimplexizerelative)
  
  # Possible values: TRUE, FALSE
  # Note: If the simplex method converges, the simplex size is near zero.
  if (!missing(tolsimplexizemethod))
    newobj <- neldermead.set(this=newobj,key='tolsimplexizemethod',value=tolsimplexizemethod)
  
  # The tolerance for the function value delta
  if (!missing(toldeltafv))
    newobj <- neldermead.set(this=newobj,key='toldeltafv',value=toldeltafv)
  
  # Possible values: TRUE, FALSE
  if (!missing(tolssizedeltafvmethod))
    newobj <- neldermead.set(this=newobj,key='tolssizedeltafvmethod',value=tolssizedeltafvmethod)
  
  # The value used in Pfeffer method initial simplex computation for non-zero parameters
  if (!missing(simplex0deltausual))
    newobj <- neldermead.set(this=newobj,key='simplex0deltausual',value=simplex0deltausual)
  
  # The value used in Pfeffer method initial simplex computation for zero parameters
  if (!missing(simplex0deltazero))
    newobj <- neldermead.set(this=newobj,key='simplex0deltazero',value=simplex0deltazero)
  
  # Possible values: "oriented", "axes", "spendley", "pfeffer"
  if (!missing(restartsimplexmethod))
    newobj <- neldermead.set(this=newobj,key='restartsimplexmethod',value=restartsimplexmethod)
  
  # The maximum number of restarts
  if (!missing(restartmax))
    newobj <- neldermead.set(this=newobj,key='restartmax',value=restartmax)
  
  # The epsilon value for O'Neill restart detection
  if (!missing(restarteps))
    newobj <- neldermead.set(this=newobj,key='restarteps',value=restarteps)
  
  # The step length for O'Neill restart detection
  if (!missing(restartstep))
    newobj <- neldermead.set(this=newobj,key='restartstep',value=restartstep)
  
  # Number of restarts performed
  if (!missing(restartnb))
    newobj <- neldermead.set(this=newobj,key='restartnb',value=restartnb)
  
  # Possible values: TRUE, FALSE
  if (!missing(restartflag))
    newobj <- neldermead.set(this=newobj,key='restartflag',value=restartflag)
  
  # Type of restart detection method: "kelley", "oneill"
  if (!missing(restartdetection))
    newobj <- neldermead.set(this=newobj,key='restartdetection',value=restartdetection)
  
  # The Kelley stagnation detection in termination criteria: TRUE/FALSE
  # (i.e. sufficient decrease of function value)
  if (!missing(kelleystagnationflag))
    newobj <- neldermead.set(this=newobj,key='kelleystagnationflag',value=kelleystagnationflag)
  
  # The Kelley stagnation detection can be normalized or not.
  # Note:
  # * in the 1997 paper "Detection and Remediation of Stagnation in Nelder-Mead
  #   algorithm", Kelley uses the constant value of 1.e-4.
  # * in the 1999 book "Iterative Methods for Optimization", Kelley uses normalization.
  # Results are slightly changed, as indicated in the book/paper (the modification is
  # not mentioned, but the iteration number when the restart is performed
  # is modified).
  if (!missing(kelleynormalizationflag))
    newobj <- neldermead.set(this=newobj,key='kelleynormalizationflag',value=kelleynormalizationflag)
  
  # The Kelley stagnation detection parameter
  if (!missing(kelleystagnationalpha0))
    newobj <- neldermead.set(this=newobj,key='kelleystagnationalpha0',value=kelleystagnationalpha0)
  
  # The current value of Kelley's alpha, after normalization, if required
  if (!missing(kelleyalpha))
    newobj <- neldermead.set(this=newobj,key='kelleyalpha',value=kelleyalpha)
  
  # Set to TRUE when the startup has been performed
  if (!missing(startupflag))
    newobj <- neldermead.set(this=newobj,key='startupflag',value=startupflag)
  
  # Number of points required in the simplex (for Box method)
  if (!missing(boxnbpoints))
    newobj <- neldermead.set(this=newobj,key='boxnbpoints',value=boxnbpoints)
  
  # Effective number of points required in the simplex (for Box method)
  if (!missing(boxnbpointseff))
    newobj <- neldermead.set(this=newobj,key='boxnbpointseff',value=boxnbpointseff)
  
  # The scaling coefficient in nonlinear inequality constraints
  # in Box method, in (0,1) range
  if (!missing(boxineqscaling))
    newobj <- neldermead.set(this=newobj,key='boxineqscaling',value=boxineqscaling)
  
  # Set to FALSE to disable the checking of the connection of the cost function
  if (!missing(checkcostfunction))
    newobj <- neldermead.set(this=newobj,key='checkcostfunction',value=checkcostfunction)
  
  # The scaling algorithm: "tox0", "tocentroid"
  if (!missing(scalingsimplex0))
    newobj <- neldermead.set(this=newobj,key='scalingsimplex0',value=scalingsimplex0)
  
  # Minimum alpha for constraints scaling
  if (!missing(guinalphamin))
    newobj <- neldermead.set(this=newobj,key='guinalphamin',value=guinalphamin)
  
  # Box's alpha coefficient for bounds constraints.
  # The value used in Box's paper was 1.e-6 (delta in
  # Richardson and Kuester's algorithm 454)
  if (!missing(boxboundsalpha))
    newobj <- neldermead.set(this=newobj,key='boxboundsalpha',value=boxboundsalpha)
  
  # Set to TRUE to enable Box termination criteria
  if (!missing(boxtermination))
    newobj <- neldermead.set(this=newobj,key='boxtermination',value=boxtermination)
  
  # The absolute tolerance on function value in Box termination criteria (beta in
  # Richardson and Kuester's algorithm 454)
  if (!missing(boxtolf))
    newobj <- neldermead.set(this=newobj,key='boxtolf',value=boxtolf)
  
  # The number of consecutive match in Box termination criteria (gamma in
  # Richardson and Kuester's algorithm 454)
  if (!missing(boxnbmatch))
    newobj <- neldermead.set(this=newobj,key='boxnbmatch',value=boxnbmatch)
  
  # Current number of consecutive match
  if (!missing(boxkount))
    newobj <- neldermead.set(this=newobj,key='boxkount',value=boxkount)
  
  # Box reflection/expansion factor
  if (!missing(boxreflect))
    newobj <- neldermead.set(this=newobj,key='boxreflect',value=boxreflect)
  
  # Set to TRUE to enable tolerance on variance
  if (!missing(tolvarianceflag))
    newobj <- neldermead.set(this=newobj,key='tolvarianceflag',value=tolvarianceflag)
  
  # Absolute tolerance on variance
  if (!missing(tolabsolutevariance))
    newobj <- neldermead.set(this=newobj,key='tolabsolutevariance',value=tolabsolutevariance)
  
  # Relative tolerance on variance
  if (!missing(tolrelativevariance))
    newobj <- neldermead.set(this=newobj,key='tolrelativevariance',value=tolrelativevariance)
  
  # The variance of the initial simplex
  if (!missing(variancesimplex0))
    newobj <- neldermead.set(this=newobj,key='variancesimplex0',value=variancesimplex0)
  
  # User-defined algorithm
  if (!missing(mymethod))
    newobj <- neldermead.set(this=newobj,key='mymethod',value=mymethod)
  
  # User-defined termination criteria
  if (!missing(myterminate))
    newobj <- neldermead.set(this=newobj,key='myterminate',value=myterminate)
  
  # Flag to enable the user-defined termination criteria
  if (!missing(myterminateflag))
    newobj <- neldermead.set(this=newobj,key='myterminateflag',value=myterminateflag)
  
  # Set to TRUE to enable greedy Nelder-Mead
  if (!missing(greedy))
    newobj <- neldermead.set(this=newobj,key='greedy',value=greedy)
  
  # Customizable output for specialized function
  if (!missing(output))
    newobj <- neldermead.set(this=newobj,key='output',value=output)
  
  # Customizable exit flag for specialized function
  if (!missing(exitflag))
    newobj <- neldermead.set(this=newobj,key='exitflag',value=exitflag)
  
  class(newobj) <- 'neldermead'
  
  return(newobj)
  
}

