## This file runs the R 'sii' code on the examples provided with The
## 'Official' C Program Developer-F¢S Kit available at-A
## http://www.sii.to/html/programs.html to ensure that it yields the
## correct results.

library(SII)

## OCTAVE.TST	0.491	Octave Procedure
OCTAVE.TST=sii(
  speech=c(50.00,40,40,30,20,0.0),
  noise=c(70,65,45,25,1.00,-15.00),
  threshold=c(0.00,0,0,0,0,7.10),
  method="octave"
  )
stopifnot(round(OCTAVE.TST$sii,3)==0.491)

## OCTAVE_1.TST	0.323	Octave Procedure With Alternative Band Importance
OCTAVE_1.TST=sii(
  speech=c(50.00,40,40,30,20,0.0),
  noise=c(70,65,45,25,1.00,-15.00),
  threshold=c(0.00,0,0,0,0,7.10),
  importance=c(0,0,1,0,0,0),
  method="octave"
  )
stopifnot(round(OCTAVE_1.TST$sii,3)==0.323)

## TO.TST	0.445	1/3-Octave Procedure
TO.TST=sii(
  speech=c(90,5,40,40,40,40,40,40,40,40,40,40,40,40,-10,-10,-10,-10),
  noise=c(10,-10,-10,75,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10,10,10),
  threshold=c(90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  method="one-third octave"
  )
stopifnot(round(TO.TST$sii,3)==0.445)

## TO_1.TST	0.438	1/3-Octave Procedure With Alternative Band Importance Function 
TO_1.TST=sii(
  speech=c(90,5,40,40,40,40,40,40,40,40,40,40,40,40,-10,-10,-10,-10),
  noise=c(10,-10,-10,75,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10,10,10),
  threshold=c(90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  importance=c(0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0),
  method="one-third octave"
  )
stopifnot(round(TO_1.TST$sii,3)==0.438)


## CB.TST		0.273	Critical Band Procedure
CB.TST=sii(
  speech=c(-10,-10,-10,90,5,40,40,40,40,40,40,40,40,40,40,40,40,-10,-10,-10,-10),
  noise=c(10,10,10,10,-10,-10,70,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10,10,10),
  threshold=c(0,0,0,90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  )
stopifnot(round(CB.TST$sii,3)==0.273)

## CB_1.TST	0.410	Critical Band Procedure With Alternative Band Importance Function
CB_1.TST=sii(
  speech=c(-10,-10,-10,90,5,40,40,40,40,40,40,40,40,40,40,40,40,-10,-10,-10,-10),
  noise=c(10,10,10,10,-10,-10,70,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10,10,10),
  threshold=c(0,0,0,90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  importance=c(0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0,0,0)
  )
stopifnot(round(CB_1.TST$sii,3)==0.410)

# ECB.TST	0.278	Equally-Contributing Critical Band Procedure
ECB.TST=sii(
  speech=c(-10,90,5,40,40,40,40,40,40,40,40,40,40,40,40,-10,-10),
  noise=c(10,10,-10,-10,70,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10),
  threshold=c(0,90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  method="equal-contributing"
  )
stopifnot(round(ECB.TST$sii,3)==0.278)

# ECB_1.TST	0.410	Equally-Contribution Critical Band Procedure With Alternative Band Importance Function
ECB_1.TST=sii(
  speech=c(-10,90,5,40,40,40,40,40,40,40,40,40,40,40,40,-10,-10),
  noise=c(10,10,-10,-10,70,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10),
  threshold=c(0,90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  importance=c(0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0),
  method="equal-contributing"
  )
stopifnot(round(ECB_1.TST$sii,3)==0.410)

