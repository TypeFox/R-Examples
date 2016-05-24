library(pcalg)

xx <- TRUE
##################################################
## DAG / CPDAG
##################################################
## CPDAG 1: Paper Fig 1
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1, 0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6
z <- c(2,4); xx <- xx & gac(mFig1,x,y,z,type)$gac
z <- c(4,5); xx <- xx & gac(mFig1,x,y,z,type)$gac
z <- c(4,2,1); xx <- xx & gac(mFig1,x,y,z,type)$gac
z <- c(4,5,1); xx <- xx & gac(mFig1,x,y,z,type)$gac
z <- c(4,2,5); xx <- xx & gac(mFig1,x,y,z,type)$gac
z <- c(4,2,5,1); xx <- xx & gac(mFig1,x,y,z,type)$gac
z <- 2; xx <- xx & !gac(mFig1,x,y,z,type)$gac
z <- NULL; xx <- xx & !gac(mFig1,x,y,z,type)$gac

## CPDAG 2: Paper Fig 5a
mFig5a <- matrix(c(0,1,0,0,0, 1,0,1,0,0, 0,0,0,0,1, 0,0,1,0,0, 0,0,0,0,0), 5,5)
type <- "cpdag"
x <- c(1,5); y <- 4
z <- c(2,3);xx <- xx & gac(mFig5a,x,y,z,type)$gac
z <- 2;xx <- xx & !gac(mFig5a,x,y,z,type)$gac

## DAG 1 from Marloes' Talk
m <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6)
type <- "dag"
x <- 1; y <- 3
z <- NULL;xx <- xx & gac(m,x,y,z,type)$gac
z <- 2;xx <- xx & !gac(m,x,y,z,type)$gac
z <- 4;xx <- xx & gac(m,x,y,z,type)$gac
z <- 5;xx <- xx & !gac(m,x,y,z,type)$gac
z <- 6;xx <- xx & !gac(m,x,y,z,type)$gac
z <- c(4,5);xx <- xx & !gac(m,x,y,z,type)$gac

## DAG 2 from Marloes' Talk
m <- matrix(c(0,1,0,1,0,0, 0,0,0,0,0,0, 0,1,0,0,1,0, 0,0,0,0,1,0, 0,0,0,0,0,1, 0,0,0,0,0,0), 6,6)
type <- "dag"
x <- 4; y <- 6
z <- 1;xx <- xx & gac(m,x,y,z,type)$gac
z <- 3;xx <- xx & gac(m,x,y,z,type)$gac
z <- 5;xx <- xx & !gac(m,x,y,z,type)$gac
z <- c(1,5);xx <- xx & !gac(m,x,y,z,type)$gac
z <- c(1,2);xx <- xx & gac(m,x,y,z,type)$gac
z <- c(1,3);xx <- xx & gac(m,x,y,z,type)$gac
z <- 2;xx <- xx & !gac(m,x,y,z,type)$gac

##################################################
## PAG
##################################################
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
mFig3b <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
mFig3c <- matrix(c(0,3,0,0, 2,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
mFig4a <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,3,3,2, 0,0,2,0,2,2, 0,0,2,1,0,2, 0,0,1,3,3,0), 6,6)
mFig4b <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,0,3,2, 0,0,0,0,2,2, 0,0,2,3,0,2, 0,0,2,3,2,0), 6,6)
mFig5b <- matrix(c(0,1,0,0,0,0,0, 2,0,2,3,0,3,0, 0,1,0,0,0,0,0, 0,2,0,0,3,0,0, 0,0,0,2,0,2,3, 0,2,0,0,2,0,0, 0,0,0,0,2,0,0), 7,7)
type <- "pag"
x<-2; y<-4; z<-NULL; xx <- xx & !gac(mFig3a,x,y,z,type)$gac
x<-2; y<-4; z<-NULL; xx <- xx & !gac(mFig3b,x,y,z,type)$gac
x<-2; y<-4; z<-NULL; xx <- xx & gac(mFig3c,x,y,z,type)$gac
x<-3; y<-4; z<-NULL; xx <- xx & !gac(mFig4a,x,y,z,type)$gac
x<-3; y<-4; z<-6; xx <- xx & gac(mFig4a,x,y,z,type)$gac
x<-3; y<-4; z<-c(1,6); xx <- xx & gac(mFig4a,x,y,z,type)$gac
x<-3; y<-4; z<-c(2,6); xx <- xx & gac(mFig4a,x,y,z,type)$gac
x<-3; y<-4; z<-c(1,2,6); xx <- xx & gac(mFig4a,x,y,z,type)$gac
x<-3; y<-4; z<-NULL; xx <- xx & !gac(mFig4b,x,y,z,type)$gac
x<-3; y<-4; z<-6; xx <- xx & !gac(mFig4b,x,y,z,type)$gac
x<-3; y<-4; z<-c(5,6); xx <- xx & !gac(mFig4b,x,y,z,type)$gac
x<-c(2,7); y<-6; z<-NULL; xx <- xx & !gac(mFig5b,x,y,z,type)$gac
x<-c(2,7); y<-6; z<-c(4,5); xx <- xx & gac(mFig5b,x,y,z,type)$gac
x<-c(2,7); y<-6; z<-c(4,5,1); xx <- xx & gac(mFig5b,x,y,z,type)$gac
x<-c(2,7); y<-6; z<-c(4,5,3); xx <- xx & gac(mFig5b,x,y,z,type)$gac
x<-c(2,7); y<-6; z<-c(1,3,4,5); xx <- xx & gac(mFig5b,x,y,z,type)$gac

## PAG in Marloes' talk
m <- matrix(c(0,0,0,3,2,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 2,0,0,0,0,3,2, 3,2,2,0,0,0,3, 0,0,0,2,0,0,0, 0,0,0,2,2,0,0),7,7)
x <- c(5,6); y <- 7
z<-NULL; xx <- xx & !gac(m,x,y,z,type)$gac
z<-1; xx <- xx & !gac(m,x,y,z,type)$gac
z<-4; xx <- xx & !gac(m,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m,x,y,z,type)$gac
z<-3; xx <- xx & !gac(m,x,y,z,type)$gac
z<-c(2,3); xx <- xx & !gac(m,x,y,z,type)$gac
z<-c(1,4); xx <- xx & gac(m,x,y,z,type)$gac
z<-c(1,4,2); xx <- xx & gac(m,x,y,z,type)$gac
z<-c(1,4,3); xx <- xx & gac(m,x,y,z,type)$gac
z<-c(1,4,2,3); xx <- xx & gac(m,x,y,z,type)$gac

##################################################
## Tests from Ema
##################################################
load("gacData.rda")
type <- "pag"
m1 <-matrices[[1]] 
x <- 6; y <- 9
z<-NULL; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-1; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-3; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-4; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,3); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,3,8); xx <- xx & gac(m1,x,y,z,type)$gac
z<-c(2,3,7,8); xx <- xx & gac(m1,x,y,z,type)$gac
z<-c(2,3,5,8); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,3,5,7,8); xx <- xx & !gac(m1,x,y,z,type)$gac

x <- c(6,8); y <- 9
z<-NULL; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-1; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-3; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-4; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,3); xx <- xx & gac(m1,x,y,z,type)$gac
z<-c(2,3,4); xx <- xx & gac(m1,x,y,z,type)$gac
z<-c(2,3,7); xx <- xx & gac(m1,x,y,z,type)$gac
z<-c(2,3,5); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,3,5,7); xx <- xx & !gac(m1,x,y,z,type)$gac

x <- 3; y <- 1
z<-NULL; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-4; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-5; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-6; xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,6); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,8); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,7,8); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,5,8); xx <- xx & !gac(m1,x,y,z,type)$gac
z<-c(2,5,7,8); xx <- xx & !gac(m1,x,y,z,type)$gac

m2 <- matrices[[2]]
x <- 3; y <-1
z<-NULL; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-2; xx <- xx & gac(m2,x,y,z,type)$gac 
z<-4; xx <- xx & !gac(m2,x,y,z,type)$gac 
z<-c(2,8); xx <- xx & !gac(m2,x,y,z,type)$gac
z<-8; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-9; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-c(2,8,9); xx <- xx & !gac(m2,x,y,z,type)$gac
z<-c(2,5); xx <- xx & gac(m2,x,y,z,type)$gac

x <- c(3,9); y <- 1
z<-NULL; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-4; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-c(2,8); xx <- xx & !gac(m2,x,y,z,type)$gac
z<-8; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-9; xx <- xx & !gac(m2,x,y,z,type)$gac
z<-c(2,8,9); xx <- xx & !gac(m2,x,y,z,type)$gac
z<-c(2,5); xx <- xx & !gac(m2,x,y,z,type)$gac

m3 <- matrices[[3]]
x <- 1; y <- 9
z<-NULL; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-3; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-5; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-7; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-8; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-c(2,3); xx <- xx & gac(m3,x,y,z,type)$gac
z<-c(5,7); xx <- xx & gac(m3,x,y,z,type)$gac

x <- 1; y <- 8
z<-NULL; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-2; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-3; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-5; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-7; xx <- xx & gac(m3,x,y,z,type)$gac
z<-9; xx <- xx & !gac(m3,x,y,z,type)$gac
z<-c(2,3); xx <- xx & gac(m3,x,y,z,type)$gac
z<-c(5,9); xx <- xx & !gac(m3,x,y,z,type)$gac 

if (!xx) stop("Problem when testing function gac.")
