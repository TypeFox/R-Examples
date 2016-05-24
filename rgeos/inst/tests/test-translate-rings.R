library(testthat)
library(rgeos)

setScale()

context("Translate Rings")

test_that("translate linear ring", {

    lr1 = readWKT("LINEARRING (1 1, 1 2, 2 2, 2 1, 1 1)")
    lr2 = readWKT("LINEARRING (1 1, 2 1, 2 2, 1 2, 1 1)")
    gclr1 = readWKT("GEOMETRYCOLLECTION( LINEARRING (1 1, 1 2, 2 2, 2 1, 1 1), LINEARRING (1 1, 1 2, 2 2, 2 1, 1 1) )")
    gclr2 = readWKT("GEOMETRYCOLLECTION( LINEARRING (1 1, 1 2, 2 2, 2 1, 1 1), LINEARRING (1 1, 2 1, 2 2, 1 2, 1 1) )")
    gclr3 = readWKT("GEOMETRYCOLLECTION( LINEARRING (1 1, 2 1, 2 2, 1 2, 1 1), LINEARRING (1 1, 2 1, 2 2, 1 2, 1 1) )")
    
    Ring11 = Ring(cbind( x=c(1,1,2,2,1),y=c(1,2,2,1,1) ),ID="1")
    Ring12 = Ring(cbind( x=c(1,1,2,2,1),y=c(1,2,2,1,1) ),ID="2")
    Ring21 = Ring(cbind( x=c(1,2,2,1,1),y=c(1,1,2,2,1) ),ID="1")
    Ring22 = Ring(cbind( x=c(1,2,2,1,1),y=c(1,1,2,2,1) ),ID="2")
    
    splr1   = SpatialRings( list(Ring11) ); #rownames(splr1@bbox) = c("x","y")
    splr2   = SpatialRings( list(Ring21) ); #rownames(splr2@bbox) = c("x","y")
    spgclr1 = SpatialRings( list(Ring11,Ring12) ); #rownames(spgclr1@bbox) = c("x","y")
    spgclr2 = SpatialRings( list(Ring11,Ring22) ); #rownames(spgclr2@bbox) = c("x","y")
    spgclr3 = SpatialRings( list(Ring21,Ring22) ); #rownames(spgclr3@bbox) = c("x","y")
    

    expect_that( lr1  , is_identical_to(splr1) )
    expect_that( lr2  , is_identical_to(splr2) )
    expect_that( gclr1, is_identical_to(spgclr1) )
    expect_that( gclr2, is_identical_to(spgclr2) )
    expect_that( gclr3, is_identical_to(spgclr3) )
    
    expect_that( splr1  , is_identical_to( translate(splr1)))
    expect_that( splr2  , is_identical_to( translate(splr2)))
    expect_that( spgclr1, is_identical_to( translate(spgclr1)))
    expect_that( spgclr2, is_identical_to( translate(spgclr2)))
    expect_that( spgclr3, is_identical_to( translate(spgclr3)))
    
    expect_that( lr1  , is_identical_to( translate(lr1) ))
    expect_that( lr2  , is_identical_to( translate(lr2) ))
    expect_that( gclr1, is_identical_to( translate(gclr1) ))
    expect_that( gclr2, is_identical_to( translate(gclr2) ))
    expect_that( gclr3, is_identical_to( translate(gclr3) ))

})