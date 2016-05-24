library(testthat)
library(rgeos)

setScale()

context("Translate Lines")

test_that("translate lines", {

    l = readWKT("LINESTRING (1 1, 2 2, 3 3, 4 4)")

    ml1 = readWKT("MULTILINESTRING ((1 1, 2 2, 3 3, 4 4),(1 1, 2 2, 3 3, 4 4))")
    ml2 = readWKT("MULTILINESTRING ((1 1, 2 2, 3 3, 4 4),(4 1, 3 2, 2 3, 1 4))")

    gcl1 = readWKT("GEOMETRYCOLLECTION( LINESTRING (1 1, 2 2, 3 3, 4 4), LINESTRING (1 1, 2 2, 3 3, 4 4) )")
    gcl2 = readWKT("GEOMETRYCOLLECTION( LINESTRING (1 1, 2 2, 3 3, 4 4), MULTILINESTRING ((1 1, 2 2, 3 3, 4 4),(4 1, 3 2, 2 3, 1 4)), LINESTRING (1 1, 2 2, 3 3, 4 4) )")


    Line1 = Line(cbind( x=1:4,y=1:4 ))
    Line2 = Line(cbind( x=4:1,y=1:4 ))
    
    Linesl = Lines( list(Line1), ID = "1" )
    Linesl2 = Lines( list(Line1), ID = "2" )
    
    Linesml1 = Lines( list(Line1, Line1), ID = "1" )
    Linesml2 = Lines( list(Line1, Line2), ID = "1" )
    
    #FIXME - weirdness with rownames in the bbox
    spl    = SpatialLines( list(Linesl) ); rownames(spl@bbox) = c("x","y")
    spml1  = SpatialLines( list(Linesml1) ); rownames(spml1@bbox) = c("x","y")
    spml2  = SpatialLines( list(Linesml2) ); rownames(spml2@bbox) = c("x","y")
    
    spgcl1 = SpatialLines( list(Linesl,Linesl2) ); rownames(spgcl1@bbox) = c("x","y")
    Linesml2@ID = "2"
    Linesl2@ID = "3"
    spgcl2 = SpatialLines( list(Linesl,Linesml2,Linesl2) ); rownames(spgcl2@bbox) = c("x","y")


    expect_that( l   , is_identical_to(spl) )
    expect_that( ml1 , is_identical_to(spml1) )
    expect_that( ml2 , is_identical_to(spml2) )
    expect_that( gcl1, is_identical_to(spgcl1) )
    expect_that( gcl2, is_identical_to(spgcl2) )
    
    expect_that( spl   , is_identical_to( translate(spl)))
    expect_that( spml1 , is_identical_to( translate(spml1)))
    expect_that( spml2 , is_identical_to( translate(spml2)))
    expect_that( spgcl1, is_identical_to( translate(spgcl1)))
    expect_that( spgcl2, is_identical_to( translate(spgcl2)))
    
    expect_that( l   , is_identical_to( translate(l) ))
    expect_that( ml1 , is_identical_to( translate(ml1) ))
    expect_that( ml2 , is_identical_to( translate(ml2) ))
    expect_that( gcl1, is_identical_to( translate(gcl1) ))
    expect_that( gcl2, is_identical_to( translate(gcl2) ))

})