library(testthat)
library(rgeos)

setScale()

context("Translate Points")

test_that("translate points", {
    
    p = readWKT("POINT(1 1)")
    mp = readWKT("MULTIPOINT(1 1, 2 2, 3 3, 4 4, 5 5)") 
    gcp1 = readWKT("GEOMETRYCOLLECTION( POINT(1 1), POINT(2 2), POINT(3 3), POINT(4 4), POINT(5 5))")
    gcp2 = readWKT("GEOMETRYCOLLECTION( POINT(1 1), POINT(2 2), MULTIPOINT(3 3, 4 4, 5 5))")
    gcp3 = readWKT("GEOMETRYCOLLECTION( POINT(1 1), POINT(2 2), MULTIPOINT(3 3, 4 4),POINT(5 5))")
    gcp4 = readWKT("GEOMETRYCOLLECTION( MULTIPOINT(1 1, 2 2), MULTIPOINT(3 3, 4 4, 5 5))")
    gcp5 = readWKT("GEOMETRYCOLLECTION( MULTIPOINT(1 1, 2 2), MULTIPOINT(3 3, 4 4),POINT(5 5))")
    
    spp = SpatialPoints(list(x=1,y=1))
    spmp = SpatialPoints(list(x=1:5,y=1:5))
    
    rownames(spp@coords) = c("1")
    expect_that( identical(p,spp), is_true())
    expect_that( identical(spp, translate(spp)), is_true())
    
    rownames(spmp@coords) = c("1","1","1","1","1")
    expect_that( identical(mp,spmp), is_true() )
    expect_that( identical(spmp,translate(spmp)), is_true())
    
    rownames(spmp@coords) = c("1","2","3","4","5")
    expect_that( identical(gcp1,spmp), is_true() )
	expect_that( identical(spmp,translate(spmp)), is_true())
    
    rownames(spmp@coords) = c("1","2","3","3","3")
    expect_that( identical(gcp2,spmp), is_true() )
    expect_that( identical(spmp,translate(spmp)), is_true())

    rownames(spmp@coords) = c("1","2","3","3","4")
    expect_that( identical(gcp3,spmp), is_true() )
    expect_that( identical(spmp,translate(spmp)), is_true())
    
    rownames(spmp@coords) = c("1","1","2","2","2")    
    expect_that( identical(gcp4,spmp), is_true() )
    expect_that( identical(spmp,translate(spmp)), is_true())
    
    rownames(spmp@coords) = c("1","1","2","2","3")
    expect_that( identical(gcp5,spmp), is_true() )
    expect_that( identical(spmp,translate(spmp)), is_true())

    
    expect_that( identical(p, translate(p) ), is_true())
    expect_that( identical(mp, translate(mp) ), is_true())
    expect_that( identical(gcp1, translate(gcp1) ), is_true())
    expect_that( identical(gcp2, translate(gcp2) ), is_true())
    expect_that( identical(gcp3, translate(gcp3) ), is_true())
    expect_that( identical(gcp4, translate(gcp4) ), is_true())
    expect_that( identical(gcp5, translate(gcp5) ), is_true())
    
})