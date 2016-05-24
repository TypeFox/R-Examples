library(testthat)
library(rgeos)

setScale()

context("Translate Polygons")

test_that("translate simple polygon", {

    p=readWKT("POLYGON((1 1,5 1,5 5,1 5,1 1))")
    ph1=readWKT("POLYGON((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2))")
    ph2=readWKT("POLYGON((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2), (3 3,3 4,4 4,4 3,3 3) ) ")
    
    mp=readWKT("MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1)),((3 5,5 7, 1 7, 3 5)))")
    mph1=readWKT("MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)),((3 5,5 7, 1 7, 3 5)))")
    mph2=readWKT("MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1)),((3 5,5 7, 1 7, 3 5),(3 5.5,4 6.5, 2 6.5, 3 5.5) ))")
    mph3=readWKT("MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)),((3 5,5 7, 1 7, 3 5),(3 5.5,4 6.5, 2 6.5, 3 5.5)) )")
    
    Poly1 = Polygon(list(x=c(1,5,5,1,1),y=c(1,1,5,5,1)), hole=FALSE)
    Poly2 = Polygon(list(x=c(2,2,3,3,2),y=c(2,3,3,2,2)), hole=TRUE)
    Poly3 = Polygon(list(x=c(3,3,4,4,3),y=c(3,4,4,3,3)), hole=TRUE)
    Poly4 = Polygon(list(x=c(3,5,1,3),y=c(5,7,7,5)), hole=FALSE)
    Poly5 = Polygon(list(x=c(3,4,2,3),y=c(5.5,6.5,6.5,5.5)), hole=TRUE)
    Poly6 = Polygon(list(x=c(5,7,7,5),y=c(3,5,1,3)), hole=FALSE)
    
    Polysp = Polygons(list(Poly1), ID="1")
    Polysph1 = Polygons(list(Poly1,Poly2), ID="1")
    Polysph2 = Polygons(list(Poly1,Poly2,Poly3), ID="1")
    
    Polysmp = Polygons(list(Poly1,Poly4), ID="1" )
    Polysmph1 = Polygons(list(Poly1,Poly2,Poly4), ID="1" )
    Polysmph2 = Polygons(list(Poly1,Poly4,Poly5), ID="1" )
    Polysmph3 = Polygons(list(Poly1,Poly2,Poly4,Poly5), ID="1" )
    
    comment(Polysp) <- "0"
    comment(Polysph1) <- "0 1"
    comment(Polysph2) <- "0 1 1"
    comment(Polysmp) <- "0 0"
    comment(Polysmph1) <- "0 1 0"
    comment(Polysmph2) <- "0 0 2"
    comment(Polysmph3) <- "0 1 0 3"

    spp     = SpatialPolygons( list(Polysp) )
    spph1   = SpatialPolygons( list(Polysph1) )
    spph2   = SpatialPolygons( list(Polysph2) )
    spmp    = SpatialPolygons( list(Polysmp) )
    spmph1  = SpatialPolygons( list(Polysmph1) )
    spmph2  = SpatialPolygons( list(Polysmph2) )
    spmph3  = SpatialPolygons( list(Polysmph3) )
    
    expect_that( identical(p    , spp), is_true())
    expect_that( identical(ph1  , spph1), is_true())
    expect_that( identical(ph2  , spph2), is_true())
    expect_that( identical(mp   , spmp), is_true())
    expect_that( identical(mph1 , spmph1), is_true())
    expect_that( identical(mph2 , spmph2), is_true())
    expect_that( identical(mph3 , spmph3), is_true())
    
    expect_that( spp   , is_identical_to( translate(spp)))
    expect_that( spph1 , is_identical_to( translate(spph1)))
    expect_that( spph2 , is_identical_to( translate(spph2)))
    expect_that( spmp  , is_identical_to( translate(spmp)))
    expect_that( spmph1, is_identical_to( translate(spmph1)))
    expect_that( spmph2, is_identical_to( translate(spmph2)))
    expect_that( spmph3, is_identical_to( translate(spmph3)))
    
    expect_that( p   , is_identical_to( translate(p)))
    expect_that( ph1 , is_identical_to( translate(ph1)))
    expect_that( ph2 , is_identical_to( translate(ph2)))
    expect_that( mp  , is_identical_to( translate(mp)))
    expect_that( mph1, is_identical_to( translate(mph1)))
    expect_that( mph2, is_identical_to( translate(mph2)))
    expect_that( mph3, is_identical_to( translate(mph3)))
    
})
