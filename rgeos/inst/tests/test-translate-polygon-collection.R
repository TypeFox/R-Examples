library(testthat)
library(rgeos)

setScale()

context("Translation Polygon Collections")

test_that("translate polygon collection", {

    gcp1=readWKT("GEOMETRYCOLLECTION( POLYGON((1 1,5 1,5 5,1 5,1 1)), POLYGON((3 5,5 7, 1 7, 3 5)), POLYGON((5 3,7 5,7 1,5 3)) )")
    gcp2=readWKT("GEOMETRYCOLLECTION( MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1)), ((3 5,5 7, 1 7, 3 5)) ), POLYGON((5 3,7 5,7 1,5 3)) )")
    gcp3=readWKT("GEOMETRYCOLLECTION( POLYGON((1 1,5 1,5 5,1 5,1 1)), MULTIPOLYGON( ((3 5,5 7, 1 7, 3 5)), ((5 3,7 5,7 1,5 3)) ))")
    gcp4=readWKT("GEOMETRYCOLLECTION( MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1)), ((3 5,5 7, 1 7, 3 5)), ((5 3,7 5,7 1,5 3)) ))")
    
    gcph1=readWKT("GEOMETRYCOLLECTION( POLYGON((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)), POLYGON((3 5,5 7, 1 7, 3 5),(3 5.5,4 6.5, 2 6.5, 3 5.5)), POLYGON((5 3,7 5,7 1,5 3)) )")
    gcph2=readWKT("GEOMETRYCOLLECTION( MULTIPOLYGON( ((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)), ((3 5,5 7, 1 7, 3 5),(3 5.5,4 6.5, 2 6.5, 3 5.5)) ), POLYGON((5 3,7 5,7 1,5 3)) )")
    gcph3=readWKT("GEOMETRYCOLLECTION( POLYGON((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)), MULTIPOLYGON( ((3 5,5 7, 1 7, 3 5),(3 5.5,4 6.5, 2 6.5, 3 5.5)), ((5 3,7 5,7 1,5 3)) ))")
    
    
    Poly1 = Polygon(list(x=c(1,5,5,1,1),y=c(1,1,5,5,1)), hole=FALSE)
    Poly2 = Polygon(list(x=c(2,2,3,3,2),y=c(2,3,3,2,2)), hole=TRUE)
    Poly3 = Polygon(list(x=c(3,3,4,4,3),y=c(3,4,4,3,3)), hole=TRUE)
    Poly4 = Polygon(list(x=c(3,5,1,3),y=c(5,7,7,5)), hole=FALSE)
    Poly5 = Polygon(list(x=c(3,4,2,3),y=c(5.5,6.5,6.5,5.5)), hole=TRUE)
    Poly6 = Polygon(list(x=c(5,7,7,5),y=c(3,5,1,3)), hole=FALSE)
   
    Polygcp11 = Polygons(list(Poly1), ID="1")
    Polygcp12 = Polygons(list(Poly4), ID="2")
    Polygcp13 = Polygons(list(Poly6), ID="3")
    Polygcp21 = Polygons(list(Poly1,Poly4), ID="1")
    Polygcp22 = Polygons(list(Poly6), ID="2")
    Polygcp31 = Polygons(list(Poly1), ID="1")
    Polygcp32 = Polygons(list(Poly4,Poly6), ID="2")
    Polygcp4  = Polygons(list(Poly1,Poly4,Poly6), ID="1")
    
    Polygcph11 = Polygons(list(Poly1,Poly2), ID="1")
    Polygcph12 = Polygons(list(Poly4,Poly5), ID="2")
    Polygcph13 = Polygons(list(Poly6), ID="3")
    Polygcph21 = Polygons(list(Poly1,Poly2,Poly4,Poly5), ID="1")
    Polygcph22 = Polygons(list(Poly6), ID="2")
    Polygcph31 = Polygons(list(Poly1,Poly2), ID="1")
    Polygcph32 = Polygons(list(Poly4,Poly5,Poly6), ID="2")

    comment(Polygcp11) <- "0"
    comment(Polygcp12) <- "0"
    comment(Polygcp13) <- "0"
    comment(Polygcp21) <- "0 0"
    comment(Polygcp22) <- "0"
    comment(Polygcp31) <- "0"
    comment(Polygcp32) <- "0 0"
    comment(Polygcp4) <- "0 0 0"

    comment(Polygcph11) <- "0 1"
    comment(Polygcph12) <- "0 1"
    comment(Polygcph13) <- "0"
    comment(Polygcph21) <- "0 1 0 3"
    comment(Polygcph22) <- "0"
    comment(Polygcph31) <- "0 1"
    comment(Polygcph32) <- "0 1 0"

    spgcp1  = SpatialPolygons( list(Polygcp11,Polygcp12,Polygcp13) )
    spgcp2  = SpatialPolygons( list(Polygcp21,Polygcp22) )
    spgcp3  = SpatialPolygons( list(Polygcp31,Polygcp32) )
    spgcp4  = SpatialPolygons( list(Polygcp4) )
    spgcph1 = SpatialPolygons( list(Polygcph11,Polygcph12,Polygcph13) )
    spgcph2 = SpatialPolygons( list(Polygcph21,Polygcph22) )
    spgcph3 = SpatialPolygons( list(Polygcph31,Polygcph32) )
    
    expect_that( gcp1 , is_identical_to(spgcp1) )
    expect_that( gcp2 , is_identical_to(spgcp2) )
    expect_that( gcp3 , is_identical_to(spgcp3) )
    expect_that( gcp4 , is_identical_to(spgcp4) )
    expect_that( gcph1, is_identical_to(spgcph1) )
    expect_that( gcph2, is_identical_to(spgcph2) )
    expect_that( gcph3, is_identical_to(spgcph3) )
    
    expect_that( spgcp1 , is_identical_to( translate(spgcp1)))
    expect_that( spgcp2 , is_identical_to( translate(spgcp2)))
    expect_that( spgcp3 , is_identical_to( translate(spgcp3)))
    expect_that( spgcp4 , is_identical_to( translate(spgcp4)))
    expect_that( spgcph1, is_identical_to( translate(spgcph1)))
    expect_that( spgcph2, is_identical_to( translate(spgcph2)))
    expect_that( spgcph3, is_identical_to( translate(spgcph3)))
    
    expect_that( gcp1 , is_identical_to( translate(gcp1)))
    expect_that( gcp2 , is_identical_to( translate(gcp2)))
    expect_that( gcp3 , is_identical_to( translate(gcp3)))
    expect_that( gcp4 , is_identical_to( translate(gcp4)))
    expect_that( gcph1, is_identical_to( translate(gcph1)))
    expect_that( gcph2, is_identical_to( translate(gcph2)))
    expect_that( gcph3, is_identical_to( translate(gcph3)))
})
