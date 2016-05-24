#> nobarsMM(y~x+anything(ga|bu)) -> y ~ x OK
#> noNonSpatialbarsMM(y~x+anything(ga|bu)) -> y ~ x OK
#> findSpatial(y~x+anything(ga|bu)) -> NULL
#> findSpatial(y~x+corrMatrix(ga|bu))
# [[1]]
# corrMatrix(ga | bu)
#> findbarsMM(y~x+anything(ga|bu))
# [[1]]
# ga | bu