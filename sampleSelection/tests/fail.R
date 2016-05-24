library( "sampleSelection" )

# tobit 2

try( selection( s ~ z1, y ~ x1 ) )

try( selection( "s ~ z1", y ~ x1 ) )

try( selection( ~ z1, y ~ x1 ) )

try( selection( s ~ z1, "y ~ x1" ) )

try( selection( s ~ z1, ~ x1 ) )

# tobit 5

try( selection( s ~ z1, list( y1 ~ x1, y2 ~ x1 ) ) )

try( selection( "s ~ z1", list( y1 ~ x1, y2 ~ x1 ) ) )

try( selection( ~ z1, list( y1 ~ x1, y2 ~ x1 ) ) )

try( selection( s ~ z1, list( "y1 ~ x1", y2 ~ x1 ) ) )

try( selection( s ~ z1, list( ~ x1, y2 ~ x1 ) ) )

try( selection( s ~ z1, list( y1 ~ x1, "y2 ~ x1" ) ) )

try( selection( s ~ z1, list( y1 ~ x1, ~ x1 ) ) )

