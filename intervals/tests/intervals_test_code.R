library( "intervals" )

######## NAs and empty intervals

u <- Intervals_full( as.numeric(NA), type = "Z" )
u <- c( u, u )
v <- Intervals_full( c(1,3,1,Inf), type = "Z" )
x <- Intervals( 1, closed = FALSE, type = "Z" ) # empty
w <- c( x, u, v )
rownames(w) <- letters[ 1:nrow(w) ]

x
u
v
w

is.na(w)
empty(w)

distance_to_nearest( u, v )
distance_to_nearest( w, v )

pts <- c( -Inf, Inf, NA, NaN, 0:2 )

distance_to_nearest( pts, v )

identical(
          distance_to_nearest( pts, v ),
          distance_to_nearest( pts, open_intervals( v ) )
          )

interval_overlap( w, v )

reduce( w )

open_intervals(w)




######## Subset assignment

u <- Intervals( 1:8, type = "Z" )
rownames( u ) <- letters[ 1:nrow(u) ]
v <- Intervals( -4:-1, type = "Z" )
rownames( v ) <- letters[ nrow(u) + 1:nrow(v) ]
w <- open_intervals(u)

u
v
w

# Basic
z <- u; z[2:3,] <- v
z

# With names
z <- u; z[c("a","d"),] <- v
z

# Missing row name
result <- try( { z <- u; z[c("a","e"),] <- v }, TRUE )
result

# Closure adjustment
z <- w; z[2:3,] <- v
z

# Size mismatch error
result <- try( { z <- w; z[2:3,] <- v[1,] }, TRUE )
result

# Intervals_full method
x <- Intervals_full( 1:6, c(TRUE,FALSE) )
rownames( x ) <- letters[ 1:nrow(x) ]
y <- Intervals_full( -4:-1 )

x
y

# Missing value names
z <- x; z[2,] <- y[1,]
z

# Missing x names
z <- y; z[1,] <- x[1,]
z

# Type mismatch error
result <- try( { z <- x; z[2:3,] <- v }, TRUE )
result

# Coercion up
type(v) <- "R"
closed(v) <- c( FALSE, TRUE )
x
v
z <- x; z[2:3,] <- v
z

# With warning due to assignment
z <- v; z[1,] <- x[3,]
z

# No row name copying with matrices
A <- matrix( 0, 2, 2 )
rownames(A) <- c("1","2")
z <- x; z[1:2,] <- A
z




######## distance_to_nearest() behavior

a <- Intervals_full( c(2,5), FALSE )
b <- Intervals_full( c(1,5,2,10), matrix(c(TRUE,FALSE,FALSE,TRUE),2,2) )

a
b

distance_to_nearest(a,b)

type(a) <- "Z"
type(b) <- "Z"

distance_to_nearest(a,b)

a <- as( a, "Intervals" )
b <- as( open_intervals( b ), "Intervals" )

a
b

distance_to_nearest(a,b)

type(a) <- "R"
type(b) <- "R"

distance_to_nearest(a,b)
