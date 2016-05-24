Life.lines <- 
function( entry.date = NA,
           exit.date = NA,
          birth.date = NA,
           entry.age = NA,
            exit.age = NA,
           risk.time = NA
         )
{
# A function allowing any three of the arguments to be specified
# and yet returns enty age and -time and exit age and -time.

# Check if any variable is supplied with class
if( conv <- any( inherits( entry.date, "Date"     ),
                 inherits(  exit.date, "Date"     ),
                 inherits( birth.date, "Date"     ),
                 inherits( entry.age , "difftime" ),
                 inherits(  exit.age , "difftime" ),
                 inherits(  risk.time, "difftime" ) ) )
  {
  # Convert "Date" and "difftime" to years
  if( inherits( entry.date, "Date"     ) )  entry.date <- as.numeric( entry.date ) / 365.35 + 1970
  if( inherits(  exit.date, "Date"     ) )   exit.date <- as.numeric(  exit.date ) / 365.35 + 1970
  if( inherits( birth.date, "Date"     ) )  birth.date <- as.numeric( birth.date ) / 365.35 + 1970
  if( inherits( entry.age , "difftime" ) )  entry.age  <- as.numeric( entry.age  ) / 365.35
  if( inherits(  exit.age , "difftime" ) )   exit.age  <- as.numeric(  exit.age  ) / 365.35
  if( inherits(  risk.time, "difftime" ) )   risk.time <- as.numeric(  risk.time ) / 365.35
  # Convert to numeric
  class( entry.date ) <- "numeric"
  class(  exit.date ) <- "numeric"
  class( birth.date ) <- "numeric"
  class( entry.age  ) <- "numeric"
  class(  exit.age  ) <- "numeric"
  class(  risk.time ) <- "numeric"
  }
  
# Find out which three items are supplied.
#
wh <- (1:6)[!is.na( list( entry.date,
  		          entry.age,
  		           exit.date,
                           exit.age,
  		   	  birth.date,
                           risk.time ) )]

# Matrix of relevant quantities.
#
LL <- rbind( entry.date,
             entry.age,
              exit.date,
              exit.age,
             birth.date,
              risk.time )

# Matrix giving the three constraints among the six quantities:
#
M <- rbind( c( -1, 1,  0,  0, 1, 0 ),
            c(  0, 0, -1,  1, 1, 0 ),
            c(  0, 1,  0, -1, 0, 1 ) )

# Now in principle we have that M %*% LL = 0.
# Partitioning M=(A1|A2), t(LL)=(t(x1),t(x2)) 
# this gives A1 %*% x1 = -A2 %*% x2

# Check if there is sufficient information  
#
if( qr( M[,-wh[1:3]] )$rank < 3 ) 
   cat( "Insufficient information to display life lines" )

# Then do the calculation
#
A1 <- M[, wh[1:3]]
A2 <- M[,-wh[1:3]]
x1 <- LL[wh[1:3],]

x2 <- -solve( A2 ) %*% A1 %*% x1
LL[-wh[1:3],] <- x2
LL <- data.frame( t(LL) )
attr( LL, "Date" ) <- conv

# Convert to dates and difftimes
if( conv )
  {
  LL[,c(1,3,5)] <- ( LL[,c(1,3,5)] - 1970 ) * 365.25
  LL[,c(2,4,6)] <-   LL[,c(2,4,6)]          * 365.25
  class( LL[,1] ) <-
  class( LL[,3] ) <-
  class( LL[,5] ) <- "Date"
  class( LL[,2] ) <-
  class( LL[,4] ) <-
  class( LL[,6] ) <- "difftime"  
  }

LL
}
