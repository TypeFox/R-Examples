
# TEST FILE FOR hash METHODS set

library(hash)
h <- hash()

# SET: key-value pairs
.set( h, "a", 1:2 )
.set( h, letters, 1:26  )
.set( h, 1:5, 1:5 ) 
.set( h, letters, 12  )

# SET: key-hash pair added in version 1.0.4
.set( h, "ha", hash( a=1, b=2 ) )
class( h[["ha"]] ) == "hash"


# SET: data.frame
.set( h, "df", data.frame( a=1:10, b=11:20) )
class( h[["df"]] ) == "data.frame" 

# SET: list
.set( h, "li", list( a=1, b=1:5, c=letters[1:3] ) )
class( h[["li"]] ) == "list" 

# SET: environment
