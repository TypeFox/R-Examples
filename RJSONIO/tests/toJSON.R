library(RJSONIO)

toJSON(list(1, 2, list(1, NA)))

toJSON(list(1, 2, list(NA)))  # collapses the sub-list into the main vector.


e = new.env(); e$a = 1:10; e$bc = letters[1:3]
cat(toJSON(e, pretty = TRUE))


a = list(x=1, y=character(0))
fromJSON( toJSON( a ) )
a = list(x=1, y=character(0), b = 1)
fromJSON( toJSON( a ) )


a = list(x = vector(), y = 123, z = "allo")
fromJSON(toJSON(a))


