library(datamap)

set.seed(1)
newMapper(
	type="EXAMPLE",
	init=function(map,symbols=c('foo','bar','baz'),len=3){

		# Install symbols that the users passes in from newMap().
		lapply(symbols,install,map)

		# Now let's add some state to the internal portion of our map.
		map$len <- len

		# Returning FALSE means failure
		return(TRUE)
	},
	get = function(x) {
		cat("I'll just get",x,"for you.\n")

		# len is pulled from the internal portion of the map
		# by lexical scoping rules. Anything can be returned here, but we 
		# default to a numeric value
		rnorm(len)
	},
	assign = function(x,val){
		cat("Calling assign",val,"to",x,".\n")
	},
	finalize = function(map){
		cat("Finalization can clear any state, like shutting down database\n")
		cat("connections, socket connections, etc.\n")
	},

	# The rest of the arguments are copied to the internal portion of the map.
	foo = 'bar'
)
m <- newMap('EXAMPLE')
m
m[['bar']]
with(m,bar)
m$get('bar')
with(m,x <- 'buzzle')
m
mapAttach(m)
baz
rm(m)
foo
detach('datamap:EXAMPLE')
