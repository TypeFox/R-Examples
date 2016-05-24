ktsolve <- function(yfunc, known=list(), guess, tool=c('BB', 'nleqslv'), show=TRUE, ...){
#revised 23 Sept 2013 to fix search/replace strings in gsub
# length(body(y)) is one greater than number of lines due to"{}"
if( !(is(yfunc,'function')) ) stop('yfunc type must be "function" ')
# if lengths are equal, we should be OK.
nullfoo <- switch( sign( length(unlist(guess)) - (length(body(yfunc)) - 3)) +2  ,
	warning('Fewer guesses than equations.  System may be underdefined.'),
	NULL,
	warning('More guesses than equations.  System may be overdefined.')
	)
# verify no duplicated names in known vs. guess
if( length( intersect(names(known), names(guess)) ) ) {
	stop(paste('The name(s) ',paste(intersect(names(known), names(guess)),collapse=' '),' are in both "known" and "guess" ') )
	}
# pick the solver - used switch() so easy to add more options
 tool<-tool[1] #get rid of the extras!
switch(tool,
	'BB' = 	require(BB, quietly=TRUE, warn.conflicts=FALSE ),
	'nleqslv' = require(nleqslv, quietly=TRUE, warn.conflicts=FALSE ),
	stop('Unknown solver package specified.') 
	)
# 	# replace variable name with value of variable,i.e. "known" in function body, if any
if (length(known)>0 ) {
	for (i in 1:length(known)) {
		if(length(grep(names(known)[i], body(yfunc)[-1])) < 1) warning("Input '", names(known)[i], '" not found in function body')
		lookfor<-paste("\\b",names(known)[i],"\\b",sep="",collapse="")
		parse(text=gsub(lookfor,known[i],body(yfunc)[-1])) -> body(yfunc)[-1]

		}
	}
# Replace each "guess" name with x[n] in yfunc's body 
for (i in 1:length(guess)) {
	if(length(grep(names(guess)[i], body(yfunc)[-1]))<1) warning("Guess '", names(guess)[i],"' not found in function body")
	subpat <- paste('x[', i, ']', sep='')
	lookfive<-paste("\\b",names(guess)[i],"\\b",sep="",collapse="")
	parse(text=gsub(lookfive, subpat, body(yfunc)[-1])) -> body(yfunc)[-1]
	}
# call solver tool with the values in 'guess' 
toolset <- cbind(c('BB','nleqslv'),c('BBsolve','nleqslv'))
do.call(toolset[toolset[,1]==tool,2], list(unlist(guess), yfunc, ...) ) -> solution
if(show) {
	cat('solution is:\n')
	switch(tool,
		'BB' = print(solution$par),
		'nleqslv' = print(solution$x)
		)
	if(length(known) > 0 ) {
		cat('"known" inputs were:\n')
		print(rbind(known))
		}
	} # end of if(show)

outs<-list(results=solution, tool=tool, yfunc=yfunc)
return(invisible(outs))
}
