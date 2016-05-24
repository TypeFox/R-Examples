# Copy the functionality of  unix  bash commands 'pushd' and 'popd'
# may 2015 -- this uses a localized environment to store the stack,
# thus meeting CRAN requirements. [i.e. not to alter GlobalEnv]
# Using approach suggested by MrFlick on SO  http://stackoverflow.com/questions/30416292

# TODO
#	1) add a little  "getstack()" function to let user know what's in the stack
#
#

.pushpop <- (function() {
 #  I need to initialize the existence  of my stack, but keep it empty.
 	.dirhist <- vector()
# pushpop is set by the calling pushd or popd function
    f <- function(path, dn=FALSE, rot=0, pull=0, pushpop)  {
 # go to case construction
	switch(pushpop+1,
	{
			#if (pushpop==0) # internal func, can't be anything else but0	
		# do popd
		if( !(exists('.dirhist')) || length(.dirhist)==0 ) {
#		if( !(exists('.dirhist')) ) {
			print('Nothing to pop.')
			return(invisible(1)) # an error value
		} 
		if(!dn & !pull ) {
		setwd(.dirhist[1])
		.dirhist<<- .dirhist[-1]
		}
#  
# In bash-land, looks like any "pull" value suppresses changing directory,
# not quite what the documentation claims. Could change this by
# adding an else here, with !dn allowing a setwd() to happen
		if(pull) {
			# "fix" negative values
			if(pull<0) pull <- length(.dirhist)+pull
			.dirhist<<-.dirhist[-pull]
			}
		return(invisible(0)) #status value
		}  ,
	 {
	 	# switch input was 1+1
			# do pushd
			newhist<-getwd() #the full path
			if(missing(path) && !dn &&!rot ) {
			# just swap the top two directories
				if(exists('.dirhist') && length(.dirhist)>1) .dirhist[1:2]<-.dirhist[2:1]
				return(invisible(0))
				}
			# dn is equivalent to the first, bold "-n" arg to pushd in bash
			if (!dn )  {
				setwd(path) }
	# bash adds directory to stack when dn is true or false
	# but NOT when path is missing...need to check for nonexistence of .dirhist
	# when trying to rotate
	#  only add newhist if intending to change dir, i.e. not rot ;
	# but can't have both a path and rot at the same time (but could have neither)
	# so added that if(!rot) to the else
				if(exists('.dirhist') && !missing(path) ) {
					ddirhist <- c(newhist,.dirhist)
					} else {
		# looks odd, but if no path and no rot, basically we exited two "if"s ago 
						if(!rot) ddirhist <- newhist	
						}		
			# looks like bash does change dir when rot!=0 so long as dn is false
			# control how "far" to rotate stack
			if (rot && exists('.dirhist') && length(.dirhist)>1 ){
				ddirhist <- .dirhist  # since it wasn't created yet
				rot <- rot%%length(ddirhist)
# safety check to make sure there was no path entered.
# later, put an input validation in pushd() so can't have rot and path together.
				if(missing(path) ) {
					ddirhist[] <- ddirhist[c( (rot+1):length(ddirhist),1:(rot-1) )]
				} 
			}
	# it's possible to get here w/o ever creating ddirhist			
			if(exists('ddirhist')) .dirhist<<-ddirhist
			return(invisible(0))
			}  ,
	{
		#switch input was 2+1
				return(.dirhist)
		}   
		) #end of switch
    }  #end of "f" function
    return(environment())
})()


pushd <- function(path, dn=FALSE,rot=0) { 
	.pushpop$f(path=path,dn=dn, rot=rot, pushpop=1 )
	}

# arg "pull" is like +/-n in bash version, to remove items from stack
popd <- function(dn=FALSE, pull=0) { 
	.pushpop$f(dn=dn, pull=pull, pushpop=0)
	}

# little helper function to return the current stack
getstack <- function(){
	.pushpop$f(pushpop=2)
}
############ bash documentation (from cygwin)
#  pushd [-n] [+n] [-n]
  # Adds  a  directory to the top of the directory stack, or rotates
              # the stack, making the new top of the stack the  current  working
              # directory.  With no arguments, exchanges the top two directories
              # and returns 0, unless the directory stack is empty.   Arguments,
              # if supplied, have the following meanings:
              # -n     Suppresses  the  normal  change  of directory when adding
                     # directories to the stack,  so  that  only  the  stack  is
                     # manipulated.
              # +n     Rotates  the  stack  so  that the nth directory (counting
                     # from the left of the list shown by  dirs,  starting  with
                     # zero) is at the top.
              # -n     Rotates  the  stack  so  that the nth directory (counting
                     # from the right of the list shown by dirs,  starting  with
                     # zero) is at the top.
              # dir    Adds dir to the directory stack at the top, making it the
                     # new current working directory.
# popd [-n] [+n] [-n]
              # Removes entries from the directory stack.   With  no  arguments,
              # removes  the  top directory from the stack, and performs a cd to
              # the new top directory.  Arguments, if supplied, have the 
              #following meanings:
              # -n     Suppresses  the  normal change of directory when removing
                     # directories from the stack, so that  only  the  stack  is
                     # manipulated.
              # +n     Removes  the nth entry counting from the left of the list
                     # shown by dirs, starting with zero.  For  example:  ``popd
                     # +0'' removes the first directory, ``popd +1'' the second.
              # -n     Removes the nth entry counting from the right of the list
                     # shown by dirs, starting with zero.  For  example:  ``popd
                     # -0''  removes the last directory, ``popd -1'' the next to
                     # last.
