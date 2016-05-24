################################
# G+C content
#################################

GC <- function(seq, forceToLower = TRUE, exact = FALSE, NA.GC = NA, oldGC = FALSE )
{
  #
  # NA propagation:
  #
  if(length(seq) == 1 && is.na(seq)) return(NA)
  #
  # Check that sequence is a vector of chars:
  #
  if(nchar(seq[1]) > 1) stop("sequence is not a vector of chars")
  #
  # Force to lower-case letters if requested:
  #
  if(forceToLower) seq <- tolower(seq)
  #
  # Compute the count of each base:
  #
  nc <- sum( seq == "c" )
  ng <- sum( seq == "g" )			
  na <- sum( seq == "a" )
  nt <- sum( seq == "t" )
  
  #
  # oldGC case:
  #
  if(oldGC){
    warning("argument oldGC is deprecated")
    return( (nc + ng)/length(seq) )
  }
  
  #
  # General case:
  #
  if(! exact){
     if(na + nc + ng + nt == 0){
       result <- NA.GC
     } else {
      result <- (nc + ng)/(na + nc + ng + nt)
    }
  } else {
		#
		# We have our first estimate of GC vs. AT base counts:
		#
		ngc <- ng + nc
		nat <- na + nt

		#
		# weak and strong bases are 100% informative with respect
		# to the GC content, we just add them:
		#
		# s : Strong (g or c)
		# w : Weak (a or t)
		#
		ngc <- ngc + sum( seq == "s" )
		nat <- nat + sum( seq == "w" )

		##########################
		# Ambiguous base section #
		##########################
		
		#
		# m : Amino (a or c)
		#
		if(na + nc != 0){
			nm <- sum( seq == "m")
			ngc <- ngc + nm*nc/(na + nc)
		nat <- nat + nm*na/(na + nc)
		}
               
		#
		# k : Keto (g or t)
		#
		if(ng + nt != 0){
			nk <- sum( seq == "k" )
			ngc <- ngc + nk*ng/(ng + nt)
			nat <- nat + nk*nt/(ng + nt)
		}
		
		#
		# r : Purine (a or g)
		#
		if(ng + na != 0){
			nr <- sum( seq == "r" )
			ngc <- ngc + nr*ng/(ng + na)
			nat <- nat + nr*na/(ng + na)
                      }      
               
		#
		# y : Pyrimidine (c or t)
		#
		if(nc + nt != 0){
			ny <- sum( seq == "y" )
			ngc <- ngc + ny*nc/(nc + nt)
			nat <- nat + ny*nt/(nc + nt)
		}
                
		#
		# v : not t (a, c or g)
		#
		if(na + nc + ng != 0){
			nv <- sum( seq == "v" )
			ngc <- ngc + nv*(nc + ng)/(na + nc + ng)
			nat <- nat + nv*na/(na + nc + ng)
		}
		#
		# h : not g (a, c or t)
		#
		if(na + nc + nt != 0){
			nh <- sum( seq == "h" )
			ngc <- ngc + nh*nc/(na + nc + nt)
			nat <- nat + nh*(na + nt)/(na + nc + nt)
		}
		#
		# d : not c (a, g or t)
		#
		if(na + ng + nt != 0){
			nd <- sum( seq == "d" )
			ngc <- ngc + nd*ng/(na + ng + nt)
			nat <- nat + nd*(na + nt)/(na + ng + nt)
		}
		#
		# b : not a (c, g or t)
		#
		if(nc + ng + nt != 0){
			nb <- sum( seq == "b" )
			ngc <- ngc + nb*(nc + ng)/(nc + ng + nt)
			nat <- nat + nb*nt/(nc + ng + nt)
		}
		#
		# n : any (a, c, g or t) is not informative, so
		# we compute the G+C content as:
		#
                if( ngc + nat == 0){
                  result <- NA.GC
                } else {
  		  result <- ngc/(ngc + nat)
                }
	}
  return(result)
}

######################
# GCpos		     #
######################

GCpos <- function(seq, pos, frame = 0, ...){
  if(nchar(seq[1]) > 1){
    warning("sequence is not a vector of chars, I'm trying to cast it into one")
    seq <- s2c(seq[1])
  }
  #
  # Take frame into account:
  #
  if(frame != 0) seq <- seq[(1 + frame):length(seq)]
  #
  # Return result:
  #
  GC(seq[seq(pos, length(seq), by = 3)], ...)
}

######################
# GC1		     #
######################

GC1 <- function(seq, frame = 0, ...) GCpos(seq = seq, pos = 1, frame = frame, ...)

######################
# GC2		     #
######################

GC2 <- function(seq, frame = 0, ...) GCpos(seq = seq, pos = 2, frame = frame, ...)

######################
# GC3		     #
######################

GC3 <- function(seq, frame = 0, ...) GCpos(seq = seq, pos = 3, frame = frame, ...)

