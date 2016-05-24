### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### triples.r v0.2
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### 
### 2005-11-22:
###   * changed setTriple warning() to stop(), added more text
### 2005-11-08:
###   * tried to add with(environment(M), { ... })
###     to getTriple, setTriple and delTriple (does not work)
###   * added error handling to setTriple (produces warning
###     if no environment exists (by checking the class,
###     throws an error if input parameters are wrong)
### 
### 2005-08-26:
###   * added garbage collection to delTriple
###   * removed useless create environment from setTriple

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# convert input subject (column names or 
# column positions) to column position
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
getSubjectId <- function(M, subject) {
    if (is.character(subject)) {
        return( which(match(colnames(M),subject)>0) );
    } else if (is.numeric(subject) && !any(subject>ncol(M)) ) {
        return( subject );
    } else return( NULL );
}

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# getTriple: return the value (=object) of the 
# requested triple(s) from the environment variables 
# "triples$S/P/O" the given matrix M. Leave out
# predicate to get all triples of the specified 
# subject. Leave out subject, to get all triples.
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
getTriple <- function(M, subject, predicate) {
    
    if ( exists("triples$S",envir=environment(M)) ) {
        
        if ( ! missing(subject) ) {
            
            spos = which( get("triples$S", envir=environment(M)) == getSubjectId(M,subject) );
            
            if ( ! missing(predicate)) {
                ppos = which( get("triples$P", envir=environment(M))[spos] == tolower(predicate) );
                objects = as.vector(get("triples$O", envir=environment(M))[spos][ppos]);
            } else {
                objects = list( as.vector(get("triples$P", envir=environment(M))[spos]), as.vector(get("triples$O", envir=environment(M))[spos]) );
            }
            
        } else {
            if ( length(get("triples$S",envir=environment(M))) == 0) {
                return( NULL );
            } else {
                return ( list( as.vector(get("triples$S", envir=environment(M))), as.vector(get("triples$P", envir=environment(M))), as.vector(get("triples$O", envir=environment(M))) ) );
            }
        }
        
        if ( length(objects)==0 ) {
            return( NULL );
        } else return( objects );
            
    } else return( NULL ); # if no triples exist -> return NULL
    
}

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# setTriple: enter a new triple into the environment 
# variables "triples$S/P/O" of the given matrix M.
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
setTriple <- function(M, subject, predicate, object) {
    
    # be aware: the environment must already be
    # created outside this function! If you use
    # your own object (not created by textmatrix, 
    # please do so by:
    #    environment(M) = new.env();
    #    class(M) = "textmatrix";
    
    if ( ! inherits(M, "textmatrix") ) {
        stop("[setTriple] - You are using a matrix which has not been generated with the \nlsa package. Therefore, you have to manually add an environment to the \nmatrix you use and set the class to 'textmatrix' to be able to save triples.\nAlternatively, you can use as.textmatrix() to convert a matrix\nto a textmatrix (be aware that it temporarilty\nneeds twice the amount of memory of the input matrix).\n");
    }
    
    # if input is vectors, check if they 
    # have the same number of elements (else break)
    if (length(subject) != length(predicate) || length(predicate)!=length(object) || length(subject)!=length(object) ) {
        stop("[setTriple] - Input vectors are not of the same length!");
    }
    
    if ( ! exists("triples$S",envir=environment(M)) ) {
        
        # if not yet existing, add 'triples$S/P/O' to 
        # environment of M and insert first triple
        assign("triples$S", factor(getSubjectId(M,subject)), envir=environment(M));
        assign("triples$P", factor(tolower(predicate)), envir=environment(M));
        assign("triples$O", factor(object), envir=environment(M));
        
    } else {
        
        if ( !any( is.na( match(getTriple(M, subject, predicate), object)) == FALSE)  ) {
            # triple does not exist, so append
            striples = get("triples$S", envir=environment(M));
            levels(striples) = unique(c(levels(striples), getSubjectId(M,subject)));
            striples[(length(striples)+1):(length(striples)+length(subject))] = getSubjectId(M,subject);
            assign("triples$S", striples, envir=environment(M));
            
            ptriples = get("triples$P", envir=environment(M));
            levels(ptriples) = unique(c(levels(ptriples), tolower(predicate)));
            ptriples[(length(ptriples)+1):(length(ptriples)+length(predicate))] = tolower(predicate);
            assign("triples$P", ptriples, envir=environment(M));
            
            otriples = get("triples$O", envir=environment(M));
            levels(otriples) = unique(c(levels(otriples), object));
            otriples[(length(otriples)+1):(length(otriples)+length(object))] = object;
            assign("triples$O", otriples, envir=environment(M));
        }
        
    } # insert triple(s)
    
} # # setTriple

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# delTriple: remove specific triple(s) from
# environment of M. Currently not very memory sensitive ;)
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# changes:
# 
# 2006-07-31:
#    * removed gc() output to make delTriple more silent
# 
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

delTriple <- function(M, subject, predicate, object) {
    
    # find position
    spos = which( get("triples$S",envir=environment(M)) == getSubjectId(M,subject) )
    ppos = which( get("triples$P",envir=environment(M))[spos] == tolower(predicate) )
    opos = which( get("triples$O",envir=environment(M))[spos][ppos] == object )
    origppos = ppos[opos]
    origspos = spos[origppos]
    
    # retract
    assign("triples$S", get("triples$S",envir=environment(M))[-origspos], envir=environment(M))
    assign("triples$P", get("triples$P",envir=environment(M))[-origspos], envir=environment(M))
    assign("triples$O", get("triples$O", envir=environment(M))[-origspos], envir=environment(M))
    
    # garbage collection
    a = gc()
        
} # delTriple

