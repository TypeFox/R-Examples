### Internal functions.

.dbg.internal<-0

mySeq <- function
### Sequence generation resembling matlab ':' operator.
### The difference from R `seq' is the behaviour in case when
### from=2, to=1 and b=1 - in that case R seq raises an error, while
### matlab returns an empty sequence.
### This function returns an empty sequence as well in this case.
##keyword<<internal
(from,##<< starting value
to,##<< ending value
by=1##<< step (defaults to 1)
){
    if ((from<to) == (by>0)) return(seq(from,to,by))
    else if (from==to) return(from)
    else return(vector(class(from),0))
}

printVar<-function
### Print the name, type and value of a variable.
##keyword<<internal
(x##<< variable to print
) {
    #cat(paste(deparse(substitute(x)),' [',class(x),', ',length(x),']: ',paste(x,collapse=' '),'\n',sep=''))
    cat(paste(deparse(substitute(x)),
        ' [',paste(class(x),collapse=','),', ',
        ifelse(!is.null(dim(x)),paste(dim(x),collapse=','),length(x)),']: ',sep=''))
    if (is.environment(x)) {
        cat(paste(paste(ls(x),collapse=','),'\n',sep=''))
    } else if (length(x)>0) {
        for (i in 1:length(x)) {
            if (is.environment(x[i])) {
                tmp<-as.data.frame(x[i])
            } else {
                tmp<-x[i]
            }
            if (!is.null(names(x[i]))) {
                cat(paste(names(x[i]),'=',tmp,ifelse(i<length(x),', ','\n'),sep=''))
            } else {
                cat(paste(tmp,ifelse(i<length(x),', ','\n'),sep=''))
            }
        }
    } else {
        cat('\n')
    }
    #print(x)
}

.findVar<-function(idfName) {
    found<-FALSE
    value<-NA
    #print(sys.nframe())
    for (i in sys.nframe():1) {
        #print(i)
        if (exists(idfName,envir=sys.frame(i),inherits=FALSE)) {
            value<-get(idfName,envir=sys.frame(i),inherits=FALSE)
            found<-TRUE
            break
        }
    }
    return(list(value=value,found=found))
}

.gfc<-function
### GFC: "get from caller" function makes given variable that appears
### in the environment of some caller of the function calling the
### `.gfc' function usable directly by the caller of the `.gfc'
### function.
###
### Example:
###  f1() { a<-1; f2()}
###  f2() { a<-.gfc(a); <<now a copy of `a' appears in f2>> }
(nm,##<< identifier (not a character name) of a variable to get
required=TRUE##<< Require this variable and err if not found? If FALSE,
## NULL return value indicates either the variable exists and has the
## NULL value assigned, or the variable has not been found.
) {
    if (.dbg.internal) cat(paste('gcf: looking for \'',deparse(substitute(nm)),'\'\n',sep=''))
    rv<-.findVar(deparse(substitute(nm)))
    if (rv$found) {
        return(rv$value)
    } else {
        if (required) {
            print(traceback())
            stop(paste('\'',deparse(substitute(nm)),'\' not found in caller stack',sep=''))
        }
        # obsoleted fallback:
        #return(attr(rv$value,'.sharedEnv')[deparse(substitute(nm))])
        return(NULL)
    }
}
