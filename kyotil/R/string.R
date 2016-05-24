# paste two strings together
# e.g. "a" %+% "b"

"%+%" <- function (a, b) {
    out=paste(a,b,sep="")
    out
}



#### file name manipulation

getFileStem=function(file.name){
    substr(file.name, 1, lastIndex(file.name, ".")-1 )
}
getStem=getFileStem
fileStem=getFileStem

getExt = function(file.name){
    substr(file.name, lastIndex(file.name, ".")+1, nchar(file.name) )
}

# s="prefix_name" is changed to "name"
remove.prefix=function(s,sep="_"){
    tmp=strsplit(s,sep)    
    sapply(tmp, function (x) {
        if (length(x)==1) return (x)
        concatList(x[-1],sep)
    })
}


remove.postfix=function(s,sep="_"){
    tmp=strsplit(s,sep)    
    sapply(tmp, function (x) {
        if (length(x)==1) return (x)
        concatList(x[-length(x)],sep)
    })
}



#### misc

escapeUnderline=function (name) {
    gsub("_", "\\_", name)
}



#### string search

firstIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
            break;
        }
    }
    ret
}

lastIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
        }
    }
    ret
}

# return TRUE if s1 starts with s2? 
startsWith=function(s1, s2){
    sapply (s1, function (s) {
        if ( substring (s, 1, nchar(s2)) == s2 ) {
            return (TRUE);
        } else {
            return (FALSE);
        }
    })
}

# return TRUE if s1 ends with s2, s1 can be a vector
endsWith=function(s1, s2){
    sapply (s1, function (s) {
        if ( substring (s, nchar(s)-nchar(s2)+1, nchar(s)) == s2 ) {
            return (TRUE);
        } else {
            return (FALSE);
        }
    })
}



# return TRUE if s1 contains s2
contain =function (s1, s2) {
    sapply (s1, function (s) {
        k=nchar (s2)
        matched=0
        for (i in 1:(nchar(s)-k+1) ) {
            if (substr(s, i, i+k-1)==s2) {
                matched=i
                break
            }
        }
        matched!=0
    })
}
