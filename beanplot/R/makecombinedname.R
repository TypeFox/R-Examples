`makecombinedname` <-
function(string1,string2) {
    if (is.na(string2))
        return(string1)
    if (string1==string2)
        return(string1)
    s1 <- unlist(strsplit(string1,NULL))
    s2 <- unlist(strsplit(string2,NULL))
    fromleft <- 0
    gfromleft <- 0
    while ((fromleft<length(s1)) && (fromleft<length(s2)) && (s1[fromleft+1]==s2[fromleft+1])) {
        if (any(grep("[. \t]",s1[fromleft+1])))
            gfromleft<-fromleft
        fromleft<-fromleft+1
    }
    fromright <- 0
    gfromright <- 0
    while ((fromright<length(s1)) && (fromright<length(s2)) && (s1[length(s1)-fromright]==s2[length(s2)-fromright])) {
         if (any(grep("[. \t]",s1[length(s1)-fromright])))
              gfromright<-fromright
         fromright<-fromright+1
    }
    if (gfromleft>gfromright)
        result<-substr(string1,1,gfromleft)
    else if (gfromright==0)
        result<-paste(string1,string2,sep="+")
    else
        result<-substr(string1,nchar(string1)-gfromright+1,nchar(string1))
    result
}