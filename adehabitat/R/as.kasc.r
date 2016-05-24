"as.kasc" <- function(l)
{

### 1. Verification that all "asc" attributes are similar
    clobj<-unlist(lapply(l,class))
    if (!all(clobj=="asc")) stop("input should be a list of \"asc\" objects")
    u<-TRUE
    la<-list()
    for (i in 1:length(l)) la[[i]]<-attributes(l[[i]])
    o<-la[[1]]
    if (o$type=="factor") {
          o<-o[names(o)!="levels"]
        }
    o<-o[names(o)!="type"]
    o<-o[names(o)!="dimnames"]

### 2. storage of attributes, but we delete the variable type and an
### eventual "levels" attribute
    if (length(l)>1) {
      for (i in 2:length(l))
        {
          tmp<-la[[i]]
          if (tmp$type=="factor") {
            tmp<-tmp[names(tmp)!="levels"]
          }
          tmp<-tmp[names(tmp)!="type"]
          tmp<-tmp[names(tmp)!="dimnames"]

          u[i]<-all(sort(unlist(tmp))==sort(unlist(o)))
        }
      if (!all(u)) stop("all the objects should have the same attributes")
    }

### 3. Computation of the kasc
    u<-as.vector(l[[1]])
    if (attr(l[[1]], "type")=="factor") {
      ct<-levels(l[[1]])
      lab<-list()
      for (j in 1:length(ct)) {
        lab[[j]]<-ct[j]
      }
      lab<-unlist(lab)
      u<-factor(u, levels=1:length(lab), labels=lab)
    }
    output<-data.frame(u)

    if (length(l)>1) {
      for (i in 2:length(l)) {
        u<-as.vector(l[[i]])
        if (attr(l[[i]], "type")=="factor") {
          ct<-levels(l[[i]])
          lab<-list()
          for (j in 1:length(ct)) {
            lab[[j]]<-ct[j]
          }
          lab<-unlist(lab)
          u<-factor(u, levels=1:length(lab), labels=lab)
        }
        output<-cbind.data.frame(output, u)
      }
    }

### 5. The attributes
    attr(output, "cellsize")<-attr(l[[1]], "cellsize")
    attr(output, "xll")<-attr(l[[1]], "xll")
    attr(output, "yll")<-attr(l[[1]], "yll")
    attr(output, "ncol")<-nrow(l[[1]])
    attr(output, "nrow")<-ncol(l[[1]])
    attr(output, "type")<-unlist(lapply(l, function(x) attr(x, "type")))
    names(output)<-names(l)
    class(output)<-c("kasc","data.frame")
    return(output)
  }

