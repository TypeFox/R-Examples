splitData3<-function(dat,nobj,ENV){

    dat<-as.matrix(dat[[1]])
    # calculates NAgroup membership
    if (any(is.na(dat))) {
      dichX <- ifelse(is.na(dat),1,0)
      strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
      gmemb <- as.vector(data.matrix(data.frame(strdata)))
    } else {
      gmemb <- rep(1,nrow(dat))
    }

    # list of data submatrices according to NA pattern
    listX<-split(data.frame(dat),gmemb) # list with data, splitted into NA groups
    # generates aggregate information for each NA block in actual cov group
    listA<-lapply(listX,aggreg3,ENV)          # list with counts, logical not NA pattern, CL vector s

    listA  ###das RESULTAT
}
