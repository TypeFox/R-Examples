`cmlprep` <-
function(X01,mt_vek,mpoints,Groups,W,gmemb)
{

  levs <- (gmemb-1)*max(Groups)+Groups              #merge Groups and gmemb vector into level vector
  
  if (length(Groups)==1) {                          #if no group contrast
     x_mt <- colSums(X01,na.rm=TRUE)                #item category raw scores as vector
     #eventuell x_mtlist auf NA gruppen aufbrechen
     x_mtlist <- list(x_mt) 
     ngroups <- 1
    
  } else {                                            #if groups defined
    ngroups <- max(Groups)                            #number of groups
    x_mtlist <- by(X01,levs,colSums,na.rm=TRUE)       #item-category raw scores for each group (as list)
    x_mtlist.G <- by(X01,Groups,colSums,na.rm=TRUE)   #item-category raw scores for each group (as list)
    #FIXME!!! use x_mtlist??
    x_mt <- as.vector(unlist(x_mtlist.G))             #as vector: g1|g2|...
  }

  end1 <- length(mt_vek)*mpoints*ngroups
  mt_ind <- rep(1:end1,rep(mt_vek,mpoints*ngroups)) #category index vector (for converting x_mt into list)
  x_tmt <- split(x_mt,mt_ind)                       #list for likelihood: item-wise * ngroups
  rtot <- sum(mt_vek)*mpoints

  ics <-  rep(sequence(mt_vek),mpoints)                 #item category scores for each item as vector
  rv <- apply(X01,1,function(x) {                       #person raw scores of 0/1 matrix
                      ics[!is.na(x)]%*%na.exclude(x)}) 

  #--- preparing index vector for item parameters ---
  if (ngroups > 1) {                                    #groups
    seglen <- sum(mt_vek)                               #length of beta vector (segment)
    gind <- rep(rep(1:ngroups,rep(seglen,ngroups)),mpoints) #parameter index vector for group extraction
  } else {
    gind <- rep(1,dim(W)[1])
  }

  #--- preparing lists for person splits ---
  rvlist <- split(rv,levs)                    #split person raw scores due to levels (NAgroup AND treatment)
  nrlist <- lapply(rvlist,function(rvel) {    #list with item raw score frequencies for each group (transposed)
                            rvtab <- table(rvel)                            #raw score frequencies
                            dnamevek <- as.numeric(unlist(dimnames(rvtab))) #different raw scores for 0 fill up
                            nr <- rep (0,rtot+1)                            #setting 0 raw score frequencies
                            nr[dnamevek+1] <- rvtab #vector with person raw scores from 1:rtot (with 0 fill up)
                            nr <- nr[-1]
                            return(nr)
                          })
                 
  
  if ((ngroups > 1) && (length(unique(gmemb)) > 1)) {          #NA groups AND Groups
    gg <- table(Groups,gmemb)
    #gg[gg > 0] <- 1
    g_NA <- as.vector(rowSums(gg))                             #How many NA-sub groups in each Group
    #grgm <- cbind(Groups, gmemb)
    #grgmst <- apply(grgm,1,function(x) {                       #merge indexes to characters
    #            paste(x[1],x[2]) })
    #GGind <- rank(unique(grgmst))    
    #levtab <- table(levs)                                      #frequencies of levels
    #FIXME!!! following line wrong index
    #gby <- rep(GGind,levtab)                                   #ordering by NAgroups nested in Group

    #this probably does the job
    gby <- levs
  } else {
    g_NA <- 1
    gby <- gmemb
  }
  
  NAstruc <- by(!is.na(X01),gby,function(x) {                  #list of unique NA structures for each Group
                                    x.u <- unique(x)
                                    as.numeric(as.matrix(x.u))}) #NA's are coded with 0
                                    
  NAcheck <- sapply(NAstruc,sum)                         #if for certain NAgroups only 1 item was presented
                                    
list(x_mt=x_mt,mt_ind=mt_ind,x_tmt=x_tmt,rtot=rtot,nrlist=nrlist,gind=gind,x_mtlist=x_mtlist,
     NAstruc=NAstruc,g_NA=g_NA,gby=gby)
}

