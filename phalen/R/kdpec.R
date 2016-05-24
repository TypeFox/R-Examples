kdpec <-
function(id,kdim,startdate,enddate,slack = 0,restartindex = FALSE) {
  
  if (class(startdate)!="Date" || class(enddate)!="Date") {
    stop("start date or end date is not a date")
  }
  
  startdate = as.numeric(startdate)
  enddate = as.numeric(enddate) +slack
  
  if (is.null(ncol(id)) == TRUE) {
    id1 = seq.int(1,length(id))
  } else {
    id1 = seq.int(1,nrow(id))
  }
  
  
  if (is.null(ncol(kdim)) == TRUE) {
    kdim1 = kdim
  } else {
    kdim1 = kdim[,1]
    i = 2
    k = ncol(kdim)
    if (k>1) {
      for(i in 2:ncol(kdim)) {
        kdim1 = paste(kdim1,kdim[,i],sep="~")
      }
    }
  }
  
  kdim1 = match(kdim1,unique(kdim1))
  f = data.frame(id1,kdim1,startdate,enddate)
  # reorder to keep things legit
  f = f[order(f$kdim1,f$startdate,f$enddate,f$id1),]
  row.names(f) = seq.int(1,nrow(f))
  f = cbind(f,"episode" = seq.int(1,nrow(f)),"feed" = seq.int(1,nrow(f)))
  
  f = sqldf("
            SELECT   finl.id1
                    ,finl.kdim1
                    ,finl.startdate
                    ,finl.enddate
                    ,clst.episode
                    ,finl.feed
            FROM    f finl
                    INNER JOIN (SELECT    rec.id1
                                          ,min(minrec.episode) as episode
                                FROM      f minrec
                                INNER JOIN f rec
                                        ON minrec.kdim1 = rec.kdim1
                                        AND minrec.feed <= rec.feed
                                        AND minrec.enddate >= rec.startdate
                                GROUP BY rec.id1) clst
                             ON finl.id1 = clst.id1
            ORDER BY finl.feed
            ")
  
  
  repeat{
    f.it = sqldf("SELECT  f1.id1
                          ,f2.episode
                 FROM     f f1
                          INNER JOIN f f2
                                 ON f1.kdim1 = f2.kdim1
                                 AND f1.episode = f2.feed
                                 AND f1.episode <> f2.episode
                 ")
    
    if (nrow(f.it)==0) {break}
    f$episode[match(f.it$id1,f$id1)] = f.it$episode
  }
  
  if (restartindex==TRUE) {
    f.rn = unique(data.frame("kdim1" = f$kdim1,"episode" = f$episode))
    f.dm = unique(f.rn$kdim1)
    
    for(i in 1:length(f.dm)) {
      f.mt = f.rn[f.rn$kdim1 == f.dm[i],]
      f.mt = cbind(f.mt,"newcluster" = seq.int(nrow(f.mt)))
      
      f.wh = data.frame(which( outer(f$kdim1, f.mt$kdim1, "==") & 
                                 outer(f$episode, f.mt$episode, "=="), 
                               arr.ind=TRUE))
      
      f$episode[f.wh$row] = f.wh$col
    }
    
  } else {
    f$episode = match(f$episode,unique(f$episode))
  }
  
  f = f[order(f$id1),]
  kdpec = data.frame(id,"kdimidx" = f$kdim1,"episode" = f$episode[match(id1,f$id1)])
}
