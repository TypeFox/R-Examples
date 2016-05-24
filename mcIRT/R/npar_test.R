
DDF <- function(reshOBJ, wm="focal")
{
  

# ------------------------------------------------------------------

  
faclev <- lapply(1:ncol(reshOBJ$d[[1]]), function(zi){
                    sort(unique(c(unique(reshOBJ$d[[1]][,zi]),unique(reshOBJ$d[[2]][,zi]))))
                                                  })

Ps <- mapply(function(daf,d01)
  {# groups
  corcat  <- grep("cor",colnames(d01))
  rands   <- rowSums(d01[,corcat])
  rands_f <- factor(rands,levels=0:length(corcat))
  
  groupP <- mapply(function(dd,fl)
      {#items
      mycol <- factor(daf[,dd],levels=fl)
      
      prp <- tapply(mycol,rands_f,function(ee)
              {
              prop.table(table(ee))
              })
      
      # delete NULL entries
      null_table <- table(mycol)
      null_table[] <- 0
      
      prp2 <- sapply(prp,function(inn)#
              {
              if(is.null(inn)) null_table
                  else inn
              })
            },dd=1:ncol(daf),fl=faclev,SIMPLIFY=FALSE)
        
  }, daf=reshOBJ$d, d01 = reshOBJ$recm,SIMPLIFY=FALSE)
  


wm_all <- mapply(function(daf,d01)
{# groups
  corcat  <- grep("cor",colnames(d01))
  rands   <- rowSums(d01[,corcat])
  rands_f <- factor(rands,levels=0:length(corcat))
  
  groupP <- mapply(function(dd,fl)
  {#items
    mycol <- factor(daf[,dd],levels=fl)
    
    prp <- tapply(mycol,rands_f,function(ee)
    {
      sum(table(ee))
    })
      
    prp2 <- lapply(prp,function(inn)
    {
      if(is.na(inn)) 0
      else inn
    })
    unlist(prp2)
  },dd=1:ncol(daf),fl=faclev,SIMPLIFY=TRUE)
  
}, daf=reshOBJ$d, d01 = reshOBJ$recm,SIMPLIFY=FALSE)




if(wm == "focal")
{
wmc <- wm_all[[2]]
} else if(wm == "reference")
  {
    wmc <- wm_all[[1]]  
  } else if(wm == "total")
      {
      wmc <- wm_all[[1]] + wm_all[[2]]
      }

##### standardized p-difference #################

items <- ncol(reshOBJ$d[[1]])

# gesamtabweichungsindex

stdpdif <- lapply(1:items,function(alli)
      {
      totsum <- sum(wmc[,alli])
      # w*(reference - focal)/sum(w)
      colSums((t(Ps[[1]][[alli]]-Ps[[2]][[alli]]) * unlist(wmc[,alli]))/totsum)
    })


#####   RMWSD   #################


rmwsd <- lapply(1:items,function(alli)
{
  totsum <- sum(wmc[,alli])
  # w*(reference - focal)/sum(w)
  sqrt(colSums((t(Ps[[1]][[alli]]-Ps[[2]][[alli]])^2 * wmc[,alli])/totsum))
})


##### VAR(P_fs - P_bs) #################


varpp <- mapply(function(ro,di)
          {
          ro^2 - di^2
          },ro=rmwsd,di=stdpdif,SIMPLIFY=FALSE)


# Beschriftung
names(stdpdif) <- colnames(reshOBJ$d[[1]])
names(rmwsd) <- colnames(reshOBJ$d[[1]])
names(varpp) <- colnames(reshOBJ$d[[1]])

ergddf <- list(stdpdif=stdpdif,rmwsd=rmwsd,varpp=varpp,Ps=Ps)
class(ergddf) <- c("DDF","standardization")

return(ergddf)

}








