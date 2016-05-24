permgpu <-
function(datobj,y,event=NULL,test,B,diag=FALSE,scale=FALSE)
  {
    error=FALSE
    if(!is.element(test,c("ttest","wilcoxon","pearson","spearman","cox","npcox")))
       {
         print("This test has not been implemented yet")
         error=TRUE
       }
   
    if((test=="cox" | test=="npcox") &(!error))
      event=pData(datobj)[[event]]
    else
      event=0
    
    y=pData(datobj)[[y]]
    EXPR=exprs(datobj)
    gnames=featureNames(datobj)
    n=ncol(EXPR)
    K=nrow(EXPR)

    if (scale)
      EXPR = scale(t(EXPR), scale=FALSE)
    else
      EXPR = t(EXPR)
 
    if(!error)
      {
        output=rep(0,3*K)
        out<-.C("permgpu",
                as.single(EXPR),
                as.single(y),
                as.single(event),
                as.integer(n),
                as.integer(K),
                as.integer(B),
                as.character(test),
                as.single(output),# 3 columns: stat, pB,PB
                NAOK=TRUE,
                PACKAGE="permGPU"
                )
        names(out)=c("EXPR","y","event","n","K","B","test","output")
        res=data.frame(gnames, stat=out[[8]][((2*K)+1):(3*K)],unadjp=out[[8]][1:K], fwerp=out[[8]][(K+1):(2*K)])
        if(diag)
          res=list(RESULTS=res,out)
      }
    else
      {
        res=NULL
      }
    return(res)
  }

