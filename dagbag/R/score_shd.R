###A family of (generalized) SHD (structure hamming distance) to aggregate the DAGs learned on bootstrap resamples 
score_shd <- function(boot.adj, alpha = 1, threshold = 0, max.step = 500, blacklist = NULL, whitelist = NULL, print = FALSE){

 p=dim(boot.adj)[1]     ## number of variables
 nb=dim(boot.adj)[3]    ## number of DAGs in the emsemble 

 if(is.null(blacklist)){
   blacklist=matrix(0,p,p)
 }

 if(is.null(whitelist)){
   whitelist=matrix(0,p,p)
 }

 final.step=0
 movement=matrix(0,max.step,3) 
 adj.matrix=matrix(0,p,p)
 c.size=0
 freq.cut=(1-threshold)/2
 

  junk<-.C("score_gshd", 
           as.integer( p), 
           as.integer(boot.adj),
           as.double(alpha),
           as.integer(nb),
           as.integer(max.step),
           as.integer(print), 
           as.integer(blacklist),
           as.integer(whitelist),
           as.double(freq.cut),
           movement=as.integer(movement), 
           adj.matrix=as.integer(adj.matrix), 
           final.step=as.integer(final.step),
           c.size=as.integer(c.size)
         )

####outputing
 temp<-new.env()

 m.temp=matrix(junk$adj.matrix,p,p)
 temp$adj.matrix=m.temp
 temp$final.step=junk$final.step
 temp$movement=matrix(junk$movement,max.step,3)
 #temp$c.size=junk$c.size

 res=as.list(temp)
 return(res)
}