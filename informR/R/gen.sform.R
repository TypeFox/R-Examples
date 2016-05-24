gen.sform<-
function(a,sform,olev=NULL){
   if(!is.null(olev)){a.uid<-sort(unique(olev))} else{a.uid<-attr(a,"a.uid")[,1]}
   if(!is.null(attr(a,"a.null"))){a.null<-which(attr(a,"a.uid")[,2]%in%attr(a,"a.null"))}  
   arb<-FALSE
   if(grepl("[[:punct:]]|[[:digit:]]",sform)){arb<-TRUE}
   plach<-lapply(vector("list",length=length(a)),function(x) x<-matrix(0,1,length(a.uid),dimnames=list(1,a.uid)))
   ev.seq<-a
   #bmat<-matrix(0,nrow=length(ev.seq),ncol=length(unique(ev.seq)),dimnames=list(1:length(ev.seq),a.uid))
  
   #New Search Engine
   if(!arb){
      b1<-sfsplit(sform)
      if(!any(grepl(b1[[1]],a))){
         if(!is.null(attr(a,"a.null"))){plach<-lapply(plach, function(x) rbind(x[,-a.null]))}
         return(plach)
      }
      b1.nchar<-lapply(b1,nchar)
      b1.evnames<-names(b1)
      #tmp<-NULL

      lotus<-lapply(b1,function(x) regmat.ind(x,a)[,2])
      for(h in 1:length(b1)){
         if(!is.null(lotus[[h]])){
            for(i in 1:length(lotus[[h]])){
               if(length(plach)>=lotus[[h]][i]+1){
                  plach[[lotus[[h]][i]+1]][,b1.evnames[[h]]]<-1
               }
            }
         }  
      }
   if(!is.null(attr(a,"a.null"))){plach<-lapply(plach, function(x) rbind(x[,-a.null]))}
   return(plach)
   }
   if(arb){
   pl.true<-FALSE
   or.true<-FALSE
      #The sfdec function determines if there are any OR statements in the s-form regex and returns 
      #a vector of sequences if TRUE
   sfdec<-function(sform){
     if(grepl("\\(*\\|*\\)",sform)){
        lb.ind<-gregexpr("\\(",sform)
        or.ind<-gregexpr("\\|",sform)
        rb.ind<-gregexpr("\\)",sform)
        or.mat<-mapply(cbind,lb.ind,or.ind,rb.ind,SIMPLIFY=FALSE)[[1]]
        op.mat<-apply(or.mat,1,function(x) rbind(substr(sform,x[1]+1,x[2]-1),substr(sform,x[2]+1,x[3]-1)))
        uni.sforms<-vector("list",length=ncol(op.mat))
        for(i in 1:ncol(op.mat)){
           for(j in 1:nrow(op.mat)){
              uni.sforms[[i]][j]<-gsub(substr(sform,or.mat[i,1],or.mat[i,3]),op.mat[j,i],sform,fixed=TRUE)   
           }
        }
     return(unlist(uni.sforms))
     }
     return(sform)
   }
   arb.sforms<-sfdec(sform)
   val.sforms<-unlist(lapply(sapply(arb.sforms,strsplit,""),function(x) any(grepl(x[1],a))))
   if(!any(val.sforms)){
      if(!is.null(attr(a,"a.null"))){plach<-lapply(plach, function(x) rbind(x[,-a.null]))}
      return(plach)
      }
   arb.sforms<-unique(arb.sforms[which(val.sforms)])
   if(any(sapply(arb.sforms,grepl,pattern="\\+"))){
   pl.true<-TRUE
   pl.val<-which(sapply(arb.sforms,grepl,pattern="\\+"))
   }
   if(grepl("\\(*\\|*\\)",sform)){or.true<-TRUE}
   if(or.true){
      if(pl.true){orb.sforms<-arb.sforms[-pl.val]} else orb.sforms<-arb.sforms
      for(g in 1:length(orb.sforms)){
         b1<-sfsplit(orb.sforms[g])
         b1.nchar<-lapply(b1,nchar)
         b1.evnames<-names(b1)
         #tmp<-NULL

         lotus<-lapply(b1,function(x) regmat.ind(x,a)[,2])
         for(h in 1:length(b1)){
            if(!is.null(lotus[[h]])){
               for(i in 1:length(lotus[[h]])){
                  if(length(plach)>=lotus[[h]][i]+1){
                     plach[[lotus[[h]][i]+1]][,b1.evnames[[h]]]<-1
                  }
               }
            }  
         }
      }
   }
   if(pl.true){
      prb.sforms<-arb.sforms[pl.val]
      for(g in 1:length(prb.sforms)){
         b1<-strsplit(prb.sforms[[g]],"")[[1]]
         b1.pl<-which(b1%in%"+")
         b1.evnames<-b1
         b1.evnames[b1.pl]<-b1[b1.pl-1]
         #tmp<-NULL
         for(h in 1:(length(b1)-1)){
            if((h+1)%in%b1.pl){ind.up<-c(b1.evnames[h+1],b1.evnames[h+2])} else ind.up<-b1.evnames[h+1]
            for(i in 1:(length(a)-1)){
               if(!h%in%b1.pl){
                  if(h==1){if(a[i]==b1.evnames[h]){plach[[i+1]][,ind.up]<-1}}
                  if(a[i]==b1.evnames[h]&&plach[[i]][,b1.evnames[h]]==1){plach[[i+1]][,ind.up]<-1}
               }
            }
         }
      }
   }
   if(!is.null(attr(a,"a.null"))){plach<-lapply(plach, function(x) rbind(x[,-a.null]))}
   return(plach)
   }
}
