`cat.freq.post` <-
function(cat.ma,num.samples){
  num.cats <- ncol(cat.ma)-2
  if(num.cats<2) stop("Must be at least 2 unique categories")
  cat.names <- names(cat.ma)[3:ncol(cat.ma)]
  
  post.samples <- array(dim=c(num.samples,num.cats,nrow(cat.ma)))

  for(a in cat.names){
    cat.ma[,paste(a,".lpi",sep="")] <- NA
    cat.ma[,paste(a,".upi",sep="")] <- NA
  }

    for(i in 1:nrow(cat.ma)){

      if(cat.ma$n[i]>0){
        x <- unlist(cat.ma[i,2] * cat.ma[i,3:(num.cats+2)])

        #The posterior under uniform dirichlet prior
        samples <- rdirichlet(n=num.samples,a=x+1)

        post.samples[,,i] <- samples
        
        for(j in 1:num.cats){
          cat.ma[i,paste(cat.names[j],c(".lpi",".upi"),sep="")] <-
            quantile(samples[,j],probs=c(.025,.975))

        }
      }
    }
  
  
  list(cat.ma=cat.ma,post.samples=post.samples)
}

