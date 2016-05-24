`item.stats` <-
function(x.i, n.x, steps = 1, C=1, correct.extreme=.03, model){
 n.i <- ncol(x.i)
 responded <- ! is.na(x.i)

 items <- list()
 n.ni <- apply(x.i,2,function(XXX) sum((! is.na(XXX))*C*n.x))
  
 S.ih <- apply(x.i,2,function(XXX) sapply(1:steps,function(YYY) sum( (XXX >= YYY)*C*n.x, na.rm=TRUE ) ) )
 if(model == "PCM") S.ih[S.ih == 0] <- NA
 T.ih <- apply(x.i,2,function(XXX) sapply(1:steps,function(YYY) sum( (XXX == YYY)*C*n.x, na.rm=TRUE ) ) )
 T.ih[is.na(S.ih)] <- NA 

 if(steps == 1){ Tx <- sum(S.ih)
 } else { Tx   <- rowSums(S.ih)  }

 Si   <- colSums( x.i*C*n.x, na.rm=T )

 if(steps==1) S.ih <- matrix(S.ih,nrow=1,dimnames=list(1,names(S.ih))) 

 list(n.i = n.i, n.ni = n.ni, S.ih = S.ih, T.ih = T.ih, steps = steps, Si = Si, Tx = Tx)
} #end item.stats

