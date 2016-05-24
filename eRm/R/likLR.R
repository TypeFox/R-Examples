likLR <- function(X, W, mpoints, Groups, model, st.err, sum0, etaStart){

  if (any(is.na(X))) {
    dichX <- ifelse(is.na(X),1,0)
    strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
    gmemb <- as.vector(data.matrix(data.frame(strdata)))
  } else {
    gmemb <- rep(1,dim(X)[1])
  }

  #data preparation, design matrix generation for various models
  if (model=="RM") { Xprep <- datprep_RM(X,W,sum0)
  } else if (model=="LLTM") { Xprep <- datprep_LLTM(X,W,mpoints,Groups,sum0)
  } else if (model=="RSM") { Xprep <- datprep_RSM(X,W,sum0)
  } else if (model=="PCM") { Xprep <- datprep_PCM(X,W,sum0)
  } else if (model=="LRSM") { Xprep <- datprep_LRSM(X,W,mpoints,Groups,sum0)
  } else if (model=="LPCM")  {Xprep <- datprep_LPCM(X,W,mpoints,Groups,sum0)
  }

  if (any(is.na(etaStart))) etaStart <- rep(0,dim(Xprep$W)[2])       #check starting vector
  if (length(etaStart) != dim(Xprep$W)[2]) stop("Vector with starting values does not match number of parameters!") 
  ng <- max(Groups)
  if ((dim(Xprep$W)[1]) != ((dim(Xprep$X01)[2])*ng)) stop("Mismatch between number of rows (beta's) in W and number of items (categories) in X!")


  Lprep <- cmlprep(Xprep$X01,Xprep$mt_vek,mpoints,Groups,Xprep$W,gmemb)                   
  parest <- fitcml(Lprep$mt_ind,Lprep$nrlist,Lprep$x_mt,Lprep$rtot,Xprep$W,
                   max(Groups),gind=Lprep$gind,x_mtlist=Lprep$x_mtlist,
                   Lprep$NAstruc,g_NA=Lprep$g_NA,st.err,etaStart,gby=Lprep$gby)      

  W1 <- Xprep$W
  #rownames(W1) <- NULL
  #colnames(W1) <- paste("eta",1:dim(W1)[2],sep="")
                           
  return(list("W" = W1, "parest" = parest, "X01" = Xprep$X01))                          #returns design matrix and results

}
