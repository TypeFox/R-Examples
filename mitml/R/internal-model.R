# prepare model input by formula
.model.byFormula <- function(data, formula, group, group.original,
                      method=c("pan","jomo","jomo.matrix")){

  method <- match.arg(method)

  ft <- terms(formula)
  tl <- attr(ft,"term.labels")
  vrs <- attr(ft,"variables")[-1]

  # responses
  yvrs <- as.character(vrs)[attr(ft,"response")]
  yvrs <- gsub("[\r\n]","",yvrs)
  y.fml <- as.formula(paste0("~",yvrs))
  yvrs <- attr(terms(y.fml), "term.labels")
  # check for untransformed yvrs
  err <- !(yvrs %in% colnames(data))
  if(any(err)) stop("Could not find: ", paste0(yvrs[err],collapse=", "), "). Target variables must be contained in the data set 'as is', and transformations must be applied beforehand.")
  # order of formula terms in data: needed for indicator matrix
  # yord <- sapply(yvrs, function(x) which(colnames(data)==x))

  # cluster id
  clt <- tl[grep("\\|",tl)]
  if(length(clt)==0) stop("Cluster indicator not found in formula\n\n",formula,"\n\nPlease specify the cluster indicator and at least one random term using the '|' operator.")
  clt <- strsplit( clt, split="[[:blank:]]*\\|[[:blank:]]*" )[[1]]
  clname <- clt[2]

  # order data
  data <- data[ order(group,data[,clname]), ]
  group.original <- group.original[ order(group) ]
  group <- group[ order(group) ]

  # predictors: fixed
  pvrs <- c(if(attr(ft,"intercept")){"(Intercept)"}, tl[-grep("\\|",tl)])
  fe.fml <- c(if(attr(ft,"intercept")){"1"}else{"0"}, tl[-grep("\\|",tl)])
  fe.fml <- as.formula(paste0("~", paste0(fe.fml,collapse="+")))
  # predictors: random
  cl.fml <- as.formula(paste("~",clt[1])) 
  cl.ft <- terms(cl.fml)
  qvrs <- c(if(attr(cl.ft,"intercept")){"(Intercept)"}, attr(cl.ft,"term.labels"))

  # model matrix for fe and cl
  attr(data,"na.action") <- identity
  mmp <- suppressWarnings( model.matrix(fe.fml, data=data) )
  mmq <- suppressWarnings( model.matrix(cl.fml, data=data) )
  pnames <- colnames(mmp)
  qnames <- colnames(mmq)
  psave <- setdiff( c(pnames,qnames), c("(Intercept)",colnames(data)) )
 
  switch( method ,
    pan={ # panImpute (matrix input)
      y <- data.matrix(data[yvrs])
      ycat <- NULL
    },
    jomo={ # jomoImpute, newer versions (data input)
      y <- data[yvrs]
      ycat <- NULL
    },
    jomo.matrix={ # jomoImpute, older versions (matrix input)
      y <- data.matrix(data[yvrs])
      cvrs <- sapply(data[,yvrs,drop=F], is.factor)
      ycat <- y[,cvrs,drop=F]
      y <- y[,!cvrs,drop=F]
    }
  )

  clus <- data[,clname]
  pred <- cbind(mmp, mmq[,!(qnames%in%pnames),drop=F])
  xcol <- which(colnames(pred)%in%pnames)
  zcol <- which(colnames(pred)%in%qnames)

  inp <- list(
    y=y, ycat=ycat, clus=clus, pred=pred, xcol=xcol, zcol=zcol, data=data,
    group=group, group.original=group.original, psave=psave, clname=clname,
    yvrs=yvrs, pvrs=pvrs, qvrs=qvrs, pnames=pnames, qnames=qnames
  )
  for(i in names(inp)) assign(i, inp[[i]], pos=parent.frame())
  invisible(NULL)

}
 

# prepare model input by type
.model.byType <- function(data, type, group, group.original,
                   method=c("pan","jomo","jomo.matrix")){

  if(ncol(data)!=length(type)) stop("Length of 'type' must be equal to the number of colums in 'data'.")
  if(sum(type==-2)<1) stop("Cluster indicator not found.")
  if(sum(type==-2)>1) stop("Only one cluster indicator may be specified.")

  data <- data[ order(group,data[,type==-2]), ]
  group.original <- group.original[ order(group) ]
  group <- group[ order(group) ]

  clname <- colnames(data)[type==-2]
  clus <- data[,clname]
  yvrs <- colnames(data)[type==1]

  switch( method ,
    pan={ # panImpute (matrix input)
      y <- data.matrix(data[yvrs])
      ycat <- NULL
    },
    jomo={ # jomoImpute, newer versions (data input)
      y <- data[yvrs]
      ycat <- NULL
    },
    jomo.matrix={ # jomoImpute, older versions (matrix input)
      y <- data.matrix(data[yvrs])
      cvrs <- sapply(data[,yvrs,drop=F], is.factor)
      ycat <- y[,cvrs,drop=F]
      y <- y[,!cvrs,drop=F]
    }
  )

  pred <- cbind(1,as.matrix(data[type%in%c(2,3)]))
  pnames <- colnames(data)[type%in%c(2,3)]
  pvrs <- c("(Intercept)",pnames)
  qnames <- colnames(data)[type==3]
  qvrs <- c("(Intercept)",qnames)
  colnames(pred) <- pvrs
  
  xcol <- 1:length(pvrs)
  zcol <- xcol[pvrs%in%qvrs]

  inp <- list(
    y=y, ycat=ycat, clus=clus, pred=pred, xcol=xcol, zcol=zcol, data=data,
    group=group, group.original=group.original, clname=clname, yvrs=yvrs,
    pvrs=pvrs, qvrs=qvrs, pnames=pnames, qnames=qnames
  )
  for(i in names(inp)) assign(i, inp[[i]], pos=parent.frame())
  invisible(NULL)

}

