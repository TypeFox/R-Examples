#-#####################################################################################
#
# PREPROCESSING FUNCTIONS
#
#-#####################################################################################

######################################################################################
#tdmPreNAroughfix
#
#' Replace <NA> values with suitable non <NA> values
#'
#' This function replaces <NA> values in a list entry or data frame column with the median 
#' (for numeric columns) or the most frequent mode (for factor columns). 
#' It does the same as \code{na.roughfix} in package randomForest, but does so faster.
#' 
#' @param  object      list or data frame
#' @param  ...         additional arguments
#'
#' @return \code{object}, the list or data frame with <NA> values replaced
#'     
#' @export
######################################################################################
#--- code borrowed from https://stat.ethz.ch/pipermail/r-help/2010-July/244390.html  ---------------
tdmPreNAroughfix <- function (object, ...) {
  res <- lapply(object, roughfix)
  structure(res, class = "data.frame", row.names = seq_len(nrow(object)))
}

roughfix <- function(x) {
  missing <- is.na(x)
  if (!any(missing)) return(x)

  if (is.numeric(x)) {
    x[missing] <- median.default(x[!missing])
  } else if (is.factor(x)) {
    freq <- table(x)
    x[missing] <- names(freq)[which.max(freq)]
  } else {
    stop("tdmPreNAroughfix only works for numeric or factor")
  }
  x
}
#---------------------------------------------------------------------------------------------------

######################################################################################
#
#' Find constant columns. 
#'
#' Find all those columns in data frame \code{dset} which are completely constant or
#' completely NA and return a vector with their names.
#' @param dset  data frame
#' @return name vector of constant columns
#' @export
#
tdmPreFindConstVar <- function(dset)
{
    const.variables = NULL
    for (n in names(dset)) {
      if (is.factor(dset[,n])) {
        if (length(levels(dset[,n]))==1)  const.variables = c(const.variables,n);
      }
      if (is.numeric(dset[,n])) {
        if (max(dset[,n],na.rm=T)==min(dset[,n],na.rm=T)) const.variables = c(const.variables,n);
      }
    }
    return(const.variables);
}

######################################################################################
#tdmPreGroupLevels
#
#' Group the levels of factor variable in \code{dset[,colname]}.
#'
#' This function reduces the number of levels for factor variables with too many levels.
#' It counts the cases in each level and orders them decreasingly. It binds the least
#' frequent levels together in a new level "OTHER" such that the remaining untouched
#' levels have more than opts$PRE.Xpgroup percent of all cases. OR it binds the levels with 
#' least cases together in "OTHER" such that the total number of new levels
#' is opts$PRE.MaxLevel. From these two choices for "OTHER" take the one which binds more 
#' variables in column "OTHER".
#' 
#' @param  dset      data frame
#' @param  colname   name of column to be re-grouped
#' @param  opts      list, here we need \itemize{
#'      \item  PRE.Xpgroup   [0.99]
#'      \item  PRE.MaxLevel  [32]  (32 is the maximum number of levels allowed for \code{\link[randomForest]{randomForest}})
#'      }
#' @return \code{dset}, a data frame with \code{dset[,colname]} re-grouped
#'     
#' @export
######################################################################################
tdmPreGroupLevels <- function(dset,colname,opts)
{
    if (is.null(opts$PRE.Xpgroup)) opts$PRE.Xpgroup=0.99;
    if (is.null(opts$PRE.MaxLevel)) opts$PRE.MaxLevel=32;

    thisCol =  dset[,colname];
    z <- split(thisCol,thisCol);       # be aware that empty values ("") are dropped
    bz <- sapply(z,length);
    bz <- bz[order(bz,decreasing=T)];
    perc <- cumsum(bz)/sum(bz);        # cumulative percentage
    w1 <- which(perc>opts$PRE.Xpgroup);
    L <- length(bz); w2 <- NULL;
    if(L>opts$PRE.MaxLevel) w2 <- (opts$PRE.MaxLevel:L);
    if(length(w1)>length(w2)) {
      othernames <- names(perc)[w1];
    } else {
      othernames <- names(perc)[w2];
    }
    othernames <- setdiff(othernames,"");

    # some 'factor-arithmetic' to give newCol exactly the right level names:
    newCol <- factor(thisCol,levels=c(levels(factor(thisCol)),"OTHER"));
    newCol[which(thisCol %in% othernames)] <- "OTHER";
    newCol <- factor(newCol, levels=setdiff(levels(newCol),othernames));
    dset[,colname] <- newCol;

    dset;
}

######################################################################################
# tdmPreLevel2Target
#
#' Relate levels of a column with a target (column). 
#'
#' Print for each level of factor variable f which ratio 0 / 1 of the binary target
#' variable it contains and how many cases are in each level
#'
#' @param  dset    data frame
#' @param  target  name of target column
#' @param  f       number of column with factor variable
#' @param  opts    list, here we need \itemize{
#'    \item opts$thresh_pR
#'    \item opts$verbose
#' }
#'
#' @note SIDE EFFECTS:
#'     some printed output
#' @export
#
tdmPreLevel2Target <- function(dset,target,f,opts)
{
  if (is.null(opts$thresh_pR))
    opts$thresh_pR=0 #23.08      # 0.0: print every level, >0: print only levels with %reorderer > thresh_pR

  z = split(as.numeric(dset[,target])-1,dset[,f]);        # why "-1"? - because as.numeric converts factors 0/1 to numbers 1/2
  bz = as.data.frame(sapply(z,length));
  bz = cbind(bz,sapply(z,mean)*100);                      # sapply(..,mean) gives the fraction of "1"'s (%reorderer)
  colnames(bz) <- c("cases","%reorderer");
  ind <- which(bz[,"%reorderer"]>=opts$thresh_pR);
  if (length(ind)>0) {
    cat1(opts,"column",names(dset)[f],":\n")
    print(bz[ind,]);
  }
}

######################################################################################
# tdmPrePCA.train
#
#'     PCA (Principal Component Analysis) for numeric columns in a data frame. 
#'
#'     tdmPrePCA.train is capable of linear PCA, based on prcomp (which uses SVD), and 
#'     of kernel PCA (either KPCA, KHA or KFA).
#'
#' @param   dset    the data frame with training (and test) data. 
#' @param   opts    a list from which we need here the following entries: \itemize{
#'     \item    PRE.PCA:   ["linear" | "kernel" | "none" ]
#'     \item    PRE.knum:  if >0 and if PRE.PCA="kernel", take only a subset of PRE.knum records from dset 
#'     \item    PRE.PCA.REPLACE:  [T] =T: replace the original numerical columns with the PCA columns; =F: add the PCA columns
#'     \item    PRE.PCA.npc:   if >0, then add for the first PRE.PCA.npc PCs the monomials of
#'                   degree 2 (see tdmPreAddMonomials)
#'     \item    PRE.PCA.numericV   vector with all column names in dset for which PCA is performed.
#'                   These columns may contain *numeric* values only.
#'     }
#' @return    \code{pca},     a list with entries: 
#'     \item{dset}{  the input data frame dset with columns numeric.variables replaced or extended (depending on \code{opts$PRE.PCA.REPLACE})
#'                    by the PCs with names PC1, PC2, ... (in case PRE.PCA=="linear")
#'                    or with names KP1, KP2, ... (in case PRE.PCA=="kernel")
#'                    and optional with monomial columns added, if PRE.PCA.npc>0. 
#'                    The number of PCs is min(nrows(dset),length(numeric.variables)).  }
#'     \item{numeric.variables}{  the new numeric column names (PCs, monomials, and optionally old numericV, if \code{opts$PRE.PCA.REPLACE==F}) }
#'     \item{pcaList}{  a list with the items \code{sdev, rotation, center, scale, x} as returned from \code{\link{prcomp}} 
#'                    plus \code{eigval}, the eigenvalues for the PCs }
#' 
#' @note CAUTION: Kernel PCA (opts$PRE.PCA=="kernel") is currently disabled in code, it *crashes*
#'       for large number of records or large number of columns.
#'
#' @seealso  \code{\link{tdmPrePCA.apply}}
#' @author Wolfgang Konen, FHK, Mar'2011 - Jan'2012
#' @export
#
tdmPrePCA.train <- function(dset,opts)  {
    if (is.null(opts$PRE.PCA.numericV))
      stop("Please define opts$PRE.PCA.numericV, the vector of numeric column names, in order to run PCA");

    if (nrow(dset) < length(opts$PRE.PCA.numericV))
      stop("There are fewer rows in dset than numeric variables. Please increase the rows in dset.
            Consider opts$PRE.allNonVali=TRUE (use all non-validation data for PCA training)");      

    numeric.variables = opts$PRE.PCA.numericV;
    ptm <- proc.time()
    fname <-  opts$filename;
    other.vars <- setdiff(names(dset),numeric.variables);
    
    if (opts$PRE.PCA=="linear") {
      cat1(opts,fname,": linear PCA with",length(numeric.variables),"numeric variables on",nrow(dset),"records ...\n");
      x <- dset[,numeric.variables];

      pcaList <- prcomp(x);
      # pcaList$rotation[,i]: ith eigenvector, i.e. its components from original space (x)
      # pcaList$sdev[i]^2: the corresponding ith eigenvalue
      pcaList$eigval = pcaList$sdev^2;
      rx = as.data.frame( as.matrix(dset[,numeric.variables]) %*% pcaList$rotation);         
      # rx contains in its rows the PC-rotated x-vectors for each case: The 1st component
      # in each row is the projection of this case on PC1, the 2nd component in each row is
      # the projection on PC2 and so on. Thus column rx[,1] contains the PC1-values of all cases,
      # column rx[,2] the PC2-values and so on.
      # The names of rx are PC1, PC2, ... (this is inherited from colnames(pcaList$rotation))

      dset <- tdmAdjustDSet(dset,numeric.variables,rx,"PC",opts$PRE.PCA.REPLACE);
    }
    
    if (opts$PRE.PCA=="kernel") {
      stop("kernel PCA currently disabled, please check later versions");
#      cat1(opts,fname,": Kernel PCA (KHA) with",length(numeric.variables),"numeric variables...\n");
#      x <- dset[,numeric.variables];
#      
#      require(kernlab);
#      if (opts$PRE.knum>0) x <- x[1:opts$PRE.knum, ];
#      # *** big mem/other problems with kpca: R crashes when PRE.knum = 3000; error code from kpca when PRE.knum >190;
#      # *** runs for PRE.knum<=150 on appAcid, but with (of course) terrible result (rgain=45%)
#      #kpc <- kpca(~.,data=x,kernel="rbfdot",kpar=list(sigma=0.2),features=40);
#      # *** use instead kfa, much faster than kpca and no crash
#      #kpc <- kfa(~.,data=x,kernel="rbfdot",kpar=list(sigma=0.2),features=0);
#      # *** another option is kha, takes longer than kfa, but perhaps more reliable numerical results
#      kpc <- kha(~.,data=x,kernel="rbfdot",kpar=list(sigma=0.2),features=20, maxiter=100);
#      # pcv(kpc)[,i]: the ith eigenvector, i.e. its components from original space (x)  (not for kfa)
#      # eig(kpc)[i]: the corresponding ith eigenvalue (not for kfa)
#      eigval <- NULL # eig(kpc);
#      rx <- as.data.frame(predict(kpc,dset[,numeric.variables]));     
#      # rx contains in its rows the PC-rotated x-vectors for each case: The 1st component
#      # in each row is the projection of this case on PC1, the 2nd component in each row is
#      # the projection on PC2 and so on. Thus column rx[,1] contains the PC1-values of all cases,
#      # column rx[,2] the PC2-values and so on.
#      names(rx) <- sub("V","KP",names(rx));                               # names are KP1, KP2, ...
#
#      dset <- tdmAdjustDSet(dset,numeric.variables,rx,"KP",opts$PRE.PCA.REPLACE);
    }
    
    if (opts$PRE.PCA=="none") {
      rx = dset[,numeric.variables];
      pcaList = NULL;
    }
  
    # adding nonlinear input combinations (monomials of degree 2 for the first opts$PRE.PCA.npc PCs)
    dset <- tdmPreAddMonomials(dset,rx,opts$PRE.PCA.npc,opts);
    numeric.variables = setdiff(names(dset),other.vars);  
    #  = union{ PC-names , monomial names , original numeric.variables (if opts$PRE.PCA.REPLACE=FALSE)}
                      
    cat1(opts,"Proc time for PCA: ",(proc.time()-ptm)[1],"\n");
    
    pca <- list(dset=dset
               ,numeric.variables=numeric.variables
               ,pcaList=pcaList
               );
}

######################################################################################
# tdmPrePCA.apply
#
#'     Apply PCA (Principal Component Analysis) to new data. 
#'
#'     The PCA rotation is taken from \code{pcaList}, a value returned from a prior call to  \code{\link{tdmPrePCA.train}}.
#'
#' @param   dset    the data frame with the new data
#' @param   pcaList a value returned from a prior call to  \code{\link{tdmPrePCA.train}} 
#' @param   opts    a list from which we need here the following entries: \itemize{
#'     \item    PRE.knum:  if >0 and if PRE.PCA="kernel", take only a subset of PRE.knum records from dset 
#'     \item    PRE.PCA.npc:   if >0, then add for the first PRE.PCA.npc PCs the monomials of
#'                   degree 2 (see tdmPreAddMonomials)
#'     \item    PRE.PCA.numericV   vector with all column names in dset for which PCA is performed.
#'                   These columns may contain *numeric* values only.
#'     }
#' @param   dtrain  [NULL] optional, only needed in case that dset is a 0-row-data frame: then we 'borrow' the columns from dtrain,
#'                  the data set returned from \code{\link{tdmPrePCA.train}} in \code{pca$dset}.
#' @return    \code{pca},     a list with entries: 
#'     \item{dset}{  the input data frame dset with columns numeric.variables replaced
#'                    by the PCs with names PC1, PC2, ... (in case PRE.PCA=="linear")
#'                    or with names KP1, KP2, ... (in case PRE.PCA=="kernel")
#'                    and optional with monomial columns added, if PRE.PCA.npc>0  }
#'     \item{numeric.variables}{  the new column names for PCs and for the monomials }
#'     
#'
#' @seealso  \code{\link{tdmPrePCA.train}}
#' @author Wolfgang Konen, FHK, Mar'2011 - Jan'2012
#' @export
######################################################################################
tdmPrePCA.apply <- function(dset,pcaList,opts,dtrain=NULL)  {
    numeric.variables = opts$PRE.PCA.numericV;
    fname <-  opts$filename;
    other.vars <- setdiff(names(dset),numeric.variables);
    
    if (nrow(dset)>0) {
      if (opts$PRE.PCA=="linear") {
        rx = as.data.frame( as.matrix(dset[,numeric.variables]) %*% pcaList$rotation);         
        # The names of rx are PC1, PC2, ... (this is inherited from colnames(pcaList$rotation))
  
        dset <- tdmAdjustDSet(dset,numeric.variables,rx,"PC",opts$PRE.PCA.REPLACE);
      }
      
      if (opts$PRE.PCA=="kernel") {
        stop("kernel PCA currently disabled, please check later versions");
      }
      
      if (opts$PRE.PCA=="none") {
        rx = dset[,numeric.variables];
      }
    
      # adding nonlinear input combinations (monomials of degree 2 for the first opts$PRE.PCA.npc PCs)
      dset <- tdmPreAddMonomials(dset,rx,opts$PRE.PCA.npc,opts);
      numeric.variables = setdiff(names(dset),other.vars);   
      #  = union{ PC-names , monomial names , original numeric.variables (if opts$PRE.PCA.REPLACE=FALSE)}
    } 
    else { # i.e. nrow(dset)==0
      if (is.null(dtrain))  stop("Can only fill a 0-row-data-frame if dtrain is present!");        
			dset<-dtrain[1,];
			dset<-dset[-1,];         # now dset is a 0-row-data-frame with the same columns as dtrain
      numeric.variables = setdiff(names(dtrain),other.vars);   
      #  = union{ PC-names , monomial names , original numeric.variables (if opts$PRE.PCA.REPLACE=FALSE)}
    }
                      
    pca <- list(dset=dset, numeric.variables=numeric.variables);
}

######################################################################################
# tdmPreSFA.train
#
#'     SFA (Slow Feature Analysis) for numeric columns in a data frame. 
#'
#'     tdmPreSFA.train uses package \code{\link[rSFA]{rSFA}}. It is assumed that classification for the variable 
#'     contained in column \code{response.var} is done. SFA seeks features in an expanded function space for which  
#'     the intra-class variation w.r.t.  \code{response.var} is as low as possible.
#'
#' @param   dset    the data frame with training (and test) data. 
#' @param   response.var  the response variable for classification. 
#' @param   opts    a list from which we need here the following entries: \itemize{
#'     \item    PRE.SFA:   [ "linear" | "2nd" | "none" ]  which stands for [ 1st | 2nd degree monomial SFA | no SFA ]
#'     \item    PRE.SFA.REPLACE:  [T] =T: replace the original numerical columns with the SFA columns; =F: add the SFA columns
#'     \item    PRE.SFA.npc:   if >0, then add for the first PRE.SFA.npc PCs the monomials of
#'                   degree 2 (see tdmPreAddMonomials)
#' 		 \item    PRE.SFA.PPRANGE:  [11] number of inputs after preprocessing, they enter into expansion 
#' 		 \item    PRE.SFA.ODIM:    [5] number of SFA output dimensions (slowest signals) to return 
#'     \item    PRE.SFA.numericV   vector with all column names in dset which are input for SFA.
#'                   These columns may contain *numeric* values only.
#'     }
#' @return    \code{sfa},     a list with entries: 
#'     \item{dset}{   the input data frame dset with columns numeric.variables replaced or extended (depending on \code{opts$PRE.SFA.REPLACE})
#'                    by the SFA components with names SF1, SF2, ... 
#'                    and with optional monomial columns added, if PRE.SFA.npc>0  }
#'     \item{numeric.variables}{    the new numeric column names of \code{dset}, i.e. SFA components, monomials (and optionally  
#'                    PRE.SFA.numericV, if \code{opts$PRE.SFA.REPLACE==F}) }
#'     \item{sfaList}{  a list with the items \code{opts (sfaOpts)}, matrices DSF and SF and many others, as returned from 
#'                    \code{\link[rSFA]{sfaStep}}  }
#'
#' @seealso  \code{\link{tdmPreSFA.apply}}
#' @author Wolfgang Konen, Martin Zaefferer, FHK, Jan'2012 - Feb'2012
#' @export
#
tdmPreSFA.train <- function(dset,response.var,opts)  {
    #require(rSFA);
    
    if (is.null(opts$PRE.SFA.numericV))
      stop("Please define opts$PRE.SFA.numericV, the vector of numeric column names, in order to run SFA");
      
    numeric.variables = opts$PRE.SFA.numericV;
    ptm <- proc.time()
    fname <-  opts$filename;
    other.vars <- setdiff(names(dset),numeric.variables);
    
  	sfaOpts<-list(epsC=0)		
  	#idim=ceiling(ncol(dset[,numeric.variables])*opts$PRE.SFA.IDIM);
  	#ppRange=max(c(idim*sfaOpts$PRE.SFA.PPRANGE,2))#sfaOpts$PRE.SFA.PPRANGE #TODO should be between 2:idim.. but not always easy
  	idim <- length(numeric.variables);
  	ppRange <- min(opts$PRE.SFA.PPRANGE,idim);        # clip ppRange if larger than idim
  	odim <- opts$PRE.SFA.ODIM;
      	
    if (opts$PRE.SFA=="linear") {
  		stop("linear expansion functions are not yet implemented, please check later versions");
    }    
    else if (opts$PRE.SFA=="2nd") {		
      cat1(opts,fname,": SFA 2nd degree with",length(numeric.variables),"numeric variables and PPRANGE =",ppRange,"on ",nrow(dset)," records ...\n");
      x <- dset[,numeric.variables];
    	realclass=dset[,response.var];
  		sfaList = rSFA::sfa2Create(ppRange, rSFA::xpDim(ppRange), "PCA2", "ORD1", 0, sfaOpts)
    	sfaList$deg=2;
    	sfaList$classes=levels(dset[,response.var]);
  	  sfaList$nclass=length(sfaList$classes);
  		if (opts$PRE.SFA.doPB) {
    		#  perform bootstrap, if desired and necessary
    		sfaList$doPB = opts$PRE.SFA.doPB; #only needed in PB
    		sfaList$ppRange = ppRange;
    		tmpResult = opts$PRE.SFA.fctPB(realclass,as.matrix(x),sfaList);
    		x=tmpResult$x
    		realclass=tmpResult$realclass
  		}
  		# perform the preprocessing step
  		sfaList=rSFA::sfaStep(sfaList, x, "preprocessing");
  		# perform the expansion step 
  		# IMPORTANT: add patterns one class at a time (!)
  		for (nC in 1:sfaList$nclass){
  		  idx = which(realclass==sfaList$classes[nC]);
  		  if (length(idx)<2)               # this check should be later in rSFA also
  		    stop(sprintf("There must be at least 2 training records for each class, but class %s has only %3d record(s)", sfaList$classes[nC],length(idx)));
  		  #cat1(opts,"--- starting expansion for class",sfaList$classes[nC], "with",length(idx),"train records ...\n");  
  			sfaList=rSFA::sfaStep(sfaList, x[idx,], "expansion","CLASSIF");
  		}
  		sfaList=rSFA::sfaStep(sfaList, x, "sfa","SVDSFA");     # close the algorithm
  		y = rSFA::sfaExecute(sfaList, x);	       # y: matrix with same numbers of rows as x and ODIM columns with column names "V1","V2",...
      sx <- as.data.frame(y[,1:odim]);
      sfaPrefix <- "SF";
      names(sx) <- sub("V",sfaPrefix,names(sx));  
      
      dset <- tdmAdjustDSet(dset,numeric.variables,sx,sfaPrefix,opts$PRE.SFA.REPLACE);

    }    
  	else if (opts$PRE.SFA=="custom" ) {
  		stop("custom expansion functions are not yet implemented, please check later versions");
	    #TODO: Implement an option for custom expansion functions!!! xpDimfun sfaExpandFun
  	}
  	
    if (opts$PRE.SFA=="none") {
      sx = dset[,numeric.variables];
      sfaList = NULL;
    }
  
    # adding nonlinear input combinations (monomials of degree 2 for the first opts$PRE.SFA.npc SFA-components)
    dset <- tdmPreAddMonomials(dset,sx,opts$PRE.SFA.npc,opts);
    numeric.variables = setdiff(names(dset),other.vars);   
    #  = union{ SFA-component-names , monomial names , original numeric.variables (if opts$PRE.SFA.REPLACE=FALSE)}
                      
    cat1(opts,"Proc time for SFA: ",(proc.time()-ptm)[1],"\n");
    
    sfa <- list(dset=dset
               ,numeric.variables=numeric.variables
               ,sfaList=sfaList
               );
}

######################################################################################
# tdmPreSFA.apply
#
#' Apply SFA (Slow Feature Analysis) to new data. 
#'
#'     The SFA projection is taken from \code{sfaList}, a value returned from a prior call to  \code{\link{tdmPreSFA.train}}.
#'
#' @param   dset    the data frame with the new data
#' @param   sfaList a value returned from a prior call to  \code{\link{tdmPreSFA.train}} 
#' @param   opts    a list from which we need here the following entries: \itemize{
#'     \item    PRE.SFA.REPLACE:  [T] =T: replace the original numerical columns with the SFA columns; =F: add the SFA columns
#'     \item    PRE.SFA.npc:   if >0, then add for the first PRE.SFA.npc PCs the monomials of
#'                   degree 2 (see tdmPreAddMonomials)
#' 		 \item    PRE.SFA.ODIM:    [5] number of SFA output dimensions (slowest signals) to return 
#'     \item    PRE.SFA.numericV   vector with all column names in dset for which SFA is performed.
#'                   These columns may contain *numeric* values only.
#'     }
#' @param   dtrain  [NULL] optional, only needed in case that dset is a 0-row-data frame: then we 'borrow' the columns from dtrain,
#'                  the data set returned from \code{\link{tdmPreSFA.train}} in \code{sfa$dset}.
#' @return    \code{sfa},     a list with entries: 
#'     \item{dset}{  the input data frame dset with columns numeric.variables replaced
#'                    by the PCs with names PC1, PC2, ... (in case PRE.SFA=="linear")
#'                    or with names KP1, KP2, ... (in case PRE.SFA=="kernel")
#'                    and optional with monomial columns added, if PRE.SFA.npc>0  }
#'     \item{numeric.variables}{  the new column names for PCs and for the monomials }
#'
#' @seealso  \code{\link{tdmPreSFA.train}}
#' @author Wolfgang Konen, Martin Zaefferer, FHK, Jan'2012 - Feb'2012
#' @export
#
tdmPreSFA.apply <- function(dset,sfaList,opts,dtrain=NULL)  {
    numeric.variables = opts$PRE.SFA.numericV;
    fname <-  opts$filename;
    other.vars <- setdiff(names(dset),numeric.variables);
    
    if (nrow(dset)>0) {
      if (opts$PRE.SFA=="2nd") {
        x <- dset[,numeric.variables];
    		y <- rSFA::sfaExecute(sfaList, x);	          # column names are "V1", "V2", ...
    		#odim=ceiling(ncol(y)*opts$PRE.SFA.ODIM)		
      	odim <- opts$PRE.SFA.ODIM;
        sx <- as.data.frame(y[,1:odim]);
        sfaPrefix <- "SF";
        names(sx) <- sub("V",sfaPrefix,names(sx));  
        # The names of sx are now SF1, SF2, ... 
  
        dset <- tdmAdjustDSet(dset,numeric.variables,sx,sfaPrefix,opts$PRE.SFA.REPLACE);
      }     
      else if (opts$PRE.SFA=="linear") {
        stop("linear expansion functions are not yet implemented, please check later versions");
      }      
      else if (opts$PRE.SFA=="none") {
        sx = dset[,numeric.variables];
      }
    
      # adding nonlinear input combinations (monomials of degree 2 for the first opts$PRE.SFA.npc PCs)
      dset <- tdmPreAddMonomials(dset,sx,opts$PRE.SFA.npc,opts);
      numeric.variables = setdiff(names(dset),other.vars);   
      #  = union{ SFA-component-names , monomial names , original numeric.variables (if opts$PRE.SFA.REPLACE=FALSE)}
    } 
    else { # i.e. nrow(dset)==0
      if (is.null(dtrain))  stop("Can only fill a 0-row-data-frame if dtrain is present!");        
			dset<-dtrain[1,];
			dset<-dset[-1,];         # now dset is a 0-row-data-frame with the same columns as dtrain
      numeric.variables = setdiff(names(dtrain),other.vars);   
      #  = union{ SFA-component-names , monomial names , original numeric.variables (if opts$PRE.SFA.REPLACE=FALSE)}
    }
                      
    sfa <- list(dset=dset, numeric.variables=numeric.variables);
}

#-#####################################################################################
# tdmAdjustDSet:
#     helper function for tdmPrePCA.* and tdmPreSFA.*
#     a) adjust the names of the columns in rx in case they are conflicting with names in dset
#     b) replace the numeric variables in dset with columns from rx or add them, depending on PRE.REPLACE
tdmAdjustDSet <- function(dset,numeric.variables,rx,prefix,PRE.REPLACE) {
      if (length(intersect(names(dset),names(rx)))>0) {
        # try "PC_1" instead of "PC1" (if prefix=="PC")                         
        names(rx) <- sub(prefix,paste(prefix,"_",sep=""),names(rx));  
      }
      if (length(intersect(names(dset),names(rx)))>0) {
        cat("Conflicting names are:", intersect(names(dset),names(rx)), "\n");
        stop("dset contains some columns with the same name as names(rx)");
      }  
      other.vars <- setdiff(names(dset),numeric.variables);
    	if (PRE.REPLACE){
        # replace columns numeric.variables with the columns from rx (PCA or SFA columns)
    		dset <- as.data.frame(dset[,other.vars]);    # these statements are needed to get the right column names 
    		names(dset) <- other.vars;                   # also in the case where length(other.vars)==1  (!!)
    		dset <- cbind(dset, rx);
    	} else {
        # don't replace columns from dset, simply add columns from rx (PCA or SFA columns)
        dset <- cbind(dset, rx);
      }
}

######################################################################################
# tdmPreAddMonomials
#
#' Add monomials of degree 2 to a data frame.
#'
#' Given the data frame \code{dset} and a data frame \code{rx} with the same number of rows,
#' add monomials of degree 2 to dset for all quadratic combinations of the first PRE.npc
#' columns of \code{rx}. The naming of these new columns is "R1x2" for the combination of cols
#' 1 and 2 and so on (if prefix="R").
#'
#' @param   dset  the target data frame
#' @param   rx    a data frame where to draw the monomials from
#' @param   PRE.npc the number of columns from \code{rx} to use (clipped to \code{ncol(rx)} if necessary)
#' @param   opts  a list from which we need here the following entries: \itemize{
#'              \item   filename
#'              \item   VERBOSE
#'              }
#' @param   degree  [2] (currently only 2 is supported)
#' @param   prefix  ["R"] character prefix for the monomial column names
#'
#' @return  data frame \code{dset} with the new monomial columns appended. If 
#'          PRE.npc==0, the data frame is returned unchanged.
#'
#' @note CAVEAT: The double for-loop costs some time (e.g. 2-4 sec for ncol(rx)=8 or 10)
#' How to fix: make a version w/o for-loop and w/o frequent assigns to dset (**TODO**)
#' @export
#
tdmPreAddMonomials <- function(dset,rx,PRE.npc,opts,degree=2,prefix="R") {
    if (PRE.npc>0) {
      if (degree!=2) stop("other degrees than 2 not yet supported");
      npc = min(PRE.npc,ncol(rx));
      for (i in 1:npc) {
        for (j in i:npc) {
          dset <- cbind(dset, rx[,i]*rx[,j]/mean(rx[,i])/mean(rx[,j]));
          names(dset)[ncol(dset)] <- paste(prefix,i,"x",j,sep="");
        }
      }
      cat1(opts,opts$filename,": added", npc*(npc+1)/2,"monomials of degree 2 for the first",npc,"PCs ...\n")
    }
    dset;
}

#---- obsolete, now integrated into tdmPreAddMonomials --------------------------------
## tdmAddMonomials:
##     helper function for tdmPrePCA.*  and tdmPreSFA.*
##     adds nonlinear input combinations (monomials of degree 2 for the first PRE.npc PCs)
#tdmAddMonomials <- function(dset,rx,PRE.npc,opts) {
#    old.names <- names(dset);
#    if (PRE.npc>0) {
#      npc = min(PRE.npc,ncol(rx));
#      dset <- tdmPreAddMonomials(dset,rx[,1:npc],degree=2);
#      Nmono <- length(setdiff(names(dset),old.names));
#      cat1(opts,opts$filename,": added", Nmono,"monomials of degree 2 for the first",npc,"PCs ...\n")
#    }
#    dset;
#}

### obsolete
#    makeX <- function(dset,numeric.variables,fname,strg,opts) {
#      cat1(opts,fname,strg,length(numeric.variables),"numeric variables...\n")
#      if(opts$READ.TST) x <- dset[dset[,opts$TST.COL]==0,numeric.variables]
#      else x <- dset[,numeric.variables];
#      x;
#    }
