######################################################################################
# tdmParaBootstrap
#
#'       Parametric bootstrap: add 'noisy copies' to a data frame (training data).
#'
#'   A normal distribution is approximated from the data given in \code{dset[,input.variables]} and new  
#'   data are drawn from this distribution for the columns \code{input.variables}. The column \code{resp} 
#'   is filled at random with levels with the same relative frequency as in \code{dset[,resp]}. 
#'   Other columns of dset are filled by copying the entries from the first row of dset.
#'
#'   @param dset     data frame with training set
#'   @param resp     name of column in dset which carries the target variable
#'   @param input.variables     vector with names of input columns 
#'   @param opts        additional parameters [defaults in brackets]
#'     \describe{
#'     \item{\code{ncopies}}{ how many noisy copies to add }
#'     \item{\code{ncsigma}}{ [1.0] multiplicative factor for each std.dev. }
#'     \item{\code{ncmethod}}{ [3] which method to use for parametric bootstrap\cr
#'          =1: each 'old' record from X in turn is the centroid for a new pattern;\cr
#'          =2: the centroid is the average of all records from the same class, 
#'              the std.dev. is the same for all classes;\cr
#'          =3: centroid as in '2', the std.dev. is the std.dev. of
#'              all records from the same class  (*recommended*)
#'         }
#'     \item{\code{TST.COL}}{ (optional) name of column in dset where each PB record is marked with a 0}
#'     }
#'         
#'   @return  data frame \code{dset} with the new parametric bootstrap records added as last rows.
#'
#' @seealso   \code{\link{tdmClassify}}
#' @author Wolfgang Konen, FHK, Nov'2011-Dec'2011
#'
#' @export
######################################################################################
tdmParaBootstrap <- function(dset,resp,input.variables,opts) {
  if(is.null(opts$ncsigma)) opts$ncsigma=1.0;                           
  if(is.null(opts$ncmethod)) opts$ncmethod=3;                           

  new_dset=NULL;
  if (opts$ncopies>0){
      # put input variables in matrix form, because matrix calculation is faster:
      #
      x = as.matrix(dset[,input.variables]);
      r0 = dset[,resp];
      new_x = NULL;
      new_realclass = NULL;
      c0 = getCentroid(x,r0,opts$ncmethod);
      s0 = getSigma(x,r0,opts$ncmethod);

      sz1 = dim(x)[1];           # number of training records in x
      sz2 = dim(x)[2];
      nn1 = floor(opts$ncopies/sz1)+1;     # example: ncopies=220, sz1=100 --> nn1=3
      # We divide the generation of the ncopies bootstrap patterns
      # into nn1 loop steps. This division is only needed for opts.ncmethod=1
      # (each training record in turn serves as centroid). If
      # opts.ncmethod=2 or 3, we could as well generate the bootstrap patterns
      # in one step.
      for (k in 1:nn1){
          x2 = c0 + opts$ncsigma * s0 * matrix(rnorm(sz1*sz2),sz1,sz2);          # /WK/ runif -> rnorm
          # i.e., calculate for each column of x the standard deviation,
          # repeat it sz1 times, multiply each entry by a random number
          # with standard normal distribution, multiply by opts.ncsigma
          # --> as a result we have
          #       std(x2-c0) = opts.ncsigma * std(x)
          # (approximately, up to statistical fluctuations).
          # The tacit assumption is here that the standard deviations in
          # each column are the same for all classes. 
          
          if (k==nn1){
              # in the last pass, add only that many randomly selected
              # records from x2 that the total number of added records 
              # amounts exactly to opts.ncopies
              idx = sample(sz1);
              idx = idx[1:(opts$ncopies%%sz1)];
              x2 = x2[idx,];
              r0 = r0[idx];
          }
          new_x = rbind(new_x, x2);
          new_realclass = factor(c(new_realclass, r0),labels=levels(r0));
      }  #for (k)
      
      # put results from matrix form back in data-frame form:
      #
      rownam = 1:opts$ncopies +  max(as.numeric(row.names(dset)));
      new_dset=data.frame(dset[rep(1,opts$ncopies),],row.names=rownam);  # make a fresh data frame with same # rows as new_x,
      names(new_dset) <- names(dset);                                    # same columns as dset and same entries as first line of dset
      if(!is.null(opts$TST.COL))                                           
        if (any(names(dset)==opts$TST.COL)) new_dset[,opts$TST.COL]=0;   # mark PB records as training data
      #new_dset[,input.variables] <- data.frame(new_x);
      new_dset[,input.variables] <- new_x;
      new_dset[,resp] <- factor(new_realclass,levels=levels(dset[,resp]));
  } #if (opts$ncopies>0)

  dset <- rbind(dset,new_dset);
}



###################################################################################
# Helper functions for tdmParaBootstrap
###################################################################################

repmat2 <- function(a,n,m) {
  if (is.vector(a)) a=as.matrix(t(a));
  kronecker(matrix(1,n,m),a);
}

# Get centroid for adding noisy copies
#
# @param x0 			matrix with training data
# @param r0			  vector with true class of training data
# @param ncmethod		see opts$ncmethod:\cr
#						=1: c0 = x0 (each pattern is its own centroid). Faster, but not exactly 'parametric bootstrap'\cr
#           >1: c0 = centroid(m(x0)) where m(x0) is the class to which x0 belongs. We recommend this method since
#						    it adds new patterns which are statistically independent from x0 (true parametric bootstrap)
#
# @return matrix \code{c0} \cr
#  - same size as x0, each row is the centroid of one noisy copy
#
###################################################################################
getCentroid <- function(x0,r0,ncmethod) {
    if (ncmethod==1){c0=x0;} 
    else{
        classes=unique(r0);
        c0 = x0*0; 
        for(i in 1:length(classes)){
            ind = which(r0==classes[i]);
            c0[ind,] = repmat2(apply(x0[ind,],2,mean), length(ind), 1); # column mean
        }
    }
	return(c0)
}

###################################################################################
# Get sigma for adding noisy copies
#
# @param x0 			training data
# @param r0			true class of training data
# @param ncmethod		see opts$ncmethod:\cr
#						<3: s0 = std(x0) (each column has a common sigma, std over all classes) faster, but perhaps not precise modeling\cr
#           =3: s0 = std(m(x0)) where m(x0) is the class to which x0 belongs  (*recommended*)
#
# @return matrix \code{s0} \cr
#  - same size as x0, each row has the std.deviations of one noisy copy
#
###################################################################################
getSigma <- function(x0,r0,ncmethod) {
    sz1 = dim(x0)[1];           #% number of training records in x0
    if (ncmethod<3){ 
        s0=repmat2(apply(x0,2,sd),sz1,1);   # column std      
    }else{
        classes=unique(r0);
        s0 = x0*0;
        for(i in 1:length(classes)){
            ind = which(r0==classes[i]);
            s0[ind,] = repmat2(apply(x0[ind,],2,sd), length(ind), 1);  # column std
        }
    }
	return(s0)
}

