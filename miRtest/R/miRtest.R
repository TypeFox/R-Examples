


#' Contingency table. Necessary for Fisher test.
#' @param gene.set Vector of gene sets.
#' @param p.val Vector with p-values.
#' @param sign Significane threshold.
#' @author Stephan Artmann
contingency.table = function (gene.set,p.val,sign=0.05) {
 sign.gs = length(which(p.val[gene.set] <= 0.05));
 sign.not.gs = length(which(p.val[-gene.set] <= 0.05));
 cont.tab = c(sign.gs,length(gene.set)-sign.gs,sign.not.gs,length(p.val)-length(gene.set)-sign.not.gs);
 dim(cont.tab) = c(2,2)
 cont.tab; 
}


#' Turn a data.frame indicating gene sets into the allocation matrix. 
#' @param df data.frame with mRNAs in its first and miRNAs in its second column.
#' @param X Expression matrix of miRNAs whose row names will be used to generate the list of miRNAs.
#' @param Y Expression matrix of mRNAs whose row names will be used to generate the list of mRNAs.
#' @param verbose Logical. Shall progress be printed?
#' @return Allocation matrix A necessary for "miR.test" function.
#' @author Stephan Artmann
#' @examples
##MAINEXAMPLE
generate.A = function (df,X=NULL,Y=NULL,verbose=TRUE) {
 mRNA = unique(df[,1]);
 miRNA = unique(df[,2]);
 if (!is.null(X)) {
  miRNA = unique(rownames(X));
  if (length(miRNA) != nrow(X)) {
   print("Some miRNAs in X occur more often than once!")
   warning("Some miRNAs in X occur more often than once!")
  }
 }
 if (!is.null(Y)) {
  mRNA = unique(rownames(Y));
  if (length(mRNA) != nrow(Y)) {
   print("Some mRNAs in Y occur more often than once!")
   warning("Some mRNAs in Y occur more often than once!")
  }
 }
 colnames(df) = c("mRNA","miRNA")
 if (verbose) {
  progress = seq(1,length(miRNA),length.out=20);
 }
 A = rep(0,length(miRNA)*length(mRNA));
 dim(A) = c(length(mRNA),length(miRNA));
 colnames(A) = miRNA;
 rownames(A) = mRNA;

 for (i in 1:length(miRNA)) {
  if (verbose && i %in% round(progress)) print(paste("miRNA",i,"of",length(miRNA),"at",Sys.time()));
  x = df[df$miRNA == miRNA[i],];
  A[match(x$mRNA,rownames(A)),i] = 1;
 }
 A
}

#' Internal algorithm: Make limma test one-sided
#' @param fit Result of "lmFit" and "eBayes" functions in "limma" package.
#' @param lower Shall one-sided p-value indicated down-regultation?
limma.one.sided = function (fit,lower=FALSE) {
  se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total <- fit$df.prior + fit$df.residual
  pt(fit$t, df=df.total, lower.tail=lower)[,2]
}

#' internal algorithm for author's convenience. Create a linear model with the limma package.
#' @param X Expression matrix.
#' @param group Group membership of replicates.
#' @param design Design as specified in limma (design matrix, see model.matrix).
#' @author Stephan Artmann
limma.test = function (X,group=NULL,design=NULL) {
  if (!is.null(group) & !is.null(design)) stop ("Just specify group *or* design in limma.test()")
  if (!is.null(group)) {
   design = model.matrix(~group);
  }
  fit = lmFit(X, design)
  fit = eBayes(fit)
  fit;
}

#' Internal function for gene set testing.
#' @param A Allocation matrix as in "miR.test" function.
#' @param X miRNA expression matrix as in `miR.test' function. Only necessary when allocation.matrix=TRUE.
#' @param Y mRNA expression matrix as in "miR.test" function.
#' @param group group as in `miR.test' function
#' @param tests Test applied, sie gene.set.tests  
#' @param permutation Shall permutation procedure for global tests be applied? Put 'FALSE' to use approximate results or give a number for the number of permutations.
#' @param nrot Number of rotations of rotation tests. Defaults to 1000 to be able to show p-values as low as 10^-3.
#' @param design If specified, group will be ignored. Design matrix as used in `limma' package. Cannot be used with global tests.
#' @param allocation.matrix Logical, is A an allocation matrix with mRNAs in its columns and miRNAs in its rows, or is it an allocation data.frame?
#' @param verbose Defaults to FALSE. If TRUE, progress is printed.
#' @return List of the following, for up- and for down-regulation: Matrix with testing results for every gene set in its rows and the applied gene set test in its columns.
#' @references 
#' Artmann, Stephan and Jung, Klaus and Bleckmann, Annalen and Beissbarth, Tim (submitted).
#' Detection of simultaneous group effects in microRNA expression and 
#' related functional gene sets.
#'
#' Brunner, E. (2009) Repeated measures under non-sphericity.
#' Proceedings of the 6th St. Petersburg Workshop on Simulation,
#' 605-609.
#' 
#' Jelle J. Goeman, Sara A. van de Geer, Floor de Kort, Hans C. van
#' Houwelingen (2004) A global test for groups of genes: testing
#' association with a clinical outcome. Bioinformatics 20, 93-99.
#'
#' Jung, Klaus and Becker, Benjamin and Brunner, Edgar and Beissbarth, Tim (submitted).
#' Comparison of Global Tests for Functinoal Gene Sets in
#' Two-Group Designs and Selection of Potentially
#' Effect-causing Genes.
#' 
#' Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M,
#' Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander,
#' WS, Hilton, DJ, and Blewitt, ME (2010). Opposing roles of polycomb
#' repressive complexes in hematopoietic stem and progenitor cells.
#' _Blood_, published online 5 May 2010.
#'
#' Mansmann, U. and Meister, R., 2005, Testing differential gene
#' expression in functional groups, _Methods Inf Med_ 44 (3).
#' 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for
#' assessing differential expression in microarray experiments.
#' _Statistical Applications in Genetics and Molecular Biology_,
#' Volume *3*, Article 3.
#'
#' Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader,
#' JE, and Smyth, GK (2010). ROAST: rotation gene set tests for
#' complex microarray experiments. _Bioinformatics_, published online
#' 7 July 2010.
#' 
#' @author Stephan Artmann
gs.test = function(A,X=NULL,Y,group=NULL,tests,permutation=FALSE,nrot=1000,design=NULL,allocation.matrix=FALSE,verbose=FALSE) {
 # Load required libraries
 #library(limma)
 #library(globaltest)
 #library(GlobalAncova)
 #library(RepeatedHighDim)

 ga.method = "approx";
 ga.perm = 0;
 gt.perm = 0;
 if (permutation) {
  ga.method="perm";
  gt.perm = permutation;
  ga.perm = permutation;
 }
 if (length(tests) == 0) stop ("No gene set tests specified in gs.test!");
 if (!is.null(design) && !is.null(group)) stop ("Group and design specified, aborting")
 if (!is.null(design)) tests = tests[tests != "GA" & tests != "globaltest" & tests != "RHD"]
 if (length(tests) == 0) stop ("Only competitive gene set tests and `ROAST' can be applied to design matrices! This is not yet implemented for the remaining tests.");
 # Prepare the p-value matrix depending on which tests have been chosen
 testNo = length(tests);
 miRno = ncol(A);
 if (!allocation.matrix) {
  miRs = rownames(X);
  miRno = length(miRs);
 }
 P.l = rep(NA,miRno*(testNo));
 dim(P.l) = c(miRno,(testNo));
 P.h = P.l;

 # Do the testing for the gene set tests chosen
 gt.options(transpose=TRUE);
 if (is.null(design)) design = model.matrix(~group);
 
 if (is.null(design)) {
  fit = limma.test(Y,group=group); 
 } else {
  fit = limma.test(Y,design=design); 
 }
 fit.l = limma.one.sided(fit,lower=TRUE);
 fit.h = 1 - fit.l;
 fit.l.adj = p.adjust(fit.l,method="BH");
 fit.h.adj = p.adjust(fit.h,method="BH");

 rank.low = rank(fit.l,ties.method="random");
 rank.high = rank(1-fit.l,ties.method="random");
 M = fit$coef[,2];
 L.romer = list();
 L.roast = list();

 if (verbose) {
  progress = seq(1,miRno,length.out=20);
 }


 for (j in 1:miRno) {
  if (verbose && j %in% round(progress)) print (paste("Gene Set",j,"of",miRno,"Gene Sets at",Sys.time()));
  if (allocation.matrix) {
   index =(A[,j] == 1);
  } else {
   if (is.null(rownames(Y))) stop ("Please specify the gene names as the row names of Y. Otherwise it is impossible to match the data in A to the genes in Y.")
   index = match(A[A[,2] == miRs[j],1],rownames(Y));
   ind = rep(FALSE,nrow(Y));
   ind[index] = TRUE;
   index = ind;
  }
  if ("KS" %in% tests) {
   if (length(which(index) > 2)) {
    ks.rank.low = ks.test(rank.low[index],"punif",min=1,max=max(rank.low),alternative="greater")$p.value;
    ks.rank.high = ks.test(rank.high[index],"punif",min=1,max=max(rank.high),alternative="greater")$p.value;
   } else {
    ks.rank.low = NA;
    ks.rank.high = NA;
   }
   P.l [j,match("KS",tests)] = ks.rank.low;
   P.h [j,match("KS",tests)] = ks.rank.high;
  }
  if ("W" %in% tests) {
    if (length(which(index) > 2)) {
     w.rank.low = wilcox.test(rank.low[index],rank.low[-which(index)],mu=0,paired=FALSE,alternative="less")$p.value;
     w.rank.high = wilcox.test(rank.high[index],rank.high[-which(index)],mu=0,paired=FALSE,alternative="less")$p.value;
    } else {
     w.rank.low = NA;
     w.rank.high = NA;
    }
    P.l [j,match("W",tests)] = w.rank.low;
    P.h [j,match("W",tests)] = w.rank.high;
  }
  if ("Fisher" %in% tests) {
   if (length(which(index) > 2)) {
    f.adj.low = fisher.test(contingency.table(gene.set=which(index),p.val=fit.l.adj),alternative="greater")$p.value;
    f.adj.high = fisher.test(contingency.table(gene.set=which(index),p.val=fit.h.adj),alternative="greater")$p.value;
   } else {
    f.adj.low = NA;
    f.adj.high = NA;
   }
    P.l [j,match("Fisher",tests)] = f.adj.low;
    P.h [j,match("Fisher",tests)] = f.adj.high;
  }
  if ("globaltest" %in% tests) {
   if (length(which(index & (M<0))) > 2) {
    gt.low = gt(group,Y[index & M<0,],permutations=gt.perm)@result[,1];
   } else {
    gt.low = NA;
   }
   P.l [j,match("globaltest",tests)] = gt.low;
   if (length(which(index & (M>0))) > 2) {
    gt.high = gt(group,Y[index & M>0,],permutations=gt.perm)@result[,1];
   } else {
    gt.high = NA;
   }
   P.h [j,match("globaltest",tests)] = gt.high;
  }
  if ("GA" %in% tests) {
   if (length(which(index & (M<0))) > 2) {
    sink("ga.out")
    ga.low = GlobalAncova(Y[index & M<0,],group=group,method=ga.method,perm=ga.perm)$test.result[2];
    sink();
   } else {
    ga.low = NA;
   }
   P.l [j,match("GA",tests)] = ga.low;
   if (length(which(index & (M>0))) > 2) {
    sink("ga.out")
    ga.high = GlobalAncova(Y[index & M>0,],group=group,method=ga.method,perm=ga.perm)$test.result[2];  
    sink()
   } else {
    ga.high = NA;
   }
   P.h [j,match("GA",tests)] = ga.high;
  }
  if ("RHD" %in% tests) {
   if (is.factor(group) && length(levels(group)) == 2) {
    l = levels(group);
    if (length(which(index & (M<0))) > 2) {
     B.low = RepeatedHighDim(Y[index & M<0,which(group==l[1])],Y[index & M<0,which(group == l[2])],paired=FALSE)$p;
    } else {
     B.low = NA;
    }
    if (length(which(index & (M>0))) > 2) {
     B.high = RepeatedHighDim(Y[index & M>0,which(group==l[1])],Y[index & M>0,which(group == l[2])],paired=FALSE)$p;
    } else {
     B.high = NA;
    }
    P.l [j,match("RHD",tests)] = B.low;
    P.h [j,match("RHD",tests)] = B.high;
   } else {
    P.l [j,match("RHD",tests)] = NA;
    P.h [j,match("RHD",tests)] = NA;
    warning("RepeatedHighDim only supports two-group comparisons!")
   }
  }
  if ("roast" %in% tests) {
   if (length(which(index) > 2)) {
    L.roast[[length(L.roast) + 1]] = which(index);
   } else {
    L.roast[[length(L.roast) + 1]] = rep(FALSE,nrow(Y));
   }
  }
  if ("romer" %in% tests) {
   if (length(which(index) > 2)) {
    L.romer[[length(L.romer) + 1]] = which(index);
   } else {
    L.romer[[length(L.romer) + 1]] = rep(NA,nrow(Y));
   }
  }
 }
 if ("roast" %in% tests) {
   if(verbose) print("Starting ROAST procedure...")
   Roast = mroast(L.roast,y=Y,design=design,nrot=nrot,adjust.method="none")$P.Value;
   P.l[,match("roast",tests)] = Roast[,3];
   P.h[,match("roast",tests)] = Roast[,2];
   if(verbose) print("Finished ROAST procedure...")
 }
 if ("romer" %in% tests) {
   if(verbose) print("Starting romer procedure...")
  Romer = romer(index=L.romer,y=Y,design=design,nrot=nrot);
  P.l[,match("romer",tests)] = Romer[,3];
  P.h[,match("romer",tests)] = Romer[,2];
   if(verbose) print("Finished romer procedure...")
 }
 list(low=P.l,high=P.h);
}

#' Fisher method of p value combination.
#' @param p1,p2 one-sided p-values that shall be combined.
#' @param check.range If set to "TRUE" values above 1 will be set to 1.
#' @return Combined p-value.
#' @author Stephan Artmann
fisher.combination = function (p1,p2,check.range=FALSE) {
 if (check.range) {
  if (p1 > 1) {
	  p1 = 1;
  }
  if (p2 > 1) {
	  p2 = 1;
  } 
 }
	s = -2*(log(p1) + log(p2));
	(1-pchisq(s,df=4));
}


#' Inverse-normal method for p value combination.
#' @param p1,p2 one-sided p-values that shall be combined.
#' @return Two-sided combined p-value.
#' @author Stephan Artmann
inverse.normal.combination = function(p1,p2) {
 S = (qnorm(p1) + qnorm(p2))/sqrt(2);
 2*(1-pnorm(abs(S)));
}


#' Internal function for author's convenience and more legible code. Applies a function to every column vector of a matrix and a vector.
#' @param M The matrix for whose column vectors mapply shall be used.
#' @param v The vector.
#' @param FUN The function.
#' @param ... Further arguments to be given to FUN.
#' @author Stephan Artmann
m.combine = function(M,v,FUN,...) {
 E = mapply(FUN,M,v,...);
 dim(E) = dim(M);
 E;
}

#' Main Function of miRtest package.
#' @author Stephan Artmann
#' @param X miRNA expression matrix with genes in rows and replicates in columns
#' @param Y mRNA expression matrix with genes in rows and replicates in columns
#' @param A Allocation data.frame or Allocation matrix. An allocation data.frame contains the mRNAs in its first column and the miRNAs in its second column. See vignette `miRtest' for information on Allocation matrices.
#' @param group.miRNA Vector of miRNA group membership, being either numeric or a factor (**this makes a difference**). E. g. if you have four replicates in a control group and three replicates in a treated group, you may choose c(1,1,1,1,2,2,2)
#' @param design.miRNA If specified, group.miRNA will be ignored. Here you can specify a design matrix as it is returned from the model.matrix `limma' function.
#' @param design.mRNA If specified, group.mRNA will be ignored. Here you can specify a design matrix as it is returned from the model.matrix `limma' function.
#' @param group.mRNA Vector of mRNA group membership, being either numeric or a factor (**this makes a difference**).E. g. if you have four replicates in a control group and three replicates in a treated group, you may choose c(1,1,1,1,2,2,2)
#' @param gene.set.tests Test to be applied for gene set testing. Can be one or more of the following: `globaltest', `GA', `RHD', `KS', `W', `Fisher', `roast', `romer', or `all' if you want to do all tests.
#' @param adjust Muliple hypothesis testing adjustment. Same options as in "p.adjust" function.
#' @param permutation Number of permutations for `globaltest' or `GlobalAncova' gene set tests. Put to "FALSE" to use the approximate p-values instead of permutation ones.
#' @param nrot Number of rotations for rotation tests `ROAST' and `romer'
#' @param allocation.matrix Logical, is A an allocation matrix with mRNAs in its columns and miRNAs in its rows, or is it an allocation data.frame?
#' @param verbose Defaults to FALSE. If TRUE, output on progress is printed.
#' @param errors Defaults to TRUE. If set to FALSE, some errors checking correct sizes of matrices are turned into warning messages.
#' @return Matrix with testing results for every miRNA in its rows and the applied gene set test in its columns. Note that result will depend on whether multiple hypothesis testing correction was applied or not.
#' @references 
#' Artmann, Stephan and Jung, Klaus and Bleckmann, Annalen and Beissbarth, Tim (submitted).
#' Detection of simultaneous group effects in microRNA expression and 
#' related functional gene sets.
#'
#' Brunner, E. (2009) Repeated measures under non-sphericity.
#' Proceedings of the 6th St. Petersburg Workshop on Simulation,
#' 605-609.
#' 
#' Jelle J. Goeman, Sara A. van de Geer, Floor de Kort, Hans C. van
#' Houwelingen (2004) A global test for groups of genes: testing
#' association with a clinical outcome. Bioinformatics 20, 93-99.
#'
#' Jung, Klaus and Becker, Benjamin and Brunner, Edgar and Beissbarth, Tim (submitted).
#' Comparison of Global Tests for Functinoal Gene Sets in
#' Two-Group Designs and Selection of Potentially
#' Effect-causing Genes.
#' 
#' Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M,
#' Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander,
#' WS, Hilton, DJ, and Blewitt, ME (2010). Opposing roles of polycomb
#' repressive complexes in hematopoietic stem and progenitor cells.
#' _Blood_, published online 5 May 2010.
#'
#' Mansmann, U. and Meister, R., 2005, Testing differential gene
#' expression in functional groups, _Methods Inf Med_ 44 (3).
#' 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for
#' assessing differential expression in microarray experiments.
#' _Statistical Applications in Genetics and Molecular Biology_,
#' Volume *3*, Article 3.
#'
#' Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader,
#' JE, and Smyth, GK (2010). ROAST: rotation gene set tests for
#' complex microarray experiments. _Bioinformatics_, published online
#' 7 July 2010.
#' 
#' @examples 
##MAINEXAMPLE
miR.test = function (X,Y,A,group.miRNA=NULL,group.mRNA=NULL,gene.set.tests="romer",design.miRNA=NULL,design.mRNA=NULL,adjust="none",permutation=FALSE,nrot=1000,allocation.matrix=FALSE,verbose=FALSE,errors=TRUE) {

# library(limma)

 ### Convert data.frames into matrices
 X = as.matrix(X);
 Y = as.matrix(Y);

 ### Check Data Input ### 
 if (length(gene.set.tests) == 0) stop("Please provide gene.set.tests")
 if (allocation.matrix) {
  if (!(ncol(A) == nrow(X))) stop("Number of columns of A must equal number of rows of X")
  if (errors && !(nrow(A) == nrow(Y))) stop("Number of rows of A must equal number of rows of Y. Check that Y does not have duplicate row names. You can disable this error message with errors=FALSE. Note that then genes ocurring more often than once will have a larger weight in the gene set test.")
 }
 if (length(gene.set.tests) == 1 && gene.set.tests == "all") gene.set.tests = c("globaltest","GA","RHD","KS","W","Fisher","roast","romer")
 if (!(all(gene.set.tests %in% c("globaltest","GA","RHD","KS","W","Fisher","roast","romer")))) stop("Check gene.set.tests and enter only one or more of the following tests: globaltest, GA, KS, RHD, W, Fisher, roast, romer or all if you want to do all tests")
 if (!is.null(group.miRNA) & !is.null(design.miRNA)) warning("group.miRNA will be ignored as design.miRNA is specified")
 if (!is.null(group.miRNA) & !is.null(group.mRNA)) {
  if (!all(levels(group.miRNA) == levels(group.mRNA))) stop ("Group names of miRNA samples must be the same as of the mRNA samples. Aborting");
 print("Assuming that group names of miRNA samples are the same as of mRNA samples!");
 }
  if("roast" %in% gene.set.tests) {
   print("Note: For compatibility reasons ROAST is not available in this version.");
   print(" It will be added to the next version of miRtest.");
   print(" To use miRtest with ROAST refer to older versions.");
   gene.set.tests = gene.set.tests[gene.set.tests != "roast"];
 }
 if (length(gene.set.tests) == 0) stop("Please provide gene.set.tests")
 if (!is.null(design.mRNA)) gene.set.tests = gene.set.tests[gene.set.tests != "GA" & gene.set.tests != "globaltest" & gene.set.tests != "RHD"]

 ### Order the allocation matrix ###
 if(allocation.matrix) {
  if (!any(is.null(colnames(A)),is.null(rownames(X)))) {
   A[,order(colnames(A))];
   X[order(rownames(X)),];
   if (!all(colnames(A) == rownames(X))) warning("Column names of A are not equal to rownames of X")
  } else {
   warning("Allocation matrix A has no colnames and/or miRNA matrix X has no rownames. Assuming that columns in A are in the same order as rows in X.")
  }
  if (!any(is.null(colnames(A)),is.null(rownames(Y)))) {
   A[order(rownames(A)),];
   Y[order(rownames(Y)),];
   if(!all(rownames(A) == rownames(Y))) warning("Row names of A are not equal to rownames of Y")
  } else {
   warning("Allocation matrix A and/or miRNA matrix X has no rownames. Assuming that rows in A are in the same order as rows in Y.")
  }
 } else {
  if (is.null(rownames(X))) stop("Please specify row names of X");
  if (is.null(rownames(Y))) stop("Please specify row names of Y");
 }
 if (!is.null(design.mRNA)) gene.set.tests = gene.set.tests[gene.set.tests != "globaltest" & gene.set.tests != "GA"];

 ### Do miRNA-wise testing ###
 if (!is.null(design.miRNA)) {
  miR = limma.test(X,design=design.miRNA);
 } else {
  miR = limma.test(X,group=group.miRNA);
 }
 miR.l = limma.one.sided(miR,lower=TRUE);
 miR.h = 1-miR.l;

 ### Do gene set testing ###
 tests = gene.set.tests;
 if(!is.null(design.mRNA)) {
  GS = gs.test(A,X,Y,group=NULL,tests,permutation=permutation,nrot=nrot,design=design.mRNA,allocation.matrix=allocation.matrix,verbose=verbose);
 } else {
  GS = gs.test(A,X,Y,group.mRNA,tests,permutation=permutation,nrot=nrot,allocation.matrix=allocation.matrix,verbose=verbose);
 }

 ### Combine the results ###
 rot.tests = rep(FALSE,length(tests));
 if ("roast" %in% tests) {
  rot.tests[match("roast",tests)] = TRUE;
 }
 if ("romer" %in% tests) {
  rot.tests[match("romer",tests)] = TRUE;
 }
 if ("W" %in% tests) {
  rot.tests[match("W",tests)] = TRUE;
 }
 P.l = GS$low;
 P.h = GS$high;

 P = rep(NA,nrow(P.l)*ncol(P.l));
 dim(P) = dim(P.l); # two-sided p-values

 if (length(which(!rot.tests)) > 0) {
  P.l[,!rot.tests] = m.combine(P.l[,!rot.tests],miR.h,fisher.combination)            ### p-value for up-regulation in miRNA
  P.h[,!rot.tests] = m.combine(P.h[,!rot.tests],miR.l,fisher.combination)            ### p-value for down-regulation in miRNA
  P[,!rot.tests] = apply(2*pmin(P.l[,!rot.tests,drop=FALSE],P.h[,!rot.tests,drop=FALSE]),c(1,2),min,1); 
 }
 if (length(which(rot.tests)) > 0) {
  P.l[,rot.tests] = m.combine(P.l[,rot.tests],miR.h,inverse.normal.combination);
  P.l[,rot.tests] [abs(P.l[,rot.tests] - miR.h) == 1] = 1;                           ### remove NaN that results from combining 0 and 1
  P.h[,rot.tests] = m.combine(P.h[,rot.tests],miR.l,inverse.normal.combination);
  P.l[,rot.tests] [abs(P.h[,rot.tests] - miR.l) == 1] = 1;                           ### remove NaN that results from combining 0 and 1
  P[,rot.tests] = pmax(P.l[,rot.tests],P.h[,rot.tests]);                             ### as rotation tests can return 1 we have to take the maximum here
 }
 P = apply(P,2,p.adjust,method=adjust);
 colnames(P) = tests;
 if (!is.null(colnames(A))) {
  if(allocation.matrix) {
   rownames(P) = colnames(A)
  } else {
   rownames(P) = rownames(X)
  }
 }
 if (is.null(dim(P)) || ncol(P) == 1) colnames(P) = "miRtest"
 P;
}


