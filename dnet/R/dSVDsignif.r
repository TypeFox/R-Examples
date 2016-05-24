#' Function to obtain SVD-based gene significance from the input gene-sample matrix
#'
#' \code{dSVDsignif} is supposed to obtain gene signficance from the given gene-sample matrix according to singular value decomposition (SVD)-based method. The method includes: 1) singular value decomposition of the input matrix; 2) determination of the eigens in consideration (if not given); 3) construction of the gene-specific project vector based on the considered eigens; 4) calculation of the distance statistic from the projection vector to zero point vector; and 5) based on distance statistic to obtain the gene significance.
#'
#' @param data an input gene-sample data matrix used for singular value decomposition
#' @param num.eigen an integer specifying the number of eigens in consideration. If NULL, this number will be automatically decided on based on the observed relative eigenexpression against randomised relative eigenexpression calculated from a list (here 100) of permutated input matrix
#' @param pval.eigen p-value used to call those eigens as dominant. This parameter is used only when parameter 'num.eigen' is NULL. Here, p-value is calcualted to assess how likely the observed relative eigenexpression are more than the maximum relative eigenexpression calculated from permutated matrix
#' @param signif the singificance to return. It can be either "pval" for using the p-value as the gene significance, or "fdr" for using the fdr as the gene significance
#' @param orient.permutation the orientation of matrix being permutated. It can be either "row" to permutate values within each row, or "column" to permutate values within each column, or "both" to permutate values both within rows and columns. Notably, when using the p-value as the gene significance, it is always to permutate values within each row.
#' @param num.permutation an integer specifying how many permutations are used
#' @param fdr.procedure the procedure to adjust the fdr. To ensure that the high distance statistic the more significance, the fdr should be adjusted either using "stepup" for step-up procedure (from the most significant to the least significant) or using "stepdown" for step-down procedure (from the least significant to the most significant)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return a vector storing gene significance
#' @note none
#' @export
#' @seealso \code{\link{dFDRscore}}
#' @include dSVDsignif.r
#' @examples
#' \dontrun{
#' # 1) generate data with an iid matrix of 1000 x 9
#' data <- cbind(matrix(rnorm(1000*3,mean=0,sd=1), nrow=1000, ncol=3), 
#' matrix(rnorm(1000*3,mean=0.5,sd=1), nrow=1000, ncol=3), 
#' matrix(rnorm(1000*3,mean=-0.5,sd=1), nrow=1000, ncol=3))
#'
#' # 2) calculate the significance according to SVD
#' # using "fdr" significance
#' fdr <- dSVDsignif(data, signif="fdr", num.permutation=10)
#' # using "pval" significance
#' pval <- dSVDsignif(data, signif="pval", num.permutation=10)
#' }

dSVDsignif <- function(data, num.eigen=NULL, pval.eigen=1e-2, signif=c("fdr","pval"), orient.permutation=c("row","column","both"), num.permutation=100, fdr.procedure=c("stepup","stepdown"), verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
    }
    
    orient.permutation <- match.arg(orient.permutation)
    signif <- match.arg(signif)
    fdr.procedure <- match.arg(fdr.procedure)
    
    # A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            msg <- sprintf("%d (out of %d)", i, B)
            time <- as.character(Sys.time())
            if(flag) message(paste(c(msg, " at ", time), collapse=""), appendLF=T)
        }
    }
    
    if(verbose){
        message(sprintf("First, singular value decomposition of the input matrix (with %d rows and %d columns)...", nrow(data), ncol(data)), appendLF=T)
    }
    ## Singular Value Decomposition of a Matrix
    res <- svd(data)
    
    if(verbose){
        message(sprintf("Second, determinate the eigens..."), appendLF=T)
    }
    if(is.null(num.eigen)){
        if(verbose){
            message(sprintf("\tvia automatically deciding on the number of dominant eigens under the cutoff of %1.2e pvalue", pval.eigen), appendLF=T)
        }
        ## automatically decide on the number of dominant eigens
        
        ## relative eigenexpression
        re_obs <- res$d*res$d/sum(res$d*res$d)
        
        ## decide on the significant eigenexpression. That is, how likely (e.g. at least 95%) they are more than the maximum relative eigenexpression calculated from permutated matrix
        B <- 100
        re_max <- sapply(1:B, function(i){

            if(verbose){
                progress_indicate(i, B, 10)
            }
            
            ## permuataion
            if(orient.permutation=="row"){
                data_perm <- apply(data, 1, sample)
                if(nrow(data) > ncol(data)) data_perm <- t(data_perm)
            }else if(orient.permutation=="column"){
                data_perm <- apply(data, 2, sample)
            }else if(orient.permutation=="both"){
                # row-wise
                data_perm <- apply(data, 1, sample)
                if(nrow(data) > ncol(data)) data_perm <- t(data_perm)
                # column-wise
                data_perm <- apply(data_perm, 2, sample)
            }
            ## maximum of relative eigenexpression from the permuated matrix
            re_perm <- svd(data_perm)
            max(re_perm$d*re_perm$d/sum(re_perm$d*re_perm$d))
        })
        ## pval for each eigenexpression
        re_pval <- sapply(re_obs, function(x){
            sum(x<re_max)/B
        })
        
        if(is.null(pval.eigen) | is.na(pval.eigen)){
            pval.eigen <- 1e-2
        }
        num.eigen <- sum(re_pval < pval.eigen)   
    }
    
    if(num.eigen > ncol(data) | num.eigen <=0 | num.eigen!=as.integer(num.eigen)){
        stop("The number of the significant eigens is overflown/valid.\n")
    }
    
    if(verbose){
        message(sprintf("\tnumber of the eigens in consideration: %d", num.eigen), appendLF=T)
    }
    
    if(verbose){
        message(sprintf("Third, construct the gene-specific projection vector,and calculate distance statistics..."), appendLF=T)
    }
    ## Construct the observed gene-specific projection vector,and calculate distance statistics DS from the projection vector to zero point vector
    DV <- data %*% res$v[,1:num.eigen] # Projection of all prototype vectors from codebook matrix into SVD subspaces
    DS <- sqrt(apply(DV*DV,1,sum)) # Calculation of the Euclidian distance from the projection vector to zero point vector
    
    if(signif=="pval"){
    
        if(verbose){
            message(sprintf("Finally, obtain gene significance (p-value) based on %d of row-wise permutations...", num.permutation), appendLF=T)
        }
        
        DS_pval <- sapply(1:nrow(data), function(i){
            if(verbose){
                progress_indicate(i, nrow(data), 10)
            }
            obs <- as.matrix(data[i,],ncol=1)
            exp <- apply(obs[,rep(1,num.permutation)],2,sample)
            DV_perm <- t(exp) %*% res$v[,1:num.eigen]
            DS_perm <- sqrt(apply(DV_perm*DV_perm,1,sum))
            sum(DS[i]<DS_perm)/num.permutation
        })
        if(!is.null(rownames(data))){
            names(DS_pval) <- rownames(data)
        }
    }else if(signif=="fdr"){
    
        if(verbose){
            message(sprintf("Finally, obtain gene significance (fdr) based on %d permutations...", num.permutation), appendLF=T)
        }
    
        if(verbose){
            message <- paste(c("\tdoing ", orient.permutation, "-wise permutations..."), collapse="")
            message(message, appendLF=T)
        }
        
        DS_perm <- sapply(1:num.permutation, function(i){
            if(verbose){
                progress_indicate(i, num.permutation, 10)
            }

            ## permuataion
            if(orient.permutation=="row"){
                data_perm <- apply(data, 1, sample)
                if(nrow(data) > ncol(data)) data_perm <- t(data_perm)
            }else if(orient.permutation=="column"){
                data_perm <- apply(data, 2, sample)
            }else if(orient.permutation=="both"){
                # row-wise
                data_perm <- apply(data, 1, sample)
                if(nrow(data) > ncol(data)) data_perm <- t(data_perm)
                # column-wise
                data_perm <- apply(data_perm, 2, sample)
            }
            
            DV_perm <- data_perm %*% res$v[,1:num.eigen]
            sqrt(apply(DV_perm*DV_perm,1,sum))
        })

        if(verbose){
            message <- paste(c("\testimating fdr..."), collapse="")
            message(message, appendLF=T)
        }
        
        DS_sorted <- sort.int(DS, decreasing=T, index.return=T)
        fdr_sorted <- vector()
        for(i in 1:length(DS_sorted$x)){
        
            if(verbose){
                progress_indicate(i, length(DS_sorted$x), 10, flag=T)
            }
            fdr_sorted[i] <- stats::median(apply(DS_perm > DS_sorted$x[i], 2, sum))/i
        }

        if(verbose){
            message <- paste(c("\tusing ", fdr.procedure, " procedure..."), collapse="")
            message(message, appendLF=T)
        }
        DS_fdr <- rep(0, length(fdr_sorted))
        if(fdr.procedure=="stepup"){
            ## step-up procedure (from the most significant to the least significant)
            for(i in 1:length(fdr_sorted)){
                if(i>=2 && fdr_sorted[i] < fdr_sorted[i-1]){
                    fdr_sorted[i] <- fdr_sorted[i-1]
                }
                DS_fdr[DS_sorted$ix[i]] <- fdr_sorted[i]
            }
        }else if(fdr.procedure=="stepdown"){
            ## step-down procedure (from the least significant to the most significant)
            for(i in length(fdr_sorted):1){
                if(i<length(fdr_sorted) && fdr_sorted[i] > fdr_sorted[i+1]){
                    fdr_sorted[i] <- fdr_sorted[i+1]
                }
                DS_fdr[DS_sorted$ix[i]] <- fdr_sorted[i]
            }
        }
        if(!is.null(rownames(data))){
            names(DS_fdr) <- rownames(data)
        }
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    
    if(signif=="pval"){
        return(DS_pval)
    }else if(signif=="fdr"){
        return(DS_fdr)
    }
}
