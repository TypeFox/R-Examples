
#' Extended InterVA method for non-standard input
#'
#' @param data  A matrix input, or data read from csv files. Sample input is included as \code{data(RandomVA3)}.
#' @param train  A matrix input, or data read from csv files in the same format
#' as \code{data}, but with an additional column specifying cause-of-death. Sample input is included as \code{data(RandomVA3)}.
#' @param causes.train the column name of the cause-of-death assignment label in training data.
#' @param causes.table list of causes to consider in the training data. Default to be NULL, which uses all the causes present in the training data.
#' @param thre numerical number between 0 and 1. Symptoms with missing rate higher than \code{thre} in the training data will be dropped from both training and testing data.
#' @param type type of data conversion when calculating conditional probability (probability of each symptom given each cause of death) for InterVA and InSilicoVA models. Both ``quantile'' and ``fixed'' usually give similar results empirically. 
#' \itemize{
#' \item{\code{quantile}: }{the rankings of the P(S|C) are obtained by matching the same quantile distributions in the default InterVA P(S|C)}\item{\code{fixed}: }{P(S|C) are matched to the closest values in the default InterVA P(S|C) table.} \item{\code{empirical}: }{no ranking is calculated, but use the empirical conditional probabilities directly.}
#'}
#' @param prior The prior distribution of CSMF. ``uniform'' uses no prior information, i.e., 1/C for all C causes and ``train'' uses the CSMF in the training data as prior distribution of CSMF.
#' @param ... not used
#'
#' @return fitted \code{interVA} object
#' @importFrom InSilicoVA extract.prob
#' @export interVA.train
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark (2016) \emph{Probabilistic
#' cause-of-death assignment using verbal autopsies.}
#' \url{http://arxiv.org/abs/1411.3042}, \emph{To appear, Journal of the American Statistical Association}
#' @references Zehang R. Li, Tyler H. McCormick, Samuel J. Clark (2014) \emph{InterVA4: An R package to analyze verbal autopsy data.}, \emph{Center for Statistics and the Social Sciences Working Paper, No.146}
#' @references http://www.interva.net/
#' @keywords InterVA4
#' @examples
#' \donttest{
#' data(RandomVA3)
#' test <- RandomVA3[1:200, ]
#' train <- RandomVA3[201:400, ]
#' out <- interVA.train(data = test, train = train, causes.train = "cause", 
#'                      prior = "train", type = "quantile")
#' }
interVA.train <- function(data, train, causes.train, causes.table = NULL, thre = 0.95, type = c("quantile", "fixed", "empirical")[1], prior = c("uniform", "train")[1],  ...){

	# extract conditional probabilities as in InSilicoVA
	prob.learn <- extract.prob(train = train, 
							  gs = causes.train, 
							  gstable = causes.table, 
							  thre = thre, 
							  type = type, 
							  isNumeric = FALSE)
	# remove unused symptoms
	col.exist <- c(colnames(data)[1], causes.train, colnames(prob.learn$symps.train))
	remove <- which(colnames(data) %in% col.exist == FALSE)
	if(length(remove) > 0){
		warning(paste(length(remove), "symptoms deleted from testing data to match training data:", 
			paste(colnames(data)[remove], collapse = ", ")),
			immediate. = TRUE)	
		data <- data[, -remove]
	}

    if(causes.train %in% colnames(data)){
        data <- data[, -which(colnames(data) == causes.train)]
    }

	if(type == "empirical"){
		probbase <- prob.learn$cond.prob
	}else{
		probbase.alpha <- prob.learn$cond.prob.alpha
		table.dev <- prob.learn$table.alpha
		table.num.dev <- prob.learn$table.num

		probbase <- matrix(0, dim(probbase.alpha)[1], dim(probbase.alpha)[2])
		colnames(probbase) <- colnames(probbase.alpha)
		rownames(probbase) <- rownames(probbase.alpha)
		for(i in 1:length(table.dev)){
			probbase[which(probbase.alpha == table.dev[i])] <- table.num.dev[i]
		}
	}

	causetext <- cbind(colnames(probbase),
                       colnames(probbase))
    if(prior == "uniform"){
        Sys_Prior <- rep(1/dim(probbase)[2], dim(probbase)[2])
    }else if(prior == "train"){
        raw <- table(train[, causes.train]) 
        raw <- raw / sum(raw)
        Sys_Prior <- rep(0, dim(probbase)[2])
        Sys_Prior[match(names(raw), colnames(probbase))] <- as.numeric(raw)
    }
	############################
    ## define mid-step functions
    ############################
    
    va <- function(ID , MALPREV, HIVPREV , PREGSTAT, PREGLIK , PRMAT , INDET , CAUSE1, LIK1, CAUSE2 , LIK2 , CAUSE3 , LIK3 , wholeprob, ...){
        ## ID
        ID <- ID
        ## The prevalence of Malaria
        MALPREV <- as.character(MALPREV)
        ## The prevalence of HIV
        HIVPREV <- as.character(HIVPREV)
        ## Make PregStat a character string of length 5
        PREGSTAT <- paste(PREGSTAT,paste(rep(" ",5),collapse=""),collapse="")
        ## Likelihood of PregStat
        PREGLIK <- PREGLIK
        ## Likelihood of Maternal Death
        PRMAT <- PRMAT
        ## Indicator of indeterminate outcome
        INDET <- as.character(INDET)
        ## The full distribution of probability on CODs
        wholeprob <- wholeprob
        va.out <- list(ID = ID, MALPREV = MALPREV, HIVPREV = HIVPREV, PREGSTAT = PREGSTAT, PREGLIK = PREGLIK, PRMAT = PRMAT, INDET = INDET, CAUSE1 = CAUSE1, LIK1 = LIK1, CAUSE2 =CAUSE2, LIK2 = LIK2, CAUSE3 = CAUSE3, LIK3 = LIK3, wholeprob = wholeprob)
        va.out
    }

    Input <- as.matrix(data)
    ## Check if there is any data at all
    if(dim(Input)[1] < 1){
        stop("error: no data input")
    }
    N <- dim(Input)[1]  ## Number of data
    S <- dim(Input)[2] - 1 ## Length of symptoms
    D <- length(Sys_Prior)
    ID.list <- Input[, 1]
    Input <- Input[, -1]
    VAresult <- vector("list",N)
   

    nd <- max(1, round(N / 100))
    np <- max(1, round(N / 10))
    
    ## Calculate the InterVA result one by one
    for(i in 1:N){
        ## print out progress
        if(i %% nd == 0){cat(".")}
        if(i %% np == 0){cat(paste(round(i/N * 100), "% completed\n", sep = ""))}
      
        ## Change input Y/NA into binary value
        Input[i, which(is.na(Input[i, ]))] <- "0"
        Input[i, which(toupper(Input[i, ]) != "Y")] <- "0"
        Input[i, which(toupper(Input[i, ]) == "Y")] <- "1"
        ## Change input as a numerical vactor
        input.current <- as.numeric(Input[i, ])       
        prob <- Sys_Prior
        temp <- which(input.current == 1)
 
        for(jj in 1:length(temp)){
            temp_sub <- temp[jj]
            for(j in 1:D){
                prob[j] <- prob[j] * as.numeric(probbase[temp_sub, j])
            }
        }
        if(sum(prob) > 0){
            prob <- prob / sum(prob)            
        }

        names(prob) <- causetext[,2]
        wholeprob <- prob
        if(max(prob) <= 0.4){
            indet <- "Indet"
            cause1<-lik1<-cause2<-lik2<-cause3<-lik3<-" "
        }
        ## Determine the output of InterVA
        if(max(prob) > 0.4){
            ## Find max likelihood
            indet <- " "
            lik1 <- round(max(prob)*100)
            cause1 <- names(prob)[which.max(prob)]
            ## Delete the max and find the second max
            prob <- prob[-which.max(prob)]
            lik2 <- round(max(prob)*100)
            cause2 <- names(prob)[which.max(prob)]
            ## Not show the second if it is too small
            if(max(prob) < 0.5 * max(prob)) lik2 <- cause2 <- " "
            
            ## Delete the second max and find the third max
            prob <- prob[-which.max(prob)]
            lik3 <- round(max(prob)*100)
            cause3 <- names(prob)[which.max(prob)]
            ## Not show the third if it is too small
            if(max(prob) < 0.5 * max(prob)) lik3 <- cause3 <- " "
        }

        Malaria = NULL
        HIV = NULL
        preg_state = NULL
        lik.preg = NULL
        lik_mat = NULL
        VAresult[[i]] <- va(ID = ID.list[i], MALPREV = Malaria, HIVPREV = HIV, PREGSTAT = preg_state, PREGLIK = lik.preg, PRMAT = lik_mat, INDET = indet, CAUSE1 = cause1, LIK1 = lik1, CAUSE2 =cause2, LIK2 = lik2, CAUSE3 = cause3, LIK3 = lik3, wholeprob = wholeprob)
       
    }
    out <- list(ID = ID.list, VA = VAresult, prior = prior, type = type, dev = TRUE)
    class(out) <- "interVA"
    return(out)
 
}