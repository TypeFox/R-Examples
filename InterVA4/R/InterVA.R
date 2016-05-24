#' Provide InterVA4 analysis on the data input.
#' 
#' This function implements the algorithm in the InterVA4 software.  It
#' produces individual cause of death and population cause-specific mortality
#' fractions.
#' 
#' InterVA performs the same tasks as the InterVA4. The output is saved in a
#' .csv file specified by user. The calculation is based on the conditional and
#' prior distribution of 68 CODs. The function also could save the full
#' probability distibution of each individual to file. All information about
#' each individual is saved to a va class object.
#' 
#' Be careful if the input file does not match InterVA input format strictly.
#' The function will run normally as long as the number of symptoms are
#' correct. Any inconsistent symptom names will be printed in console as
#' warning. If there's wrong match of symptom from warning, please change in
#' the input to correct orders.
#' 
#' @param Input A matrix input, or data read from csv files in the same format
#' as required by InterVA4. Sample input is included as data(SampleInput).
#' @param HIV An indicator of the level of prevalence of HIV. The input should
#' be one of the following: "h"(high),"l"(low), or "v"(very low).
#' @param Malaria An indicator of the level of prevalence of Malaria. The input
#' should be one of the following: "h"(high),"l"(low), or "v"(very low).
#' @param directory The directory to store the output from InterVA4. It should
#' either be an existing valid directory, or a new folder to be created. If no
#' path is given, the current working directory will be used.
#' @param filename The filename the user wish to save the output. No extension
#' needed. The output is in .csv format by default.
#' @param output "classic": The same deliminated output format as InterVA4; or
#' "extended": deliminated output followed by full distribution of cause of
#' death proability.
#' @param append A logical value indicating whether or not the new output
#' should be appended to the existing file.
#' @param replicate A logical value indicating whether or not the calculation
#' should replicate original InterVA4 software exactly.If replicate = F, causes
#' with small probability are not dropped out of calculation in intermediate
#' steps, and a possible bug in original InterVA4 implementation is fixed.  If
#' replicate=T, then the output values will be exactly as they would be from
#' calling the InterVA4 program. Since version 1.6, two control variables are added
#' to control the two bugs respectively. Setting this to TRUE will overwrite both to
#' TRUE.
#' @param replicate.bug1 This logical indicator controls whether or not the bug
#' in InterVA4.2 involving the symptom "skin_les" will be replicated or not. It
#' is suggested to set to FALSE.
#' @param replicate.bug2 This logical indicator controls whether the causes
#' with small probability are dropped out of calculation in intermediate
#' steps or not. It is suggested to set to FALSE.
#' @param groupcode A logical value indicating whether or not the group code
#' will be included in the output causes.
#' @param write A logical value indicating whether or not the output (including 
#' errors and warnings) will be saved to file.
#' @param ... not used
#' @return \item{ID }{identifier from batch (input) file} \item{MALPREV
#' }{selected malaria prevalence} \item{HIVPREV }{selected HIV prevalence}
#' \item{PREGSTAT }{most likely pregnancy status} \item{PREGLIK }{likelihood of
#' PREGSTAT} \item{PRMAT }{ likelihood of maternal death} \item{INDET
#' }{indeterminate outcome} \item{CAUSE1 }{ most likely cause} \item{LIK1 }{
#' likelihood of 1st cause} \item{CAUSE2 }{ second likely cause} \item{LIK2 }{
#' likelihood of 2nd cause} \item{CAUSE3 }{ third likely cause} \item{LIK3 }{
#' likelihood of 3rd cause} \item{wholeprob}{ full distribution of causes of
#' death}
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' @seealso \code{\link{InterVA.plot}}
#' @references http://www.interva.net/
#' @keywords InterVA
#' @examples
#' 
#' data(SampleInput)
#' ## to get easy-to-read version of causes of death make sure the column
#' ## orders match interVA4 standard input this can be monitored by checking
#' ## the warnings of column names
#' 
#' sample.output1 <- InterVA(SampleInput, HIV = "h", Malaria = "l", directory = "VA test", 
#'     filename = "VA_result", output = "extended", append = FALSE, replicate = FALSE)
#' 
#' ## to get causes of death with group code for further usage
#' sample.output2 <- InterVA(SampleInput, HIV = "h", Malaria = "l", directory = "VA test", 
#'     filename = "VA_result_wt_code", output = "classic", append = FALSE, 
#'     replicate = FALSE, groupcode = TRUE)
#' 
InterVA<-function(Input, HIV, Malaria, directory = NULL, filename = "VA_result", output="classic", append=FALSE, groupcode = FALSE, replicate = FALSE, replicate.bug1 = FALSE, replicate.bug2 = FALSE, write = TRUE, ...){
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
    PREGSTAT <- paste(PREGSTAT,paste(rep(" ",5-nchar(PREGSTAT)),collapse=""),collapse="")
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

save.va <- function(x, filename, write){
    ## This function saves va object to file in the deliminated format of InterVA4.
    ## The input is a va object and a filename (without extension).
    ## The output is a .csv file.
    ##
    ## Delete the full probability distribution.
    if(!write){return}

    x <- x[-14]
    x <- as.matrix(x)
    filename <- paste(filename, ".csv", sep = "") 
    write.table(t(x), file=filename, sep = ",", append = TRUE,row.names = FALSE,col.names = FALSE)    
}
save.va.prob <- function(x, filename, write){
    ## This function saves va object to file in the deliminated format of InterVA4
    ## followed by a full probability distribution on CODs.
    ## The input is a va object and a filename (without extension).
    ## The output is a .csv file.
    ##
    ## Extract the full probability distribution.
    if(!write){return}

    prob <- unlist(x[14])
    x <- x[-14]
    ## Reformat the matrix with probability distribution.
    x <- unlist(c(as.matrix(x),as.matrix(prob)))
    filename <- paste(filename, ".csv", sep = "") 
    write.table(t(x), file=filename, sep = ",", append = TRUE,row.names = FALSE,col.names = FALSE)    
}
    ## overwrite replication options if needed
    if(replicate){
      warning("option 'replicate' is turned on, all bugs in InterVA-4 is replicated\n", immediate. = TRUE)
      replicate.bug1 <- TRUE
      replicate.bug2 <- TRUE
    }
    ########################
    ## Read in data files
    ########################
    # if no directory is provided, set to default working directory
    if(is.null(directory)) directory = getwd()
    # create and set the directory if it does not exist
    # if it does exist, then fine, do not need to print warning
    dir.create(directory, showWarnings = FALSE)
    globle.dir <- getwd()
    setwd(directory)
    
    # data(probbase)
    data("probbase", envir = environment())
    probbase <- get("probbase", envir  = environment())
    
    probbase <- as.matrix(probbase)
    # data(causetext) 
    data("causetext", envir = environment())
    causetext <- get("causetext", envir  = environment())
    # decide whether to use group code
    if(groupcode){
    		causetext <- causetext[,-2]
    }else{
    		causetext <- causetext[,-3]
    	}
    
    ## Build the skeleton of the error log.
    if(write){
         cat(paste("Error log built for InterVA", Sys.time(), "\n"),file="errorlog.txt",append = FALSE)
         cat(paste("Warning log built for InterVA", Sys.time(), "\n"),file="warnings.txt",append = FALSE)     
    }
    ######################################################
    ## Input should be a matrix with each rows containing:
    ## Field 1: ID number
    ## Field 2-22: descriptors of death
    ## Field 23-246: Indicators
    ## Input should have proper Column names!
    #######################################################
    Input <- as.matrix(Input)
    ## Check if there is any data at all
    if(dim(Input)[1] < 1){
        stop("error: no data input")
    }
    N <- dim(Input)[1]  ## Number of data
    S <- dim(Input)[2]  ## Length of individial field
    ##  Check if the length of input variable matches the probbase dataset
    if(S != dim(probbase)[1] ){
        stop("error: invalid data input format. Number of values incorrect")
    }
    ## Check if the last field is the correct one
    if(tolower(colnames(Input)[S]) != "scosts"){
        stop("error: the last variable should be 'scosts'")
    }
    ## check the column names and give warning
    data("SampleInput", envir = environment())
    SampleInput <- get("SampleInput", envir  = environment())
    valabels = colnames(SampleInput)
    count.changelabel = 0
    for(i in 1:S){
        if(tolower(colnames(Input)[i]) != tolower(valabels)[i]){
            warning(paste("Input columne '", colnames(Input)[i], "' does not match InterVA standard: '", 
                    valabels[i], "'", sep = ""),
                    call. = FALSE, immediate. = TRUE)
            count.changelabel = count.changelabel + 1
        }         
    }
    if(count.changelabel > 0){
        warning(paste(count.changelabel, "column names changed in input. \n If the change in undesirable, please change in the input to match standard InterVA4 input format."), call. = FALSE, immediate. = TRUE)
        colnames(Input) <- valabels
    }
    
    ## Change conditional probability labels into values
    probbase[probbase=="I"]<-1
    probbase[probbase=="A+"]<-0.8
    probbase[probbase=="A"]<-0.5
    probbase[probbase=="A-"]<-0.2
    probbase[probbase=="B+"]<-0.1
    probbase[probbase=="B"]<-0.05
    probbase[probbase=="B-"]<-0.02
    probbase[probbase=="B -"]<-0.02
    probbase[probbase=="C+"]<-0.01
    probbase[probbase=="C"]<-0.005
    probbase[probbase=="C-"]<-0.002
    probbase[probbase=="D+"]<-0.001
    probbase[probbase=="D"]<-0.0005
    probbase[probbase=="D-"]<-0.0001
    probbase[probbase=="E"]<-0.00001
    probbase[probbase=="N"]<-0
    probbase[probbase==""]<-0
    
    ## Extract Prior distribution from the dataset
    ## The first 13 values are not CODs
    probbase[1,1:13]<-rep(0,13)
    ## The first row in the dataset is the expected value of probs, i.e. priors
    Sys_Prior <- as.numeric(probbase[1,])
    # Number of indicators + 13 description variables. A_group:14-16;B_group:17:76;D_group:77:81
    D <- length(Sys_Prior)
    ## Modify the prior based on HIV and Malaria prevalence
    ## 19 = B_HIVAIDS; 21 = B_MALAR; 39 = B_SICKLE
    HIV <- tolower(HIV)
    Malaria <- tolower(Malaria)
    if(!(HIV %in% c("h", "l", "v"))  || !(Malaria %in% c("h", "l", "v"))){
      stop("error: the HIV and Malaria indicator should be one of the three: 'h', 'l', and 'v'")
    }
    if(HIV == "h") Sys_Prior[19] <- 0.05
    if(HIV == "l") Sys_Prior[19] <- 0.005
    if(HIV == "v") Sys_Prior[19] <- 0.00001
    if(Malaria == "h"){
    	Sys_Prior[21] <- 0.05
    	Sys_Prior[39] <- 0.05
    }
    if(Malaria == "l"){
    	Sys_Prior[21] <- 0.005
    	Sys_Prior[39] <- 0.00001
    }
    if(Malaria == "v"){
    	Sys_Prior[21] <- 0.00001
    	Sys_Prior[39] <- 0.00001
    }
    ## Prepare the output
    ID.list <- rep(NA, N)
    VAresult <- vector("list",N)
    ## If append is FALSE, build the skeleton of the new file for output
    if(write && append == FALSE) {
    	header=c("ID","MALPREV","HIVPREV","PREGSTAT","PREGLIK","PRMAT","INDET",
    	"CAUSE1","LIK1","CAUSE2","LIK2","CAUSE3","LIK3")
    	if(output == "extended") header=c(header,as.character(causetext[,2]))
        write.table(t(header),file=paste(filename,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")

    }
    ## add progress indicators now
    nd <- max(1, round(N / 100))
    np <- max(1, round(N / 10))
    
    ## Calculate the InterVA result one by one
    for(i in 1:N){
        ## print out progress
        if(i %% nd == 0){cat(".")}
        if(i %% np == 0){cat(paste(round(i/N * 100), "% completed\n", sep = ""))}
      
        ## Save the current death ID
        index.current <- as.character(Input[i, 1])
        ## Change input Y/NA into binary value
        Input[i, which(is.na(Input[i, ]))] <- "0"
        Input[i, which(toupper(Input[i, ]) != "Y")] <- "0"
        Input[i, which(toupper(Input[i, ]) == "Y")] <- "1"
        ## Change input as a numerical vactor
        input.current <- as.numeric(Input[i,])
        input.current[1] <- 0
        ## Check if age is specified in the input
        ## If not specified, mark as error and skip the case
        if(sum(input.current[2:8]) < 1 ){
            if(write){
              cat(paste(index.current," Error in age indicator: Not Specified ","\n"), file="errorlog.txt", append=TRUE)
            }
            next
        }

        ## Check if sex is specified in the input
        ## If not, mark as error and skip the case
        if(sum(input.current[9:10]) < 1){
            if(write){
              cat(paste(index.current," Error in sex indicator: Not Specified ","\n"), file="errorlog.txt", append=TRUE)
            }
            next
        }
        ## Check if there is any symptoms
        ## 2-22 & 224-246 are not symptoms, but personal profile, or life style
        ## This range is set in the InterVA file
        if(sum(input.current[23:223]) < 1 ){
            if(write){
              cat(paste(index.current," Error in indicators: No symptoms specified ","\n"), file="errorlog.txt", append=TRUE)
            }
            next
        }
        
        ## Repeat twice the check of "ask if" and "don't ask".
        ## If there is contradictory with "ask if" or "don't ask", follow the following rules:
        ## If B is the "don't ask" for A but B has value 1 --> make sure A has value 0;
        ## If B is the "ask if" for A but B has value 0 --> change B into value 1
        for(k in 1:2){
            for(j in 1:(S-1)){
                if(input.current[j + 1] == 1 ){
                    # Note here the first element in input is index; the first element in probbase is expected.
                    Dont.ask <- probbase[j + 1, 4:11]
                    Dont.ask.list <- input.current[match(toupper(Dont.ask), toupper(colnames(Input)))]
                    Dont.ask.list[ is.na(Dont.ask.list)] <- 0
                    
                    if( sum( Dont.ask.list ) > 0 ) {
                    	input.current[j + 1] <- 0
                        	if(write){
                                cat(index.current, "   ", paste(probbase[j+1, 2], "  value inconsistent with ", Dont.ask[which(Dont.ask.list > 0)], " - cleared in working file \n"), file="warnings.txt", append=TRUE)
                            }
                    	}
                 }
                 # Note input.current[j+1] might be changed in the step above!
                     if(input.current[j + 1] == 1 ){
                    # Note here the first element in input is index; the first element in probbase is expected.
                    Ask.if <- probbase[j + 1, 12]
                    if( !is.na(match(toupper(Ask.if), toupper(colnames(Input))))  ){
                        if(input.current[match(toupper(Ask.if), toupper(colnames(Input)) )] == 0){
                            input.current[match(toupper(Ask.if), toupper(colnames(Input)) )] <- 1
                            if(write){
                                cat(index.current, "   ", paste(probbase[j+1, 2], "  not flagged in category ", Ask.if, " - updated in working file \n"), file="warnings.txt", append=TRUE)
                            }
                        }
                    }
                }   
            }
        }
        
        
        ## This seems to be a bug in InterVA
        ## So if the user wishes to replicate entirely as InterVA
        ## the replicate option should be set to TRUE
        ## effect: whenever skin = 1 --> skin_les = 1
        if(replicate.bug1 == TRUE && input.current[84] == 1){
        	input.current[85] <- 1
        }
        
        
        
        ## Initialize ReproductiveAge, Preg_State and Likelihood of Preg
        reproductiveAge <- 0
        preg_state <- " "
        lik.preg <- " "
        ## Determine if at ReproductiveAge
        if(input.current[10] == 1 && (input.current[4] == 1 || input.current[5]==1) ) reproductiveAge <- 1
        ## Find the indicator of Symptoms
        prob <- Sys_Prior[14:D] #The first 13 fields are not indicators
        temp <- which(input.current[2:length(input.current)] == 1)
        
        # Calculate likelihood for each CODs
        # loop through each indicator
        for(jj in 1:length(temp)){
        	temp_sub <- temp[jj]
        	for(j in 14:D){
            prob[j-13] <- prob[j-13] * as.numeric(probbase[temp_sub + 1, j])
        }
        # Normalize A group
        if(sum(prob[1:3]) > 0) prob[1:3] <- prob[1:3]/sum(prob[1:3])
        # Normalize B group 
		if(sum(prob[4:63]) > 0) prob[4:63] <- prob[4:63]/sum(prob[4:63])
        # delete too small probs
        if(replicate.bug2){prob[prob < 0.000001] <- 0}
        }
              
        names(prob) <- causetext[,2]
        prob_A <- prob[1:3] # Extracting only A_group
        prob_B <- prob[4:63] # Extracting only COD
        
        ## Determine Preg_State and Likelihood
        if(sum(prob_A) == 0 || reproductiveAge == 0){
            preg_state <- "Indet"
            lik.preg <- 0
        }
        if(which.max(prob_A) == 1 && prob_A[1] != 0 && reproductiveAge == 1){
            preg_state <- "nrp"
            lik.preg <- round(prob_A[1]/sum(prob_A)*100)
        }
        if(which.max(prob_A) == 2 && prob_A[2] != 0 && reproductiveAge == 1){
            preg_state <- "pr6w"
            lik.preg <- round(prob_A[2]/sum(prob_A)*100)
        }
        if(which.max(prob_A) == 3 && prob_A[3] != 0 && reproductiveAge == 1){
            preg_state <- "preg"
            lik.preg <- round(prob_A[3]/sum(prob_A)*100)
        }
        
        ## Calculate likelihood of marternal death
        lik_mat <- " "
        if(reproductiveAge == 1 && sum(prob_A) != 0) lik_mat <- round((prob_A[2]+prob_A[3])/sum(prob_A)*100)
        
        ## Normalize the probability of CODs
        if(sum(prob_B) != 0)  prob_B<-prob_B/sum(prob_B)
        prob.temp <- prob_B
        if(max(prob.temp) <= 0.4){
            indet <- "Indet"
            cause1<-lik1<-cause2<-lik2<-cause3<-lik3<-" "
        }
        ## Determine the output of InterVA
        if(max(prob.temp) > 0.4){
            ## Find max likelihood
            indet <- " "
            lik1 <- round(max(prob.temp)*100)
            cause1 <- names(prob.temp)[which.max(prob.temp)]
            ## Delete the max and find the second max
            prob.temp <- prob.temp[-which.max(prob.temp)]
            lik2 <- round(max(prob.temp)*100)
            cause2 <- names(prob.temp)[which.max(prob.temp)]
            ## Not show the second if it is too small
            if(max(prob.temp) < 0.5 * max(prob_B)) lik2 <- cause2 <- " "
            
            ## Delete the second max and find the third max
            prob.temp <- prob.temp[-which.max(prob.temp)]
            lik3 <- round(max(prob.temp)*100)
            cause3 <- names(prob.temp)[which.max(prob.temp)]
            ## Not show the third if it is too small
            if(max(prob.temp) < 0.5 * max(prob_B)) lik3 <- cause3 <- " "
        }
        ## Save the result as a list object
        ID.list[i] <- index.current
        VAresult[[i]] <- va(ID = index.current, MALPREV = Malaria, HIVPREV = HIV, PREGSTAT = preg_state, PREGLIK = lik.preg, PRMAT = lik_mat, INDET = indet, CAUSE1 = cause1, LIK1 = lik1, CAUSE2 =cause2, LIK2 = lik2, CAUSE3 = cause3, LIK3 = lik3, wholeprob = c(prob_A,prob_B))
        ## Determine the form of file saved
        if(output=="classic") save.va(VAresult[[i]],filename=filename, write)
        if(output=="extended") save.va.prob(VAresult[[i]],filename=filename, write)
    }
    setwd(globle.dir)
    out <- list(ID = ID.list[which(!is.na(ID.list))], VA = VAresult[which(!is.na(ID.list))], Malaria = Malaria, HIV = HIV)
    class(out) <- "interVA"
    return(out)
}



