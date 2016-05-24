##' print22 function
##'
##' A function to print details of the 2 by 2 table for use with the function twotwoinfo.
##'
##' @return prints the names of the arguments of twotwofunction info to screen in their correct place in the 2 by 2 table
##' @seealso \link{twotwoinfo}
##' @export


print22 <- function(){

    cat("            | Outcome 1 | Outcome 2 | TOTALS \n")
    cat("---------------------------------------------\n")
    cat("Exposed     |    e1     |    e2     |   et   \n")
    cat("Not Exposed |    u1     |    u2     |   ut   \n")
    cat("---------------------------------------------\n")
    cat("TOTALS      |    o1t    |    o2t    |   T    \n")
}


##' latexformat function
##'
##' A function to format text or numeric variables using scientific notation for LaTeX documents.
##'
##' @param x a numeric, or character
##' @param digits see ?format
##' @param scientific see ?format
##' @param ... other arguments to pass to the function format
##' @return ...
##' @export

latexformat <- function(x,digits=3,scientific=-3,...){
    if(is.character(x)){
        x <- as.numeric(x)
    }
    xtxt <- format(x,digits=digits,scientific=scientific,...)
    if(length(grep("e",xtxt))>0){
        spl <- unlist(strsplit(xtxt,"e"))
        xtxt <- paste(spl[1],"$\\times10^{",as.character(as.numeric(spl[2])),"}$",sep="")
    }    
    return(xtxt)
}

##' twotwoinfo function
##'
##' A function to compute and diplay information about 2 by 2 tables for copying into LaTeX documents. Computes odds ratios and relative 
##' risks together with confidence intervals for 2 by 2 table and prints to screen in LaTeX format. The funciton will try to fill in any 
##' missing values from the 2 by 2 table. Type print22() at the console to see what each argument refers to.
##'
##' @param e1 type print22() at the console 
##' @param u1 type print22() at the console
##' @param o1t type print22() at the console 
##' @param e2 type print22() at the console 
##' @param u2 type print22() at the console 
##' @param o2t type print22() at the console 
##' @param et type print22() at the console 
##' @param ut type print22() at the console 
##' @param T type print22() at the console 
##' @param lev significance level for confidence intervals. Default is 0.95 
##' @param LaTeX whether to print the 2 by 2 information as LaTeX text to the screen, including the table, odds ratio, relative risk and confidence intervals 
##' @param digits see ?format
##' @param scientific see ?format
##' @param ... other arguments passed to function format
##' @return Computes odds ratios and relative risks together with confidence intervals for 2 by 2 table and prints to screen in LaTeX format.
##' @seealso \link{print22}
##' @export

twotwoinfo <- function(e1=NA,u1=NA,o1t=NA,e2=NA,u2=NA,o2t=NA,et=NA,ut=NA,T=NA,lev=0.95,LaTeX=TRUE,digits=3,scientific=-3,...){
 
    M <- matrix(c(e1,u1,o1t,e2,u2,o2t,et,ut,T),3,3)
    rownames(M) <- c("Exposed","Not Exposed","TOTALS")
    colnames(M) <- c("Outcome 1","Outcome 2","TOTALS")
    
    for(i in 1:10){
        rowNA <- apply(M,1,function(x){sum(is.na(x))})
        if (any(rowNA==1)){
            idx <- which(rowNA==1)[1]
            vec <- M[idx,]
            vecidx <- which(is.na(vec))
            if(vecidx==1){
                vec[1] <- vec[3] - vec[2]
            }
            else if(vecidx==2){
                vec[2] <- vec[3] - vec[1]
            }
            else{
                vec[3] <- vec[1] + vec[2]
            }
            M[idx,] <- vec
        }
        
        colNA <- apply(M,2,function(x){sum(is.na(x))})
        if (any(colNA==1)){
            idx <- which(colNA==1)[1]
            vec <- M[,idx]
            vecidx <- which(is.na(vec))
            if(vecidx==1){
                vec[1] <- vec[3] - vec[2]
            }
            else if(vecidx==2){
                vec[2] <- vec[3] - vec[1]
            }
            else{
                vec[3] <- vec[1] + vec[2]
            }
            M[,idx] <- vec
        }
    }
    
    if (any(is.na(M))){
        print(M)
        stop("Cannot complete 2x2 table")
    }
    
    if(LaTeX){
        latextable(M,rownames=rownames(M),colnames=c("",colnames(M)),digits=10,scientific=0)
    }
    
    RR <- (M[1,1]/M[1,3])/(M[2,1]/M[2,3]) # rel risk of outcome 1 in the exposed group compared with the unexposed group
    logRR <- log(RR)
    
    OR <- (M[1,1]/M[1,2])/(M[2,1]/M[2,2]) # odds of outcome 1 in the exposed group compared with the unexposed group
    logOR <- log(OR)
    
    selogRR <- sqrt(1/M[1,1] + 1/M[2,1] - 1/M[1,3] - 1/M[2,3])
    selogOR <- sqrt(1/M[1,1] + 1/M[1,2] + 1/M[2,1] + 1/M[2,2])
    
    z <- round(qnorm(1-(1-lev)/2),2)
    
    ciRR <- exp(c(logRR-z*selogRR,logRR+z*selogRR))
    ciOR <- exp(c(logOR-z*selogOR,logOR+z*selogOR))
    
    if(LaTeX){  
        write("Relative risk/Odds ratio of Outcome 1 in the exposed compared with unexposed group:","")
    
        write("","")

    
        write(paste("$$\\text{Relative Risk} = \\frac{",M[1,1],"/",M[1,3],"}{",M[2,1],"/",M[2,3],"}=",latexformat(RR),"$$",sep=""),"")
        write(paste("$$\\text{Odds Ratio} = \\frac{",M[1,1],"/",M[1,2],"}{",M[2,1],"/",M[2,2],"}=",latexformat(OR),"$$",sep=""),"")          
 
        write("","")

        write(paste("$$\\se(\\text{log Relative Risk}) = \\sqrt{\\frac1{",M[1,1],"} + \\frac1{",M[2,1],"} - \\frac1{",M[1,3],"} - \\frac1{",M[2,3],"}}=",latexformat(selogRR),"$$",sep=""),"")
        write(paste("$$\\se(\\text{log Odds Ratio}) = \\sqrt{\\frac1{",M[1,1],"} + \\frac1{",M[1,2],"} + \\frac1{",M[2,1],"} + \\frac1{",M[2,2],"}}=",latexformat(selogOR),"$$",sep=""),"")
        
        write("","")
        
        write(paste("For the relative risk, a ",floor(100*lev),"\\% CI is $\\exp[\\log RR\\pm",z,"\\times\\se(\\log RR)]$.",sep=""),"")
        write(paste("For the odds ratio, a ",floor(100*lev),"\\% CI is $\\exp[\\log OR\\pm",z,"\\times\\se(\\log OR)]$.",sep=""),"")
        
        write("","")
        
        write(paste("Relative risk of disease in exposed group compared with unexposed group is ",latexformat(RR),", ",floor(100*lev),"\\% confidence interval ",latexformat(ciRR[1])," to ",latexformat(ciRR[2]),".",sep=""),"")
        write(paste("Odds of disease in exposed group compared with unexposed group is ",latexformat(OR),", ",floor(100*lev),"\\% confidence interval ",latexformat(ciOR[1])," to ",latexformat(ciOR[2]),".",sep=""),"")
    }

    return(list(RR=RR,selogRR=selogRR,ciRR=ciRR,OR=OR,selogOR=selogOR,ciOR=ciOR))    

}