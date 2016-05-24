#' @title Configural Frequencies Analysis for two Samples. 
#' @export S2CFA
#' @exportClass S2CFA
#' @description Calculates coefficients for the two-sample CFA. Instead of differentiating between types and antitypes, two-sample CFA looks for discrimination types, that is configurations with significant differences in frequencies between two subsamples.
#' @details no details at the moment ...
#' 
#' @param patternfreq an object of class \code{"Pfreq"}, which is data in pattern frequencies representation - see function \code{\link{dat2fre}}. The variable defining the two subsaples (a variable with max. two categories) must be located in the last but one column of the object of class \code{"Pfreq"} 
#' 
#' @param alpha a numeric giving the alpha level for testing (default set to \code{alpha=.05})
#' 
#' @param ccor a logical (TRUE / FALSE) determining wether to apply a continuity correction or not. When set to \code{ccor=TRUE} continuity correction is applied. For \code{ccor=FALSE} no continuity correction is applied.
#'   
#' @param ... additional parameters passed through to other functions.
#' @return an object of class \code{S2CFA} with results.
#' @references Stemmler, M. (2014). \emph{Person-Centered Methods – Configural Frequency Analysis (CFA) and Other Methods for the Analysis of Contingency Tables.} Cham Heidelberg New York Dordrecht London: Springer.
#' @references Stemmler, M., & Hammond, S. (1997). Configural frequency analysis of dependent samples for intra-patient treatment comparisons. \emph{Studia Psychologica, 39}, 167–175.
#' 
#' @examples #######################################
#' ############### some examples #########
#' ######### example from Marks Textbook
#' data(Lienert1978)
#' res1 <- S2CFA(Lienert1978)
#' summary(res1)
#' res2 <- S2CFA(Lienert1978, ccor=TRUE) # with continuity correction
#' summary(res2)
#' ######### example with biger numbers
#' data(suicide)
#' ftab(suicide) # 'Epoche' may divide the sample into 2 subsamples  
#' suicide_2s <- suicide[, c(1,3,2) ] # reorder data that 'Epoche' is the last column
#' ftab(suicide_2s) # check reordering
#' suicide_2s_fre <- dat2fre(suicide_2s)
#' res3 <- S2CFA(suicide_2s_fre)
#' summary(res3)
#' res4 <- S2CFA(suicide_2s_fre, ccor=TRUE) # with continuity correction
#' summary(res4)

############### start of function definition ##################
S2CFA<-function(patternfreq, alpha=.05, ccor=FALSE, ...){
if(any(class(patternfreq)=="Pfreq") != TRUE){stop("patternfreq must be an object of class 'Pfreq'","\n","see func. dat2fre()", call. = TRUE) }
#---------- Abfragen von steuerargumenten -----------------
# bisher keine
#---------- Aufbereiten von patternfreq -----------------
subsamples <- names(table(patternfreq[,(ncol(patternfreq)-1)]))
if(length(subsamples)!=2){stop("wrong definition of subsamples! - check argument 'patternfreq'")}
pattern_1 <- apply(patternfreq[ (patternfreq[,(ncol(patternfreq)-1)])==(subsamples[1])  ,1:(ncol(patternfreq)-2)], 1,paste, collapse=" ")
pattern_2 <- apply(patternfreq[ (patternfreq[,(ncol(patternfreq)-1)])==(subsamples[2])  ,1:(ncol(patternfreq)-2)], 1,paste, collapse=" ")
observed_1 <- (patternfreq[ (patternfreq[,(ncol(patternfreq)-1)])==(subsamples[1])  ,(ncol(patternfreq))] )
observed_2 <- (patternfreq[ (patternfreq[,(ncol(patternfreq)-1)])==(subsamples[2])  ,(ncol(patternfreq))] )
oplist <- mapply(FUN=function(o1,o2){a=o1;b=o2; c=sum(observed_1)-o1; d=sum(observed_2)-o2; A=a+b; B=c+d;C=a+c; D=b+d; n=A+B; data.frame(a,b,c,d,A,B,C,D,n)  },o1=observed_1 ,o2=observed_2, SIMPLIFY = FALSE)
names(oplist) <- pattern_1
# oplist
n <- sum(c(observed_1,observed_2))
with.out<-function(x, y) {x[!x %in% y]} # hilfsfunktion
#----------------- calculate expected counts ------------------
expected_1 <- unlist(lapply(oplist, function(x){(x["A"]*x["C"])/x["n"]}))
expected_2 <- unlist(lapply(oplist, function(x){(x["A"]*x["D"])/x["n"]}))
names(expected_1) <- pattern_1
names(expected_2) <- pattern_2
#----------------- calculate exact fisher Test ------------------
#### factorial mit kürzen und package{gmp}  weniger limitiert
oplist_1 <- lapply(oplist, function(x){ list(oben= tabulate(sort(with.out(c((1:x[1,"A"]),(1:x[1,"B"]),(1:x[1,"C"]),(1:x[1,"D"])),0)),n)      , unten= tabulate(sort(with.out(c((1:x[1,"n"]),(1:x[1,"a"]),(1:x[1,"b"]),(1:x[1,"c"]),(1:x[1,"d"])),0)),n) ) } )
# oplist_1
#### kürzen
oplist_2 <- lapply(oplist_1, function(x){x$oben - x$unten}) #  alle mit "-" stehen dann unten
# oplist_2
oplist_3_ind <- lapply(oplist_2, function(x){list(ind_oben = x>0, ind_unten = x<0)}) # index was bleibt oben und unten
# oplist_3_ind
oplist_3 <- lapply(oplist_2, function(x){xv <- abs(x); names(xv) <- 1:n; xv}) # absolut werte haufigkeiten der produkttherme oben und unten
# oplist_3
oplist_4 <- mapply(function(ind, wert){ list(oben=wert[ind$ind_oben],unten=wert[ind$ind_unten]) }, ind=oplist_3_ind, wert=oplist_3, SIMPLIFY = FALSE) # names sind die Produktherme; vector ist die anzahl des jew. Produktherms 
# oplist_4
###### berechnungen mit gekürzter Produktherm liste 'oplist_4'
# library(gmp) # für funkktion 'as.bigz'
oplist_5 <- lapply(oplist_4, function(x){ list(oben=as.bigz(rep(x=as.numeric(names(x$oben)), times=x$oben)), unten=as.bigz(rep(x=as.numeric(names(x$unten)), times=x$unten))      )} )
# oplist_5
oplist_6 <- lapply(oplist_5, function(x){ list(oben= prod(x$oben), unten= prod(x$unten)) } )
# oplist_6
fisher_res <- unlist(lapply(oplist_6, function(x){as.double(div.bigq (x$oben , x$unten)) }))
# fisher_res

#----------------- calculate traditional Chi-square Test ------------------
# oplist
if(ccor==FALSE){
  # x <- oplist[[5]]
oplist_Chi_1 <- lapply(oplist, function(x){ list(oben=(prod(as.bigz(c(x[1,"n"], ((x[1,"a"]*x[1,"d"]) - (x[1,"b"]*x[1,"c"])) , ((x[1,"a"]*x[1,"d"]) - (x[1,"b"]*x[1,"c"])))))), unten=(prod(as.bigz(c(x[1,"A"], x[1,"B"], x[1,"C"], x[1,"D"]))))  ) } )
# x <- oplist_Chi_1[[5]]
  Chi <- unlist(lapply(oplist_Chi_1, function(x){ if(as.double(x$unten)!=0){as.double(div.bigq( x$oben, x$unten))}else{NaN} })) 
}
if(ccor==TRUE){
  # x <- oplist[[2]]
  # 36* (abs(3*4-14*15)  -   36 * 0.5)    ^2  / (17*19*18*18)
  # Chi <- unlist(lapply(oplist, function(x){ ( x[1,"n"] * (abs((x[1,"a"]*x[1,"d"]) - (x[1,"b"]*x[1,"c"])) - (x[1,"n"]) * 0.5)^2 )  / (x[1,"A"]*x[1,"B"]*x[1,"C"]*x[1,"D"])   } )) #  as.bigz fehlt noch
  oplist_Chi_1_cont <- lapply(oplist, function(x){ list(oben=(prod(as.bigz(c(x[1,"n"],(abs((x[1,"a"]*x[1,"d"])-(x[1,"b"]*x[1,"c"]))-x[1,"n"]*0.5),(abs((x[1,"a"]*x[1,"d"])-(x[1,"b"]*x[1,"c"]))-x[1,"n"]*0.5))))) , unten=(prod(as.bigz(c(x[1,"A"], x[1,"B"], x[1,"C"], x[1,"D"]))))  ) } )
  Chi <- unlist(lapply(oplist_Chi_1_cont, function(x){ if(as.double(x$unten)!=0){as.double(div.bigq( x$oben, x$unten))}else{NaN} }))  
}

df <- rep(1, times=length(Chi))

pChi <- (1-pchisq(Chi,df))

bonferroni <- alpha/length(expected_1)

erg <- data.frame(pattern_1, expected_1, observed_1, pattern_2, expected_2, observed_2, ex.fisher.test=fisher_res, Chi=Chi, df=df, pChi=pChi)
names(erg)[1:6] <- paste(c( rep(subsamples[1],3),rep(subsamples[2],3) ), rep(c("pat.","exp.","obs."),2)    ,sep=".")

result <- list( local.test = erg, bonferroni.alpha=bonferroni, global.test = NULL) 

class(result)<-c("S2CFA","list")

return(result)
}
# End of function ------ 