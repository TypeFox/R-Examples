#' @export andersentest.pers
#' @title Andersen's Likelihood Ratio Test for Object of class "pers"
#' @description The Andersen likelihood ratio test is based on splitting the dataset into subgroups of persons. One can argue that it is a significance testable version of the more descriptive graphical model check - see \code{\link{grm}}.
#' @details Andersen (1973) proposed to split the dataset by [raw] score groups, which can be achieved setting the argument \code{split = "score"}. However as pointed out by Rost (2004) there might be several different splitting criteria for testing subsample invariance of the raschmodel. Thus the argument \code{split} provides some other options for splitting the data - see description of arguments.   
#' 
#' @param pers_obj an object of class\code{"pers"} - see function \code{\link{pers}}.
#' 
#' @param split Specifies the splitting criterion. Basically there are three different options available - each with several modes - which are controlled by passing the corresponding character expression to the argument. 
#' 
#' 1) Using the rawscore for splitting into subsamples with the following modes: \code{split = "median"} median raw score split - high score group and low score group; \code{split = "mean"} mean raw score split - high score group and low score group.
#' Finaly \code{split = "score"} that is splitting \code{daten} into as many subsamples as there are raw score groups - discarding min and max (theoretical) score group - which matches the concept proposed by Andersen (1973).
#' 
#' 2) Dividing the persons in \code{daten} into subsamples with equal size by random allocation with the following modes: \code{split = "random"} (which is equivalent to \code{split = "random.2"}) divides persons into two subsamples with equal size. In general the number of desired subsamples must be expressed after the dot in the character expression - e.g. \code{split = "random.6"} divides persons into 6 subsamples (with equal size) by random allocation etc. 
#' 
#' 3) The third option is using a manifest variable as a splitting criterion. In this case a vector with the same length as number of cases in \code{daten} must be passed to the argument grouping the data into subsamples. This vector should be coded as \code{"factor"} or a \code{"numeric"} integer vector with min = 1.
#' 
#' @param splitseed numeric, used for \code{set.seed(splitseed)} for random splitting - see argument \code{split}.
#' 
#' @param pot optional argument, at default (\code{pot=NULL}) setting is read from \code{pers_obj} - see description for \code{\link{pair}}.
#' @param zerocor optional argument, at default (\code{zerocor=NULL}) setting is read from \code{pers_obj} - see description for \code{\link{pair}}.
#' 
#' pot=pers_obj$pair$fuargs$pot, zerocor=pers_obj$pair$fuargs$zerocor
#' 
#' @return A (list) object of class \code{"andersentest.pers"} ...
#' @exportClass andersentest.pers
#'            
#'@references Andersen, E. B. (1973). A goodness of fit test for the rasch model. \emph{Psychometrika, 38}(1), 123–140. 
#'@references Rost, J. (2004). \emph{Lehrbuch Testtheorie - Testkonstruktion} (2 nd Ed.) Huber: Bern.
#'
#' @examples data(bfiN) # loading example data set
#' 
#' data(bfi_cov) # loading covariates to bfiN data set
#' model <- pers(pair(bfiN,m=6))
#' andersentest.pers(model, split = bfi_cov$gender)
#' andersentest.pers(model, split = "random")
#' andersentest.pers(model, split = "median")
#' ### unsing simulated data:
#' data("sim200x3")
#' model2 <- pers(pair(sim200x3))
#' andersentest.pers(model2, split = "median")

############## funktions beginn ########################################################

andersentest.pers<-function(pers_obj, split="median", splitseed="no", pot=NULL, zerocor=NULL){
  if(!is(pot, "NULL")){pers_obj$pair$fuargs$pot <- pot}
  if(!is(zerocor, "NULL")){pers_obj$pair$fuargs$pot <- pot}
  ### auslesen von pers_obj
  daten <- pers_obj$pair$resp  
  
  #### abfragen der teilungskriterien und teiler vorbereiten
    teil <- split  # übergabe an internes argument
    if(!(length(teil) > 1)) {  
      if(teil=="no"){
        teiler<-rep(1,dim(daten)[1])
        #OK
      }
      if(teil=="random"){
        if (class(splitseed)=="numeric"){set.seed(splitseed)}
        teiler<-as.numeric(cut(sample(1:(dim(daten)[1])),2))
        #OK
      }
      if(nchar(teil)>6){
        nteil<-as.numeric(unlist(strsplit(teil,".",TRUE))[2]) 
        if (class(splitseed)=="numeric"){set.seed(splitseed)}
        teiler<-as.numeric(cut(sample(1:(dim(daten)[1])),nteil))
        #OK
      }     
      if(teil=="mean"){
        daten<-as.matrix(daten)
        rscore<-rowSums(daten,na.rm = TRUE)
        teiler<-factor(ifelse(rscore > round(mean(rscore)),"above mean" ,"mean and below" ))
        #OK
      }
      if(teil=="median"){
        daten<-as.matrix(daten)
        rscore<-rowSums(daten,na.rm = TRUE)
        teiler<-factor(ifelse(rscore > median(rscore),"above median" ,"median and below" ))
        #OK
      }
      if(teil=="score"){
        daten<-as.matrix(daten)
        rscore<-rowSums(daten,na.rm = TRUE)
        rscore_ne <- rscore
        rscore_ne[rscore_ne==0] <- 1
        maxscore <- sum(apply(pers_obj$pair$threshold,1,function(x){length(na.omit(x))})) # OK
        rscore_ne[rscore_ne==maxscore] <-maxscore-1
#         nam_scoregrps <- as.numeric(names(table(rscore)))
#         fre_scoregrps <- as.vector(table(rscore))
#         maxscore <- sum(apply(pers_obj$pair$threshold,1,function(x){length(na.omit(x))})) # OK
#         minscore <-  0
#         out <- c(which(nam_scoregrps==minscore),which(nam_scoregrps==maxscore) )
#         if(length(out)>0){empscrg <- nam_scoregrps[-out]}else{empscrg <- nam_scoregrps}
#         empscrg
#         
        teiler<-factor(rscore_ne)
        #OK
      }
      
    }
    
    if((class(teil)=="integer") | (class(teil)=="numeric") | (class(teil)=="factor")){
      #teiler<-daten[,teil]
      if( (dim(daten)[1])!=length(teil) ){stop("length of argument 'split' dose not match with 'data'")}
      teiler<-teil
      #if (class(teiler)=="factor"){teiler<-(as.numeric(teiler))}
      #if (min(teiler!=1)){stop("argument teil is not valid specified")}
      # daten<-daten[,-teil]
      #OK
    }
    #### ENDE abfragen der teilungskriterien und teiler vorbereiten  
    # vorbereiten des objektes datalist anhand des vectors teiler
    subsamp <- names(table(teiler))
    
    datalist<-vector("list",length=length(subsamp)) #vorber. leere datalist   
    for (i in 1:length(datalist)){
      datalist[[i]]<-daten[which(teiler==subsamp[i]),]  #hier die zuordnung der subsamples aus daten
    }
    names(datalist) <- paste(subsamp,"sample")
     
    L_pair_erg <- lapply(datalist, pair, m=pers_obj$pair$m, pot=pers_obj$pair$fuargs$pot, zerocor=pers_obj$pair$fuargs$zerocor) 
  
  if(!is(pers_obj$fuargs$incidenz,"NULL")){ # wenn eine incidenz matrix übergeben wurde muss auch die entsprechend geteilt werden
    incilist<-vector("list",length=length(subsamp)) #vorber. leere incidenz list
    for (i in 1:length(incilist)){
      incilist[[i]]<-pers_obj$fuargs$incidenz[which(teiler==subsamp[i]),]  #hier die zuordnung der subsamples aus daten
    }
    names(incilist) <- paste(subsamp,"sample")
    
    L_pers_erg <- mapply(FUN=function(x,y,...){pers(itempar=x, incidenz=y) }, x=L_pair_erg, y=incilist, MoreArgs=(pers_obj$fuargs[c("daten","na_treat","limit","iter","Nrel","tecout")]) , SIMPLIFY = FALSE)
    
  }else{L_pers_erg <- lapply(L_pair_erg, pers, daten=pers_obj$fuargs$daten, incidenz=pers_obj$fuargs$incidenz, na_treat=pers_obj$fuargs$na_treat, limit=pers_obj$fuargs$limit, iter=pers_obj$fuargs$iter, Nrel=pers_obj$fuargs$Nrel, tecout=pers_obj$fuargs$tecout)}
  
  logLik_subsamples  <- sapply(L_pers_erg,function(x){logLik.pers(x)[[1]]})
  df_subsamples <- sapply(L_pers_erg,function(x){attr(x = logLik.pers(x),which = "df")})
  # df_subsamples <- sapply(logLik_subsamples, attr, which = "df")
  logLik_totalsamples  <- logLik.pers(pers_obj)[[1]]
  df_totalsamples <- attr(x=logLik.pers(pers_obj),which = "df")
  
    #chi_test <- 2*((sum(unlist(logLik_subsamples))) - logLik_totalsamples[1])
    chi_test <-max( 0,(-2*( (logLik_totalsamples[1]) - (sum(unlist(logLik_subsamples))) )))# new 21.11.2015
    df_test <- (sum(df_subsamples) - (df_totalsamples)  ) 
    p_test <- 1-pchisq(q=chi_test, df=df_test)
  result <- list("Chi^2"=chi_test, df = df_test, p = p_test)
class(result) <- c("andersentest.pers","list")
  return(result) 
}  
