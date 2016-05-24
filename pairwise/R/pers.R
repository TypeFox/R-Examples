#' @title WLE - Rasch Person Parameter
#' @export pers
#' @description This is the (new) main function for calculation of person estimates based on answering dichotomous or polytomous items according theRasch Model (Rasch, 1960) and Partial Credit Model (Masters, 1982), given the item parameters (object of class \code{"pair"} - as a result of \code{\link{pair}()}) and and the datamatrix (argument \code{daten}) containing the person respose vectors (rows), using an WL approach, introduced by Warm (1989).
#' @details no detail in the moment.
#' @param itempar The item parameter prior calculated or estimated. A list object of class \code{"pair"} as a result of applying the function \code{\link{pair}()} to the data. Or an 'ordinary' \code{"matrix"} with \code{nrow = k} (number of items) and \code{ncol = m} (maximum number of thresholds), holding  the 'thurstonian' thresholds of the respective item. Some matrix entries may be \code{NA}, depending on the number of categories of the respective item.
#' @param daten A \code{"matrix"} (or \code{"data.frame"}) optionaly with named colums (names of items) and named rows (person IDs). This argument can be left empty when the argument itempar (above) is of class \code{"pair"}. \code{daten} holds polytomous or dichotomous (or mixted category numbers) responses of \code{n} respondents (rows) on \code{k} items (colums) coded starting with 0 for lowest category to \code{m-1} for highest category, with m beeing a vector (with length k) with the number of categories for the respective item. Responses in \code{daten} must be stored as \code{"integers"} (not \code{"factors"} !) and may have missing values.
#' @param incidenz This argument is only relevant when items are assigned to different booklets. For such a booklet-design a \code{"matrix"} should be assigned to this argument, with the same dimensions like \code{daten}, containig 0 and 1 integer codes, giving the information (for every person) if the respective item was in the respective booklet (coded 1) given to the person or not (coded 0).    
#' @param na_treat optionaly an integer (vector) defining the type of treatment to missing responses in the argument \code{daten}. If set to \code{na_treat=NULL} (default) missing responses are treated as missings and the respective person is assigned to an corresponding missing group for estimation. An option is to set \code{na_treat} to any integer value between 0 (lowest category) and the numeric code for the maximum ctaegory of the respective item.
#' @param limit numeric giving the limit at which accuracy the WL-algorithm stops.
#' @param iter numeric giving the maximum numer of iteration to perform.
#' @param Nrel logical with default set to \code{Nrel=FALSE} to include persons with perfect response vectors for calculating WLE reliability. If set to \code{Nrel=TRUE} persons with perfect response vectors are excluded for calculating WLE reliability.
#' @param tecout logical default set to \code{FALSE}. If set to \code{TRUE} the result will be a (very) long list with estimation details for every case in \code{daten}. In case of a booklet-design the list entries will be divided by "booklet".
#' @return An object of class \code{c("pers", "data.frame")} or a (very long) \code{"list"} (when setting on \code{techout=TRUE}) containing the person parameters.
#' @exportClass pers
#' @references Masters, G. (1982). A Rasch model for partial credit scoring. \emph{Psychometrika, 47}(2), 149–174.
#' @references Rasch, G. (1960). \emph{Probabilistic models for some intelligence and attainment tests.} Copenhagen: Danmarks pædagogiske Institut.
#' @references Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory. \emph{Psychometrika, 54}(3), 427–450.
#' @examples ############
#' data(sim200x3)
#' result <- pers(itempar=pair(sim200x3))
#' summary(result)
#' plot(result)
#' logLik(result) # Log-Likelihood for 'estimated' model
#' logLik(result, sat=TRUE) # Log-Likelihood for saturated model
#' AIC(logLik(result)) # AIC for 'estimated' model
#' AIC(logLik(result, sat=TRUE)) # AIC for saturated model
#' BIC(logLik(result)) # BIC for 'estimated' model
#' BIC(logLik(result, sat=TRUE)) # BIC for saturated model
#' ###### following example requires package eRm ######
#' # require(eRm)
#' # # itemparameter with eRm:
#' # itempar_eRm <- thresholds(PCM(sim200x3))$ threshtable[[1]][,2:3]
#' # # pairwise personparameter with eRm-itemparameter and data:
#' # summary(pers(itempar=itempar_eRm,daten=sim200x3))
#' # # eRm personparameter:
#' # person.parameter(PCM(sim200x3))
#' # # personparameter with pairwise:
#' # summary(pers(pair(sim200x3))) 

pers<-function(itempar, daten=NULL, incidenz=NULL,na_treat=NULL,limit=0.00001,iter=50,Nrel=FALSE,tecout=FALSE){
fuargs <- list(itempar=itempar, daten=daten, incidenz=incidenz, na_treat=na_treat, limit=limit, iter=iter, Nrel=Nrel, tecout=tecout)

#  daten=NULL; incidenz=NULL;na_treat=NULL;limit=0.00001;iter=50;tecout=FALSE
##################################################################
# na_treat=NULL;limit=0.00001;iter=50;tecout=FALSE
# m wird aus der struktur des arguments itempar berechnet
#### needs internal functions: ------
# dataprep1<-function(X)
# missing_group<-function(X,s=NULL,all=FALSE)
# GewLL = function(np,sb,fw,m,  k=length(m),mscs=sum(m-1))
# PersPar = function(np,sb,m, imax=30,limit=0.0001, k=length(m),mscs=sum(m-1) )

##################################################################
#################### Begin der funktion  #########################
##################################################################

###################### check der Argumente #######################
##### infos aus itempar class matrix---------------------
if(any(class(itempar)=="matrix")){
  if(length(daten)==0){stop("no data assigned to argument daten")}
  
  threshold <-lapply(1:nrow(itempar), function(i) {na.omit(itempar[i,])})

  ### berechnung der sb; sb als liste!!!!
  sb<-(lapply(threshold,function(x){cumsum(x)})) # richtig: umrechnung threshold (tau) in sb !!!!!!!!! OK
  sb_return <- sb
  m <- apply(itempar,1,function(x){length(na.omit(x))+1}) # OK
  k <- length(m)
  if(length(unique(m))!=1){
    maxLen <- max(sapply(sb, length))
    # create a new list with elements padded out with 0s
    newsb <- lapply(sb, function(.ele){c(.ele, rep(0, maxLen))[1:maxLen]})
    sb <- do.call(rbind, newsb)
  }
  if(length(unique(m))==1){
    sb<-do.call(rbind, sb)
  }
  ### some minor consitency checks between daten and itempar
  if(dim(daten)[2]!=k){stop("missmatch: number of items in daten and itempar")}
  # aufbereiten für ausgabe mit den separat eingegbenen daten aus argument daten
  itempar_pair <- list(threshold=itempar,sigma=rowMeans(itempar,na.rm=TRUE),sb=sb_return ,resp=dataprep1(daten)) 
  class(itempar_pair) <- c("pair", "list")
}

##### infos aus itempar: class(itempar)=="pair" & daten != NULL ---------------------
if(   (any(class(itempar)=="pair")) & (length(daten)!=0)    ){
  if(dim(daten)[2]!=length(itempar$sb)) stop("missmatch of number of items in data and itempar")
  sb <- itempar$sb
  m <- sapply(sb,length)+1
  k <- length(m)
  if(length(unique(m))!=1){
    maxLen <- max(sapply(sb, length))
    # create a new list with elements padded out with 0s
    newsb <- lapply(sb, function(.ele){c(.ele, rep(0, maxLen))[1:maxLen]})
    sb <- do.call(rbind, newsb)
  }
  if(length(unique(m))==1){
    sb<-do.call(rbind, sb)
  }
  # aufbereiten für ausgabe mit den separat eingegbenen daten aus argument daten
  if(any(colnames(daten) != names(sb))){stop("mismatch: item names in daten and itempar")}
  itempar_pair <- list(threshold=itempar[["threshold"]],sigma=itempar[["sigma"]],sb=itempar[["sb"]] ,resp=dataprep1(daten)) 
  class(itempar_pair) <- c("pair", "list")
}

##### infos aus itempar:  class(itempar)=="pair" & daten == NULL---------------------
if(   (any(class(itempar)=="pair")) & (length(daten)==0)    ){
  daten <- itempar$resp
  sb <- itempar$sb
  m <- sapply(sb,length)+1
  k <- length(m)
  if(length(unique(m))!=1){
    maxLen <- max(sapply(sb, length))
    # create a new list with elements padded out with 0s
    newsb <- lapply(sb, function(.ele){c(.ele, rep(0, maxLen))[1:maxLen]})
    sb <- do.call(rbind, newsb)
  }
  if(length(unique(m))==1){
    sb<-do.call(rbind, sb)
  }
  itempar_pair <- itempar # übernahme wie im argument eingegeben
} 


######################################################################################
##### sortierbare zeilennamen für daten; item namen neu nur wenn fehlend -------------
daten<-dataprep1(daten)

##### some category checks with m and k ---------------------------------------------
#if (all(apply(daten,2,function(x){min(x,na.rm=TRUE)})== 0) == FALSE){stop("item categories must start with 0") }
#if(length(m)==0){m<-apply(daten,2,function(x){max(x,na.rm=TRUE)+1})}
#if(length(m)==1){m<-rep(m,dim(daten)[2])} 
#if(any (m < apply(daten,2,function(x){max(x,na.rm=TRUE)+1}))){stop("some items in data have more categories than defined in m","\n","max item categories in data are: ",paste(apply(daten,2,function(x){max(x,na.rm=TRUE)+1}),collapse=" , "),"\n", "but m was defined: ",paste(m,collapse=" , ")  )}
#if(dim(daten)[2]!=k){stop("number of items dose not match with argument itempar") }
  

##### herstellen bzw. übergabe der incidenz matrix -----------------------------------
if(is(incidenz, "NULL")){incidenz <- matrix(1,nrow = dim(daten)[1],ncol = dim(daten)[2])}
if(!is(incidenz, "NULL")){
  incidenz <- incidenz
  all_incpat<-(do.call("paste",c(as.data.frame(incidenz), sep = ""))) # zeilen der incidenz als pattern
  uni_incpat<-unique(all_incpat)# unique davon
}

##### check der incidenz matrix - bookletdesign --------------------------------------
uni_nullbook<-which(rowSums(do.call(rbind,lapply(strsplit(uni_incpat,"",fixed = TRUE),as.numeric)))==0) # ident. von personen mit booklets ohne ein Item
if(length(uni_nullbook)!=0){stop("some persons got booklets with no items and must be removed ", "\n", "-> check your design!")
# folgendes ist erst mal auskommentiert:
# wenn die Incidenz matrix zeilen vectoren mit sum =0 enthält stimmt eh was nicht.      
#       uni_incpat <- uni_incpat[-uni_nullbook]
#       all_nullbook <- which(rowSums(do.call(rbind,lapply(strsplit(all_incpat,"",fixed = TRUE),as.numeric)))==0)
#       nobP <- row.names(d)[all_nullbook]
#       # test <- d[all_nullbook,]
#       # test1 <- PBA_COG_S[all_nullbook,] 
#       cat(length(nobP), "Persons got booklets with no items and were removed for estimation", "\n", "-> check your design")
    }

### aufteilen der daten nach "booklets" in daten_booklet_list mit dem entspr. itempars ----
# ist eine Liste der Länge eins wenn kein bookletdesign.
    daten_booklet_list<-list()
    for (i in 1:length(uni_incpat)){
      item_index<-which(as.numeric(unlist(strsplit(uni_incpat[i],"",fixed = TRUE)))==1)
      person_index<-which(all_incpat == uni_incpat[i]) # all_incpat %in% uni_incpat[i]
      sb_tmp<-as.matrix(sb[item_index, ]) # ergänzung (daten nur 1 Item)
      if(  (length(person_index)>1) & (length(item_index)>1)   ){
        dat_tmp<-as.matrix(daten[person_index,item_index])}# ergänzung (daten nur 1 Item)
      if(length(person_index)==1){
        dat_tmp<-as.data.frame(t(as.matrix(daten[person_index,item_index])))# ergänzung (daten nur 1 Item)
        rownames(dat_tmp)<-rownames(daten)[person_index] 
      }
      if(length(item_index)==1){
        dat_tmp<-as.data.frame((as.matrix(daten[person_index,item_index])))# ergänzung (daten nur 1 Item)
        colnames(dat_tmp)<-colnames(daten)[item_index] 
      }
      ##### anwenden von na_treat auf daten_booklet_list
      if (length(na_treat)!=0){dat_tmp[is.na(dat_tmp)]<-na_treat}
      #if(dim(dat_tmp)[2]==1){colnames(dat_tmp)<-names(sb_tmp) }# ergänzung (daten nur 1 Item)
      m_tmp<-m[which( (unlist(strsplit(uni_incpat[i],"",fixed = TRUE)))=="1")]
      # print(dim(dat_tmp));print(dim(sb_tmp)); print(length(m_tmp))
      daten_booklet_list[[i]] <- list("sb"=sb_tmp ,"daten"=dat_tmp, "m"= m_tmp)# änderung (daten nur 1 Item)
    }
    names(daten_booklet_list)<-uni_incpat
# ENDE aufteilen der daten nach "booklets" in liste
  # lapply(daten_booklet_list, function(x){ dim(x$daten)})  # Check   
  # lapply(daten_booklet_list, function(x){ ftab(x$daten)}) # Check
  # sink(file="test.txt");lapply(daten_booklet_list, function(x){ print(x$daten)}) ;sink() # Check

#daten_booklet_list ok; na_treat ok

#### hier mit schleife über alle "booklets" -------------------------------------------------
ergL<-list()
for (i in 1:length(daten_booklet_list)){
  sbi <- daten_booklet_list[[i]]$sb
  di <- daten_booklet_list[[i]]$daten
  mi <- daten_booklet_list[[i]]$m
  
  X <- missing_group(as.data.frame(di))
  X_gr <- X$Pin
  d_gr <- X_gr$X.nonmis.group # daten ohne die missing !!!!!!!!!!!!!!
  m_gr <- lapply(X_gr$use.index,function(x){mi[x]}) # Anz. kategorie vectoren
  #print(m_gr)
  k_gr <- lapply(m_gr,length) # Anz. items integers
  #print(k_gr)
  mscs_gr <- lapply(m_gr,function(x){sum(x-1)}) # Max. Scoresumme integers
  #print(mscs_gr)
  np_gr <- mapply(function(x,y){tabulate(rowSums(x)+1,nbins=sum(y)+1)},d_gr, mscs_gr, SIMPLIFY=FALSE)
  #print(np_gr)
  sb_gr <- mapply(function(x){as.matrix(sbi[x,])},X_gr$use.index, SIMPLIFY=FALSE)
  
  est.details <- mapply(PersPar,np_gr,sb_gr,m_gr,iter,limit,SIMPLIFY=FALSE)
  scoL <- lapply(d_gr,rowSums)
  e1 <- mapply(function(x,y){(data.frame( PERS=names(x),raw=x, WLE=y$fw[x+1], SE.WLE=y$se.fw[x+1], ITER=rep(y$itmax,length(x)), WLL=rep(y$wLL,length(x)) ,row.names=names(x),stringsAsFactors=FALSE))},scoL, est.details, SIMPLIFY = FALSE)
  
  nag <- as.vector(unlist(mapply(rep, names(e1), sapply(e1,function(x){dim(x)[1]}))))# as.vector ergänzt
  e2 <- do.call(rbind,e1)
  e3 <- data.frame(NA.group=nag, e2) ;row.names(e3)<-NULL
  
  ## personen mit WLE = NA:
  PeNA <- X$Pout$person.index
  if(length(PeNA)!=0){ # weglassen falls länge null 
  mi_PeNA <- rep(names(PeNA),sapply(PeNA,length))
  wh_PeNA <- unlist(lapply(PeNA,names))
  e4 <- data.frame(NA.group=mi_PeNA, PERS=wh_PeNA,raw=NA, WLE=NA, SE.WLE=NA, ITER=NA, WLL=NA) ;row.names(e4)<-NULL
  ## ende personen mit WLE = NA:
  e5 <- rbind(e3,e4)
  ergL[[i]] <- list(est.result=e5, est.details=est.details) 
  }
  if(length(PeNA)==0){
    ergL[[i]] <- list(est.result=e3, est.details=est.details)  
  }
  
}
names(ergL) <- uni_incpat

############### AUSGABE ---------------------------
if(tecout==TRUE){
#   if(exists("nobP")==TRUE){
#     result <- list(ergL, nobP)  
#   }
 # else {result <- list(ergL)} 
{result <- list(ergL)} 
}
#------------------------------------
if(tecout==FALSE){
  
  ee1<-lapply(ergL,function(x){ x[["est.result"]]} )
  WLE <-  unlist(lapply(ee1,function(x){ x[["WLE"]]} ),T)
  PERS <-  unlist(lapply(ee1,function(x){ x[["PERS"]]} ),T)
  raw <-  unlist(lapply(ee1,function(x){ x[["raw"]]} ),T)
  SE.WLE <-  unlist(lapply(ee1,function(x){ x[["SE.WLE"]]} ),T)
  ITER <-  unlist(lapply(ee1,function(x){ x[["ITER"]]} ),T)
  NA.group <-  unlist(lapply(ee1,function(x){ x[["NA.group"]]} ),T)
  WLL <- unlist(lapply(ee1,function(x){ x[["WLL"]]} ),T) #added 24.11.2015

  res1 <- data.frame(persID=PERS,NA.group,raw,WLE,SE.WLE,ITER,WLL,row.names = PERS)
  res2 <- res1[order(res1$persID) ,]
  res3 <- data.frame(book=all_incpat,res2)
   
  
  # WLE Reliability rost 2004 p 381 see also http://www.rasch.org/erp7.htm
  reldat1 <- data.frame(WLE=res3$WLE,SE.WLE=res3$SE.WLE)
  N=dim(reldat1)[1]
  if (Nrel==TRUE){
    reldat1 <- reldat1[!(res3$raw==sum(m-1) | res3$raw==0),  ]
    N2=dim(reldat1)[1]
  } else{N2=N}  
  
  reldat2 <- reldat1[complete.cases(reldat1),]
  # dim(reldat1) 
  
  r.WLE.rel <- var(reldat2$WLE) / (mean((reldat2$SE.WLE)^2) + var(reldat2$WLE))
  
  n.WLE.rel <- dim(reldat2)[1]
  WLE.rel <- list(r.WLE.rel=r.WLE.rel,n.WLE.rel=n.WLE.rel,N.perf=N-N2)
  # END WLE Reliability  
  res3 <- as.list(res3) # changed 27-3-2015 --> pers is now a list
  res3$resp <- daten_booklet_list # changed 27-3-2015 --> pers is now a list
  result <- list(pers=res3 , pair=itempar_pair, WLE.rel=WLE.rel, fuargs=fuargs)
  
  class(result) <- c("pers","list")
  
}

  return(result)
} 
