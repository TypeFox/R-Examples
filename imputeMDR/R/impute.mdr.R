.PACKAGE <- "imputeMDR" 

impute.mdr <-
function(dataset, colresp, cs, combi, cv.fold = 10, na.method=0, max_iter=30, randomize = FALSE) {

      vec.matches<- function (vec,mat){
          # return line number of matrix where the values of the line matche with given vector 
          len <-length(vec)
          if (len==1){      # for a number and a vector input
              matched <- which(mat==as.character(vec))
          } else {
              seqs<- seq(nrow(mat))
              matched <- seqs
              for (i in 1: len){
                  val <- which(mat[,i]==as.character(vec[i]))
                  matched <- val[which(val %in% matched)]
              }
          }
          return( matched)
      }

      errRate2.C <- function(comb, train, test, threshold, na_method, loc_na, loc_target, target_num, max_iter) {
        z <- .C("err_rate", as.integer(comb), as.integer(nrow(comb)), as.integer(ncol(comb)),
                as.integer(unlist(train)), as.integer(nrow(train)), as.integer(ncol(train)),
                as.integer(unlist(test)), as.integer(nrow(test)),
                as.double(threshold), 
                err.train = double(ncol(comb)), 
                err.test = double(ncol(comb)), 
                as.integer(na_method) ,
                as.integer(loc_na-1), as.integer(loc_target-1), as.integer(target_num), length(loc_na), as.integer( max_iter),NAOK = FALSE,PACKAGE="imputeMDR")
        # na_method = 0:Complete, 1:Category, 2:Available, 3:EM
        cbind(train = z$err.train, test = z$err.test)
      }

    resp <- dataset[, colresp]
	
	### check data integrity ###
    # response values must be 0 or 1
    if(length(which(!unique(resp)%in%c(0,1))) > 0)
         stop("Input data error : Response values must be 0 or 1 !")
    
    if(length(which(is.na(dataset)))  >0)
          stop("Input data error : Missing values must be coded as 3 !")
   	
	vals<-function(x){	#count the number of values that is not one of 0,1,2, and 3
           length(which(!unique(x) %in% c(0:3)))
    }
       apply(dataset,1,vals)->a
   if (sum(a)>0)
            stop("Input data error : Only 0,1,2 and 3 are allowed in input dataset !")

    if(na.method== 0 && length(which(dataset==3)) > 0)
           stop("Selected missing data handling method are only allowed for complete dataset !")
    
            
    case <- which(resp == cs)
    ctl <- which(resp != cs)
    if (randomize) {
        case <- sample(case)
        ctl <- sample(ctl)
    }
    resp <- as.integer(resp == cs)
    snp <- dataset[, -colresp]
	snp.names <- colnames(dataset[, -colresp])

    cv.case <- matrix(0, nrow = length(case), ncol = 2)
    cv.case[, 1] <- 1:length(case) %% cv.fold + 1
    cv.case[, 2] <- case
    cv.ctl <- matrix(0, nrow = length(ctl), ncol = 2)
    cv.ctl[, 1] <- 1:length(ctl) %% cv.fold + 1
    cv.ctl[, 2] <- ctl
    cv <- rbind(cv.case, cv.ctl)
    d <- cbind(resp, snp)

    ## result
    z <- list()

    ## combinations
    comb <- list()
    k <-combi
    comb <- combn(ncol(snp), k)
    comb <- comb + 1
    z <- list()
    z$min.comb <- matrix(0, nrow = k, ncol = cv.fold)
    z$train.erate <- numeric(cv.fold)
    z$test.erate <- numeric(cv.fold)

      ## EM variable### 
      # loc_na : line with missing 
      # loc_target : line to fill out
      # target_num : number of lines to fill out
      int=combi
      tab <- expand.grid(rep(list(0:3),int))
      if(combi>1){
            loc_na <- which(apply(tab==3,1,any))
            #loc_na <- loc_na[-which(loc_na==which(apply(tab==3,1,all)))]
            for (i in  1:(length(loc_na)-1)){
                  cond <- which(tab[loc_na[i],]!=3)
                  if(i==1){
                        loc_target <- vec.matches(tab[loc_na[i],cond],tab[,cond])
                        loc_target <- loc_target[-which(loc_target%in%loc_na)]
                        target_num <- length(loc_target)
                  } else if(length(cond)>0){
                        loc_tmp <- vec.matches(tab[loc_na[i],cond],tab[,cond])
                        loc_tmp <- loc_tmp[-which(loc_tmp%in%loc_na)]
                        loc_target <- c(loc_target, loc_tmp)
                        target_num <- c(target_num,length(loc_tmp))
                  }
            }
            # cases of missing for all values
            loc_tmp<-seq(1:nrow(tab))[-loc_na]
            loc_target <- c( loc_target, loc_tmp)
            target_num <- c( target_num, length(loc_tmp))
      }else if(combi==1){
            loc_na<-which(tab==3)
            loc_target <- seq(1:nrow(tab))[-loc_na]
            target_num <-  length(loc_target)
      }


    for (i in 1:cv.fold) {
        train <- d[-cv[cv[, 1] == i, 2], ]
        test <- d[cv[cv[, 1] == i, 2], ]
        threshold <- length(which(train[, 1] == cs)) / 
        length(which(train[, 1] != cs)) 

      ## train and test error rate for all combination of SNPs
      errs <- errRate2.C(comb, train, test, threshold, na_method=na.method, loc_na, loc_target, target_num, max_iter)
        min.loc <- which.min(errs[, 1])
        z$min.comb[, i] <- comb[, min.loc] - 1
        z$train.erate[i] <- errs[min.loc, 1]
        z$test.erate[i] <- errs[min.loc, 2]
    }
	# 5.15.2010 : print SNP names of input data in the results
    # z$min.comb[z$min.comb >= colresp] <- z$min.comb[z$min.comb >= colresp] + 1
	z$min.comb <- matrix(snp.names[z$min.comb],combi,)

    ## maximum count of repeatedly selecting the same model across CV
  
   
    name <- apply(t(z$min.comb),1,function(x) {
        t<-NULL;
        for(i in 1:length(x)) 
            t<-paste(t,as.character(x[i]),sep=";") ;return(t)})
	z$best.combi<-unlist(strsplit(names(sort(table(name),decreasing=TRUE))[1],";"))[-1]
	
    cv.result <- cbind(t(z$min.comb),z$train.erate,z$test.erate)
    colnames(cv.result)<- c(paste("SNP",1:combi,sep=""),"train.err","test.err")
    return(list(cv.result=cv.result,best=z$best.combi))
}

.onLoad <- function(lib,pkg){
  library.dynam("imputeMDR",pkg,lib)
}
