#    Do not delete!
#  File name   	HotDeckImputation.R
#  Part of:		   HotDeckImputation (GNU R contributed package)
#  Author:			Dieter William Joenssen
#  Copyright:		Dieter William Joenssen
#  Email:			Dieter.Joenssen@googlemail.com
#  Created:		   14 May 2013
#  Last Update: 	21 October 2014
#  Description:	R code for Package HotDeckImputation. Implemented functionions include following:
#                 Topics:
#                 -impute.NN_HD           ~ An extremely comprehensive implementation Nearest-Neighbor hot deck algorithms
#                 -impute.mean            ~ Imputes the mean value of the complete cases in a variable for the missing values
#                 -reweight.data          ~ Reweights data with given wieghts, preprocessing for distance calculation
#				          -match.d_r_vam          ~ performs the Vogel's approximation method matching for donors-recipients
#                 -match.d_r_odd          ~ performs the ODD matching for donors-recipients
#                 -impute.fast_NN_HD      ~ performs fast version of nearest neighbor hot deck, not precomputing the distance matrix
#                 Hidden functions:
#                 -.create.weightsvector  ~ create the weights vector by interpreting the user input from impute.NN_HD
#                 -.create.dr_list        ~ create the lists of donors and recipients
#                 -.calculate.distancematrix  ~ calculate the distance matrix interpreting the user input from impute.NN_HD
#                 -.create.donorlimit     ~ interprets the user input from impute.NN_HD to create a single value donor limit
#                 -.match.donor_recipient ~ Matches the donors to the recipients interpreting the user input from impute.NN_HD
#                                           essentially one big switch for some C functions
#                 -.dump_diagnostics      ~ dumps the diagnostics from impute.NN_HD to a given environment
#                 -.impute                ~ performs actual imputation (copying of values)
#                 -.DATA_recode           ~ recodes data for distance calculation
#                 -.normalized.principal_components~ hmm just what it says... useful for calculating the Mahalanobis distance
impute.NN_HD<-function(DATA=NULL,distance="man",weights="range",attributes="sim", comp="rw_dist",donor_limit=Inf,optimal_donor="no",
                list_donors_recipients = NULL, diagnose = NULL)
{
   arguements<-list(distance=distance,weights=weights,attributes=attributes, comp=comp,donor_limit=donor_limit,optimal_donor=optimal_donor,
                    list_donors_recipients = list_donors_recipients, diagnose = diagnose)
   #DATA<-data.matrix(DATA)
   original_data<-DATA
   if(class(DATA)=="data.frame")
   {
     DATA<-.DATA_recode(DATA)
     recoded<-TRUE
     #test what type of weights and create appropriate vector
     weights <- .create.weightsvector(DATA = original_data, weights = weights)
     weights<-weights[DATA$weightsindex]
     DATA<-DATA$DATA_recode
   }
   else
   {
     DATA<-data.matrix(DATA)
     #test what type of weights and create appropriate vector
     weights <- .create.weightsvector(DATA = DATA, weights = weights)
     recoded<-FALSE
   }

   n<-dim(DATA)[1]
   m<-dim(DATA)[2]
   
   #currently only the simultaneous imputation of all missing values in a row is supported
   if(attributes=="sim")
   {
      #if no distance matrix is specified, then donors and recipients need to be found
      if(is.null(list_donors_recipients))
      {
        list_donors_recipients <- .create.dr_list(DATA = DATA, attributes = attributes)
      }
      
      #if distance is not matrix, matrix must be computed
      if(!is.matrix(distance))
      {
        distance<- .calculate.distancematrix(DATA = DATA, distance = distance, weights = weights, comp = comp, list_donors_recipients)
      }
  
      #make donor_limit static
      donor_limit <- .create.donorlimit(donor_limit = donor_limit, list_donors_recipients = list_donors_recipients)
  
      ##match.donor_recipient
      
      list_recip_donor<- .match.donor_recipient(list_donors_recipients = list_donors_recipients, distance = distance, 
                                                optimal_donor = optimal_donor, donor_limit = donor_limit)
      
      envir <- as.environment(".GlobalEnv")
      
      .dump_diagnostics(diagnose = diagnose, diagnostics = 
                          list(arguements = arguements, distances=distance, list_donors_recipients=list_recip_donor), 
                        envir = envir)
      
      #IMPUTE!
      DATA<- .impute(DATA = original_data, list_recip_donor = list_recip_donor, recoded = recoded)   
   }
   else
   {stop("Unimplemented way of handling attributes used: \"seq \"")}
   
   return(DATA)
}

#this function imputes the complete cases' mean by variable
impute.mean <-function(DATA=NULL)
{
   out <- .C("c_impute_mean",
             data = as.double(DATA),
             n = as.integer(dim(DATA)[1]),
             m = as.integer(dim(DATA)[2]),
             NAOK=TRUE
             )
   return(matrix(out[[1]],ncol=dim(DATA)[2]))
}

#this function preprocess the datamatrix to apply the given weights
reweight.data<-function(DATA=NULL,weights=NULL,minkovski_factor=1)
{
   weights <- (weights)^(1/minkovski_factor)
   out <- .C("reweight_data",
             data = as.double(DATA),
             n = as.integer(dim(DATA)[1]),
             m = as.integer(dim(DATA)[2]),
             weights = as.double(weights),
             NAOK=TRUE
             )
   return(matrix(out[[1]],ncol=dim(DATA)[2]))
}

#this function performs the vogel approximation heuristic matching
match.d_r_vam<-function(distance = NULL,recipients=NULL,donors=NULL,donor_limit=NULL)
{
   
   demand_fictive_recipient <- sum(donor_limit) - length(recipients);
   if(demand_fictive_recipient == 0)
   {
      distribution <- .C("c_VAM",
                         distribution = as.integer(rep(0,length(distance))),
                         distance = as.double(distance),
                         supply = as.integer(donor_limit),
                         demand = as.integer(rep(1,length(recipients))),
                         nrow_s =as.integer(length(donors)),
                         ncol_d =as.integer(length(recipients)),
                         NAOK=TRUE)$distribution
      
      distribution<-matrix(distribution,nrow=length(donors),byrow=FALSE)
   }
   if(demand_fictive_recipient>0)
   {
      distance<-cbind(distance, rep(0,length(donors)))
      
      distribution <- .C("c_VAM",
                         distribution = as.integer(rep(0,length(distance))),
                         distance = as.double(distance),
                         supply = as.integer(donor_limit),
                         demand = as.integer(c(rep(1,length(recipients)),demand_fictive_recipient)),
                         nrow_s =as.integer(length(donors)),
                         ncol_d =as.integer(length(recipients)+1),
                         NAOK=TRUE)$distribution
      
      distribution<-matrix(distribution,nrow=length(donors),byrow=FALSE)
      distribution<-distribution[,-(ncol(distribution))]
      
   }
   if(demand_fictive_recipient<0)
   {
      distance<-rbind(distance, rep(0,length(recipients)))
      
      distribution <- .C("c_VAM",
                         distribution = as.integer(rep(0,length(distance))),
                         distance = as.double(distance),
                         supply = as.integer(c(donor_limit,demand_fictive_recipient*-1)),
                         demand = as.integer(rep(1,length(recipients))),
                         nrow_s =as.integer(length(donors)+1),
                         ncol_d =as.integer(length(recipients)),
                         NAOK=TRUE)$distribution
      
      distribution<-matrix(distribution,nrow=(length(donors)+1),byrow=FALSE)
      distribution<-distribution[-(nrow(distribution)),]
      
   }
   donor_index<-apply(distribution,2,which.max)
   return(cbind(recipient=recipients,donor=donors[donor_index]))
}

#this function performs an optimal matching (by using the Rglpk package)
match.d_r_odd<-function(distance = NULL,recipients=NULL,donors=NULL,donor_limit=NULL)
{
   #Get the number of donors and recipients
   nd<-nrow(distance)
   nr<-ncol(distance)
   #setup the coefficient matrix A
   A<-matrix(0,nrow=nr+nd,ncol=nr*nd)
   #recipient requirements; colsums of distribution matrix must = 1
   for(i in 1:nr)
   {A[i,]<-c(rep(0,(i-1)*nd),rep(1,nd),rep(0,(nr-i)*nd))}
   #donor limit; rowsums of distribution matrix must <= donor_limit
   for(i in (nr+1):(nr+nd))
   {A[i,]<-c(rep(0,(i-nr-1)),rep(c(1,rep(0,nd-1)),nr-1),1,rep(0,nd - (i-nr)))}
   
   #coefficients of objective function correspond to the distances
   obj=as.vector(distance);

   #recipient requirements and donor limit as explained for A
   dir=c(rep("==",nr), rep("<=",nd));
   
   #recipient requirements and donor limit as explained for A
   rhs=c(rep(1,nr),donor_limit);

   #all variables are of integer type (maybe sometimes binary but never continuous)
   types = rep("I",nd*nr);
   
   #pass LP that was setup above to through the interface to glpk
   res<-Rglpk_solve_LP(obj,mat=A, dir, rhs, types, max = FALSE)
   
   #make a distribution matrix back out of the values of the base variables
   distribution<-matrix(res$solution,nrow=nd)
   
   #get the corresponding donors
   donor_index<-apply(distribution,2,which.max)
   return(cbind(recipient=recipients,donor=donors[donor_index]))
}

#invisible functions
#create the weights vector by interpreting the user input from impute.NN_HD
.create.weightsvector<-function(DATA = NULL, weights = "range")
{
  n<-dim(DATA)[1]
  m<-dim(DATA)[2]
  
  if(length(weights)==1)
  {
    if(is.character(weights))
    {
      weights = match.arg(arg=weights,choices=c("none","range","sd","var"),several.ok=FALSE)
      if(weights=="none")
      {weights <-rep(1,m)}
      else
      {
        if(weights=="range")
        {weights <-1/apply(apply(DATA,MARGIN=2,range,finite=TRUE),MARGIN=2,diff)}
        else
        {
          if(weights =="sd")
          {weights <-1/apply(DATA,MARGIN=2,sd,na.rm=TRUE)}
          else
          {#must be == "var"
            weights <-1/apply(DATA,MARGIN=2,var,na.rm=TRUE)
          }
        }
      }
    }
    else
    {
      ## weights must be 1 numeric
      if(is.numeric(weights))
      {weights<-rep(weights,m)}
      else
      {stop(paste("unimplemented weights used: ",weights))}
    }
  }
  else
  {
    #more than one value for weights is given
    if(!is.numeric(weights))
    {#if it is not numeric it must be either a character vector or list
      if(is.character(weights))
      {
        if(length(weights)==m)
        {stop("unimplemented feature, option 3 for weights\n")}
        else
        {
          stop("improper usage of option 3 for weights, number of weights must equal number of attributes\n")
        }
        weights<-rep(1,m)
      }
      else
      {
        #must be list
        if(length(weights)==m)
        {stop("unimplemented feature, option 5 for weights\n")}
        else
        {
          stop("improper usage of option 5 for weights, number of weights must equal number of attributes\n")
        }
        weights<-rep(1,m)
      }
    }
  }
  weights[weights==0]<-1
  
  return(weights)
}

#create the lists of donors and recipients
.create.dr_list<-function(DATA = NULL, attributes = "sim")
{
  MV_Indicator_M<-matrix(.C("c_test_NA",data=as.double(DATA),
                            min=as.integer(ncol(DATA)),nin=as.integer(nrow(DATA)),
                            NAOK=TRUE)$data,nrow=nrow(DATA))
  NA_counts <- rowSums(MV_Indicator_M)
  if(attributes == "sim")
  {
    donors <- which(NA_counts==0)
    list_donors_recipients<-list(donors=donors, recipients=(1:dim(MV_Indicator_M)[1])[-donors])
    rm(list=c("donors","MV_Indicator_M","NA_counts"))
  }
  else
  {
    #dont know how you would get here :) but lets make sure
    stop("Unimplemented way of handling attributes used: \"seq \"")
    donors <- which(NA_counts<dim(MV_Indicator_M)[2])
    recipients <- which(rowSums(MV_Indicator_M)>0)
  }
  
  return(list_donors_recipients)
}

#calculate the distance matrix interpreting the user input from impute.NN_HD
.calculate.distancematrix<- function(DATA = NULL, distance = "man", weights = rep(1,dim(DATA)[2]), comp = "rw_dist", list_donors_recipients = NULL)
{
  if(is.numeric(distance))
  {p<-distance}
  if(is.character(distance))
  {
    #perform check, allowing for partial matches
    distance = match.arg(arg=distance,choices=c("man","eukl","tscheb","mahal"),several.ok=FALSE)
    if(distance == "man"){p <- 1}
    if(distance == "eukl"){p <- 2}
    if(distance == "tscheb"){p <- 1}
    if(distance == "mahal"){p <- 2}
    
  }
  
  #reweight data
  if(!all(weights==1))
  {
    DATA<- reweight.data(DATA,weights,p)
  } 
  #compensate for missing data
  #perform check, allowing for partial matches
    comp = match.arg(arg=comp,choices=c("rw_dist","mean","rseq","rsim"),several.ok=FALSE)
    #if(comp=="rw_dist")
    #{#no action needed}
    if(comp=="mean")
    {DATA<-impute.mean(DATA)}
    if(comp=="rseq")
    {
      #will work with real and integer data matrix
      #random if order is random
      DATA<-impute.SEQ_HD(DATA=DATA,initialvalues=0, navalues=NA, modifyinplace = FALSE);
      #stop("unimplemented feature, comp = \"rseq\" \n")
    }
    if(comp=="rsim")
    {stop("unimplemented feature, comp = \"rsim\" \n")}
  
  #special case for mahalanobis distance,
  if(distance == "mahal")
  {
    #yea if there is a way to combine these two ways, feel free to contact me ;-)
    #Theoretically this should be a bad choice anyway... available case for principal components, yuck
    if(comp=="rw_dist")
    {
      stop("distance = \"mahal\" may not be combined with comp==\"rw_dist\" \n")
    }
    #if you think this is cool, feel free to contact me
    DATA<-.normalized.principal_components(DATA)
  }
  
  if(distance == "tscheb")
  {
    #quick solution, extra c function for the tschebeycheff distance... do people read this?
    out <- .C("c_dist_hot_deck_tscheb",
              distances = as.double(rep(0,length(list_donors_recipients$donors)*length(list_donors_recipients$recipients))),
              data = as.double(DATA),
              n = as.integer(dim(DATA)[1]),
              m = as.integer(dim(DATA)[2]),
              donors = as.integer(list_donors_recipients$donors-1),
              recipients = as.integer(list_donors_recipients$recipients-1),
              n_donors = as.integer(length(list_donors_recipients$donors)),
              n_recipients = as.integer(length(list_donors_recipients$recipients)),
              NAOK=TRUE)
    
    distances<-matrix(out[[1]],nrow=length(list_donors_recipients$donors),byrow=FALSE)
    rownames(distances)<-list_donors_recipients$donors
    colnames(distances)<-list_donors_recipients$recipients
  }
  else
  {
    out <- .C("c_dist_hot_deck",
              distances = as.double(rep(0,length(list_donors_recipients$donors)*length(list_donors_recipients$recipients))),
              data = as.double(DATA),
              n = as.integer(dim(DATA)[1]),
              m = as.integer(dim(DATA)[2]),
              donors = as.integer(list_donors_recipients$donors-1),
              recipients = as.integer(list_donors_recipients$recipients-1),
              n_donors = as.integer(length(list_donors_recipients$donors)),
              n_recipients = as.integer(length(list_donors_recipients$recipients)),
              p = as.double(p),
              NAOK=TRUE)
    
    distances<-matrix(out[[1]],nrow=length(list_donors_recipients$donors),byrow=FALSE)
    rownames(distances)<-list_donors_recipients$donors
    colnames(distances)<-list_donors_recipients$recipients
  }
  return(distances)
}

#interprets the user input from impute.NN_HD to create a single value donor limit
.create.donorlimit<-function(donor_limit = Inf, list_donors_recipients = NULL)
{
  if(donor_limit >0)
  {
    if(donor_limit< 1)
    {
      #dynamic donor limit used
      donor_limit <- ceiling(donor_limit * length(list_donors_recipients$recipients))
    }else{
      #static donor limit used
      donor_limit<-ceiling(donor_limit)
    }
  }else{
    stop(paste("improper usage of donor_limit: ",donor_limit))}
  if(is.infinite(donor_limit))
  {donor_limit<-length(list_donors_recipients$recipients)}
  
  if((donor_limit*length(list_donors_recipients$donors))<length(list_donors_recipients$recipients))
  {stop(paste("Donor limit set too stringent. Number of donors * donor-limit must be >= Number of recipients.\n Compare:",
              donor_limit*length(list_donors_recipients$donors),"vs.",length(list_donors_recipients$recipients)))
  }
  return(donor_limit)
}

#Matches the donors to the recipients interpreting the user input from impute.NN_HD
#essentially one big switch
.match.donor_recipient<- function(list_donors_recipients = NULL, distance = NULL, optimal_donor = "no", donor_limit = 1)
{
  #partial matching for optimal donor
  optimal_donor = match.arg(arg=optimal_donor,choices=c("no", "rand", "mmin","modifvam","vam", "odd", "modi"),several.ok=FALSE)
  
  if(optimal_donor == "no")
  {
    out <- .C("c_match_d_r",
              recipient_donor_list = as.integer(rep(-1,length(list_donors_recipients$recipients))),
              distance = as.double(distance),
              recipient_order = as.integer( (1:dim(distance)[2])-1),
              donor_limit = as.double(rep(donor_limit,dim(distance)[1])),
              n_donors = as.integer(dim(distance)[1]),
              n_recipients = as.integer(dim(distance)[2]),
              NAOK=TRUE)
    list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[out[[1]]])
    
  }
  if(optimal_donor == "rand")
  {
    out <- .C("c_match_d_r",
              recipient_donor_list = as.integer(rep(-1,length(list_donors_recipients$recipients))),
              distance = as.double(distance),
              recipient_order = as.integer(sample(1:dim(distance)[2])-1),
              donor_limit = as.double(rep(donor_limit,dim(distance)[1])),
              n_donors = as.integer(dim(distance)[1]),
              n_recipients = as.integer(dim(distance)[2]),
              NAOK=TRUE)
    list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[out[[1]]])
    
  }
  if(optimal_donor == "mmin")
  {
    
    distribution <- .C("c_matmin",
                       distribution = as.integer(rep(0,length(distance))),
                       distance = as.double(distance),
                       supply = as.integer(rep(donor_limit,dim(distance)[1])),
                       demand = as.integer(rep(1,length(list_donors_recipients$recipients))),
                       nrow_s =as.integer(length(list_donors_recipients$donors)),
                       ncol_d =as.integer(length(list_donors_recipients$recipients)),
                       NAOK=TRUE)$distribution
    
    distribution<-matrix(distribution,nrow=length(list_donors_recipients$donors),byrow=FALSE)
    donor_index<-apply(distribution,2,which.max)
    list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[donor_index])
    
  }
  if(optimal_donor == "modifvam")
  {
    distribution <- .C("c_modifVAM",
                       distribution = as.integer(rep(0,length(distance))),
                       distance = as.double(distance),
                       supply = as.integer(rep(donor_limit,dim(distance)[1])),
                       demand = as.integer(rep(1,length(list_donors_recipients$recipients))),
                       nrow_s =as.integer(length(list_donors_recipients$donors)),
                       ncol_d =as.integer(length(list_donors_recipients$recipients)),
                       NAOK=TRUE)$distribution
    
    distribution<-matrix(distribution,nrow=length(list_donors_recipients$donors),byrow=FALSE)
    donor_index<-apply(distribution,2,which.max)
    list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[donor_index])
    
  }
  if((optimal_donor == "odd")|(optimal_donor == "modi"))
  {
    list_recip_donor<-match.d_r_odd(distance = distance,recipients=list_donors_recipients$recipients,
                                    donors=list_donors_recipients$donors,donor_limit=rep(donor_limit,dim(distance)[1]))
    
  }
  if(optimal_donor == "vam")
  {
    list_recip_donor<-match.d_r_vam(distance = distance,recipients=list_donors_recipients$recipients,
                                    donors=list_donors_recipients$donors,donor_limit=rep(donor_limit,dim(distance)[1]))
    
  }
  #if(optimal_donor == "modi")
  #{stop("unimplemented feature, optimal_donor = \"modi\" \n try optimal_donor = \"odd\" for same results using a different algorithm")}
  return(list_recip_donor)
}

#dumps the diagnostics from impute.NN_HD to a given environment
.dump_diagnostics<-function(diagnose = NULL, diagnostics = list(), envir = NULL)
{
  if(!is.null(diagnose))
  { 
    if(is.character(diagnose))
    {
      reserved_words<-c("if","else","repeat","while","function","for","in","next",
                        "break","TRUE","FALSE","NULL","Inf","NaN","NA","NA_integer_",
                        "NA_real_","NA_complex_","NA_character_","c","q","s","t","C","D","F","I","T")
      if(diagnose %in% reserved_words)
      {
        warning("diagnose is character string representing a reserved word. No diagnostics saved.")
      }else
      {
        assign(diagnose, diagnostics , envir = envir)
      }
      
    }else
    {
      warning("diagnose is not NULL but also not a character string. No diagnostics saved.")
    }
  }
  
  invisible(NULL)
}

#performs actual imputation (copying of values)
.impute<-function(DATA = NULL, list_recip_donor= NULL, recoded = FALSE)
{
  n <- dim(DATA)[1]
  m <- dim(DATA)[2]
  if(recoded)
  {
    for(i in 1:dim(list_recip_donor)[1])
    {
      current_recip<-list_recip_donor[i,1]
      current_donor<-list_recip_donor[i,2]
      missing_data_recip<-is.na(DATA[current_recip,])
      #replace
      DATA[current_recip,missing_data_recip]<-DATA[current_donor,missing_data_recip]
    }
  }
  else
  {
    DATA <- .C("c_impute_NN_HD",
               data = as.double(DATA),
               n = as.integer(n),
               m = as.integer(m),
               recipients = as.integer(list_recip_donor[,1]-1),
               n_recipients = as.integer(dim(list_recip_donor)[1]),
               donors = as.integer(list_recip_donor[,2]-1),
               NAOK=TRUE)$data
    DATA<-matrix(DATA,nrow=n,ncol=m)
  }
  return(DATA)
}

#recode data from data.frame to make into data.matrix
.DATA_recode<-function(DATA=NULL,conversionmode = "auto",varstodummy=c())
{
  #check if DATA is of class data.frame
  if(class(DATA)!="data.frame")
  {stop("DATA is not a: \"data.frame \"")}
  DATA_recode<-c()
  options(na.action='na.pass')
  
  weightsindex<-c()
  cnames<-c()
  
  conversionmode = match.arg(arg=conversionmode,choices=c("auto","manual"),several.ok=FALSE)
  #automatically detect variables which will be converted to dummy codeing
  if(conversionmode=="auto")
  {
    #find out what class each column is and then itterate over the columns
    DATA.classes<-sapply(DATA, class)
    
    for(i in 1:length(DATA.classes))
    {
      #if the column contains factors or characters (I hope this covers all non-numerics)
      if(DATA.classes[i]%in%c("factor","character"))
      {
        if(DATA.classes[i]=="character")
        {DATAtoconvert<-as.factor(DATA[,i])}
        else
        {DATAtoconvert<-DATA[,i]}
        
        if(length(levels(DATAtoconvert))==1)
        {temp<-rep(1,length(DATAtoconvert));times<-1}
        else
        {temp<-model.matrix(object = ~ DATAtoconvert-1);times<-dim(temp)[2]}
        
        cnames<-c(cnames,paste(names(DATA)[i],levels(DATAtoconvert),sep="."))
        DATA_recode<-cbind(DATA_recode,temp)
        
        weightsindex<-c(weightsindex,rep(i,times))
        
        next
      }
      else  
      {
        temp<-DATA[,i]
        DATA_recode<-cbind(DATA_recode,temp)
        
        cnames<-c(cnames,colnames(DATA)[i])
        
        weightsindex<-c(weightsindex,i)
      }
    }
  }
  if(conversionmode=="manual") #do not automatically detect variables but convert those given
  {
    for(i in 1:dim(DATA)[2])
    {
      if(i %in% varstodummy)
      {
        
        DATAtoconvert<-as.factor(DATA[,i])
        
        if(length(levels(DATAtoconvert))==1)
        {temp<-rep(1,length(DATAtoconvert));times<-1}
        else
        {temp<-model.matrix(object = ~ DATAtoconvert-1);times<-dim(temp)[2]}
        
        cnames<-c(cnames,paste(names(DATA)[i],levels(DATAtoconvert),sep="."))
        DATA_recode<-cbind(DATA_recode,temp)
        
        weightsindex<-c(weightsindex,rep(i,times))
        next
      }
      else
      {
        temp<-DATA[,i]
        DATA_recode<-cbind(DATA_recode,temp)
        cnames<-c(cnames,colnames(DATA)[i])
        
        weightsindex<-c(weightsindex,i)
        
      }
    }
  }
  DATA_recode<-as.matrix(DATA_recode)
  colnames(DATA_recode)<-cnames
  return(list(DATA_recode=DATA_recode,weightsindex=weightsindex))
}

#calculate the normalized principal components
.normalized.principal_components<-function(DATA=NULL)
{
  covarianceM<-cov(DATA)
  eigen_val_vec<-eigen(covarianceM,symmetric = TRUE)
  DATA<-DATA%*%(eigen_val_vec$vectors)
  DATA<-sweep(DATA,MARGIN=2,sqrt(eigen_val_vec$values),"/")
  return(DATA)
}