"varPed" <-
function(x, gender=NULL, lag=c(0,0), relational=FALSE, lag_relational=c(0,0), restrict=NULL, keep=FALSE, USvar=NULL, merge=FALSE, NAvar=NULL){

      if(relational!=FALSE){
        if(relational!="OFFSPRING" & relational!="MATE" & relational!="OFFSPRINGV" & relational!="MATEV"){stop("relational must either be 'OFFSPRING', 'OFFSPRINGV', 'MATE' or 'MATEV'")}
      }

      if((relational=="OFFSPRING" | relational=="MATE" | relational=="OFFSPRINGV" | relational=="MATEV") & is.null(restrict)==FALSE){
        if(restrict==TRUE){restrict<-"=="}
        if(restrict==FALSE){restrict<-"!="}
      }  # fix for making code from versions <2.1 back compatible
     
      sex<-with(parent.frame(), sex)               
      id<-with(parent.frame(), id)
      off_record<-with(parent.frame(), off_record)         
      data<-with(parent.frame(), data)                 # these are objects taken from the environment in which
      keepDam<-with(parent.frame(), keepDam)           # varPed is called - typically MCMCped
      keepSire<-with(parent.frame(), keepSire)  
      time_var<-with(parent.frame(), timevar) 
      namevar<-x

#      if(length(namevar)==1){
#        if(class(data[,x])%in%c("factor", "numeric", "integer")==FALSE){
#        stop("variables must be numeric or factors")
#        }
#      }else{
#        if(any(unlist(lapply(data[,x], class))%in%c("factor", "numeric", "integer")==FALSE)){
#          stop("variables must be numeric or factors")
#        }
#      }

      x<-data[match(x, names(data))]                                   # gets variable(s): Jarrod this used NOT to be a matrix

#      if(length(USvar)>0 & relational==FALSE){
#        x[which(is.na(x)==TRUE)]<-USvar                # Fills missing values if specified in USvar
#      }     
      hermaphrodite<-length(sex)==0                    # is this an hermaphroditic system
      sex_specific<-length(gender)==1                  # is the variable sex-specific
      ud<-with(parent.frame(), USdam)
      if(length(ud)!=1 | ud[1]!=FALSE){
        ud<-TRUE
      }
      us<-with(parent.frame(), USsire)
      if(length(us)!=1 | us[1]!=FALSE){
        us<-TRUE
      }

#####################################################################################################################
###################################  restricting variables ##########################################################
#####################################################################################################################

if(length(restrict)!=0){  

  PedDesMatrix<-list(Dam=list(id=NULL), Sire=list(id=NULL), Dam_restrict=list(id=NULL), Sire_restrict=list(id=NULL)) 

  if(hermaphrodite==FALSE){
     PedDesMatrix$Dam$id<-unique(id[which(sex=="Female")])
     PedDesMatrix$Sire$id<-unique(id[which(sex=="Male")])
     PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female")])
     PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male")])
  }else{
     PedDesMatrix$Dam$id<-unique(id)
     PedDesMatrix$Sire$id<-unique(id)
     PedDesMatrix$Dam_restrict$id<-PedDesMatrix$Dam$id
     PedDesMatrix$Sire_restrict$id<-PedDesMatrix$Sire$id
  }

  not_after_off<-c((time_var>=(time_var[off_record]+lag[1]) & time_var<=(time_var[off_record]+lag[2])) | is.na(time_var))

  if(relational==FALSE){
    if(hermaphrodite==FALSE){
      if("Female"%in%gender | sex_specific==FALSE){ 
        PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (x==restrict | is.na(x)==TRUE) & not_after_off)])
        if(keep==FALSE){
          PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
        }
      }
      if("Male"%in%gender | sex_specific==FALSE){
        PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (x==restrict | is.na(x)==TRUE) & not_after_off)])
        if(keep==FALSE){
           PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
        }
      }
    }else{
      if(sex_specific==TRUE){
        if("Male"%in%gender){
           PedDesMatrix$Sire_restrict$id<-unique(id[which((x==restrict | is.na(x)==TRUE) & not_after_off)])
        }
        if("Female"%in%gender){
          PedDesMatrix$Dam_restrict$id<-unique(id[which((x==restrict | is.na(x)==TRUE) & not_after_off)])
        }
      }else{
        PedDesMatrix$Dam_restrict$id<-unique(id[which((x==restrict | is.na(x)==TRUE) & not_after_off)])
        PedDesMatrix$Sire_restrict$id<-unique(id[which((x==restrict | is.na(x)==TRUE) & not_after_off)])
      }
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
  }

  if(relational=="OFFSPRING" | relational=="OFFSPRINGV"){
    restrict.comp<-substr(restrict, nchar(restrict), nchar(restrict))%in%c("=", "<", ">")
    # if true then the restriction is based on comparing parent and offspring phenotye
    # if false then the restriction is based on comparing teh distance between parent and offspring phenotye with what ever appears to the right of the inequality

    if(lag_relational[1]!=0 | lag_relational[2]!=0){
       off_time<-time_var[off_record]
       off_records<-which(id%in%id[off_record] & ((time_var>=(lag_relational[1]+off_time) & time_var<=(lag_relational[2]+off_time)) | is.na(time_var)))
    }
    if(hermaphrodite==FALSE){
      if("Female"%in%gender | sex_specific==FALSE){
        if(lag_relational[1]==0 & lag_relational[2]==0){
          if(restrict.comp){
           PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (eval(parse(text=paste("x[1:dim(x)[1],]", restrict, "x[off_record,]"))) | rowSums(is.na(x))>0) & not_after_off)])
          }else{
            if(relational=="OFFSPRING"){
              PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (eval(parse(text=paste("sqrt(colSums((t(as.matrix(x))-as.matrix(x)[off_record,])^2))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (eval(parse(text=paste("colSums((t(as.matrix(x))-as.matrix(x)[off_record,]))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
            }
          }
        }else{
          if(restrict.comp){
            PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (eval(parse(text=paste("rowSums(outer(x[1:dim(x)[1],], x[off_records,], {restrict}))>0"))) | rowSums(is.na(x))>0) & not_after_off)])
          }else{
            if(relational=="OFFSPRING"){
              PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){sqrt(colSums((X-t(x[off_records,]))^2))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){colSums((X-t(x[off_records,])))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
            }
          }
        }
        if(keep==FALSE){
          PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
        }
      }
      if("Male"%in%gender | sex_specific==FALSE){
        if(lag_relational[1]==0 & lag_relational[2]==0){
          if(restrict.comp){
            PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (eval(parse(text=paste("x[1:dim(x)[1],]", restrict, "x[off_record,]"))) | rowSums(is.na(x))>0) & not_after_off)])
          }else{
            if(relational=="OFFSPRING"){
              PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (eval(parse(text=paste("sqrt(colSums((t(as.matrix(x))-as.matrix(x)[off_record,])^2))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (eval(parse(text=paste("colSums((t(as.matrix(x))-as.matrix(x)[off_record,]))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
            }
          }
        }else{
          if(restrict.comp){
            PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (eval(parse(text=paste("rowSums(outer(x[1:dim(x)[1],], x[off_records,], {restrict}))>0"))) | rowSums(is.na(x))>0) & not_after_off)])
          }else{
            if(relational=="OFFSPRING"){
              PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){sqrt(colSums((X-t(x[off_records,]))^2))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){colSums((X-t(x[off_records,])))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
            }
          }
        }
        if(keep==FALSE){
           PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
        }
      }
    }else{
      if(sex_specific==TRUE){
        if("Male"%in%gender){
          if(lag_relational[1]==0 & lag_relational[2]==0){
            if(restrict.comp){
              PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("x[1:dim(x)[1],]", restrict, "x[off_record,]"))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              if(relational=="OFFSPRING"){
                PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("sqrt(colSums((t(as.matrix(x))-as.matrix(x)[off_record,])^2))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
              }else{
                PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("colSums((t(as.matrix(x))-as.matrix(x)[off_record,]))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
              }
            }
          }else{
            if(restrict.comp){
              PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(outer(x[1:dim(x)[1],], x[off_records,], {restrict}))>0"))) | rowSums(is.na(x))>0) & not_after_off)])
             }else{
               if(relational=="OFFSPRING"){
                 PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){sqrt(colSums((X-t(x[off_records,]))^2))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
               }else{
                 PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){colSums((X-t(x[off_records,])))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
               }
             }
          }
        }
        if("Female"%in%gender){
          if(lag_relational[1]==0 & lag_relational[2]==0){
            if(restrict.comp){
              PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("x[1:dim(x)[1],]", restrict, "x[off_record,]"))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              if(relational=="OFFSPRING"){
                PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("sqrt(colSums((t(as.matrix(x))-as.matrix(x)[off_record,])^2))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
              }else{
                PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("colSums((t(as.matrix(x))-as.matrix(x)[off_record,]))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
              }
            }
          }else{
            if(restrict.comp){
              PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(outer(x[1:dim(x)[1],], x[off_records,], {restrict}))>0"))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
              if(relational=="OFFSPRING"){
                PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){sqrt(colSums((X-t(x[off_records,]))^2))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
              }else{
                PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){colSums((X-t(x[off_records,])))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
              }
            }
          }
        }
      }else{
        if(lag_relational[1]==0 & lag_relational[2]==0){
          if(restrict.comp){
            PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("x[1:dim(x)[1],]", restrict, "x[off_record,]"))) | rowSums(is.na(x))>0) & not_after_off)])
            PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("x[1:dim(x)[1],]", restrict, "x[off_record,]"))) | rowSums(is.na(x))>0) & not_after_off)])
          }else{
            if(relational=="OFFSPRING"){
              PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("sqrt(colSums((t(as.matrix(x))-as.matrix(x)[off_record,])^2))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
              PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("sqrt(colSums((t(as.matrix(x))-as.matrix(x)[off_record,])^2))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
            }else{
             PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("colSums((t(as.matrix(x))-as.matrix(x)[off_record,]))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
              PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("colSums((t(as.matrix(x))-as.matrix(x)[off_record,]))", restrict))) | rowSums(is.na(x))>0) & not_after_off)])
            }
          }
        }else{
          if(restrict.comp){
            PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(outer(x[1:dim(x)[1],], x[off_records,], {restrict}))>0"))) | rowSums(is.na(x))>0) & not_after_off)])
            PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(outer(x[1:dim(x)[1],], x[off_records,], {restrict}))>0"))) | rowSums(is.na(x))>0) & not_after_off)])
           }else{
             if(relational=="OFFSPRING"){
               PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){sqrt(colSums((X-t(x[off_records,]))^2))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
               PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){sqrt(colSums((X-t(x[off_records,]))^2))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
             }else{
               PedDesMatrix$Dam_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){colSums((X-t(x[off_records,])))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
               PedDesMatrix$Sire_restrict$id<-unique(id[which((eval(parse(text=paste("rowSums(t(apply(X=x, 1, function(X){colSums((X-t(x[off_records,])))}))", restrict, ")>0"))) | rowSums(is.na(x))>0) & not_after_off)])
             }
           }
        } 
      }
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
  }

  if(hermaphrodite==TRUE){
    if(us==TRUE){
      remove.extra<-match(max(id), PedDesMatrix$Dam$id)
      if(length(remove.extra)>0){
        PedDesMatrix$Dam$id<-PedDesMatrix$Dam$id[-match(max(id), PedDesMatrix$Dam$id)]
        PedDesMatrix$Dam_restrict$id<-PedDesMatrix$Dam_restrict$id[-match(max(id), PedDesMatrix$Dam_restrict$id)]
      }
    }
    if(ud==TRUE){
      if(us==FALSE){
        remove.extra<-match(max(id), PedDesMatrix$Sire$id)
        if(length(remove.extra)>0){
          PedDesMatrix$Sire$id<-PedDesMatrix$Sire$id[-match(max(id), PedDesMatrix$Sire$id)]
          PedDesMatrix$Sire_restrict$id<-PedDesMatrix$Sire_restrict$id[-match(max(id), PedDesMatrix$Sire_restrict$id)]
        }
      }else{
        remove.extra<-match(id[which(id==c(max(id)-1))], PedDesMatrix$Sire$id)
        if(length(remove.extra)>0){
          PedDesMatrix$Sire$id<-PedDesMatrix$Sire$id[-match(id[which(id==c(max(id)-1))], PedDesMatrix$Sire$id)]
          PedDesMatrix$Sire_restrict$id<-PedDesMatrix$Sire_restrict$id[-match(id[which(id==c(max(id)-1))], PedDesMatrix$Sire_restrict$id)]
        }
      }
    }
  }
}
#####################################################################################################################
##########################################  true variables ##########################################################
#####################################################################################################################


if(length(restrict)==0){

    PedDesMatrix<-list(Dam=list(var_name=NULL, var_type=NULL, X=NULL, merge=FALSE), Sire=list(var_name=NULL, var_type=NULL, X=NULL, merge=FALSE), DamSire=list(var_name=NULL, var_type=NULL, X=NULL, merge=FALSE))    # output; design matrix+info

    if(dim(x)[2]==1){                          # gets class of variable(s)
      facnum<-class(x[,1])
      if(facnum=="integer"){
        facnum<-"numeric"
      }
    }else{
      facnum<-apply(x, 2, class)
      facnum[which(facnum=="integer")]<-"numeric"
    }        

    predict_ped=NULL

    off_time<-time_var[off_record]

    if(relational=="OFFSPRING" & lag_relational[1]!=0 & lag_relational[2]!=0){
      off_var<-as.matrix(x)[off_record,]
    }else{
      off_var<-as.matrix(x)[which(id%in%id[off_record] & ((time_var>=(lag_relational[1]+off_time) & time_var<=(lag_relational[2]+off_time)) | is.na(time_var))),]
    }

############### Covariates of fecundity, or distance from offspring ###############################

    if(relational!="MATE"){

      for(g in 1:c(2-length(gender))){

        if(sex_specific==FALSE & g==1){gender<-"Female"}
        if(sex_specific==FALSE & g==2){gender<-"Male"}
        if("Female"%in%gender){ 
          var_tmp<-subset(x, id%in%keepDam==TRUE)
          time_tmp<-subset(time_var, id%in%keepDam==TRUE)
          id_tmp<-subset(id, id%in%keepDam==TRUE)
        }
        if("Male"%in%gender){
          var_tmp<-subset(x, id%in%keepSire==TRUE)
          time_tmp<-subset(time_var, id%in%keepSire==TRUE)
          id_tmp<-subset(id, id%in%keepSire==TRUE)
        }

        time_for_P<-which((time_tmp>=(lag[1]+off_time) & time_tmp<=(lag[2]+off_time)) | is.na(time_tmp))
  
        # OFFSPRING RELATIONAL - NUMERIC

        if((relational=="OFFSPRING" | relational=="OFFSPRINGV") & "numeric"%in%facnum){   
          var_tmp<-as.matrix(as.matrix(var_tmp)[time_for_P,])      
          id_tmp<-id_tmp[time_for_P]
           if(is.null(dim(var_tmp))){
            var_tmp<-t(as.matrix(var_tmp)) 
          }
          dup_off_var<-rep(1, length(var_tmp[,1]))%*%t(off_var)
          if(relational=="OFFSPRING"){
            predict_ped<-rowSums((var_tmp-dup_off_var)^2)^0.5
          }else{
            predict_ped<-rowSums(var_tmp-dup_off_var)
          }
          if(length(NAvar)>0){
            predict_ped[which(is.na(predict_ped)==T)]<-NAvar          # Fills missing values if specified in USvar
          }
          namex<-namevar
          predict_ped<-tapply(predict_ped, id_tmp, mean, na.rm=T)
          if("Male"%in%gender){
            predict_ped<-predict_ped[match(keepSire,names(predict_ped))]
            if(us & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          if("Female"%in%gender){
            predict_ped<-predict_ped[match(keepDam,names(predict_ped))]
            if(ud & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          predict_ped<-matrix(predict_ped, length(predict_ped),1)
        }

        # OFFSPRING RELATIONAL - FACTOR

        if(relational=="OFFSPRING" & "factor"%in%facnum){
          var_tmp<-var_tmp[time_for_P,]
          id_tmp<-id_tmp[time_for_P]
          NAvec<-which(is.na(var_tmp)==TRUE)
          var_tmp<-var_tmp%in%off_var
          if(length(NAvar)>0){
            var_tmp[NAvec]<-NAvar      # Fills missing values if specified in USvar
          }else{
            var_tmp[NAvec]<-NA
          }    
          predict_ped<-tapply(var_tmp, id_tmp, mean, na.rm=T)>0
          namex<-paste(namevar, c(TRUE, FALSE), sep=".")
          if("Male"%in%gender){
            predict_ped<-predict_ped[match(keepSire,names(predict_ped))]
            if(us & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          if("Female"%in%gender){
            predict_ped<-predict_ped[match(keepDam,names(predict_ped))]
            if(ud & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          predict_ped<-matrix(predict_ped, length(predict_ped),1)      
        }

# FECUNDITY - NUMERIC

        if(relational==FALSE & "numeric"%in%facnum){
          var_tmp<-var_tmp[time_for_P,]
          id_tmp<-id_tmp[time_for_P]
          var_tmp<-tapply(var_tmp, id_tmp, mean, na.rm=T)
          id_tmp<-names(var_tmp)
          predict_ped<-var_tmp
          if(length(NAvar)>0){
            predict_ped[which(is.na(predict_ped)==T)]<-NAvar                # Fills missing values if specified in USvar
          }
          if("Male"%in%gender){
            predict_ped<-predict_ped[match(keepSire,names(predict_ped))]
            if(us & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          if("Female"%in%gender){
            predict_ped<-predict_ped[match(keepDam,names(predict_ped))]
            if(ud & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          predict_ped<-matrix(predict_ped, length(predict_ped),1)
          namex<-namevar
        }

# FECUNDITY - FACTOR

        if(relational==FALSE & "factor"%in%facnum){
          var_tmp<-var_tmp[time_for_P,]
          id_tmp<-id_tmp[time_for_P]
          if(any(tapply(var_tmp, id_tmp, function(x){sum(duplicated(x)==FALSE)>1}))){
            stop(paste(namevar, "varies over time and is a factor"))
          }
          var_tmp<-unlist(tapply(var_tmp, id_tmp, function(x){x[1]}, simplify=FALSE))
          NAvec<-which(is.na(var_tmp)==TRUE)
          if(length(NAvar)>0){
            var_tmp[NAvec]<-NAvar                # Fills missing values if specified in USvar
          }else{
            var_tmp[NAvec]<-levels(var_tmp)[1]
          }
          predict_ped<-as.matrix(model.matrix(~var_tmp)[,-1])
          rownames(predict_ped)<-names(var_tmp)
          if(length(NAvar)==0){
            predict_ped[NAvec,]<-NA
          }
          if("Male"%in%gender){
            predict_ped<-as.matrix(predict_ped[match(keepSire,rownames(predict_ped)),])
            if(us & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          if("Female"%in%gender){
            predict_ped<-as.matrix(predict_ped[match(keepDam,rownames(predict_ped)),])
            if(ud & is.null(USvar)==FALSE){
              predict_ped[length(predict_ped)]<-USvar
            }
          }
          namex<-paste(namevar, levels(var_tmp)[-1], sep=".")
          colnames(predict_ped)<-levels(var_tmp)[-1]
        }

        if(gender=="Female"){
          PedDesMatrix$Dam$X<-predict_ped
          if(sex_specific==TRUE){
            PedDesMatrix$Dam$var_name<-namex
          }else{
            PedDesMatrix$Dam$var_name<-paste(namex, "linked", sep=".")
          }
          PedDesMatrix$Dam$var_type<-facnum[1]
          if(merge==TRUE){
            PedDesMatrix$Dam$merge<-TRUE
          }
        }

        if(gender=="Male"){
          PedDesMatrix$Sire$X<-predict_ped
          if(sex_specific==TRUE){
            PedDesMatrix$Sire$var_name<-namex
          }else{
            PedDesMatrix$Sire$var_name<-paste(namex, "linked", sep=".")
          }
          PedDesMatrix$Sire$var_type<-facnum[1]
          if(merge==TRUE){
            PedDesMatrix$Sire$merge<-TRUE
          }
        }
      }
    }


############### Covariates of distance from mate ###############################

    if(relational=="MATE" | relational=="MATEV"){

      var_tmpF<-subset(x, id%in%keepDam==TRUE)
      time_tmpF<-subset(time_var, id%in%keepDam==TRUE)
      id_tmpF<-subset(id, id%in%keepDam==TRUE)
      var_tmpM<-subset(x, id%in%keepSire==TRUE)
      time_tmpM<-subset(time_var, id%in%keepSire==TRUE)
      id_tmpM<-subset(id, id%in%keepSire==TRUE)

      if(sex_specific==FALSE){gender="Female"}

      if(gender=="Female"){
        lagF<-lag
        lagM<-lag_relational
      }else{
        lagM<-lag
        lagF<-lag_relational
      }

      timePM<-c((time_tmpM>=(lagM[1]+off_time) & time_tmpM<=(lagM[2]+off_time)) | is.na(time_tmpM))
      timePF<-c((time_tmpF>=(lagF[1]+off_time) & time_tmpF<=(lagF[2]+off_time)) | is.na(time_tmpF))

      var_tmpM<-as.matrix(subset(var_tmpM, timePM))
      var_tmpF<-as.matrix(subset(var_tmpF, timePF))

      id_tmpM<-subset(id_tmpM, timePM)
      id_tmpF<-subset(id_tmpF, timePF)

      distmat<-matrix(0, nrow=length(id_tmpM), ncol=length(id_tmpF))
      id<-paste(rep(id_tmpF, each=length(id_tmpM)), rep(id_tmpM, length(id_tmpF)))


      if("numeric"%in%facnum){
        if(length(id)>0){
          if(relational=="MATE"){
            for(d in 1:length(var_tmpM[1,])){
              distmat<-distmat+(outer(c(var_tmpM[,d]), c(var_tmpF[,d]), "-")^2)
            }
            predict_ped<-tapply(c(distmat^0.5), id, mean, na.rm=T)
          }else{
            if(gender=="Female"){
              for(d in 1:length(var_tmpM[1,])){
                distmat<-distmat+(outer(c(var_tmpF[,d]), c(var_tmpM[,d]), "-"))
              }
            }else{
              for(d in 1:length(var_tmpM[1,])){
                distmat<-distmat+(outer(c(var_tmpM[,d]), c(var_tmpF[,d]), "-"))
              }
            }
            predict_ped<-tapply(distmat, id, mean, na.rm=T)       
          }
          predict_ped<-predict_ped[match(paste(rep(keepDam, each=length(keepSire)), rep(keepSire, length(keepDam))), names(predict_ped))]
        }else{
          predict_ped<-matrix(NA, length(keepDam)*length(keepSire), 1)
        }
        if(length(NAvar)>0){
          predict_ped[which(is.na(predict_ped)==TRUE)]<-NAvar                # Fills missing values if specified in USvar
        }
        if((us | ud) & is.null(USvar)==FALSE){
           us.samples<-c()
           if(ud){
             us.samples<-c(us.samples, (length(keepDam)-1)*length(keepSire)+1:length(keepSire))
           }
           if(us){
             us.samples<-c(us.samples, 0:length(keepDam)*length(keepSire))
           }
           predict_ped[us.samples]<-USvar 
        }
        namex<-namevar
      }

      if("factor"%in%facnum){
        if(length(id)>0){     
          distmat<-outer(c(var_tmpM), c(var_tmpF), "==")
          predict_ped<-tapply(distmat, id, mean, na.rm=T)>0
          predict_ped<-predict_ped[match(paste(rep(keepDam, each=length(keepSire)), rep(keepSire, length(keepDam))), names(predict_ped))]
        }else{
          predict_ped<-matrix(NA, length(keepDam)*length(keepSire), 1)
        }
        if(length(NAvar)>0){
          predict_ped[which(is.na(predict_ped)==TRUE)]<-NAvar                # Fills missing values if specified in USvar
        }
        if((us | ud) & is.null(USvar)==FALSE){
           us.samples<-c()
           if(ud){
             us.samples<-c(us.samples, (length(keepDam)-1)*length(keepSire)+1:length(keepSire))
           }
           if(us){
             us.samples<-c(us.samples, 0:length(keepDam)*length(keepSire))
           }
           predict_ped[us.samples]<-USvar 
        }
        namex<-paste(namevar, c(TRUE, FALSE), sep=".")
      }

      predict_ped<-matrix(predict_ped, length(predict_ped),1)
      PedDesMatrix$DamSire$X<-predict_ped
      PedDesMatrix$DamSire$var_name<-namex
      PedDesMatrix$DamSire$var_type<-facnum[1]

      if(merge==TRUE){
        PedDesMatrix$DamSire$merge<-TRUE
      }
    }
  }
PedDesMatrix
}
