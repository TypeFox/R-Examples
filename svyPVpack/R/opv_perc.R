opv_perc <-
function(by, svydat, pvs=NULL)
  {
    
if(!any(is.null(pvs)))
  {
  avar <- c(pvs, all.vars(by))
  dum1 <- apply(svydat$variables[,avar], 1, function(x) any(is.na(x)))
  svydat.sub <- subset(svydat, dum1 == FALSE)
  
  } else
      {
      svydat.sub  <- svydat
      }
  
    tabnamsplit <- all.vars(by)
   
  ISFACTOR<- sapply(data.frame(svydat.sub$variables[,names(svydat.sub$variables) %in% tabnamsplit]),is.factor)

  if(!all(ISFACTOR))
      {# wenn zumindest eine der by variablen KEIN Factor ist!
      
      tabnamsplitF <- paste0("factor(",tabnamsplit,")")
      crosstinp    <- paste(tabnamsplitF,collapse=":")
      perc1        <- svymean(x=  ~ eval(parse(text=crosstinp)) , design=svydat.sub, na.rm=TRUE)
      
      perc1 <- as.data.frame(perc1)
      rownames(perc1) <- gsub("^eval\\(parse\\(text = crosstinp\\)\\)","",rownames(perc1))
      colnames(perc1) <- c("Proportion", "SE Proportion")
      

    } else 
        {# wenn alle by variablen factors sind
          crosstinp   <- paste(tabnamsplit,collapse=":")
          perc1       <- svymean(x=  ~ eval(parse(text=crosstinp)) , design=svydat.sub, na.rm=TRUE)
          
          perc1 <- as.data.frame(perc1)
          rownames(perc1) <- gsub("^eval\\(parse\\(text = crosstinp\\)\\)","",rownames(perc1))
          colnames(perc1) <- c("Proportion", "SE Proportion")
      
        }
   
splitnames <- strsplit(rownames(perc1),":")
splitndf   <- data.frame(matrix(unlist(splitnames),ncol=length(splitnames[[1]]), byrow=TRUE))
colnames(splitndf) <- paste0("Group",1:length(all.vars(by)))
Percdat            <- data.frame(splitndf,perc1)
rownames(Percdat)  <- NULL

Percdat
  }
