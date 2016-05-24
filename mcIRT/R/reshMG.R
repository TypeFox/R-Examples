reshMG <-
function(da,items=0,groups=NA,correct,design="nodif",echo=TRUE,TYPE="NRM", paraM="bock")
{
  # d = data frame (including factors of reponses --> responses as integer beginning with 1 ending with m. contains the responses of the examinees in one of the m possible categories of item i
  #items = numerical vector which refers to the columns of the items
  # groups = numerical vector (of length = 1) which refers to the grouping variable
  # correct = vector of length(correct) == ncol(d), which contains the number of categories for each item which are counted as true - starts with 0!
  # design = this argument refers to the function which builds the 'group design matrix'. 
  
  #---------------------------
  # find items and groups ----
  #---------------------------

  if(length(groups)>1) stop("Only one group variable is permitted!")
  if(!TYPE %in% c("NRM","NLM")) stop("TYPE must be NRM or NLM!")
  if(!paraM %in% c("bock","01")) stop("paraM must be bock or 01!")
  
  if(is.na(groups) & any(items == 0))
    {
      gr <- factor(rep("A",nrow(da)))
      items <- 1:ncol(da)
    } else if(is.na(groups)){
        gr <- factor(rep("A",nrow(da)))
        } else {
            gr <- da[,groups] 
          }
        

  # doing controls - its necessary - trust me
  ctrlERG <- dctrl(da[,items],correct,items)
  
  if(ctrlERG$probl)
  {
    print(ctrlERG$problemlist)
    stop("we've got a problem")
  }
  
  if(TYPE == "NLM") 
    {
      
    da1 <- data.frame(mapply(function(rst,co) # shift categories depending on the correct option
    {
     if(co != 0)
       {
        cor0  <- rst == co
        chang <- rst == 0
        
        rst[cor0]  <- 0
        rst[chang] <- co
        return(rst)
       } else 
         {
          return(rst) 
         }
    }, rst = da[,items], co=correct, SIMPLIFY=FALSE))
    #
    d_new <- split(da1,gr) ## categories start with 0
    d_new <- lapply(d_new,as.matrix)

    ###  more controls
    egal <- dctrl2(d_new,da,items)
    
    recm <- lapply(levels(gr), function(LE)
    {
      
      dummycod <- mapply(function(well,item)
      {
        
        LEV <- sort(unique(da1[,1]))
        LAB <- LEV
        LAB[LEV == well] <- 0
        LAB <- c("cor",paste0("dt_",LAB[-1]))
        COL2 <- factor(da1[gr == LE,item],levels=LEV, labels=LAB)
        #COL2 <- factor(ifelse(da1[gr == LE,item] == well,"cor",paste("dt_",da1[gr == LE,item],sep="")))
        #COL2 <- factor(ifelse(da1[gr == LE,item] == 0,"cor",paste("dt_",da1[gr == LE,item],sep="")))
        caa  <- model.matrix(~ -1 + COL2, data = model.frame(~ -1 + COL2, na.action=na.pass))
        colnames(caa) <- gsub("COL2",paste("I",item,sep=""),colnames(caa))
        return(caa)
      },well=correct,item = 1:length(items),SIMPLIFY = FALSE)  
      #return(dummycod)
      
      do.call("cbind",dummycod) # new creates one matrix
    })
    
  aDD <- mapply(function(well,item) # fast original
    {
      LEV <- sort(unique(da1[,1]))
      LAB <- LEV
      LAB[LEV == well] <- 0
      LAB <- c("cor",paste0("dt_",LAB[-1]))
      COL2 <- factor(da1[,item],levels=LEV, labels=LAB)
      tabcat  <- table(COL2)
      categ   <- levels(COL2)
      anz_cat <- length(categ)
      addit   <- list(tabcat=tabcat,categ=categ,anz_cat=anz_cat)
      return(addit)
    },well=correct,item = 1:length(items),SIMPLIFY = FALSE) 
    
    
    
    
  } else {
    #
    d_new <- split(da[,items],gr) ## categories start with 0
    d_new <- lapply(d_new,as.matrix)

    ###  more controls
    egal <- dctrl2(d_new,da,items)
  # create dummy coded data
  recm <- lapply(levels(gr), function(LE)
  {
    
    dummycod <- mapply(function(well,item)
    {
      
      LEV <- sort(unique(da[,item]))
      #LEV <- c(sort(uq[uq != well]),uq[uq == well])
      LAB <- paste0("dt_",LEV)
      LAB[LEV == well] <- "cor"
      COL2 <- factor(da[gr == LE,item],levels=LEV,labels=LAB)
#       COL2 <- factor(ifelse(da[gr == LE,item] == well,"cor",paste("dt_",da[gr == LE,item],sep="")))
       caa  <- model.matrix(~ -1 + COL2, data = model.frame(~ -1 + COL2, na.action=na.pass))
      colnames(caa) <- gsub("COL2",paste("I",item,sep=""),colnames(caa))
      return(caa)
    },well=correct,item = items,SIMPLIFY = FALSE)  
    #return(dummycod)
    
   do.call("cbind",dummycod) # new creates one matrix
  })
  

    aDD <- mapply(function(well,item)
      {
      LEV <- sort(unique(da[,item]))
      LAB <- paste0("dt_",LEV)
      LAB[LEV == well] <- "cor"
      COL2 <- factor(da[,item],levels=LEV,labels=LAB)
        tabcat  <- table(COL2)
        categ   <- levels(COL2)
        anz_cat <- length(categ)
        addit   <- list(tabcat=tabcat,categ=categ,anz_cat=anz_cat)
        #
        return(addit)
      },well=correct,item = items,SIMPLIFY = FALSE)
      


    }
  ### d1uc new::::::
  ## anything which is.na() is converted to zero
 d1uc <-  lapply(recm, function(x){
  
            grerg <- apply(x, 2, function(spal){
              
              spal[is.na(spal)] <- 0
              spal
              
              })
            
            grerg
                    })
          
# -----------------------------  
  
  
  cat("data = reshaped\n")

  # descriptives
   coluN <- max(sapply(da[,items],function(x)length(table(x))))

  
#   if(is.na(groups))
#   {
#     absF1 <- sapply(da[,items],function(TA)table(TA,useNA="ifany"))
#     rownames(absF1)[1:coluN] <- paste("category",1:coluN,sep="")  
#   } else {
#     absF1 <- lapply(da[,items],function(TA)
#     {
#       tt1 <- table(TA,gr,useNA="ifany")
#       rownames(tt1)[1:coluN] <- paste("category",1:coluN,sep="")
#       tt1
#     })
#   }
  

    absF1 <- lapply(da[,items],function(TA)
    {
      tt1 <- table(TA,gr,useNA="ifany")
      rownames(tt1)[1:nrow(tt1)] <- paste0("category",1:nrow(tt1))
      tt1
    })





  # controlling the input design !!
  if(is.list(design))
  {
    
    if(!all(sapply(design,ncol) == length(items)))(stop("The number of items in the design != number of items in the dataset.\n"))
    
    design <- ctrl_design(design=design,aDD=aDD,gr=gr,TYPE=TYPE)  
  }
  
  
  # create Q matrix
  if(paraM == "bock")
  {
    gdema <- grDMb(aDD,gr,design,TYPE=TYPE)
    
  } else if(paraM == "01")
  {
    
    gdema <- grDM(aDD,gr,design,TYPE=TYPE)
    
  }
  
  attr(gdema,"paraM") <- paraM
  attr(d_new,"correct") <- correct
  
  if(echo){print(absF1)}
  
  cat("group information added \n")
  reshret <- list(recm=recm,aDD=aDD,d=d_new,gr=gr,Qmat=gdema,d1uc=d1uc,design=design)

 
  # classify
  if(TYPE=="NRM")
      {
        class(reshret) <- "reshNRM" 
      } else if(TYPE=="NLM") {
                              class(reshret) <- "reshNLM"  
                             }
      

##### WARNING ##############

testN <- lapply(d_new,function(x1)
{
  apply(x1,2,function(aa4)
  {
    ct <- table(aa4)[table(aa4) < 10]
    ctn <- paste("cat",names(ct),"with",ct,"obs")
    if(length(ct) == 0) NA
    else ctn
  })
})

if(!all(sapply(testN,function(dk) all(is.na(dk)))))
{
  pit <- lapply(testN,function(dk) dk[!is.na(dk)])
  warning("There are categories with a small number of observations.\n",immediate. = TRUE)
  print(pit)
  cat("\n")
}


  return(reshret)
}
