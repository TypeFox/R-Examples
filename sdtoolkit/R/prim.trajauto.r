
`prim.trajauto` <-
function(x, y, box.init=NULL, peel.alpha=0.05, paste.alpha=0.01,
                 mass.min=0.05, threshold=0, pasting=TRUE, verbose=FALSE,
                 threshold.type=1, paste.all=FALSE,coverage=TRUE,showbounds=TRUE,style="ineq",npts=NA,ninter=NA,nbump=10,repro=TRUE,dfrac=.5){

options(digits=4)

  d <- ncol(x)
  n <- nrow(x)
  k.max <- ceiling(1/mass.min)
  num.boxes <- k.max

  freqmat <- NULL
  ##if (is.vector(x)) x <- as.matrix(t(x))
  y.mean <- mean(y)
  mass.init <- length(y)/n

  if (is.null(box.init))
  {
    box.init <- apply(x, 2, range)
    box.diff <- box.init[2,] - box.init[1,]
    box.init[1,] <- box.init[1,] - 10*paste.alpha*box.diff
    box.init[2,] <- box.init[2,] + 10*paste.alpha*box.diff
  }

  ## find first box
  #k <- 1

  boxseq <- find.traj(x=x, y=y, box=box.init, peel.alpha=peel.alpha,
                   paste.alpha=paste.alpha, mass.min=mass.min,
                   threshold=min(y)-0.1*abs(min(y)), d=d, n=n, pasting=pasting, verbose=verbose, paste.all=paste.all)

  satisfied=FALSE

  while(!satisfied){

    trajinf <- traj.info(x=x, y=y, box.seq=boxseq, npts=npts,ninter=ninter)
    
    
#    cat("You should soon see a peeling trajectory to your right.","\n")
#    cat("Changes in the color of points represent a change in the number of","\n"
#    ,"dimensions that have been restricted.","\n")
    
    
#    flush.console()

#    margtrajs <- ranker(1:length(boxseq[[1]]),x,y,trajinf,npts,ninter)

#    checpts <- trajplot(trajinf,coverage=coverage,margtrajs=margtrajs)  #addedmargtrajs as argument requires modified trajplot
    
    #cat(trajinf$y.mean)
    checpts <- which(trajinf$y.mean==max(trajinf$y.mean))  #added for auto
    
    if(max(trajinf$y.mean)!=1){    #also added for auto
      cat("Weird, wasn't able to get 100% density","\n")
    }

    ####ADDED TO DIAGNOSE BUG - PRINTS OUT ORDER (HOPEFULLY!) OF POINT SELECTION
    
#    print("Here are the points selected on the trajectory:")
#    print(checpts)
#    cat("\n")

 #   do dimranking on these boxes
 #   if(length(checpts)==0){
           
  #    cat("(No points graphically identified, but that's ok if you meant to do that...)","\n","\n")
    
 #   }else{
    
      morestats <- ranker(checpts,x,y,trajinf,npts,ninter)
      
      pvallist <- pvallister(checpts,x,y,trajinf)
    
 #     traj.stats(trajinf,morestats,pvallist,checpts,showbounds=showbounds,style=style)

 #  }
                            
 #   boxind <- readline(cat("Which box would you like? (or enter \"done\" to just keep past box sequence)","\n",
 #   "Or enter \"more\" to inspect more boxes from this trajectory","\n"))
     
  boxind <- checpts #added for auto
  
    if(boxind=="more"){}
    else if(boxind=="done"){thebox<-boxind; satisfied=TRUE}
    else {
      
      satisfied=TRUE
      
      boxind <- as.numeric(boxind)
      
      msin <- c(1:length(checpts))[checpts==boxind]
      
      if (length(msin)==0){
        cat("Warning, using box not specifically identified on peeling trajectory.","\n","\n")
                morestats <- ranker(boxind,x,y,trajinf,npts,ninter)
                pvallist <- pvallister(boxind,x,y,trajinf)
        msin <- 1
      }
      
      #EDIT BOX:
      
  #    cat("Here is the box you selected:","\n","\n")
 
 #clumsy bug fix:
#      tempmorestats <- list(morestats[[1]][msin],morestats[[2]][msin],morestats[[3]][msin])
      
   #   curmats <- traj.stats(trajinf,list(morestats[[msin]]),pvallist,boxind,showbounds=showbounds,style=style)

      #traj.stats(trajinf,morestats,boxind,showbounds=showbounds,style=style)


#### PLOTTING THE ACTUAL BOX IF DESIRED

    #  psatisfied <- FALSE

      #CHECK FOR BOX BEING AT LEAST TWO DIMENSIONS
      
     # if(sum(trajinf$dimlist[[boxind]]$either)>1){      
        
     # plotbox <- readline(cat("Would you like to see a plot of this box over a scatter plot of the dataset?","\n",
     #   "(High value points will be filled in, low value points will be hollow)","\n","Enter 'y' or 'n'.","\n"))
        
      #cat("\n")
      
      
      #if(plotbox=="y"){
        
      #  while(!psatisfied){
        
        
          #CONVERT BOX to 'old' form so box can plot it
      #    boxy <- list(box=trajinf$box[[boxind]],dimlist=trajinf$dimlist[[boxind]])
      #    ofbox <- boxconverter(boxy)
          
        
          #EXTRACT TWO MOST IMPORTANT DIMS
     #     d1name <- as.character(curmats[[2]][1,1])
     #    d1 <- curmats[[2]][1,4]
          
     #     if(curmats[[2]][1,1]==curmats[[2]][2,1]){
     #       d2name <- as.character(curmats[[2]][3,1])
     #       d2 <- curmats[[2]][3,4]
     #     }else{
     #       d2name <- as.character(curmats[[2]][2,1])
     #       d2 <- curmats[[2]][2,4]
     #     }
          
        
     #     cat("Enter the row numbers for two dimensions (in the left most column, 
     # not actual name) you would like to use for the x and y axis separated","\n",
     #     "by a comma, or enter 'd' to accept the default, which is the two highest","\n","ranked variables.")
     #     xandy <- readline(cat("\n","  (In this case",d1name,"and",d2name,".)","\n"))
          
     #     if(xandy!="d"){
            
            #EXTRACT x and y
      #      
#            ds <- strsplit(xandy,",")[[1]]
#            
#            d1t <- as.numeric(ds[1])
#            d2t <- as.numeric(ds[2])
#            
#            d1 <- curmats[[2]][d1t,4]
#            d2 <- curmats[[2]][d2t,4]
#            
#            d1name <- as.character(curmats[[2]][d1t,1])
#            d2name <- as.character(curmats[[2]][d2t,1])
#            
#            #OLD and I think bad way:
#            #d1name <- curmats[[2]][,1][which(curmats[[2]][,1]==d1)]
#            #d2name <- curmats[[2]][,1][which(curmats[[2]][,1]==d2)]
#            
#          }
#          
#          #OPEN NEW DEVICE SO TRAJ STAYS OPEN  - but keep first box
#          #IF, stuff...
#      
#          #PLOT POINTS
#      
#          options("device")$device()
#          colptplot(cbind(x,y),xdim=d1,ydim=d2,outdim=(ncol(x)+1),lowcol="transparent",hicol="black",xname=d1name,yname=d2name)
#      
#      
#          #PLOT BOX
#          pbox(sdobj=ofbox,xdim=d1,ydim=d2,boxnum=NA,fromtype="oldbox",lwd=3,gborder="blue",mdborder="red",col=NA)
#
#        
#          bringToTop(-1)
          #PROMPT FOR PLOT OF DIFFERENT VARIABLES - loop back to TOP
#          redo <- readline(cat("Would you like to plot different variables? ('y#' or 'n')","\n"))
#          if (redo == "y"){
#            psatisfied <- FALSE
#          } else{
#            psatisfied <- TRUE
#          }
          
          #IF NOT FIRST COVERING - offer to plot an earlier box in the sequence?
         
        #}
      
      #}
      
      #}
      
#      cat("\n")
      
      if(repro){
        cat("Calculating reproducibility statistics...","\n")
        cat("You can turn this time consuming step off in the future by adding the","\n",
        "argument 'repro=FALSE' when you call sdprim.","\n")
        cat("\n")
        flush.console()      
  ######      #Section for automatically assessing reproducibility
        #IF testreproducibility
        maincov <- trajinf$marcoverage[[boxind]]
        maindens <- trajinf$y.mean[[boxind]]
        
        #print(maincov)
        #print(maindens)
        
        boxseqlist <- list()
        
        # #matrix indicating which dimensions were restricted, by resampling 
        cdimsmat <- matrix(nrow=nbump,ncol=ncol(x)) #to match coverage
        ddimsmat <- matrix(nrow=nbump,ncol=ncol(x)) #to match density
        
        for (b in 1:nbump){
          
  #        if(b==5){debug(find.traj)}
          
          curs <- sample(nrow(x),dfrac*nrow(x))  #resample ['current sample']
          
          boxseqlist[[b]] <- find.traj(x=x[curs,], y=y[curs], box=box.init, peel.alpha=peel.alpha,
                     paste.alpha=paste.alpha, mass.min=mass.min,                  
                     threshold=min(y)-0.1*abs(min(y)), d=d, n=n, pasting=pasting, verbose=verbose, paste.all=paste.all)
          
          resamptrajinf <- traj.info(x=x[curs,], y=y[curs], box.seq=boxseqlist[[b]], npts=dfrac*nrow(x),ninter=sum(y[curs])) 
        
          #FIND BOX WITH CLOSEST COVERAGE
          
          cdiffs <- abs(resamptrajinf$marcoverage-maincov)
          closestcov <- max(which(cdiffs==min(cdiffs)))
          
          cdimsmat[b,] <- resamptrajinf$dimlist[[closestcov]]$either
          
          #FIND BOX WITH CLOSEST DENSITY
          ddiffs <- abs(resamptrajinf$y.mean-maindens)
          closestcov <- max(which(ddiffs==min(ddiffs)))
          
          ddimsmat[b,] <- resamptrajinf$dimlist[[closestcov]]$either
          
          #FOR EACH, FIND WHICH DIMENSIONS WERE RESTRICTED
          
          #REPORT HISTOGRAM/TABLE
        
        }
        
        covf <- apply(cdimsmat,2,sum)/nbump
        denf <- apply(ddimsmat,2,sum)/nbump
        
        freqmat <- rbind(covf,denf)
        colnames(freqmat) <- colnames(x)
        rownames(freqmat) <- c("For eq cov","For eq dens")
        
        cat("\n")
        cat("Here are the frequencies with which dimensions were restricted","\n")
        cat("to get equivalent coverage/density on resamplings of the dataset.","\n")
        cat("This used",nbump,"resamplings with",dfrac,"of the data sampled each time.") 
        cat("\n")
        cat("\n")
        print(t(freqmat))
        
      }  #END reproducibility stat section
  
      #cat("\n")
      
      #traj.stats(trajinf,list(morestats[[msin]]),pvallist,boxind,showbounds=showbounds,style=style)
       
      #rvar <- readline(cat("Would you like to remove any variables? (\"y\" or \"n\")?","\n","\n"))

      rvar <- "n"
            
      if(rvar=="y"){
        typr <- readline(cat("Would you like to specify total number of dimensions to remove?","\n","  (will be removed beginning with least important)","\n", "OR specify individual restrictions to eliminate?","\n","(Enter \"t\" for total and \"ir\" for individual restrictions)","\n","\n")) 
         
        if(typr=="t"){
         
          totakeout <- as.numeric(readline(cat("How many?","\n")))
          #identify which dimensions it's actually talking about:
          
          dsr <- nrow(morestats[[msin]])
          
          dimsout <- morestats[[msin]][dsr:(dsr-totakeout+1),1] 
          
                                 
          newbox <- trajinf$box[[boxind]]
            
          newbox[1,dimsout] <- -Inf
          newbox[2,dimsout] <- Inf
          
          boxseq <- list(box=list(),num.class=vector(length=1))
          boxseq$box[[1]] <- newbox
          boxseq$num.class[1] <- 1
          
          #Get stats on the new box:  [should just be the same as remove var stats at present]
          trajinf <- traj.info(x=x, y=y, box.seq=boxseq, npts=npts,ninter=ninter)
          morestats <- ranker(checpts=1,x,y,trajinf,npts,ninter)
          pvallist <- pvallister(checpts=1,x,y,trajinf)
        
         
          #Write out to confirm it's ok:
          cat("You chose to remove the following variables:","\n")
          print(colnames(x)[dimsout])
          cat("Here is the box that results:","\n")
          
          traj.stats(trajinf,morestats,pvallist,cboxes=1,showbounds=showbounds,style=style)
      
          boxind <- 1
          msin   <- 1     
        
        } else if(typr=="ir"){
          
          remrs <- readline(cat("Enter one or more restrictions to remove, separated by commas","\n","(Specify using the rownumbers in the leftmost column)","\n","\n","Note that, unlike the option for removing least important dimensions,","\n","if you want to remove a dimension that is restricted from both above and below","\n","you will need to give the number for both restrictions.","\n","\n"))
         
          remrs <- as.numeric(as.vector(strsplit(remrs,",")[[1]]))
          
#          osdimsout <- curmats[[2]][remrs,4] #os = onesided - then need to figure out whether top or bottom

          #NEEDS TO BE FIXED AFTER REINTRODUCING CURMATS!
          ldimsout <- NA # curmats[[2]][remrs,4][curmats[[2]][remrs,2]==" > "]
          udimsout <- NA # curmats[[2]][remrs,4][curmats[[2]][remrs,2]==" < "]

          newbox <- trajinf$box[[boxind]]
          
          newbox[1,ldimsout] <- -Inf
          newbox[2,udimsout] <- Inf         
 
          dirrem    <- c(rep("lower",length(ldimsout)),rep("upper",length(udimsout)))
          
          uinl <- udimsout %in% ldimsout
        
          if (sum(uinl)>0){
            
            for (di in udimsout[uinl]){
              dirrem[which(ldimsout==di)] <- "both"
            }
           
            uiis <- c(1:length(udimsout))[uinl]
            udimsout <-  udimsout[-uiis]
            dirrem <- dirrem[-(length(ldimsout)+uiis)]

          } 
          
          remdnames <- colnames(x)[c(ldimsout,udimsout)]

         
        boxseq <- list(box=list(),num.class=vector(length=1))
        boxseq$box[[1]] <- newbox
        boxseq$num.class[1] <- 1
        
        #Get stats on the new box:  [should just be the same as remove var stats at present]
        trajinf <- traj.info(x=x, y=y, box.seq=boxseq, npts=npts,ninter=ninter)
        morestats <- ranker(checpts=1,x,y,trajinf,npts,ninter)
        pvallist <- pvallister(checpts=1,x,y,trajinf)
       
#        #Write out to confirm it's ok:
        cat("You chose to remove the following dimension restrictions:","\n","\n")

#        column one is name, two is "upper", "lower", "both"

        remmat <- cbind(remdnames,dirrem)
        colnames(remmat) <- c("Variable","Restriction") 
        rownames(remmat) <- c(1:nrow(remmat))

        print(remmat,quote=FALSE)

        cat("\n","Here is the box that results:","\n")
        
        traj.stats(trajinf,morestats,pvallist,cboxes=1,showbounds=showbounds,style=style)
    
        boxind <- 1
        msin   <- 1     


          #OTHER STUFF
        } else{ #SOMEWARNING THAT LOOPS BACK TO THE TOP
        }      
                           
      
      } #END remove variable section                                     
      
      
      #CONVERT THE CURRENT BOX TO LIST FORM FOR USE IN BOX PLOTTING

      boxy <- list(box=trajinf$box[[boxind]],dimlist=trajinf$dimlist[[boxind]])
      ofbox <- boxconverter(boxy)  #"old form" box
                
      
      thebox <- list(y.mean = trajinf$y.mean[[boxind]], box = trajinf$box[[boxind]],
                   mass = trajinf$mass[[boxind]], dimlist = trajinf$dimlist[[boxind]], morestats=morestats[[msin]],
                   marcoverage = trajinf$marcoverage[[boxind]], relcoverage = trajinf$relcoverage, pvallist = pvallist[[boxind]],
                   freqmat=freqmat, index=boxind, lbox=ofbox)
      }

    }

    return(thebox)

  }

