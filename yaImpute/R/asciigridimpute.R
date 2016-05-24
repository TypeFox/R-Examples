
AsciiGridPredict <- 
function(object,xfiles,outfiles,xtypes=NULL,lon=NULL,
         lat=NULL,rows=NULL,cols=NULL,
         nodata=NULL,myPredFunc=NULL,...)
{
   if (missing(xfiles)   || is.null(xfiles))   stop ("xfiles required")
   if (missing(outfiles) || is.null(outfiles)) stop ("outfiles required")
   if (!inherits(outfiles,"list")) outfiles=as.list(outfiles)
   if (is.null(names(outfiles)) && length(outfiles) == 1) names(outfiles) = "predict"
   if (is.null(names(xfiles))) stop ("xfiles elements must be named")
   if (is.null(names(outfiles))) stop ("outfiles elements must be named")

   return (
      AsciiGridImpute(object,xfiles,outfiles,xtypes=xtypes,lon=lon,
                      lat=lat,rows=rows,cols=cols,
                      nodata=nodata,myPredFunc=myPredFunc,...)
          )
}


AsciiGridImpute <-
function(object,xfiles,outfiles,xtypes=NULL,ancillaryData=NULL,
         ann=NULL,lon=NULL,lat=NULL,rows=NULL,cols=NULL,nodata=NULL,
         myPredFunc=NULL,...)
{
   if (missing(xfiles)   || is.null(xfiles))   stop ("xfiles required")
   if (missing(outfiles) || is.null(outfiles)) stop ("outfiles required")
   if (is.null(names(xfiles))) stop ("xfiles elements must be named")
   if (is.null(names(outfiles))) stop ("outfiles elements must be named")

#  outLegend is a list of factor levels and their index values for output
#  inLegend is a the same idea but for input.

   outLegend=NULL
   inLegend=NULL

#  nasum is a matrix of the number of NAs generated for each row

   nasum=NULL

#  make sure there is a type for every xfile

   if (is.null(xtypes))
   {
      xtypes=xfiles
      xtypes[]="numeric"
   }
   else
   {
      tmp=xtypes
      xtypes=xfiles
      xtypes[]=NA
      xtypes[names(tmp)]=tmp
      xtypes[is.na(xtypes)]="numeric"
   }

#  there needs to be an input file for every xvar in the imputation.

   if (is.null(object) || !inherits(object,"yai"))
   {
      have=names(xfiles)
   }
   else
   {
      have=intersect(xvars(object),names(xfiles))
      if (length(have) != length(xvars(object)))
      {
         lout = if (length(have)==0) xvars(object) else 
                setdiff(xvars(object),have)
         stop(paste("required maps are missing for variables:",
              paste(lout ,collapse=", ")))
      }
      # trim the list of xfiles to those needed.
      xfiles = xfiles[match(xvars(object),names(xfiles))]
   }

#  deal with ancillaryData and build allY

   allY = 0
   if (!is.null(ancillaryData) && inherits(object,"yai"))
   {
      if (length(intersect(class(ancillaryData),c("matrix","data.frame"))) > 0 &&
          nrow(ancillaryData[rownames(object$yRefs),,FALSE]) == nrow(object$yRefs) &&
          length(intersect(rownames(ancillaryData),rownames(object$yRefs))) == nrow(object$yRefs) )
      {
         toKeep = intersect(union(colnames(ancillaryData),colnames(object$yRefs)),names(outfiles) )
         if (length(toKeep) == 0) allY = NULL
         else
         {
            fromAn=intersect(toKeep,colnames(ancillaryData))
            fromRe=setdiff(intersect(toKeep,colnames(object$yRefs)),fromAn)
            if (length(fromAn)>0 && length(fromRe)>0)
               allY = data.frame(cbind(ancillaryData[rownames(object$yRefs),],object$yRefs)[,toKeep],
               row.names=rownames(object$yRefs))
            else if (length(fromAn)>0) allY = data.frame(ancillaryData[rownames(object$yRefs),toKeep],
               row.names=rownames(object$yRefs))
            else if (length(fromRe)>0) allY = data.frame(object$yRefs[,toKeep],
               row.names=rownames(object$yRefs))
            colnames(allY)=toKeep
         }
      }
      if (is.null(names(allY))) stop ("ancillaryData can not be used because no variables match the names in outfiles")
   }
   if (is.null(colnames(allY)) && inherits(object,"yai"))
   {
      toKeep = intersect(colnames(object$yRefs),names(outfiles))
      if (length(toKeep) == 0) allY = NULL
      else
      {
        allY = data.frame(object$yRefs[,toKeep],row.names=rownames(object$yRefs))
        colnames(allY) = toKeep
      }
   }

#  if using yai, deal with ann

   if (inherits(object,"yai") && is.null(ann)) ann=object$ann

#  set some flags used below

   predYes = length(intersect(names(outfiles),c("predict" ))) == 1
   distYes = length(intersect(names(outfiles),c("distance"))) == 1  && inherits(object,"yai")
   useidYes= length(intersect(names(outfiles),c("useid"   ))) == 1  && inherits(object,"yai")

   sumIlls=NULL

#  make a list of input file handles and open the files.

   infh = vector("list",length=length(xfiles))
   names(infh)=names(xfiles)
   for (i in 1:length(xfiles)) infh[[i]]=file(xfiles[[i]],open="rt")
   on.exit(lapply(infh,close))

#  get and check headers from each input file

   header=NULL
   postWarn=TRUE
   for (i in 1:length(infh))
   {
      newhead = readLines(infh[[i]],n=6)
      if (is.null(header)) header=newhead
      else
      {
         if (!identical(newhead,header))
         {
            cat ("Map headers are not equal\nHeader from file: ",xfiles[[i-1]],"\n")
            print (header)
            cat ("\nHeader from file: ",xfiles[[i]],"\n")
            print (newhead)
            flush.console()
            if (postWarn) warning ("map headers don't match.")
            postWarn=FALSE
         }
         header=newhead
      }
   }

#  write the "common" header to all the output files.

   getVal=function(header,tok)
   {
      for (lin in header)
      {
         l = scan (text=lin,what="character",quiet=TRUE)
         if (toupper(l[1]) == tok) return (as.numeric(l[2]))
      }
      return (NA)
   }

   nc  = getVal(header,"NCOLS")
   nr  = getVal(header,"NROWS")
   xllc= getVal(header,"XLLCORNER")
   yllc= getVal(header,"YLLCORNER")
   csz = getVal(header,"CELLSIZE")
   nodv= getVal(header,"NODATA_VALUE")
   if (any(unlist(lapply(c(nc,nr,xllc,yllc,csz,nodv),is.na)))) stop ("header error in input maps")

   if (!is.null(lon))
   {
     cols=c(0,0)
     lon=sort(lon)
     cols[1]=round((lon[1]-xllc)/csz)+1
     cols[2]=round((lon[2]-xllc)/csz)
   }
   if (!is.null(lat))
   {
     rows=c(0,0)
     lat=sort(lat)
     rows[2]=nr-round((lat[1]-yllc)/csz)
     rows[1]=nr-round((lat[2]-yllc)/csz)+1
   }

   if (is.null(rows) && is.null(cols) && is.null(nodata)) #header does not change
   {
      for (i in 1:length(outfiles)) cat (header,file=outfiles[[i]],sep="\n")
      newnr = nr
      newnc = nc
      nodata= nodv
      rows=c(1,nr)
      cols=c(1,nc)
   }
   else #header changes
   {
 
      if (is.null(nodata)) nodata=nodv

      if (is.null(rows))   rows=c(1,nr)
      if (rows[1]<1)       rows[1]=1
      if (rows[2]>nr)      rows[2]=nr
      if (rows[1]>rows[2]) rows[1]=rows[2]

      if (is.null(cols))   cols=c(1,nc)
      if (cols[1]<1)       cols[1]=1
      if (cols[2]>nc)      cols[2]=nc
      if (cols[1]>cols[2]) cols[1]=cols[2]
      
      newnr = rows[2]-rows[1]+1
      newnc = cols[2]-cols[1]+1

      if (rows[2] != nr) yllc = yllc+(csz*(nr-rows[2]))
      if (cols[1] != 1)  xllc = xllc+(csz*(cols[1]-1))

      for (i in 1:length(outfiles))
      {
         cat("NCOLS         ",as.character(newnc), "\n",file=outfiles[[i]],sep="")
         cat("NROWS         ",as.character(newnr), "\n",file=outfiles[[i]],sep="",append=TRUE)
         cat("XLLCORNER     ",as.character(xllc),  "\n",file=outfiles[[i]],sep="",append=TRUE)
         cat("YLLCORNER     ",as.character(yllc),  "\n",file=outfiles[[i]],sep="",append=TRUE)
         cat("CELLSIZE      ",as.character(csz ),  "\n",file=outfiles[[i]],sep="",append=TRUE)
         cat("NODATA_VALUE  ",as.character(nodata),"\n",file=outfiles[[i]],sep="",append=TRUE)
      } 
   }

   # set up the xlevels. In randomForest version >= 4.5-20, the xlevels
   # are stored in the forest.

   xlevels = object$xlevels
   if (is.null(xlevels) && inherits(object,"randomForest")) xlevels=object$forest$xlevels
   if (!is.null(xlevels))
   {
      # we need just the variable names...
      names(xlevels) = lapply(names(xlevels), function (x) unlist(all.vars(formula(paste("~",x)))[1]))

      if (length(xlevels)>0) for (i in names(xlevels)) if (is.numeric(xlevels[[i]]) &&
                                  length(xlevels[[i]]) == 1) xlevels[[i]] = NULL
      if (length(xlevels) == 0) xlevels=NULL
      else
      {
      	inLegend=lapply(xlevels,data.frame)
      	inLegend=unionDataJoin(inLegend)
      	names(inLegend)=names(xlevels)
      }
   }
   nskip=rows[1]-1

   dpr = max(newnr %/% 100,1)

   if (nskip > 0) cat("Rows to skip: ",nskip," ")
   cat("Rows per dot: ",dpr," Rows to do:",newnr,"\nToDo: ")

   for (ir in 1:floor(newnr/dpr)) cat (".")
   cat ("\nDone: ");flush.console()

   nodout=suppressWarnings(as.numeric(nodata))

   ircur=0
   for (ir in rows[1]:rows[2])
   {
      indata = vector("list",length=length(xfiles))
      names(indata)=names(infh)
      for (i in 1:length(infh))
      {
         indata[[i]]=scan(infh[[i]],nlines=1,what=vector(mode=xtypes[[i]],length=0),
                     skip=nskip,quiet=TRUE)
         nod = indata[[i]] == nodv
         if (any(nod)) indata[[i]][nod] = NA
      }
      nskip=0
      rowLens = unlist(lapply(indata,length))
      if (any(rowLens[1] != rowLens)) 
      {
         cat ("Row lengths for row ",ir,"\n")
         print(rowLens)
         flush.console()
         stop ("Unequal row lengths.")
      }
      if (newnc == nc) newdata=data.frame(indata)
      else             newdata=data.frame(indata)[cols[1]:cols[2],,FALSE]
      width=floor(log(nrow(newdata),10))+1
      rownames(newdata)=paste(ir,formatC(1:nrow(newdata),width=width,flag="0"),sep="x")
      origRowNames=rownames(newdata)

      newdata=na.omit(newdata)
      omitted=origRowNames[as.vector(attr(newdata,"na.action"))]
      moreOrigRowNames=NULL      
      if (!is.null(xlevels) && length(omitted)<length(origRowNames))
      {
         moreOrigRowNames=rownames(newdata)
         factorMatch = get("factorMatch",asNamespace("yaImpute"))
         newdata=factorMatch(newdata,xlevels)
         ills = attr(newdata,"illegalLevelCounts")
         if (class(ills)=="list")
         {
            if (is.null(sumIlls))
            {
               sumIlls = ills
               warning ("NA's generated due to illegal level(s).")
            }
            else
            {
               addIllegalLevels = get("addIllegalLevels",asNamespace("yaImpute"))
               sumIlls = addIllegalLevels(sumIlls,ills)
            }
            omitted=c(omitted,moreOrigRowNames[as.vector(attr(newdata,"na.action"))])
         }
      }
      # tag the vector so newtargets() will not duplicate
      # the creation of this attribute data
      else  attr(newdata,"illegalLevelCounts")=0   
                                                   

      if (length(omitted)==length(origRowNames)) # all missing.
      {
         outdata=data.frame(matrix(nodout,length(omitted),length(outfiles)),
                            row.names=omitted)
         names(outdata)=names(outfiles)
      }
      else
      {
         if (!is.null(myPredFunc))
         {
            outdata=myPredFunc(object,newdata,...)
            if (!inherits(outdata,"data.frame"))
            {
               cns = colnames(outdata)
               outdata=data.frame(predict=outdata,row.names=rownames(newdata))
               if (!is.null(cns)) colnames(outdata) = cns
            } else rownames(outdata)=rownames(newdata)               
         }
         else if (is.null(object))
         {
            outdata=newdata
         }
         else if (inherits(object,"yai"))
         {
            outdata = NULL
            saveNames=rownames(newdata)
            rownames(newdata)=paste("m",as.character(1:nrow(newdata)),sep="!")
            new = newtargets(object,newdata,k=NULL,ann=ann)
            if (!is.null(allY)) outdata = impute(new,ancillaryData=allY,observed=FALSE,...)
            rownames(outdata)=saveNames
            if (distYes)  dists = data.frame(distance=new$neiDstTrgs[,1],
                                  row.names=rownames(newdata))
            else          dists = NULL
            if (useidYes) useIds= data.frame(useid=match(new$neiIdsTrgs[,1],
                                  rownames(object$xRefs)),row.names=rownames(newdata))
            else          useIds= NULL
            if (!is.null(outdata) && !is.null(dists) ) outdata=cbind(outdata,dists)
            else if (is.null(outdata)) outdata=dists
            if (!is.null(outdata) && !is.null(useIds)) outdata=cbind(outdata,useIds)
            else if (is.null(outdata)) outdata=useIds
         }
         else
         {
            outdata=predict(object,newdata,...)
            if (!inherits(outdata,"data.frame"))
            {
               cns = colnames(outdata)
               outdata=data.frame(predict=outdata,row.names=rownames(newdata))
               if (!is.null(cns)) colnames(outdata) = cns
            } else rownames(outdata)=rownames(newdata)               
         }
         if (is.null(outLegend))
         {
           outLegend=vector("list",ncol(outdata))
           names(outLegend)=names(outdata)
           for (n in names(outLegend)) outLegend[[n]]=if (is.factor(outdata[,n])) 
             levels(outdata[,n]) else NULL
         } else {
           for (n in names(outLegend))
           {
             if (is.factor(outdata[,n]))
             { 
               for (lev in levels(outdata[,n]))
               {
                 if (length(grep(lev,outLegend[[n]],fixed=TRUE)) == 0)
                   outLegend[[n]] = c(outLegend[[n]],lev)
               }
             }
           }
         }
         #convert factors to numbers that match the outLegend
         for (n in colnames(outdata)) if (is.factor(outdata[,n]))
           outdata[,n] <- match(levels(outdata[,n])[outdata[,n]],outLegend[[n]])
         if (nrow(outdata) != nrow(newdata))
         {
            cat ("First six lines non-missing predictions for row ",ir,"\n")
            print(head(outdata))
            cat ("First six lines of non-missing xfiles data for row ",ir,"\n")
            head(head(newdata))
            flush.console()
            stop ("Unexpected results for row = ",ir)
         }
         outrs = nrow(outdata) # the predict might send NA's, change them to no data
         outtmp = na.omit(outdata)
         if (outrs > nrow(outtmp))
         {
            nasum = if (is.null(nasum)) c(ir,outrs-nrow(outtmp)) else
                            rbind(nasum,c(ir,outrs-nrow(outtmp)))
            for (icol in 1:ncol(outdata)) outdata[is.na(outdata[,icol]),icol]=nodout
         }
         if (length(omitted)>0)
         {
            # add omitted observations back into data frame in the proper location

            more = data.frame(matrix(nodout,length(omitted),length(names(outdata))),
                              row.names=omitted)
            names(more)=names(outdata)
            for (i in 1:ncol(more)) if (is.factor(outdata[,i])) more[,i]=as.factor(more[,i])
            outdata =  rbind(outdata,more)
            outdata = outdata[sort(rownames(outdata),index.return = TRUE)$ix,,FALSE]
         }
      }
      for (i in 1:length(outfiles))
      {
         vname=names(outfiles)[i]
         if (length(intersect(names(outdata),vname))==0)
         {
         	  cat ("\nFirst six lines of predicted data for map row: ",ir,"\n")
         	  print(head(outdata))
         	  flush.console()
         	  stop (vname," is not present in the predicted data")
         }
         cat (outdata[,vname],"\n",file=outfiles[[i]],append=TRUE)
      }

      ircur=ircur+1
      if (ircur>=dpr)
      {
         ircur=0
         cat (".");flush.console()
      }
   }
   cat ("\n");flush.console()

   for (n in names(outLegend))
   {
      outLegend[[n]]=as.data.frame(outLegend[[n]],stringsAsFactors=FALSE)
      names(outLegend[[n]])=n
   }
   if (! is.null(nasum))
   {
      nasum=data.frame(maprow=nasum[,1],count=nasum[,2])
      cat ("Summary of unexpected NA values generated from predict:")
      print (nasum)
      warning("Unexpected NA values generated")
   }
   if (length(outLegend)>0)
   {
     outLegend=unionDataJoin(outLegend)
     cat ("Legend of levels in output grids:\n")
     print (outLegend)
   }
   else outLegend=NULL

   if (! is.null(inLegend))
   {
     cat ("Legend of levels in input grids (assumed):\n")
     print (inLegend)
   }
   if (!is.null(sumIlls))
   {
      cat ("Factors with levels found in input maps but not present in data used to fit the model\n")
      cat ("(number of map cells; counts do not include cells coded NODATA on other maps):\n")
      if (class(sumIlls[[1]]) == "table") print (sumIlls)
      else
      {
         for (a in names(sumIlls)) 
         {
           sumIlls[[a]]=as.data.frame(sumIlls[[a]])
           colnames(sumIlls[[a]])=a
         }
         sumIlls = unionDataJoin(sumIlls)
         print (sumIlls)
         cat ("\n")
      }
   }
   invisible(list(unexpectedNAs=nasum,illegalLevels=sumIlls,outputLegend=outLegend,inputLegend=inLegend))
}
