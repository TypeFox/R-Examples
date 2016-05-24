`writeGLIMcmds` <-
function(dm,blnUndec,
   blnIntcovs,outFile,ncov,nintpars,intFile,nintcovs.out,glimCmdFile,
   covnames,covlevels,objnames,undecnames)



################################################################
# GLIM commands output
################################################################
{

    #env<-get("ENV",environment(patt.design))
    cat("exporting files for GLIM...(no output in R console)\n")
    flush.console()

    #env<-get("ENV",environment(patt.design))
    #for (i in 1:length(ls(env)))
    #  do.call("assign",list(ls(env)[i],get(ls(env)[i],env)))

    ## design file


    colnames(dm)[1]<-"!y"

    # no interaction effects
    if (!blnIntcovs) {
          write.table(dm, outFile, quote=FALSE, row.names=FALSE)
    # interaction effects
    } else {
          eint<-ncol(dm)-ncov            # last column  of interaction covs
          bint<-eint-nintpars+1          # first column of interaction covs

          write.table(dm[,-(bint:eint)], outFile, quote=FALSE, row.names=FALSE)

          # interaction design - transposed to allow for ARRAY in GLIM
          write(t(dm[,bint:eint]),intFile,ncolumns=nintcovs.out)
    }

    glim.sl<-dim(dm)[1]

    ## commands
    txt<-paste("$SL ",glim.sl,sep="")
    write(txt,file=glimCmdFile)

    if (blnIntcovs) {
         txt<-paste("$DATA y ",paste(colnames(dm)[-c(1,bint:eint)],collapse=" "),sep="")
    } else {
         txt<-paste("$DATA y ",paste(colnames(dm)[-1],collapse=" "),sep="")
    }
    write(txt,file=glimCmdFile,append=TRUE)

    txt<-paste("$DINP '",outFile,"' 132",sep="")
    write(txt,file=glimCmdFile,append=TRUE)

    if (blnIntcovs) {
         txt<-paste("$ARRAY IA ", paste(glim.sl,nintpars,sep = ",",collapse=" "),sep="")
         write(txt,file=glimCmdFile,append=TRUE)
         txt<-paste("$DATA IA ")
         write(txt,file=glimCmdFile,append=TRUE)
         txt<-paste("$DINP '",intFile,"' 132",sep="")
         write(txt,file=glimCmdFile,append=TRUE)
    }

    if (ncov>0) {
         txt<-paste("$FAC ")
         facs<-na.omit(data.frame(covnames,covlevels))  # for casewise in later version
         txt<-paste(txt, paste(facs[,1],facs[,2],sep = " ",collapse=" "),sep="")
         write(txt,file=glimCmdFile,append=TRUE)

         txt<-paste("$ELIM ",paste(covnames,collapse="*"),sep="")
         write(txt,file=glimCmdFile,append=TRUE)
    }

    txt<-paste("$LIST OBJ =",paste(objnames,collapse=","),sep=" ")
    write(txt,file=glimCmdFile,append=TRUE)


    if (blnUndec) {
         txt<-paste("$LIST U =",paste(undecnames,collapse=","),sep=" ")
         write(txt,file=glimCmdFile,append=TRUE)
    }

    write("$ERR P $YV y",file=glimCmdFile,append=TRUE)

    txt<-"$CYCLE 15 1 1e-5"
    write(txt,file=glimCmdFile,append=TRUE)

    txt<-"$FIT OBJ + U"
    write(txt,file=glimCmdFile,append=TRUE)

    write("$DISP E",file=glimCmdFile,append=TRUE)

    write("$RETURN",file=glimCmdFile,append=TRUE)
}
