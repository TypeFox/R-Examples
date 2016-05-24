"rserror" <-
function(err)
{
    bin2int<-function(x){
      r<-vector(mode="integer",length=0)
      while(x>1){
        r<-c(x%%2,r)
        x<-floor(x/2)
      }
      r<-c(1,r)
      rev(r)
    }
    errortxt<-c("\tn < 0 or n > 4096\n",
                "\tk < 0 or k > 128\n",
                "\tn_eff < 0 or n_eff > 8191\n",
                "\tburn_in < 0\n",
                "\tstep <= 0\n",
                "\tone or more entries in the input matrix are different from 1 or 0\n",
                "\tthe input matrix is of Guttman form; the sample space has only one element\n",
                "\tseed < 0 or seed > 2147483646\n"
                )

    if (err == 0) {
      cat("\nno error\n")
    } else {
      x <- bin2int(err)
      errstring <- paste("\n",paste(errortxt[(1:length(x))*x],sep="",collapse=""))
      stop(errstring, call.=FALSE)
    }
    invisible(err)
}

