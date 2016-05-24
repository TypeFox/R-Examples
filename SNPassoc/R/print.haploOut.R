`print.haploOut` <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   cat("\n")
   cat("     Haplotype using SNPs:", attr(x, "label.snp"), " adjusted by:", attr(x, 
            "varAdj"), "\n")
   cat(" Interaction \n")
   cat("-------------------------\n")
   etiq<-dimnames(x[[1]])[[2]]

   if(attr(x,"quantitative"))
    {  
      etiq[2]<-paste(etiq[2],"(dif)")
      etiq[5]<-paste(etiq[5],"(dif)")
    }
   else
    {  
      etiq[2]<-paste(etiq[2],"(OR)")
      etiq[5]<-paste(etiq[5],"(OR)")
    }

   dimnames(x[[1]])[[2]]<-etiq
   print(x[[1]])
   cat("\n")
   cat("p interaction:",x[[4]],"\n")
   cat("\n",paste(attr(x,"varInt"),"within haplotype"), "\n")
   cat("-------------------------\n")
   print(x[[2]])
   cat(paste("haplotype within",attr(x,"varInt")), "\n")
   cat("-------------------------\n")
   print(x[[3]])

}

