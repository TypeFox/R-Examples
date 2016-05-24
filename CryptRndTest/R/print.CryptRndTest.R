print.CryptRndTest <- function(x,...){    
  
    if (x$name=="Achi"){
      cat("\n","  ","Test results for Adaptive Chi-Square test")
      cat("\n","  ","-----------------------------------------")
      cat("\n","   ","Value  :",paste(round(x$statistic,digits=3)))      
      if (round(x$p.value,digits=3)==1){
        cat("\n","    ","p-Value:",">0.999")
      }else if (round(x$p.value,digits=3)==0){
        cat("\n","    ","p-Value:","<0.001")
      }else{
        cat("\n","    ","p-Value:",round(x$p.value,digits=3))
      }      
      if (x$result.acsq==0){
          cat("\n","   ","Result : Reject H0")
      }else if (x$result.acsq==1){
          cat("\n","   ","Result : Not reject H0")
      }
      cat("\n","  ","-----------------------------------------")
    } else if (x$name=="BDS"){
      cat("\n","  ","Test results for Birthday Spacings test      ")
      cat("\n","  ","----------------------------------------------")
      cat("\n","   ","with Anderson-Darling goodness of fit test:")
      cat("\n","    ","Value  :",paste(round(x$AD.statistic,digits=3)))
      if (round(x$AD.pvalue,digits=3)==1){
        cat("\n","    ","p-Value:",">0.999")
      }else if (round(x$AD.pvalue,digits=3)==0){
        cat("\n","    ","p-Value:","<0.001")
      }else{
        cat("\n","    ","p-Value:",round(x$AD.pvalue,digits=3))
      } 
      if (x$AD.result==0){
        cat("\n","    ","Result : Reject H0")
      }else if (x$AD.result==1){
        cat("\n","    ","Result : Not reject H0")
      }
      cat("\n","  ","----------------------------------------------")
      cat("\n","   ","with Kolmogorov-Smirnov goodness of fit test:")
      cat("\n","    ","Value  :",paste(round(x$KS.statistic,digits=3)))
      if (round(x$KS.pvalue,digits=3)==1){
        cat("\n","    ","p-Value:",">0.999")
      }else if (round(x$KS.pvalue,digits=3)==0){
        cat("\n","    ","p-Value:","<0.001")
      }else{
        cat("\n","   ","p-Value:",round(x$KS.pvalue,digits=3))
      }      
      if (x$KS.result==0){
        cat("\n","    ","Result : Reject H0")
      }else if (x$KS.result==1){
        cat("\n","    ","Result : Not reject H0")
      }
      cat("\n","  ","----------------------------------------------")
      cat("\n","   ","with Chi-Square goodness of fit test:")
      cat("\n","    ","Value  :",paste(round(x$CS.statistic,digits=3)))
      if (round(x$CS.pvalue,digits=3)==1){
        cat("\n","    ","p-Value:",">0.999")
      }else if (round(x$CS.pvalue,digits=3)==0){
        cat("\n","    ","p-Value:","<0.001")
      }else{
        cat("\n","   ","p-Value:",round(x$CS.pvalue,digits=3))
      }  
      if (x$CS.result==0){
        cat("\n","    ","Result : Reject H0")
      }else if (x$CS.result==1){
        cat("\n","    ","Result : Not reject H0")
      }
      cat("\n","  ","----------------------------------------------")
    }else if (x$name=="BS"){
      cat("\n","  ","Test results for Book Stack test")
      cat("\n","  ","-----------------------------------------")
      cat("\n","   ","Value  :",paste(round(x$statistic,digits=3)))      
      if (round(x$p.value,digits=3)==1){
        cat("\n","   ","p-Value:",">0.999")
      }else if (round(x$p.value,digits=3)==0){
        cat("\n","   ","p-Value:","<0.001")
      }else{
        cat("\n","   ","p-Value:",round(x$p.value,digits=3))
      }      
      if (x$BS.result==0){
        cat("\n","   ","Result : Reject H0")
      }else if (x$BS.result==1){
        cat("\n","   ","Result : Not reject H0")
      }
      cat("\n","  ","-----------------------------------------")
    }else if (x$name=="GCD"){
      cat("\n","  ","Test results for Greatest Common Divisor test")
      cat("\n","  ","-----------------------------------------------")
      if (x$test.k==TRUE){
        cat("\n","  ","Results for the tests over k:")
        cat("\n","   ","with Anderson-Darling goodness of fit test:")
        if (round(x$sig.value.k[4],digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$sig.value.k[4],digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$sig.value.k[4],digits=3))
        } 
        if (x$AD.result.k==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$AD.result.k==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Kolmogorov-Smirnov goodness of fit test:")
        if (round(x$sig.value.k[1],digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$sig.value.k[1],digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$sig.value.k[1],digits=3))
        } 
        if (x$KS.result.k==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$KS.result.k==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Chi-Square goodness of fit test:")
        if (round(x$sig.value.k[2],digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$sig.value.k[2],digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$sig.value.k[2],digits=3))
        } 
        if (x$CSQ.result.k==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$CSQ.result.k==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Jarque-Bera goodness of fit test:")
        if (round(x$sig.value.k[3],digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$sig.value.k[3],digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$sig.value.k[3],digits=3))
        } 
        if (x$JB.result.k==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$JB.result.k==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","  ","-----------------------------------------------")
      }           
      if (x$test.g==TRUE){
        cat("\n","  ","Results for the tests over k:")
        cat("\n","   ","with Kolmogorov-Smirnov goodness of fit test:")
        if (round(x$sig.value.g[1],digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$sig.value.g[1],digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$sig.value.g[1],digits=3))
        } 
        if (x$KS.result.g==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$KS.result.g==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Chi-Square goodness of fit test:")
        if (round(x$sig.value.g[2],digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$sig.value.g[2],digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$sig.value.g[2],digits=3))
        } 
        if (x$CSQ.result.g==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$CSQ.result.g==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","  ","-----------------------------------------------")
      }      
    }else if (x$name=="RW"){
      cat("\n","  ","Test results for Random Walk tests")
      cat("\n","  ","-----------------------------------------------")
      if (x$Exc==TRUE){
        cat("\n","  ","Results for Excursion tests:")
        cat("\n","   ","with Anderson-Darling goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$AD.statistic.Excursion,digits=3)))
        if (round(x$AD.pvalue.Excursion,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$AD.pvalue.Excursion,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$AD.pvalue.Excursion,digits=3))
        } 
        if (x$AD.result.Excursion==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$AD.result.Excursion==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Kolmogorov-Smirnov goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$KS.statistic.Excursion,digits=3)))
        if (round(x$KS.pvalue.Excursion,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$KS.pvalue.Excursion,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$KS.pvalue.Excursion,digits=3))
        } 
        if (x$KS.result.Excursion==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$KS.result.Excursion==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Chi-Square goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$CS.statistic.Excursion,digits=3)))
        if (round(x$CS.pvalue.Excursion,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$CS.pvalue.Excursion,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$CS.pvalue.Excursion,digits=3))
        } 
        if (x$CS.result.Excursion==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$CS.result.Excursion==1){
          cat("\n","    ","Result : Not reject H0")
        }    
        cat("\n","  ","-----------------------------------------------")
      }
      if (x$Exp==TRUE){
        cat("\n","  ","Results for Expansion tests:")
        cat("\n","   ","with Anderson-Darling goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$AD.statistic.Expansion,digits=3)))
        if (round(x$AD.pvalue.Expansion,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$AD.pvalue.Expansion,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$AD.pvalue.Expansion,digits=3))
        } 
        if (x$AD.result.Expansion==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$AD.result.Expansion==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Kolmogorov-Smirnov goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$KS.statistic.Expansion,digits=3)))
        if (round(x$KS.pvalue.Expansion,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$KS.pvalue.Expansion,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$KS.pvalue.Expansion,digits=3))
        } 
        if (x$KS.result.Expansion==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$KS.result.Expansion==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Chi-Square goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$CS.statistic.Expansion,digits=3)))
        if (round(x$CS.pvalue.Expansion,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$CS.pvalue.Expansion,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$CS.pvalue.Expansion,digits=3))
        } 
        if (x$CS.result.Expansion==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$CS.result.Expansion==1){
          cat("\n","    ","Result : Not reject H0")
        }    
        cat("\n","  ","-----------------------------------------------")
      }
      if (x$Hei==TRUE){
        cat("\n","  ","Results for Height tests:")
        cat("\n","   ","with Anderson-Darling goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$AD.statistic.Height,digits=3)))
        if (round(x$AD.pvalue.Height,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$AD.pvalue.Height,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$AD.pvalue.Height,digits=3))
        } 
        if (x$AD.result.Height==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$AD.result.Height==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Kolmogorov-Smirnov goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$KS.statistic.Height,digits=3)))
        if (round(x$KS.pvalue.Height,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$KS.pvalue.Height,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$KS.pvalue.Height,digits=3))
        } 
        if (x$KS.result.Height==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$KS.result.Height==1){
          cat("\n","    ","Result : Not reject H0")
        }
        cat("\n","   ","with Chi-Square goodness of fit test:")
        cat("\n","    ","Value  :",paste(round(x$CS.statistic.Height,digits=3)))
        if (round(x$CS.pvalue.Height,digits=3)==1){
          cat("\n","    ","p-Value:",">0.999")
        }else if (round(x$CS.pvalue.Height,digits=3)==0){
          cat("\n","    ","p-Value:","<0.001")
        }else{
          cat("\n","    ","p-Value:",round(x$CS.pvalue.Height,digits=3))
        } 
        if (x$CS.result.Height==0){
          cat("\n","    ","Result : Reject H0")
        }else if (x$CS.result.Height==1){
          cat("\n","    ","Result : Not reject H0")
        }    
        cat("\n","  ","-----------------------------------------------")
      }
    }else if (x$name=="TBT"){
      cat("\n","  ","Test results for Topological Binary test")
      cat("\n","  ","-----------------------------------------")
      cat("\n","   ","Value  :",paste(round(x$statistic,digits=3)))                    
      if (x$result.TBT==0){
        cat("\n","   ","Result : Reject H0")
      }else if (x$result.TBT==1){
        cat("\n","   ","Result : Not reject H0")
      }
      cat("\n","  ","-----------------------------------------")
    }
    
  }