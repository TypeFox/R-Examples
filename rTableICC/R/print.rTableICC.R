print.rTableICC <-
function(x,...){
  cat("Call:\n")
  print(x$call)
  if (x$ICC==TRUE){
    if (x$structure=="2x2xK"){
      cat("\n","  ","Process summary:")
      cat("\n","  ","----------------")  
      cat("\n","  ", paste(sum(x$N)),"observations in",paste(x$K),"centers were successfully generated under", paste(x$sampling),"sampling plan!")
      cat("\n","  ","Number of clusters for each center is as the following:")
      for (i in 1:x$K){
        cat("\n","    ",paste(x$M[i]),"for Center", paste(i))
      }
      if (sum(x$cluster.size==0)==0){
        cat("\n","  ", "Each cluster includes at least one individual." )
      }else{
        cat("\n","  ", paste(sum(x$cluster.size==0)),"clusters include no individual." )
      }
      cat("\n","  ", paste(sum(x$cluster.size==1)),"clusters include one individual.")  
      cat("\n","  ", paste(sum(x$cluster.size>1)),"clusters include more than one individual.")
      cat("\n","  ")
      cat("\n","The number of t sized clusters in the set of clusters in which all individuals fall in cell (j,k) for j,k=1,2:\n")
      cat("\n","  g.t = \n")
      names.1=array("",(x$K*2))
      names.2=array("",2)
      names.3=array("",(x$T-1))
      say=0
      for (j in 1:x$K){
        for (k in 1:2){
          say=say+1
          names.1[say]=paste("Center-",j,"R-",k)
        }
      }
      for (j in 1:2){
        names.2[j]=paste("C-",j)
      }
      for (i in 2:x$T){
        names.3[(i-1)]=paste("Cluster of size",i)  
      }
      names.all=list(names.1,names.2,names.3)
      dimnames(x$g.t) = names.all
      print(x$g.t)      
      
      cat("\n","The number of clusters of size t outside the set of clusters in which all individuals fall in a single cell: g.tilde = (", paste(x$g.tilde),")\n")

      cat("\n","Generated random table in two dimensions : \n")
      row.names=array("",x$K)
      col.names=c("R1C1","R1C2","R2C1","R2C2")
      for (i in 1:x$K){
        row.names[i]=paste("Center-",i)
      }
      rownames(x$rTable)=row.names
      colnames(x$rTable)=col.names
      print(x$rTable)
      
      if (x$print.regular==TRUE){
        cat("\n","  Generated random table in three dimensions : \n")
        names.1=array("",2)
        names.2=array("",2)
        names.3=array("",x$K)      
        for (j in 1:2){
          names.1[j]=paste("R-",j)
          names.2[j]=paste("C-",j)
        }
        for (i in 1:x$K){
          names.3[i]=paste("Center-",i)  
        }
        names.all=list(names.1,names.2,names.3)
        dimnames(x$rTable.regular) = names.all
        cat("    \n")
        print(x$rTable.regular)
      }
      if (x$print.raw==TRUE){
        cat("\n","  Generated random table in raw data format = \n")
        print(x$rTable.raw)
      }
    }else if (x$structure=="RxC"){
      cat("\n","  ","Process summary:")
      cat("\n","  ","----------------")
      cat("\n","  ", paste(x$N),"observations in",paste(x$M),"clusters were successfully generated under", paste(x$sampling),"sampling plan!")
      if (sum(x$cluster.size==0)==0){
        cat("\n","  ", "Each cluster includes at least one individual." )
      }else{
        cat("\n","  ", paste(sum(x$cluster.size==0)),"clusters include no individual." )
      }
      cat("\n","  ", paste(sum(x$cluster.size==1)),"clusters include one individual.")  
      cat("\n","  ", paste(sum(x$cluster.size>1)),"clusters include more than one individual.")
      cat("\n","  ")
      cat("\n","The number of t sized clusters in the set of clusters in which all individuals fall in cell (j,k) for j=1,...,R and k=1,...,C: g.t =\n")
      
      names.1=array("",x$R)
      names.2=array("",x$C)
      names.3=array("",x$T-1)
      for (j in 1:x$R){
        names.1[j]=paste("R-",j)
      }
      for (j in 1:x$C){
        names.2[j]=paste("C-",j)
      }
      for (i in 2:x$T){
        names.3[(i-1)]=paste("Cluster of size",i)  
      }
      names.all=list(names.1,names.2,names.3)
      dimnames(x$g.t) = names.all
      print(x$g.t)
      
      cat("\n","The number of clusters of size t outside the set of clusters in which all individuals fall in a single cell: g.tilde = (", paste(x$g.tilde),")\n")

      cat("\n","Generated random table in row format = (",paste(x$rTable),")\n")
      
      if (x$print.regular==TRUE){
        cat("\n","  Generated random table in RxC format = \n")
        row.namess=array("",x$R)
        col.namess=array("",x$C)
        for (i in 1:x$R){
          row.namess[i]=paste("R-",i)  
        }
        for (j in 1:x$C){
          col.namess[j]=paste("C-",j)
        }
        print(x$rTable.regular)
        rownames(x$rTable.regular) = row.namess
        colnames(x$rTable.regular) = col.namess
        print(x$rTable.regular)
      }
      if (x$print.raw==TRUE){
        cat("\n","Generated random table in raw data format = ")
        print(x$rTable.raw)
      }
    }    
  }else if (x$ICC==FALSE){
    if (x$structure=="2x2xK"){
      cat("\n","  ","Process summary:")
      cat("\n","  ","----------------")
      cat("\n","  ", paste(x$N),"observations in",paste(x$K),"centers were successfully generated under", paste(x$sampling),"sampling plan! \n")

      cat("\n","Generated random table in 2x2xK format = \n")
      
      names.1=array("",2)
      names.2=array("",2)
      names.3=array("",x$K)      
      for (j in 1:2){
        names.1[j]=paste("R-",j)
        names.2[j]=paste("C-",j)
      }
      for (i in 1:x$K){
        names.3[i]=paste("Center-",i)  
      }
      names.all=list(names.1,names.2,names.3)
      dimnames(x$rTable) = names.all
      cat("    \n")
      print(x$rTable)
      if (x$print.raw==TRUE){
        cat("\n","Generated random table in raw data format = \n")
        print(as.matrix(x$rTable.raw))
      }
    }else if (x$structure=="RxC"){
      cat("\n","  ","Process summary:")
      cat("\n","  ","----------------")
      cat("\n","  ", paste(x$N),"observations across",paste(x$R),"rows and",paste(x$C),"columns were successfully generated under", paste(x$sampling),"sampling plan! \n")

      cat("\n","Generated random table in RxC format = \n")
      row.names=array("",x$R)
      col.names=array("",x$C)
      for (i in 1:x$R){
        row.names[i]=paste("R-",i)  
      }
      for (j in 1:x$C){
        col.names[j]=paste("C-",j)
      }
      rownames(x$rTable) = row.names
      colnames(x$rTable) = col.names
      print(x$rTable)
      if (x$print.raw==TRUE){
        cat("\n","Generated random table in raw data format = \n")
        print(as.matrix(x$rTable.raw))
      }
    }
  }
  
}
