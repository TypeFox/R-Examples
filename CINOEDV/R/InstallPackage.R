InstallPackage <-
function(){
  # Install Packages that sometimes CINOEDV needs.
  # In current version, it needs "R.matlab" and "igraph".
  #
  # Junliang Shang
  # 3.11/2014
  
  # R.matlab
  if (length(grep("R.matlab",.packages(all.available=TRUE)))!=0){
    cat("    package R.matlab successfully installed.\n")
  }else
  {
    install.packages("R.matlab") 
  }
  
  # igraph
  if (length(grep("igraph",.packages(all.available=TRUE)))!=0){
    cat("    package igraph successfully installed.\n")
  }else
  {
    install.packages("igraph")
  }
  
  # ggplot2
  if (length(grep("ggplot2",.packages(all.available=TRUE)))!=0){
    cat("    package ggplot2 successfully installed.\n")
  }else
  {
    install.packages("ggplot2")
  }
  
  # reshape2
  if (length(grep("reshape2",.packages(all.available=TRUE)))!=0){
    cat("    package reshape2 successfully installed.\n")
  }else
  {
    install.packages("reshape2")
  }
}
