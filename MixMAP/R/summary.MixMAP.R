#summary Method
setMethod("summary",
    signature(object = "MixMAP"),
    function (object, ...) 
    {
		   	       
  			#Print the number of genes detected and the total number of genes considered
           cat(paste0("\nNumber of Genes Detected: ",object@num.genes.detected[1],"\n"))
           cat(paste0("Total Number of Genes: ",object@num.genes.detected[2],"\n\n"))
           cat("Top Genes:\n")
           
           #If any genes were detected, print up to ten of them
           if (object@num.genes.detected[1]>0){
           tmp<-object@detected.genes[order(object@detected.genes$predUpper),c(1,2,5,8,9,10)]
           tmp<-tmp[order(tmp$postEst),]
           print(as.data.frame(tmp[1:min(10,dim(tmp)[1]),]))
           }
           
           #If no genes are detected, display "No Genes Detected"
           if (object@num.genes.detected[1]==0){
          print("No Genes Detected")
           }
    
 
    }
)

