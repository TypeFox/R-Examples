sim.index <- function(list1, list2, FUN, index=c("sd","ji","all")){
  job1 <- FUN(list1)
  job2 <- FUN(list2)
  if(".attrs.id" %in% names(job1) & ".attrs.id" %in% names(job2)){
      if(index == "sd"){
        sd <- 2*length(intersect(job1$.attrs.id,job2$.attrs.id))/(length(unique(job1$.attrs.id))+length(unique(job2$.attrs.id))) # Sorensen-Dice index
        list("Sorensen-Dice index"=round(sd,2))
      }
      
      else if(index == "ji"){
        ji <- length(intersect(job1$.attrs.id,job2$.attrs.id))/length(union(job1$.attrs.id,job2$.attrs.id)) # Jaccard index
        list("Jaccard index"=round(ji,2))
      }
      
      else{
        sd <- 2*length(intersect(job1$.attrs.id,job2$.attrs.id))/(length(unique(job1$.attrs.id))+length(unique(job2$.attrs.id))) # Sorensen-Dice index
        ji <- length(intersect(job1$.attrs.id,job2$.attrs.id))/length(union(job1$.attrs.id,job2$.attrs.id)) # Jaccard index
        print(list("Sorensen-Dice index"=round(sd,2),"Jaccard index"=round(ji,2))) # print both in a list
      }
  }
  else{
   message("This function cannot yet handle this job data type. Please try another.")
  }
}