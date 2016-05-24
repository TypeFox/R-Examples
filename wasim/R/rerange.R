`rerange` <-
function(data,min.goal=0, max.goal=1,min.data=min(data),
max.data=max(data),center=NA){
   if(!is.na(center)){
     #split linear transformation
     #the data is split into two segments
     #one below center, one above
     #each segment undergoes its own linear transf.
     below <- data <=center
     above <- data > center
     newData <- data
     newData[below] <- rerange(c(center,data[below]), min.goal=min.goal, max.goal=max.goal/2)[-1]
     newData[above] <- rerange(c(center,data[above]), min.goal=max.goal/2, max.goal=max.goal)[-1]
     return(newData)
   } else {
     #normal linear transformation
     max.data<-max.data # Watch out! This is required because of lazzy evaluation!
     data<-data-min.data
     data<-data/(max.data-min.data)*(max.goal-min.goal)
     data<-data+min.goal
     return(data)
   }
}

