`partition` <-
function(ms){ 

             n<-length(ms)

             temp<-list(ms[[1]])

             m<-1

             for (i in 1:(n-1)){

                           temp1<-part(temp[[m]],ms[[i+1]])

                           temp<-c(temp[-m],temp1)

                           m<-length(temp)

             }

             temp1<-list()

             for (i in 1:(2*n-1)){

                           temp1<-c(temp1,temp[i][length(temp[[i]])>0])

             }

             part.list<-list(union(temp[1],temp[2]))

             for (i in 1:(n-2)){

                           temp.list1<-union(temp[2*i], temp[(2*i+1)])

                           temp.list2<-list(union(temp.list1,temp[(2*i+2)]))

                           part.list<-c(part.list,temp.list2)

             }

             part.list<-c(part.list,list(union(temp[2*(n-1)], temp[(2*(n-1)+1)])))

             list(part.list,temp,temp1)

}

