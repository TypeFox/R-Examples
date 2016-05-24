grep_split <-
function(keyword,x){
   # keyword = "th | os"
   #x =  c("The","licenses", "for","most","software")
  if(is.null(keyword)){
       extract.nums = nums = 1:length(x)
       }else{
         if(length(grep("\\|",keyword)) >0)
           keyword = str_trim(unlist(strsplit(keyword,"\\|")))
         extract.nums = c()
         for(i in keyword){
            keywords = unlist(strsplit(i," "))
            nums = grep(keywords[1],x,ignore.case = TRUE)
            if(length(keywords)==1){
                nums = nums
                }else{
                   for(j in 2:length(keywords)){
                       nums.2 = grep(keywords[j],x,ignore.case = TRUE)
                       nums = intersect(nums,nums.2)
                   }
             }
           extract.nums = unique(c(extract.nums,nums))   
         }    
   }   
     return(extract.nums)
}
