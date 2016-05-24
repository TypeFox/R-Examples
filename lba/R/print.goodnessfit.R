print.goodnessfit <- function(x,
                              digits = 2L,
                              ...){

 object <- x 

 if(class(object)[1]=='goodnessfit.lba.ls'){

 col12 <- formatC(round(object[[2]],digits),format='f',digits=digits)
 ncol12 <- nchar(col12)

 col22 <- formatC(round(object[[1]],digits),format='f',digits=digits)
 ncol22 <- nchar(col22)

 par2 <- c(ncol12,ncol22)

 col01 <- formatC(round(object[[4]],digits),format='f',digits=digits)
 ncol01 <- nchar(col01)

 col02 <- formatC(round(object[[3]],digits),format='f',digits=digits)
 ncol02 <- nchar(col02)

 col14 <- formatC(round(object[[6]],digits),format='f',digits=digits)

 col24 <- formatC(round(object[[5]],digits),format='f',digits=digits)

 col23 <- formatC(round(object[[7]],digits),format='f',digits=digits)

 col13 <- formatC(round(object[[8]],digits),format='f',digits=digits) 

 col26 <- formatC(round(object[[9]],digits),format='f',digits=digits) 

 col16 <- formatC(round(object[[10]],digits),format='f',digits=digits)

 col18 <- formatC(round(object[[12]],digits),format='f',digits=digits)
 ncol18 <- nchar(col18)
 # 
 col28 <- formatC(round(object[[11]],digits),format='f',digits=digits)
 ncol28 <- nchar(col28)

 par3 <- c(ncol18,ncol28)

 col19 <- formatC(round(object[[14]],digits),format='f',digits=digits)
 ncol19 <- nchar(col19)

 col29 <- formatC(round(object[[13]],digits),format='f',digits=digits)
 ncol29 <- nchar(col19) 

 par4 <- c(ncol19,ncol29)

 col110 <- formatC(round(object[[16]],digits),format='f',digits=digits)
 ncol110 <- nchar(col110)

 col210 <- formatC(round(object[[15]],digits),format='f',digits=digits)
 ncol210 <- nchar(col210)

 par5 <- c(ncol110,ncol210)

 aux <- rbind(par2,par3,par4,par5)
 sums <- rowSums(aux)
 aux1 <- which(sums == max(sums))

 space <- aux[aux1[1],]

 col112 <- formatC(round(object[[17]],digits),format='f',digits=digits)

 col113 <- formatC(round(object[[18]],digits),format='f',digits=digits)

 col114 <- formatC(round(object[[19]],digits),format='f',digits=digits)

 col115 <- formatC(round(object[[20]],digits),format='f',digits=digits)  

 col116 <- formatC(round(object[[21]],digits),format='f',digits=digits)   

 col117 <- formatC(round(object[[22]],digits),format='f',digits=digits)

 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')   
 cat(paste('|                        ',
           paste('STATISTICS',paste(rep(' ',space[1]+space[2]+sum(space)+6),collapse=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')    
 cat(paste('|                                 ',
           paste(paste('LBM(K)',paste(rep(' ',2*space[1]-2),collapse=''),sep=''),paste('LBM(1)',paste(rep(' ',2*space[2]-3),collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|             degree of freedom:  ',
           paste(paste(col12 ,    paste(rep(' ',2*space[1]+4-nchar(col12)),collapse=''),sep=''),paste(col22,paste(rep(' ',2*space[2]+3-nchar(col22)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')     

 cat(paste('|                      ',paste('OTHER MEASURES',paste(rep(' ',4+space[1]+space[2]+sum(space)),collapse=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')          

 cat(paste('|  Residual sum of square (RSS):  ',
           paste(paste(col01 ,    paste(rep(' ',2*space[1]+4-nchar(col01))    ,collapse=''),sep=''),paste(col02,   paste(rep(' ',2*space[2]+3-nchar(col02))    ,collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col24,paste(rep(' ',2*space[2]+3-nchar(col24)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|           required per budget:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col14,paste(rep(' ',2*space[2]+3-nchar(col14)),collapse=''),sep=''),sep=''),sep=''),'|\n')  
 cat(paste('|required per defree of freedom:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col23,paste(rep(' ',2*space[2]+3-nchar(col23)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')           

 cat(paste('|                  Weighted RSS:  ',
           paste(paste(col26 ,    paste(rep(' ',2*space[1]+4-nchar(col26))    ,collapse=''),sep=''),paste(col13,   paste(rep(' ',2*space[2]+3-nchar(col13))    ,collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col16,paste(rep(' ',2*space[2]+3-nchar(col16)),collapse=''),sep=''),sep=''),sep=''),'|\n')  
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')             
 cat(paste('|        Index of dissimilarity:  ',
           paste(paste(col18 ,    paste(rep(' ',2*space[1]+4-nchar(col18))    ,collapse=''),sep=''),paste(col28,   paste(rep(' ',2*space[2]+3-nchar(col28))    ,collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')            
 cat(paste('|  Prop. correctly classf. data:  ',
           paste(paste(col19 ,    paste(rep(' ',2*space[1]+4-nchar(col19))    ,collapse=''),sep=''),paste(col29,   paste(rep(' ',2*space[2]+3-nchar(col29))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col210,paste(rep(' ',2*space[2]+3-nchar(col210)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|           required per budget:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col110,paste(rep(' ',2*space[2]+3-nchar(col110)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|required per defree of freedom:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col112,paste(rep(' ',2*space[2]+3-nchar(col112)),collapse=''),sep=''),sep=''),sep=''),'|\n')    
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')             
 cat(paste('|        Mean angular deviation:  ',
           paste(paste(col114 ,    paste(rep(' ',2*space[1]+4-nchar(col114))    ,collapse=''),sep=''),paste(col113,   paste(rep(' ',2*space[2]+3-nchar(col113))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col115,paste(rep(' ',2*space[2]+3-nchar(col115)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|           required per budget:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col116,paste(rep(' ',2*space[2]+3-nchar(col116)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|required per defree of freedom:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col117,paste(rep(' ',2*space[2]+3-nchar(col117)),collapse=''),sep=''),sep=''),sep=''),'|\n')     
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')              
 cat('\n') 

 } else {

 col1 <- formatC(round(object[[4]],digits),format='f',digits=digits)
 ncol1 <- nchar(col1)

 col2 <- formatC(round(object[[3]],digits),format='f',digits=digits)
 ncol2 <- nchar(col2)

 par1 <- c(ncol1,ncol2)

 col12 <- formatC(round(object[[2]],digits),format='f',digits=digits)
 ncol12 <- nchar(col12)

 col22 <- formatC(round(object[[1]],digits),format='f',digits=digits)
 ncol22 <- nchar(col22)

 par2 <- c(ncol12,ncol22)

 col13 <- formatC(round(object[[8]],digits),format='f',digits=digits)

 col23 <- formatC(round(object[[7]],digits),format='f',digits=digits)

 col14 <- formatC(round(object[[6]],digits),format='f',digits=digits)

 col24 <- formatC(round(object[[5]],digits),format='f',digits=digits)

 col16 <- formatC(round(object[[10]],digits),format='f',digits=digits)

 col26 <- formatC(round(object[[9]],digits),format='f',digits=digits)

 col18 <- formatC(round(object[[12]],digits),format='f',digits=digits)
 ncol18 <- nchar(col18)
# 
 col28 <- formatC(round(object[[11]],digits),format='f',digits=digits)
 ncol28 <- nchar(col28)

 par3 <- c(ncol18,ncol28)

 col19 <- formatC(round(object[[14]],digits),format='f',digits=digits)
 ncol19 <- nchar(col19)

 col29 <- formatC(round(object[[13]],digits),format='f',digits=digits)
 ncol29 <- nchar(col19) 

 par4 <- c(ncol19,ncol29)

 col110 <- formatC(round(object[[16]],digits),format='f',digits=digits)
 ncol110 <- nchar(col110)

 col210 <- formatC(round(object[[15]],digits),format='f',digits=digits)
 ncol210 <- nchar(col210)

 par5 <- c(ncol110,ncol210)

 aux <- rbind(par1,par2,par3,par4,par5)
 sums <- rowSums(aux)
 aux1 <- which(sums == max(sums))

 space <- aux[aux1[1],]

 col112 <- formatC(round(object[[17]],digits),format='f',digits=digits)

 col113 <- formatC(round(object[[18]],digits),format='f',digits=digits)

 col114 <- formatC(round(object[[19]],digits),format='f',digits=digits)

 col115 <- formatC(round(object[[20]],digits),format='f',digits=digits)  

 col116 <- formatC(round(object[[22]],digits),format='f',digits=digits)   

 col117 <- formatC(round(object[[21]],digits),format='f',digits=digits)

 col118 <- formatC(round(object[[23]],digits),format='f',digits=digits)

 col119 <- formatC(round(object[[24]],digits),format='f',digits=digits)

 col120 <- formatC(round(object[[25]],digits),format='f',digits=digits)      

 col121 <- formatC(round(object[[27]],digits),format='f',digits=digits)

 col122 <- formatC(round(object[[26]],digits),format='f',digits=digits)        

 col123 <- formatC(round(object[[29]],digits),format='f',digits=digits)         

 col124 <- formatC(round(object[[28]],digits),format='f',digits=digits)         

 col125 <- formatC(round(object[[30]],digits),format='f',digits=digits)          

 col126 <- formatC(round(object[[31]],digits),format='f',digits=digits)          

 col127 <- formatC(round(object[[32]],digits),format='f',digits=digits)           

 col128 <- formatC(round(object[[34]],digits),format='f',digits=digits)            

 col129 <- formatC(round(object[[33]],digits),format='f',digits=digits)

 col130 <- formatC(round(object[[35]],digits),format='f',digits=digits)

 col131 <- formatC(round(object[[36]],digits),format='f',digits=digits)              

 col132 <- formatC(round(object[[37]],digits),format='f',digits=digits)              

 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')   
 cat(paste('|                        ',
           paste('STATISTICS',paste(rep(' ',space[1]+space[2]+sum(space)+6),collapse=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')    
 cat(paste('|                                 ',
           paste(paste('LBM(K)',paste(rep(' ',2*space[1]-2),collapse=''),sep=''),paste('LBM(1)',paste(rep(' ',2*space[2]-3),collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|         Likelihood ratio test:  ',
           paste(paste(col1  ,    paste(rep(' ',2*space[1]+4-nchar(col1)),collapse=''),sep=''),paste(col2, paste(rep(' ',2*space[2]+3-nchar(col2))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|             degree of freedom:  ',
           paste(paste(col12 ,    paste(rep(' ',2*space[1]+4-nchar(col12)),collapse=''),sep=''),paste(col22,paste(rep(' ',2*space[2]+3-nchar(col22)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|                   probability:  ',
           paste(paste(col13 ,    paste(rep(' ',2*space[1]+4-nchar(col13)),collapse=''),sep=''),paste(col23,paste(rep(' ',2*space[2]+3-nchar(col23)),collapse=''),sep=''),sep=''),sep=''),'|\n')  
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')     

 cat(paste('|           Pearson  chi-square:  ',
           paste(paste(col14 ,    paste(rep(' ',2*space[1]+4-nchar(col14)),collapse=''),sep=''),paste(col24,paste(rep(' ',2*space[2]+3-nchar(col24)),collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|             degree of freedom:  ',
           paste(paste(col12 ,    paste(rep(' ',2*space[1]+4-nchar(col12)),collapse=''),sep=''),paste(col22,paste(rep(' ',2*space[2]+3-nchar(col22)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|                   probability:  ',
           paste(paste(col16 ,    paste(rep(' ',2*space[1]+4-nchar(col16)),collapse=''),sep=''),paste(col26,paste(rep(' ',2*space[2]+3-nchar(col26)),collapse=''),sep=''),sep=''),sep=''),'|\n')  
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')      

 cat(paste('|           ',paste('CRITERION FOR MODEL SELECTION',paste(rep(' ',space[1]+space[2]+sum(space)),collapse=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')       
 cat(paste('|                           AIC:  ',
           paste(paste(col18 ,    paste(rep(' ',2*space[1]+4-nchar(col18))    ,collapse=''),sep=''),paste(col28,   paste(rep(' ',2*space[2]+3-nchar(col28))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|                           BIC:  ',
           paste(paste(col19 ,    paste(rep(' ',2*space[1]+4-nchar(col19)),collapse=''),sep=''),paste(col29,paste(rep(' ',2*space[2]+3-nchar(col29)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|                          CAIC:  ',
           paste(paste(col110 ,    paste(rep(' ',2*space[1]+4-nchar(col110)),collapse=''),sep=''),paste(col210,paste(rep(' ',3+space[2]),collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')       

 cat(paste('|             ',paste('INCREMENTAL FIT INDICES',paste(rep(' ',4+space[1]+space[2]+sum(space)),collapse=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')        
 cat(paste('|              Normed fit index:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col112,paste(rep(' ',2*space[2]+3-nchar(col112)),collapse=''),sep=''),sep=''),sep=''),'|\n')  
 cat(paste('|     Normed fit index modified:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col113,paste(rep(' ',2*space[2]+3-nchar(col113)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|                  Bollen index:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col114,paste(rep(' ',2*space[2]+3-nchar(col114)),collapse=''),sep=''),sep=''),sep=''),'|\n')    
 cat(paste('|            Tucker-Lewis index:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col115,paste(rep(' ',2*space[2]+3-nchar(col115)),collapse=''),sep=''),sep=''),sep=''),'|\n')     
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')         

 cat(paste('|                   ',paste('OTHER MEASURES',paste(rep(' ',7+space[1]+space[2]+sum(space)),collapse=''),sep=''),sep=''),'|\n')
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')          
 cat(paste('|                Unweighted RSS:  ',
           paste(paste(col116 ,    paste(rep(' ',2*space[1]+4-nchar(col116))    ,collapse=''),sep=''),paste(col117,   paste(rep(' ',2*space[2]+3-nchar(col117))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col118,paste(rep(' ',2*space[2]+3-nchar(col118)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|           required per budget:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col119,paste(rep(' ',2*space[2]+3-nchar(col119)),collapse=''),sep=''),sep=''),sep=''),'|\n')  
 cat(paste('|required per defree of freedom:  ',
           paste(paste(' ',    paste(rep(' ',2*space[1]+3),collapse=''),sep=''),paste(col120,paste(rep(' ',2*space[2]+3-nchar(col120)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')           

 cat(paste('|        Index of dissimilarity:  ',
           paste(paste(col121 ,    paste(rep(' ',2*space[1]+4-nchar(col121))    ,collapse=''),sep=''),paste(col122,   paste(rep(' ',2*space[2]+3-nchar(col122))    ,collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')            
 cat(paste('|  Prop. correctly classf. data:  ',
           paste(paste(col123 ,    paste(rep(' ',2*space[1]+4-nchar(col123))    ,collapse=''),sep=''),paste(col124,   paste(rep(' ',2*space[2]+3-nchar(col124))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col125,paste(rep(' ',2*space[2]+3-nchar(col125)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|           required per budget:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col126,paste(rep(' ',2*space[2]+3-nchar(col126)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|required per defree of freedom:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col127,paste(rep(' ',2*space[2]+3-nchar(col127)),collapse=''),sep=''),sep=''),sep=''),'|\n')    
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')             
 cat(paste('|        Mean angular deviation:  ',
           paste(paste(col128 ,    paste(rep(' ',2*space[1]+4-nchar(col128))    ,collapse=''),sep=''),paste(col129,   paste(rep(' ',2*space[2]+3-nchar(col129))    ,collapse=''),sep=''),sep=''),sep=''),'|\n')
 cat(paste('|                   improvement:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col130,paste(rep(' ',2*space[2]+3-nchar(col130)),collapse=''),sep=''),sep=''),sep=''),'|\n') 
 cat(paste('|           required per budget:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col131,paste(rep(' ',2*space[2]+3-nchar(col131)),collapse=''),sep=''),sep=''),sep=''),'|\n')   
 cat(paste('|required per defree of freedom:  ',
           paste(paste(' ',    paste(rep(' ', 2*space[1]+3),collapse=''),sep=''),paste(col132,paste(rep(' ',2*space[2]+3-nchar(col132)),collapse=''),sep=''),sep=''),sep=''),'|\n')     
 cat(paste('|---------------------------------',
           paste(rep('-',space[1]+space[2]+sum(space)+7),collapse=''),sep=''),'|\n')              
 cat('\n')     

 }
 }
