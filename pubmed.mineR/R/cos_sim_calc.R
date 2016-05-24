

cos_sim_calc = function(nummatrix) 
{cossimdata = NULL; for (i in 1:(nrow(nummatrix) - 1)) 
       {for (j in (i+1):nrow(nummatrix)) 
                     { tempa = sum(nummatrix[i,]*nummatrix[j,]) /(sqrt(sum(nummatrix[i,]^2))* sqrt(sum(nummatrix[j,]^2))) ;  tempb = names(nummatrix[i,1]); tempc = names(nummatrix[j,1]); tempd = cbind(tempb,tempc,tempa);  write.table(tempd, file = "cossimdata.txt", append = T, row.names = F, col.names = F, sep = "\t", quote = F  )}
        } 
        
}
