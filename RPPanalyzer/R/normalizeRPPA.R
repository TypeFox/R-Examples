normalizeRPPA <-

function(x,method="row",normalizer="housekeeping",useCol="BCA",writetable=F,vals="logged"){

        # method=c("row","extValue","proteinDye","housekeeping")
        
         if (method=="row"){
        
            xn <- norm.texas(x,writetable=writetable,vals=vals)
            }

         if (method=="extValue"){
         
            xn <- norm.BCA(x,proteinc=useCol,writetable=writetable,vals=vals)
            }
            
         if(method=="proteinDye"){

            xn <- norm.static.I(x,writetable=writetable,vals=vals)
            }
            
         if(method=="housekeeping"){

            xn <- norm.static.II(x,normalizer=normalizer,writetable=writetable,vals=vals)
            }
            
         if(method=="spotbyspot"){

            xn <- spotbyspot(x,normalizer=normalizer,vals=vals)
            }
            
            return(xn)
            
}
