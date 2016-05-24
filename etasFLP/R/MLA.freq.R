MLA.freq <-
function(x){

                f   =table(x)
                k   =length(f)

                m   =matrix(0,k,7)
                m   =as.data.frame(m)
            
                if(is.factor(x))    m[,1]   =(rownames(f)) else     m[,1]   =as.numeric(rownames(f))
                f   =as.numeric(f)
                m[,2]   =f
                n   =sum(f)
                m[,3]   =f/n
                m[,4]   =cumsum(f)
                m[,5]   =m[,4]/n
                m[,6]   =n-m[,4]+m[,2]
                m[,7]   =m[,6]/n
            colnames(m) =c("x","freq", "rel.freq.","cum.freq","rel.cum.","back.cum","back.cum.rel." )
            
return(m)   
}
