recoder <-
function(var, recode, other=NULL, na = NA, as.what = NULL)
{
    x.var = var
    z.var = var
    mentioned = rep(F, length(var))
    
    cases  = str_split(recode, ';')[[1]]
    for (onecase in  cases)
    {
        leftright = str_split(onecase, ':')[[1]] #->
        left = str_trim (leftright[1] )
        right = str_trim (leftright[2] )
        if (grepl('[=><i]' , left ,perl=TRUE)) {
            left1 = str_replace_all(left,"[&]",") & (var ")
            left2 = str_replace_all(left1,"[|]",") | (var ")
            left3 = paste('(var',left2, ')')
            left4 = str_replace_all(left3,"<-","< -")
            index  = which( eval(parse(text = left4) ) )
        } else {
            if (is.character(var) | is.factor(var) )
            index  = which( var %in% eval(parse(text = left) ) )
            else index  = which( var %in% left)
        }
        if (grepl('[$]', x = right)  ){
            rex = str_replace_all(right,"[$]","x.var[index]")
            right2 = eval(parse(text = rex ) )
        } else right2 = eval(parse(text = right ) )
        
        z.var[index] =  right2
        mentioned[index] = TRUE
    }
    z.var[is.na(var)]  = na
    
    if (!is.null(other) )
    z.var[!(mentioned | is.na(var) )] = other
    
    if (!is.null(as.what) )
    z.var = eval(parse(text = paste('as.',as.what,'(z.var)',sep="")))
    
    return(z.var)
}
