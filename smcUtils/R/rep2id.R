
rep2id = function(rep, engine="R") 
{
    engine=pmatch(engine, c("R","C"))
    sum = sum(rep)

    switch(engine,
    {
        # R implementation
        id = integer(sum)
        current.index = integer(1)
        for (i in 1:length(rep)) 
        {
            if (rep[i] != 0) 
            {
                id[current.index + 1:rep[i]] = i
                current.index = current.index + rep[i]
            }
        }
        return(id)
    },
    {
        # C implementation
        out = .C("rep2id_R", 
                 as.integer(rep),
                 as.integer(sum), 
                 id=integer(sum))
        return(out$id)
    })  
}

