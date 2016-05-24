`marking` <-
function(dataset, 
                   frame, 
                   cts = "counts", 
                   streamFrame = NULL, 
                   allowanceFrame= 2, 
                   newcolname = "wearing")
{

    cat("frame is ", frame, "\n")
    cat("streamFrame is ", streamFrame, "\n")
    cat("allowanceFrame is ", allowanceFrame, "\n")

    ct = as.vector(dataset[,names(dataset) == cts])

    if(is.null(streamFrame)){
        streamFrame = round(0.5*frame)
    }

    #all the NA's in the original counts data will be treated as 0 counts
    ct1 = ct
    ct[is.na(ct)] = 0

    size = dim(dataset)[1]
    wearing = rep("nw", size)

    ct_bool = ct > 0
    rowPos = nthOccurance (dataVct = ct_bool, value= TRUE)

    #getting section start and end positions
    startpos = rowPos[1]
    endpos = c()
    prev = TRUE
    for(q  in 2: (length(rowPos)))
    {
        if(prev)
        {
            if( rowPos[q] - rowPos[q-1]>1 )
            {
                endpos = c(endpos, rowPos[q-1])
                startpos = c(startpos, rowPos[q])
            }
        }
        else
        {
            startpos = c(startpos, rowPos[q])
            prev = TRUE
        }
        if(q == length(rowPos))
        {endpos = c(endpos, rowPos[q])}
    }#end of q

    #ele3 should be handled here on startpos/endpos level
    allowancewin = endpos-startpos
    for(r in 1:length(allowancewin))
    {
        if(allowancewin[r] < allowanceFrame)
        {
            #upstream
            usStart = startpos[r] - streamFrame
            usEnd = startpos[r] - 1
            if(usStart <=0)
            {usStart = 1}
            if(usEnd <= 0)
            {usStart = 1}
            if(usEnd-usStart == 0){
                usSignal = "nowearing"
            }else {
                if(sum(ct_bool[usStart:usEnd]) >0){
                    usSignal = "wearing"
                }else {
                    usSignal = "nowearing"    
                }
            }

            #downstream
            dsEnd = endpos[r] + streamFrame
            dsStart = endpos[r] + 1
            if(dsEnd >size)
            {dsEnd = size}
            if(dsStart > size)
            {dsStart = size}
            if(dsEnd-dsStart == 0){
                dsSignal = "nowearing"
            }else {
                if(sum(ct_bool[dsStart:dsEnd]) >0){
                    dsSignal = "wearing"
                }else {
                    dsSignal = "nowearing"    
                }
            }  

            if(usSignal == "nowearing" & dsSignal == "nowearing")
            {
                startpos[r] = -1
                endpos[r] = -1
            }      
        }#end of if/allowancewin
    }#end of for/r

    startpos = startpos[startpos != -1]
    endpos = endpos[endpos!=-1]
    #end of ele3

    #now get the non-wearing gap
    #frame is the gap allowed between time section.  ie if 90 minutes allowed
    #between two wearing sections, them frame = 90
    gap = startpos[-1] - endpos[1:length(endpos)-1]
    endgap = endpos[1:length(gap)]
    startgap = startpos[-1]
    endgap[gap<= frame] = NA
    startgap [gap <= frame] = NA
    startgap = c(startpos[1], startgap)
    endgap = c(endgap, endpos[length(gap)+1])

    newstartpos = startgap[!is.na(startgap)]
    newendpos = endgap[!is.na(endgap)]

    for(w in 1: length(newendpos)){
        wearing[newstartpos[w]:newendpos[w]] = "w"
    }

    tlen= length(wearing)
    wearing[tlen] = wearing[tlen-1]

    wearing[is.na(ct1)] = NA

    oldnames = names(dataset)
    rst = cbind(dataset, wearing = wearing)
    names(rst) = c(oldnames, newcolname)
    return(rst)
}

