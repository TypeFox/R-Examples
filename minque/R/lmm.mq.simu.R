lmm.mq.simu <-
function(formula, data = list(),v0=v,b0=b,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    result=genmod.simuold(gdata,v0,b0,SimuNum,JacNum,JacRep,ALPHA)
    return(result)
}
