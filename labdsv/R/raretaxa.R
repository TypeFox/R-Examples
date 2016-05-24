raretaxa <- function (taxa,min=1,log=FALSE,type='b',panel='all') 
{
    rare <- apply(taxa>0,2,sum) <= min
    tmp <- apply(taxa[,rare]>0,1,sum)

    if (panel == 'all' || panel == 1) {
        if (log) {
            plot(rev(sort(tmp[tmp>0])),log='y',type=type,
                 xlab='Plot',ylab='Rare Species/Plot')
        } else {
            plot(rev(sort(tmp[tmp>0])),type=type,
                xlab='Plot',ylab='Rare Species/Plot') 
        }
        if (panel == 'all') readline('Hit return')
    }

    if (panel == 'all' || panel == 2) {
        tmp <- apply(taxa[,rare],2,sum)/apply(taxa[,rare]>0,2,sum)
        if (log) {
            plot(rev(sort(tmp)),type=type,log='y',xlab='Species',ylab='Mean Abundance')
        } else {
            plot(rev(sort(tmp)),type=type,xlab='Species',ylab='Mean Abundance')
        }
        if (panel == 'all')  readline('Hit return')
    }

    if (panel == 'all' || panel == 3) {
        tmp <- apply(taxa[,rare],1,sum)
        plot(rev(sort(tmp[tmp>0])),type=type,log='y',xlab='Plot',ylab='Total Abundance')
    }
}

