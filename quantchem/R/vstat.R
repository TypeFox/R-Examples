"vstat" <-
function (x, ...) 
{
if (!is.data.frame(x)) { stop("x must be a data.frame") }

        avector = unlist(x);
        afactor = c()
        xx = list(x, ...);

if (length(xx) > 1) {
        avector = unlist(xx);
        bfactor = c(); nn = 0;
        for (i in 1:length(xx)) for (j in 1:length(xx[[i]])) {
        nn = nn + 1; afactor=c(afactor,rep(nn,nrow(xx[[i]])));
        bfactor=c(bfactor,rep(i,nrow(xx[[i]])));
        }       
        afactor = as.factor(afactor);
        bfactor = as.factor(bfactor);
        an1 = anova(aov(avector ~ afactor + bfactor));
        an2 = anova(aov(avector ~ bfactor + afactor));

        res = rbind(an2[1:2,],an1)
        res = rbind(res,c(res[3,1]+res[4,1],res[3,2]+res[4,2],0,NA,NA));
        res[5,3] = res[5,2]/res[5,1];
        res = cbind(res[,1:2],res[,2]/res[5,2]*100,res[,3],sqrt(res[,3]),sqrt(res[,3])/mean(avector)*100,res[,4:5]);
        rownames(res)=c("Groups","w/Groups","Series","w/Series","Total");
        
} else {
        for (i in 1:length(x)) { afactor=c(afactor,rep(i,nrow(x))); }
        afactor = as.factor(afactor);
        res = anova(aov(avector ~ afactor));
        res = rbind(res,c(res[1,1]+res[2,1],res[1,2]+res[2,2],0,NA,NA));
        res[3,3] = res[3,2]/res[3,1];
        res = cbind(res[,1:2],res[,2]/res[3,2]*100,res[,3],sqrt(res[,3]),sqrt(res[,3])/mean(avector)*100,res[,4:5]);
        rownames(res)=c("Series","w/Series","Total");

}
        colnames(res)=c("Df","Sum Sq","Percent","Mean Sq","SD","RSD","F","Pr(>F)");
        class(res)=c("anova","data.frame");
        attr(res,"heading")=c("\nVariability of results:\n");
        return(res);

}
