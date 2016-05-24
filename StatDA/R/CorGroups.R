"CorGroups" <-
function (dat, grouping, labels1, labels2, legend, ndigits = 4, method="pearson", ...) 
{
# Computes the correlation metrix for sub-groups of data
#
# dat ... data values (probably log10-transformed
# grouping ... factor with levels for different groups
# labels1, labels2 ... labels for groups
# legend ... plotting legend
# ndigits ... number of digits to be used for plotting the numbers
# method ... correlation method: "preason", "spearman", or "kendall"
#
    grplev=levels(grouping)
    ngroups=length(grplev)
    plot.new()
    p = ncol(dat)
    plot.window(xlim = c(-1, p-1), ylim = c(-1, p-1), xaxs = "i", yaxs = "i")
    for (i in 1:(p-1)) {
        text(i - 0.5, -0.5, labels1[i], srt = 0)
        segments(i,0,i,p-i)
    }
    for (i in 2:p) {
        text(-0.5, p - i + 0.5, labels2[i])
        segments(0,i-1,p+1-i,i-1)
    }
    segments(0,0,p-1,0)
    segments(0,0,0,p-1)
    for (i in 2:p) {
        for (j in 1:(i - 1)) {
            for (ii in 1:ngroups){
                 text(j - 0.5, p - i + (ii-1)/ngroups, 
                 round(cor(x=dat[grouping==grplev[ngroups+1-ii],c(i,j)],method=method,
                       y=NULL,use="all.obs")[2,1], ndigits), adj=c(0.5,-0.2))
            }
        }
    }
legend("topright",legend=legend)
}
