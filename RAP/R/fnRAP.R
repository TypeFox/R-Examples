fnRAP <-
function(a,nr,nc)
{

### nr = number of rows in original matrix
### nc = number of cols in original matrix
tot_tab=0

if ( (nr == 2 && nc > 2) || (nr > 2 && nc == 2) || (nr > 2 && nc > 2)) {

mat=matrix(a,nr,nc, byrow=TRUE)
outd = as.data.frame(t(c(dim(mat), toString(c(1, 2, 3, 4, 5)), toString(c(1, 2, 3, 4, 5)), round(chisq.test(mat)$p.value, 4), chisq.test(mat)$p.value < 0.05)))
ori_con = round(chisq.test(mat)$p.value, 4) # ms_code-for original conclusion: test result of given matrix
outlist=list()

for (i in c(2:nr))
{
outlist[[i]]=list()

for (j in c(2:nc))
{
outlist[[i]][[j]]=list()

rpermmat = t(combn(c(1:nr),i))
cpermmat = t(combn(c(1:nc),j))

if(i==nr) rfact = 1 else rfact = factorial(nr)/(factorial(i) * factorial(nr-i))

if(j==nc) cfact = 1 else cfact = factorial(nc)/(factorial(j) * factorial(nc-j))


if(i==nr && j==nc)  {
names(outd) = c("No. of rows", "No. of cols", "Selected rows", "Selected cols", "Pvalue", "Pvalue significant at 5%?")
outlist[[i]][[j]] = list()
outlist[[i]][[j]][[1]] = mat
outlist[[i]][[j]][[2]] = outd 
}

else {
for (ii in c(1:rfact))
for (jj in c(1:cfact))
{

tempvar = round(chisq.test(mat [ c(rpermmat[ii,]), c(cpermmat[jj,]) ])$p.value, 4)
                              if ((ori_con < 0.05 && tempvar > 0.05) || (ori_con > 0.05 && tempvar < 0.05)) #ms_code-for checking the reversal in conclusions
{
outlist[[i]][[j]][(ii-1)*cfact + jj] = tempvar
outd = rbind(outd, t(c(dim(mat [ c(rpermmat[ii,]), c(cpermmat[jj,]) ]), toString(c(rpermmat[ii,])), toString(c(cpermmat[jj,])), tempvar, tempvar < 0.05)))
                                        tot_tab=tot_tab+1
}
}

}

}


}


### outlist[[nr]][[nc]]

print(outlist[[nr]][[nc]])



} ### closing the first if condition

else { print ("Can only generate matrices of order greater than 2x2!")}

}

