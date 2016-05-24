### This tests the result of the first example of the article
print("##############################")
print("This tests the result of the first example of the article \n")
dd2=matrix(c(4,0,0,0,30,0,0,0,23),ncol=3,byrow=TRUE)
callhaplotype(dd2)
callhaplotype(dd2)/(2*57)
### This tests the cubic routine
print("##############################")
print("This tests the cubic routine")
haplotypeit(4,0,0,0,30,0,0,0,23)
### Formated of 4 digits
print("Formated of 4 digits")
round(as.numeric(Re(haplotypeit(4,0,0,0,30,0,0,0,23)$A)),digit=4)
round(as.numeric(Re(haplotypeit(4,0,0,0,30,0,0,0,23)$B)),digit=4)
round(as.numeric(Re(haplotypeit(4,0,0,0,30,0,0,0,23)$C)),digit=4)
round(as.numeric(Re(haplotypeit(4,0,0,0,30,0,0,0,23)$D)),digit=4)
### The second example
print("##############################")
print("The second example: \n")
dd=matrix(c(1212, 2, 0, 679, 0,0,75,0,0), byrow=TRUE, nrow=3)
colnames(dd)=c("CC","CT","TT")
rownames(dd)=c("CC","CT","TT")
callhaplotype(dd)
### Check the result of the cubic equation of the second example
print("##############################")
print("Check the result of the cubic equation of the second example: \n")
temp2haplo =as.numeric(t(dd));
haplotypeit(temp2haplo[1],temp2haplo[2],temp2haplo[3],temp2haplo[4],temp2haplo[5],temp2haplo[6],temp2haplo[7],temp2haplo[8],temp2haplo[9]);
rm(temp2haplo)
### Third example
print("##############################")
print("Third example : \n")
dd3=matrix(c(1030,678,123,1,1,0,0,0,0),ncol=3,byrow=TRUE)
colnames(dd3)=c("AA","AG","GG")
rownames(dd3)=c("CC","CT","TT")
callhaplotype(dd3)
###  Check for alternative solutions
print("##############################")
print("Check for alternative solutions: \n")
temp2haplo =as.numeric(t(dd3));
haplotypeit(temp2haplo[1],temp2haplo[2],temp2haplo[3],temp2haplo[4],temp2haplo[5],temp2haplo[6],temp2haplo[7],temp2haplo[8],temp2haplo[9]);
rm(temp2haplo)
print("##############################")
