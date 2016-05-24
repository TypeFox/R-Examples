HW.test <-
function(count_11,count_12,count_22){
    allele1=count_11*2+count_12
	 allele2=count_22*2+count_12
	 count=count_11+count_12+count_22
	 freq1=allele1/(allele1+allele2)
	 freq2=allele2/(allele1+allele2)
	 num1=((count_11-freq1*freq1*count)*(count_11-freq1*freq1*count))/(freq1*freq1*count)
	 num2=((count_12-2*freq1*freq2*count)*(count_12-2*freq1*freq2*count))/(2*freq1*freq2*count)
	 num3=((count_22-freq2*freq2*count)*(count_22-freq2*freq2*count))/(freq2*freq2*count)
	 X2=num1+num2+num3
	 X2p=1-pchisq(X2,1)
	 return(list('X2'=X2,'p.value'=X2p))
}
