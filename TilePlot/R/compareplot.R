`compareplot` <-
function(arrayset1, arrayset2, array1label, array2label, title)
{

data_matrix = matrix(nrow=length(arrayset1),ncol=2)

data_matrix[,1] = arrayset1
data_matrix[,2] = arrayset2

j=1

new_data_matrix1=vector()
new_data_matrix2=vector()


for(i in 1:dim(data_matrix)[1])
{
if(data_matrix[i,1]!=0)
	{
	if(data_matrix[i,2]!=0)
		{
			new_data_matrix1[j] <- data_matrix[i,1]
			new_data_matrix2[j] <- data_matrix[i,2]
			j=j+1
		}
	}
}

new_data_frame <-data.frame(cbind(x=log(new_data_matrix1),y=log(new_data_matrix2)))
rsquared = summary(lm("y ~ x",new_data_frame))$r.squared

plot(log(new_data_matrix1),log(new_data_matrix2),xlab=array1label,ylab=array2label,main=title, sub=as.character(round(rsquared, digits=3)))
abline(summary(lm("y ~ x",new_data_frame))$coefficients[1,1],summary(lm("y ~ x",new_data_frame))$coefficients[2,1])
}

