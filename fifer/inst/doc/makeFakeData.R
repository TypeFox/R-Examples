#### create fake names
pream = c("CD_", "TGF_", "IL_", "BCI_", "HemoLeptin_", "Glucose_", "TNF_", "B_regs_")
midam = 1:10
end_thing = LETTERS[1:5]

fake.names = intersperse(pream, midam, end_thing)
mns = runif(length(fake.names), 0, 100)
sd = mns*.15
fake.data = data.frame(matrix(rnorm(length(fake.names)*60), byrow=T, ncol=length(fake.names)))
names(fake.data) = fake.names
	
	#### induce a fake correlation
fake.cor = random.correlation(length(fake.names))	

	#### do choleski decomp
fake.data = data.matrix(fake.data) %*% chol(fake.cor)

	#### create a dependent variable (disease)
sign.columns = sample(1:ncol(fake.data), 5)	
coefs = runif(5, .2, 1)
disease = coefs[1]*fake.data[,sign.columns[1]] + coefs[2]*fake.data[,sign.columns[2]] + 
				coefs[3]*fake.data[,sign.columns[3]] + coefs[4]*fake.data[,sign.columns[4]] + coefs[5]*fake.data[,sign.columns[5]] + rnorm(nrow(fake.data),0, .25)


	#### convert means/variances
for (i in 1:length(fake.names)){
	fake.data[,i] = scaleIt(fake.data[,i], mean=mns[i], sd=sd[i])
}	

	#### put names back in 
fake.data = data.frame(fake.data)
names(fake.data) = sort(fake.names)
fake.data$disease = "case"
fake.data$disease[disease<median(disease)] = "control"

	#### put in fake demographics
gender = c("Male", "Female")
ethnicity = c("AA", "EA", "His", "NA")
fake.data$ID = lapply(1:nrow(fake.data),function(x){paste(sample(c(0:9, letters, LETTERS), 5, replace=TRUE),
                                 collapse="")})
fake.data$gender = sample(gender, size=nrow(fake.data), replace=T)	
fake.data$ethnicity = sample(ethnicity, size=nrow(fake.data), replace=T)	
fake.data$age = round(rnorm(nrow(fake.data), mean=40, sd=7))
names(fake.data)
fake.data = fake.data[,c("ID", "disease", "gender", "ethnicity", "age", r("B_regs_10A", "TNF_9E", names(fake.data), names=T))]
fakeMedicalData = fake.data
head(fake.data)


	#### randomly poke holes in the dataset
columnsMis = sample(1:ncol(fakeMedicalData), .05*dim(fakeMedicalData)[2])
fakeMedicalData[sample(1:nrow(fakeMedicalData), .1*nrow(fakeMedicalData)), columnsMis[1:5]] = NA
fakeMedicalData[sample(1:nrow(fakeMedicalData), .3*nrow(fakeMedicalData)), columnsMis[6:length(columnsMis)]] = NA
missing.vals(fakeMedicalData)

save(fakeMedicalData, file="research/RPackages/fifer/data/fakeMedicalData.rda")