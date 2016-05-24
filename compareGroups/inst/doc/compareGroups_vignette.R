### R code from vignette source 'compareGroups_vignette.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: compareGroups_vignette.rnw:106-107
###################################################
library(compareGroups)


###################################################
### code chunk number 2: compareGroups_vignette.rnw:152-153
###################################################
data(predimed)


###################################################
### code chunk number 3: compareGroups_vignette.rnw:160-168
###################################################
dicc<-data.frame(
"Name"=I(names(predimed)),
"Label"=I(unlist(lapply(predimed,label))),
"Codes"=I(unlist(lapply(predimed,function(x) paste(levels(x),collapse="; "))))
)
dicc$Codes<-sub(">=","$\\\\geq$",dicc$Codes)  
#print(xtable(dicc,align=rep("l",4)),include.rownames=FALSE,size="small",tabular.environment="longtable", sanitize.text.function=function(x) x)
print(xtable(dicc,align=rep("l",4)),include.rownames=FALSE,size="small", sanitize.text.function=function(x) x)


###################################################
### code chunk number 4: compareGroups_vignette.rnw:194-196
###################################################
predimed$tmain <- with(predimed, Surv(toevent, event == 'Yes'))
label(predimed$tmain) <- "AMI, stroke, or CV Death"


###################################################
### code chunk number 5: compareGroups_vignette.rnw:216-217
###################################################
compareGroups(group ~ . , data=predimed)


###################################################
### code chunk number 6: compareGroups_vignette.rnw:228-229
###################################################
compareGroups(group ~ . -toevent - event, data=predimed)


###################################################
### code chunk number 7: compareGroups_vignette.rnw:236-238
###################################################
res<-compareGroups(group ~ age + sex + smoke + waist + hormo, data=predimed)
res


###################################################
### code chunk number 8: compareGroups_vignette.rnw:257-259
###################################################
compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
              subset = sex=='Female')


###################################################
### code chunk number 9: compareGroups_vignette.rnw:266-267
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 10: compareGroups_vignette.rnw:270-272
###################################################
compareGroups(group ~ age + sex + smoke + waist + hormo, data=predimed, 
              selec = list(hormo= sex=="Female", waist = waist>20 ))


###################################################
### code chunk number 11: compareGroups_vignette.rnw:277-278
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 12: compareGroups_vignette.rnw:281-283
###################################################
compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
              selec = list(waist= !is.na(hormo)), subset = sex=="Female")


###################################################
### code chunk number 13: compareGroups_vignette.rnw:289-290
###################################################
options(width=80,keep.source=TRUE)


###################################################
### code chunk number 14: compareGroups_vignette.rnw:292-294
###################################################
compareGroups(group ~ age + sex + bmi + bmi + waist + hormo, data=predimed, 
              selec = list(bmi.1=!is.na(hormo)))


###################################################
### code chunk number 15: compareGroups_vignette.rnw:307-308
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 16: compareGroups_vignette.rnw:310-312
###################################################
compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
              method = c(waist=2))


###################################################
### code chunk number 17: compareGroups_vignette.rnw:314-315
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 18: compareGroups_vignette.rnw:330-331
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 19: compareGroups_vignette.rnw:333-335
###################################################
compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
              method = c(waist=NA), alpha= 0.01)


###################################################
### code chunk number 20: compareGroups_vignette.rnw:337-338
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 21: compareGroups_vignette.rnw:347-351
###################################################
cuts<-"lo:55=1; 56:60=2; 61:65=3; 66:70=4; 71:75=5; 76:80=6; 81:hi=7"
predimed$age7gr<-car::recode(predimed$age, cuts)
compareGroups(group ~ age7gr, data=predimed, method = c(age7gr=NA))
compareGroups(group ~ age7gr, data=predimed, method = c(age7gr=NA), min.dis=8)


###################################################
### code chunk number 22: compareGroups_vignette.rnw:367-368
###################################################
compareGroups(age7gr ~ sex + bmi + waist, data=predimed, max.ylev=7)


###################################################
### code chunk number 23: compareGroups_vignette.rnw:373-374
###################################################
compareGroups(group ~ sex + age7gr, method= (age7gr=3), data=predimed, max.xlev=5)


###################################################
### code chunk number 24: compareGroups_vignette.rnw:392-394
###################################################
compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
              include.label= FALSE)


###################################################
### code chunk number 25: compareGroups_vignette.rnw:401-402
###################################################
options(width=80, keep.source=FALSE)


###################################################
### code chunk number 26: compareGroups_vignette.rnw:405-408
###################################################
resu1<-compareGroups(group ~ age + waist, data=predimed, 
                       method = c(waist=2))
createTable(resu1)


###################################################
### code chunk number 27: compareGroups_vignette.rnw:415-418
###################################################
resu2<-compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
                       method = c(waist=2), Q1=0.025, Q3=0.975)
createTable(resu2)


###################################################
### code chunk number 28: compareGroups_vignette.rnw:425-426
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 29: compareGroups_vignette.rnw:430-431
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 30: compareGroups_vignette.rnw:433-435
###################################################
compareGroups(group ~ age + smoke + waist + hormo, data=predimed, 
                method = c(waist=2), Q1=0, Q3=1)


###################################################
### code chunk number 31: compareGroups_vignette.rnw:437-438
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 32: compareGroups_vignette.rnw:445-449
###################################################
predimed$smk<-predimed$smoke
levels(predimed$smk)<- c("Never smoker", "Current or former < 1y", "Never or former >= 1y", "Unknown")
label(predimed$smk)<-"Smoking 4 cat."
cbind(table(predimed$smk))


###################################################
### code chunk number 33: compareGroups_vignette.rnw:482-483
###################################################
compareGroups(group ~ age + smk, data=predimed, simplify=FALSE)


###################################################
### code chunk number 34: compareGroups_vignette.rnw:514-517
###################################################
res<-compareGroups(group ~ age + sex + smoke + waist + hormo, method = c(waist=2),
                   data=predimed)
summary(res[c(1, 2, 4)])


###################################################
### code chunk number 35: compareGroups_vignette.rnw:530-531
###################################################
plot(res[c(1,2)], file="./figures/univar/", type="pdf")


###################################################
### code chunk number 36: compareGroups_vignette.rnw:549-550
###################################################
plot(res[c(1,2)], bivar=TRUE, file="./figures/bivar/")


###################################################
### code chunk number 37: compareGroups_vignette.rnw:574-576
###################################################
res<-compareGroups(group ~ age + sex + smoke + waist + hormo, data=predimed)
res


###################################################
### code chunk number 38: compareGroups_vignette.rnw:581-582
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 39: compareGroups_vignette.rnw:584-586
###################################################
res<-update(res, . ~. - sex +  bmi + toevent, subset = sex=='Female', 
              method = c(waist=2, tovent=2), selec = list(bmi=!is.na(hormo)))


###################################################
### code chunk number 40: compareGroups_vignette.rnw:588-589
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 41: compareGroups_vignette.rnw:591-592
###################################################
res


###################################################
### code chunk number 42: compareGroups_vignette.rnw:606-612
###################################################
library(SNPassoc)
data(SNPs)
tab <- createTable(compareGroups(casco ~ snp10001 + snp10002 + snp10005 
                                 + snp10008 + snp10009, SNPs))
pvals <- getResults(tab, "p.overall")
p.adjust(pvals, method = "BH")


###################################################
### code chunk number 43: compareGroups_vignette.rnw:629-631
###################################################
res1<-compareGroups(htn ~ age + sex + bmi + smoke, data=predimed, ref=1)
createTable(res1, show.ratio=TRUE)


###################################################
### code chunk number 44: compareGroups_vignette.rnw:637-638
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 45: compareGroups_vignette.rnw:640-643
###################################################
res2<-compareGroups(htn ~ age + sex + bmi + smoke, data=predimed, 
                    ref=c(smoke=1, sex=2))
createTable(res2, show.ratio=TRUE)


###################################################
### code chunk number 46: compareGroups_vignette.rnw:648-649
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 47: compareGroups_vignette.rnw:656-659
###################################################
res<-compareGroups(htn ~ age + sex + bmi + hormo + hyperchol, data=predimed, 
                   ref.no='NO')
createTable(res, show.ratio=TRUE)


###################################################
### code chunk number 48: compareGroups_vignette.rnw:668-670
###################################################
res<-compareGroups(htn ~ age + bmi, data=predimed)
createTable(res, show.ratio=TRUE)


###################################################
### code chunk number 49: compareGroups_vignette.rnw:675-676
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 50: compareGroups_vignette.rnw:678-681
###################################################
res<-compareGroups(htn ~ age + bmi, data=predimed, 
                   fact.ratio= c(age=10, bmi=2))
createTable(res, show.ratio=TRUE)


###################################################
### code chunk number 51: compareGroups_vignette.rnw:683-684
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 52: compareGroups_vignette.rnw:693-695
###################################################
res<-compareGroups(htn ~ age + sex + bmi + hyperchol, data=predimed)
createTable(res, show.ratio=TRUE)


###################################################
### code chunk number 53: compareGroups_vignette.rnw:701-703
###################################################
res<-compareGroups(htn ~ age + sex + bmi + hyperchol, data=predimed, ref.y=2)
createTable(res, show.ratio=TRUE)


###################################################
### code chunk number 54: compareGroups_vignette.rnw:717-719
###################################################
plot(compareGroups(tmain ~ sex, data=predimed), bivar=TRUE, file="./figures/bivar/")
plot(compareGroups(tmain ~ age, data=predimed), bivar=TRUE, file="./figures/bivar/")


###################################################
### code chunk number 55: compareGroups_vignette.rnw:754-755
###################################################
options(width=80,keep.source=FALSE)


###################################################
### code chunk number 56: compareGroups_vignette.rnw:757-760
###################################################
res<-compareGroups(sex ~  age + tmain, timemax=c(tmain=3),
                   data=predimed)
res


###################################################
### code chunk number 57: compareGroups_vignette.rnw:762-763
###################################################
options(width=120,keep.source=FALSE)


###################################################
### code chunk number 58: compareGroups_vignette.rnw:772-774
###################################################
plot(res[2], file="./figures/univar/")
plot(res[2], bivar=TRUE, file="./figures/bivar/")


###################################################
### code chunk number 59: compareGroups_vignette.rnw:799-800
###################################################
options(width=100,keep.source=FALSE)


###################################################
### code chunk number 60: compareGroups_vignette.rnw:803-806
###################################################
res<-compareGroups(group ~ age + sex + smoke + waist + hormo, data=predimed, 
              selec = list(hormo=sex=="Female"))
restab<-createTable(res)


###################################################
### code chunk number 61: compareGroups_vignette.rnw:811-812
###################################################
print(restab,which.table='descr')


###################################################
### code chunk number 62: compareGroups_vignette.rnw:817-818
###################################################
print(restab,which.table='avail')


###################################################
### code chunk number 63: compareGroups_vignette.rnw:834-835
###################################################
update(restab, hide = c(sex="Male"))


###################################################
### code chunk number 64: compareGroups_vignette.rnw:844-846
###################################################
res<-compareGroups(group ~ age + sex + htn + diab, data=predimed)
createTable(res, hide.no='no', hide = c(sex="Male"))


###################################################
### code chunk number 65: compareGroups_vignette.rnw:855-856
###################################################
createTable(res, digits= c(age=2, sex = 3))


###################################################
### code chunk number 66: compareGroups_vignette.rnw:865-866
###################################################
createTable(res, type=1)


###################################################
### code chunk number 67: compareGroups_vignette.rnw:871-872
###################################################
createTable(res, type=3)


###################################################
### code chunk number 68: compareGroups_vignette.rnw:884-885
###################################################
createTable(res, show.n=TRUE)


###################################################
### code chunk number 69: compareGroups_vignette.rnw:892-893
###################################################
createTable(res, show.descr=FALSE)


###################################################
### code chunk number 70: compareGroups_vignette.rnw:900-901
###################################################
createTable(res, show.all=TRUE)


###################################################
### code chunk number 71: compareGroups_vignette.rnw:907-908
###################################################
createTable(res, show.p.overall=FALSE)


###################################################
### code chunk number 72: compareGroups_vignette.rnw:915-916
###################################################
createTable(res, show.p.trend=TRUE)


###################################################
### code chunk number 73: compareGroups_vignette.rnw:926-927
###################################################
createTable(res, show.p.mul=TRUE)


###################################################
### code chunk number 74: compareGroups_vignette.rnw:933-934
###################################################
createTable(update(res, subset= group!="Control diet"), show.ratio=TRUE)


###################################################
### code chunk number 75: compareGroups_vignette.rnw:939-941
###################################################
createTable(compareGroups(tmain ~  group + age + sex, data=predimed), 
            show.ratio=TRUE)


###################################################
### code chunk number 76: compareGroups_vignette.rnw:950-952
###################################################
createTable(compareGroups(tmain ~  group + age + sex, data=predimed), 
            show.ratio=TRUE, digits.ratio= 3)


###################################################
### code chunk number 77: compareGroups_vignette.rnw:959-962
###################################################
tab<-createTable(compareGroups(tmain ~  group + age + sex, data=predimed), 
                 show.all = TRUE)
print(tab, header.labels = c("p.overall" = "p-value", "all" = "All"))


###################################################
### code chunk number 78: compareGroups_vignette.rnw:974-977
###################################################
restab1 <- createTable(compareGroups(group ~ age + sex, data=predimed))
restab2 <- createTable(compareGroups(group ~ bmi + smoke, data=predimed))
rbind("Non-modifiable risk factors"=restab1, "Modifiable risk factors"=restab2)


###################################################
### code chunk number 79: compareGroups_vignette.rnw:990-991
###################################################
rbind("Non-modifiable"=restab1,"Modifiable"=restab2)[c(1,4)]


###################################################
### code chunk number 80: compareGroups_vignette.rnw:996-997
###################################################
rbind("Modifiable"=restab1,"Non-modifiable"=restab2)[c(4,3,2,1)]


###################################################
### code chunk number 81: compareGroups_vignette.rnw:1008-1013
###################################################
res<-compareGroups(group ~ age +  smoke + bmi + htn , data=predimed)
alltab <- createTable(res,  show.p.overall = FALSE)
femaletab <- createTable(update(res,subset=sex=='Female'), show.p.overall = FALSE)
maletab <- createTable(update(res,subset=sex=='Male'), show.p.overall = FALSE)
cbind("ALL"=alltab,"FEMALE"=femaletab,"MALE"=maletab)


###################################################
### code chunk number 82: compareGroups_vignette.rnw:1018-1019
###################################################
cbind(alltab,femaletab,maletab,caption=NULL)


###################################################
### code chunk number 83: compareGroups_vignette.rnw:1024-1025
###################################################
cbind(alltab,femaletab,maletab)


###################################################
### code chunk number 84: compareGroups_vignette.rnw:1038-1040
###################################################
print(createTable(compareGroups(group ~ age + sex + smoke + waist + hormo,
                                data=predimed)), which.table='both')


###################################################
### code chunk number 85: compareGroups_vignette.rnw:1043-1045
###################################################
print(createTable(compareGroups(group ~ age + sex + smoke + waist + hormo,
                                data=predimed)),  nmax=FALSE)


###################################################
### code chunk number 86: compareGroups_vignette.rnw:1051-1053
###################################################
summary(createTable(compareGroups(group ~ age + sex + smoke + waist + hormo, 
                                  data=predimed)))


###################################################
### code chunk number 87: compareGroups_vignette.rnw:1058-1062
###################################################
res<-compareGroups(group ~ age + sex + smoke + waist + hormo, data=predimed)
restab<-createTable(res, type=1, show.ratio=TRUE )
restab
update(restab, show.n=TRUE)


###################################################
### code chunk number 88: compareGroups_vignette.rnw:1068-1069
###################################################
update(restab, x = update(res, subset=c(sex=='Female')), show.n=TRUE)


###################################################
### code chunk number 89: compareGroups_vignette.rnw:1080-1081
###################################################
createTable(compareGroups(group ~ age + sex + smoke + waist + hormo, data=predimed))


###################################################
### code chunk number 90: compareGroups_vignette.rnw:1083-1084
###################################################
createTable(compareGroups(group ~ age + sex + bmi, data=predimed))[1:2, ]


###################################################
### code chunk number 91: compareGroups_vignette.rnw:1130-1133
###################################################
restab<-createTable(compareGroups(group ~ age + sex + smoke + waist + hormo, 
                                  data=predimed))
export2latex(restab)


###################################################
### code chunk number 92: compareGroups_vignette.rnw:1170-1172 (eval = FALSE)
###################################################
## ?report   # to know more about report function
## ?regicor  # info about REGICOR data set


###################################################
### code chunk number 93: compareGroups_vignette.rnw:1187-1191
###################################################
# from a compareGroups object
data(regicor)
res <- compareGroups(year ~ .-id, regicor)
missingTable(res)


###################################################
### code chunk number 94: compareGroups_vignette.rnw:1194-1197 (eval = FALSE)
###################################################
## # or from createTable objects
## restab <- createTable(res, hide.no = 'no')
## missingTable(restab)


###################################################
### code chunk number 95: compareGroups_vignette.rnw:1203-1209
###################################################
# first create time-to-cardiovascular event
regicor$tcv<-with(regicor,Surv(tocv,cv=='Yes'))
# create the table
res <- compareGroups(tcv ~ . -id-tocv-cv-todeath-death, regicor, include.miss = TRUE)
restab <- createTable(res, hide.no = 'no')
restab


###################################################
### code chunk number 96: compareGroups_vignette.rnw:1226-1228
###################################################
data(SNPs)
head(SNPs)


###################################################
### code chunk number 97: compareGroups_vignette.rnw:1233-1235
###################################################
res<-compareSNPs(casco ~ snp10001 + snp10002 + snp10003, data=SNPs)
res


###################################################
### code chunk number 98: compareGroups_vignette.rnw:1243-1245
###################################################
res<-compareSNPs(~ snp10001 + snp10002 + snp10003, data=SNPs)
res


###################################################
### code chunk number 99: compareGroups_vignette.rnw:1267-1270
###################################################
export2latex(createTable(compareGroups(group ~  age + sex + smoke + bmi + waist + 
wth + htn + diab + hyperchol + famhist + hormo + p14 + toevent + event,
                          data=predimed), hide.no="No",hide = c(sex="Male")))


###################################################
### code chunk number 100: compareGroups_vignette.rnw:1332-1337
###################################################
export2latex(createTable(compareGroups(htn ~  age + sex + smoke + bmi + waist + 
                              wth + diab + hyperchol + famhist + hormo + 
                            p14 + toevent + event,
                            data=predimed), hide.no="No",hide = c(sex="Male"),
            show.ratio=TRUE, show.descr=FALSE))


###################################################
### code chunk number 101: compareGroups_vignette.rnw:1350-1352
###################################################
export2latex(createTable(compareGroups(tmain ~  group + age + sex, data=predimed), 
            show.ratio=TRUE))


