### chunk number 1: 
oopt <- options()
options("width"=75, "prompt"=" ","continue"="  ")


### chunk number 2: 
data(math)
tstart <- proc.time()


### chunk number 3: 
m0  <- rcox(~me:ve:al+al:an:st,data=math)


### chunk number 4: 
m0  <- rcox(vcc=list(~me, ~ve, ~al, ~an, ~st), 
            ecc=list(~me:ve, ~me:al, ~ve:al, ~al:an, ~al:st, ~an:st), 
            data=math)


### chunk number 5: 
m1  <- rcox(~al:an:st, vcc=list(~me+st, ~ve+an), 
                ecc=list(~me:ve+me:al, ~ve:al+al:st),
                data=math)
#m1


### chunk number 6: 
m1c  <- rcox(~al:an:st, vcc=list(~me+st, ~ve+an), 
                ecc=list(~me:ve+me:al, ~ve:al+al:st),
                data=math, type="rcor")


### chunk number 7: 
getecc(m1)


### chunk number 8: 
summary(m1)


### chunk number 9: 
summary(m1, type="KC")


### chunk number 10: 
summary(m1c, type="KC")


### chunk number 11: 
plot(m1)


### chunk number 12: 
update(m1, joinecc=list(~an:st, ~me:ve + me:al))


### chunk number 13: 
update(m1, joinecc=getecc(m1)[c("ecc2","ecc3")])


### chunk number 14: 
update(m1, splitvcc=~ve+an)


### chunk number 15: 
update(m1, splitvcc=getvcc(m1)["vcc3"])


### chunk number 16:  eval=FALSE
## update(m1, addecc=~me:an+ve:st)
## update(m1, dropecc=~me:ve+me:al)


### chunk number 17: 
ctab <- comparecc(m1, cc1=list(~me:ve+me:al, ~ve:al+al:st), 
          cc2=list(~an:st,~al:an), type="ecc", stat="dev")


### chunk number 18: 
comparecc(m1, cc1=~me:ve+me:al, cc2=NULL, type="ecc")
comparecc(m1,  type="ecc")
comparecc(m1,  type="vcc")


### chunk number 19: 
join1(m1, scope=list(~an:st, ~me:ve + me:al, ~ve:al + al:st),type="ecc")


### chunk number 20: 
drop1(m1, scope=list(~al:an, ~an:st, ~me:ve+me:al))


### chunk number 21: 
#split1(m1, type="ecc")
split1(m1, scope=list(~ve:al + al:st, ~me:ve+me:al), type="ecc")


### chunk number 22: 
add1(m1)


### chunk number 23: 
m01<-stepjoin1(m0,type="vcc")
m02<-stepjoin1(m01, type='ecc')


### chunk number 24: 
stepdrop1(m1, criterion="test", alpha=0.01)


### chunk number 25: 
m2  <- rcox(vcc=list(~al+me+st, ~ve+an), 
              ecc=list(~me:ve+me:al+ve:al, ~al:an+al:st+an:st), 
              data=math)


### chunk number 26: 
m3 <- stepsplit1(m2, type="vcc")
m4 <- stepsplit1(m3, type="ecc")


### chunk number 27:  eval=FALSE
## stepadd1(m1, criterion="test")#, alpha=0.8,details=2)


### chunk number 28:  eval=FALSE
## rcox(vcc=list(~me+ve+al,~st), data=math)
## rcox(vcc=list(list("me","ve","al"),list("st")), data=math)


### chunk number 29:  eval=FALSE
## rcox(ecc=list(~me:ve+me:al, ~ve:al), data=math)
## rcox(ecc=list(list(c("me","ve"),c("me","al")), list(c("ve","al"))), data=math)


### chunk number 30: 
add1(m1, scope=list(c("an", "me"),c("me", "st")))
add1(m1, scope=list(~an:me, ~me:st))


### chunk number 31: 
options(oopt)


