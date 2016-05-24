### R code from vignette source 'gui_tutorial.rnw'

###################################################
### code chunk number 1: gui_tutorial.rnw:816-819
###################################################
require(sdcMicro)
sdc <- createSdcObj(testdata2, keyVars=c('urbrur','roof','walls','water','electcon','sex'), 
	numVars=c('expend','income','savings'), w='sampling_weight', hhId='ori_hid')


###################################################
### code chunk number 2: gui_tutorial.rnw:824-825
###################################################
slotNames(sdc) 


###################################################
### code chunk number 3: gui_tutorial.rnw:868-869
###################################################
print(sdc, "risk")	


