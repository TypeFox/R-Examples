classname <-'version'
funcname <-'.version.version'
standGen <- function(object) standardGeneric('.version.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.version<-'
standGen <- function(x, value) standardGeneric('.version.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'version'
funcname <-'.version.build'
standGen <- function(object) standardGeneric('.version.build')
standMethod <- function(object) object@build
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.build<-'
standGen <- function(x, value) standardGeneric('.version.build<-')
standMethod <- function(x, value) {x@build<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'version'
funcname <-'.version.update'
standGen <- function(object) standardGeneric('.version.update')
standMethod <- function(object) object@update
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.update<-'
standGen <- function(x, value) standardGeneric('.version.update<-')
standMethod <- function(x, value) {x@update<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'version'
funcname <-'.version.svnrev'
standGen <- function(object) standardGeneric('.version.svnrev')
standMethod <- function(object) object@svnrev
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.svnrev<-'
standGen <- function(x, value) standardGeneric('.version.svnrev<-')
standMethod <- function(x, value) {x@svnrev<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.path'
standGen <- function(object) standardGeneric('.experiment.path')
standMethod <- function(object) object@path
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.path<-'
standGen <- function(x, value) standardGeneric('.experiment.path<-')
standMethod <- function(x, value) {x@path<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.name'
standGen <- function(object) standardGeneric('.experiment.name')
standMethod <- function(object) object@name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.name<-'
standGen <- function(x, value) standardGeneric('.experiment.name<-')
standMethod <- function(x, value) {x@name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.subject.num'
standGen <- function(object) standardGeneric('.experiment.subject.num')
standMethod <- function(object) object@subject.num
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.subject.num<-'
standGen <- function(x, value) standardGeneric('.experiment.subject.num<-')
standMethod <- function(x, value) {x@subject.num<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.subject.names'
standGen <- function(object) standardGeneric('.experiment.subject.names')
standMethod <- function(object) object@subject.names
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.subject.names<-'
standGen <- function(x, value) standardGeneric('.experiment.subject.names<-')
standMethod <- function(x, value) {x@subject.names<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.condition.num'
standGen <- function(object) standardGeneric('.experiment.condition.num')
standMethod <- function(object) object@condition.num
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.condition.num<-'
standGen <- function(x, value) standardGeneric('.experiment.condition.num<-')
standMethod <- function(x, value) {x@condition.num<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.condition.names'
standGen <- function(object) standardGeneric('.experiment.condition.names')
standMethod <- function(object) object@condition.names
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.condition.names<-'
standGen <- function(x, value) standardGeneric('.experiment.condition.names<-')
standMethod <- function(x, value) {x@condition.names<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.expRda'
standGen <- function(object) standardGeneric('.experiment.expRda')
standMethod <- function(object) object@expRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.expRda<-'
standGen <- function(x, value) standardGeneric('.experiment.expRda<-')
standMethod <- function(x, value) {x@expRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.optionsRda'
standGen <- function(object) standardGeneric('.experiment.optionsRda')
standMethod <- function(object) object@optionsRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.optionsRda<-'
standGen <- function(x, value) standardGeneric('.experiment.optionsRda<-')
standMethod <- function(x, value) {x@optionsRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.startRda'
standGen <- function(object) standardGeneric('.experiment.startRda')
standMethod <- function(object) object@startRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.startRda<-'
standGen <- function(x, value) standardGeneric('.experiment.startRda<-')
standMethod <- function(x, value) {x@startRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.dataRda'
standGen <- function(object) standardGeneric('.experiment.dataRda')
standMethod <- function(object) object@dataRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.dataRda<-'
standGen <- function(x, value) standardGeneric('.experiment.dataRda<-')
standMethod <- function(x, value) {x@dataRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.modelRda'
standGen <- function(object) standardGeneric('.experiment.modelRda')
standMethod <- function(object) object@modelRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.modelRda<-'
standGen <- function(x, value) standardGeneric('.experiment.modelRda<-')
standMethod <- function(x, value) {x@modelRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.statsRda'
standGen <- function(object) standardGeneric('.experiment.statsRda')
standMethod <- function(object) object@statsRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.statsRda<-'
standGen <- function(x, value) standardGeneric('.experiment.statsRda<-')
standMethod <- function(x, value) {x@statsRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.regRda'
standGen <- function(object) standardGeneric('.experiment.regRda')
standMethod <- function(object) object@regRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.regRda<-'
standGen <- function(x, value) standardGeneric('.experiment.regRda<-')
standMethod <- function(x, value) {x@regRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.funcRda'
standGen <- function(object) standardGeneric('.experiment.funcRda')
standMethod <- function(object) object@funcRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.funcRda<-'
standGen <- function(x, value) standardGeneric('.experiment.funcRda<-')
standMethod <- function(x, value) {x@funcRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.subjectPrefix'
standGen <- function(object) standardGeneric('.experiment.subjectPrefix')
standMethod <- function(object) object@subjectPrefix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.subjectPrefix<-'
standGen <- function(x, value) standardGeneric('.experiment.subjectPrefix<-')
standMethod <- function(x, value) {x@subjectPrefix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.conditionPrefix'
standGen <- function(object) standardGeneric('.experiment.conditionPrefix')
standMethod <- function(object) object@conditionPrefix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.conditionPrefix<-'
standGen <- function(x, value) standardGeneric('.experiment.conditionPrefix<-')
standMethod <- function(x, value) {x@conditionPrefix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.modelPrefix'
standGen <- function(object) standardGeneric('.experiment.modelPrefix')
standMethod <- function(object) object@modelPrefix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.modelPrefix<-'
standGen <- function(x, value) standardGeneric('.experiment.modelPrefix<-')
standMethod <- function(x, value) {x@modelPrefix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.subjectDir'
standGen <- function(object) standardGeneric('.experiment.subjectDir')
standMethod <- function(object) object@subjectDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.subjectDir<-'
standGen <- function(x, value) standardGeneric('.experiment.subjectDir<-')
standMethod <- function(x, value) {x@subjectDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.conditionDir'
standGen <- function(object) standardGeneric('.experiment.conditionDir')
standMethod <- function(object) object@conditionDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.conditionDir<-'
standGen <- function(x, value) standardGeneric('.experiment.conditionDir<-')
standMethod <- function(x, value) {x@conditionDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.dataDir'
standGen <- function(object) standardGeneric('.experiment.dataDir')
standMethod <- function(object) object@dataDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.dataDir<-'
standGen <- function(x, value) standardGeneric('.experiment.dataDir<-')
standMethod <- function(x, value) {x@dataDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.weightsDir'
standGen <- function(object) standardGeneric('.experiment.weightsDir')
standMethod <- function(object) object@weightsDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.weightsDir<-'
standGen <- function(x, value) standardGeneric('.experiment.weightsDir<-')
standMethod <- function(x, value) {x@weightsDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.avgDir'
standGen <- function(object) standardGeneric('.experiment.avgDir')
standMethod <- function(object) object@avgDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.avgDir<-'
standGen <- function(x, value) standardGeneric('.experiment.avgDir<-')
standMethod <- function(x, value) {x@avgDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.regDir'
standGen <- function(object) standardGeneric('.experiment.regDir')
standMethod <- function(object) object@regDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.regDir<-'
standGen <- function(x, value) standardGeneric('.experiment.regDir<-')
standMethod <- function(x, value) {x@regDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.funcDir'
standGen <- function(object) standardGeneric('.experiment.funcDir')
standMethod <- function(object) object@funcDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.funcDir<-'
standGen <- function(x, value) standardGeneric('.experiment.funcDir<-')
standMethod <- function(x, value) {x@funcDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.betaDir'
standGen <- function(object) standardGeneric('.experiment.betaDir')
standMethod <- function(object) object@betaDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.betaDir<-'
standGen <- function(x, value) standardGeneric('.experiment.betaDir<-')
standMethod <- function(x, value) {x@betaDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.modelDir'
standGen <- function(object) standardGeneric('.experiment.modelDir')
standMethod <- function(object) object@modelDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.modelDir<-'
standGen <- function(x, value) standardGeneric('.experiment.modelDir<-')
standMethod <- function(x, value) {x@modelDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.statsDir'
standGen <- function(object) standardGeneric('.experiment.statsDir')
standMethod <- function(object) object@statsDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.statsDir<-'
standGen <- function(x, value) standardGeneric('.experiment.statsDir<-')
standMethod <- function(x, value) {x@statsDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.modeldatDir'
standGen <- function(object) standardGeneric('.experiment.modeldatDir')
standMethod <- function(object) object@modeldatDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.modeldatDir<-'
standGen <- function(x, value) standardGeneric('.experiment.modeldatDir<-')
standMethod <- function(x, value) {x@modeldatDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.avgdatFile'
standGen <- function(object) standardGeneric('.experiment.avgdatFile')
standMethod <- function(object) object@avgdatFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.avgdatFile<-'
standGen <- function(x, value) standardGeneric('.experiment.avgdatFile<-')
standMethod <- function(x, value) {x@avgdatFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.avgWFile'
standGen <- function(object) standardGeneric('.experiment.avgWFile')
standMethod <- function(object) object@avgWFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.avgWFile<-'
standGen <- function(x, value) standardGeneric('.experiment.avgWFile<-')
standMethod <- function(x, value) {x@avgWFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.avgtstatFile'
standGen <- function(object) standardGeneric('.experiment.avgtstatFile')
standMethod <- function(object) object@avgtstatFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.avgtstatFile<-'
standGen <- function(x, value) standardGeneric('.experiment.avgtstatFile<-')
standMethod <- function(x, value) {x@avgtstatFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.modelDataFile'
standGen <- function(object) standardGeneric('.experiment.modelDataFile')
standMethod <- function(object) object@modelDataFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.modelDataFile<-'
standGen <- function(x, value) standardGeneric('.experiment.modelDataFile<-')
standMethod <- function(x, value) {x@modelDataFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.modelnamesRda'
standGen <- function(object) standardGeneric('.experiment.modelnamesRda')
standMethod <- function(object) object@modelnamesRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.modelnamesRda<-'
standGen <- function(x, value) standardGeneric('.experiment.modelnamesRda<-')
standMethod <- function(x, value) {x@modelnamesRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.residualFile'
standGen <- function(object) standardGeneric('.experiment.residualFile')
standMethod <- function(object) object@residualFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.residualFile<-'
standGen <- function(x, value) standardGeneric('.experiment.residualFile<-')
standMethod <- function(x, value) {x@residualFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.derivativeFile'
standGen <- function(object) standardGeneric('.experiment.derivativeFile')
standMethod <- function(object) object@derivativeFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.derivativeFile<-'
standGen <- function(x, value) standardGeneric('.experiment.derivativeFile<-')
standMethod <- function(x, value) {x@derivativeFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.weightFile'
standGen <- function(object) standardGeneric('.experiment.weightFile')
standMethod <- function(object) object@weightFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.weightFile<-'
standGen <- function(x, value) standardGeneric('.experiment.weightFile<-')
standMethod <- function(x, value) {x@weightFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.lowresFile'
standGen <- function(object) standardGeneric('.experiment.lowresFile')
standMethod <- function(object) object@lowresFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.lowresFile<-'
standGen <- function(x, value) standardGeneric('.experiment.lowresFile<-')
standMethod <- function(x, value) {x@lowresFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.lowresAvg'
standGen <- function(object) standardGeneric('.experiment.lowresAvg')
standMethod <- function(object) object@lowresAvg
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.lowresAvg<-'
standGen <- function(x, value) standardGeneric('.experiment.lowresAvg<-')
standMethod <- function(x, value) {x@lowresAvg<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.logFile'
standGen <- function(object) standardGeneric('.experiment.logFile')
standMethod <- function(object) object@logFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.logFile<-'
standGen <- function(x, value) standardGeneric('.experiment.logFile<-')
standMethod <- function(x, value) {x@logFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'experiment'
funcname <-'.experiment.version'
standGen <- function(object) standardGeneric('.experiment.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.experiment.version<-'
standGen <- function(x, value) standardGeneric('.experiment.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.expRda'
standGen <- function(object) standardGeneric('.settings.expRda')
standMethod <- function(object) object@expRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.expRda<-'
standGen <- function(x, value) standardGeneric('.settings.expRda<-')
standMethod <- function(x, value) {x@expRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.optionsRda'
standGen <- function(object) standardGeneric('.settings.optionsRda')
standMethod <- function(object) object@optionsRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.optionsRda<-'
standGen <- function(x, value) standardGeneric('.settings.optionsRda<-')
standMethod <- function(x, value) {x@optionsRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.startRda'
standGen <- function(object) standardGeneric('.settings.startRda')
standMethod <- function(object) object@startRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.startRda<-'
standGen <- function(x, value) standardGeneric('.settings.startRda<-')
standMethod <- function(x, value) {x@startRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.dataRda'
standGen <- function(object) standardGeneric('.settings.dataRda')
standMethod <- function(object) object@dataRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.dataRda<-'
standGen <- function(x, value) standardGeneric('.settings.dataRda<-')
standMethod <- function(x, value) {x@dataRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.modelRda'
standGen <- function(object) standardGeneric('.settings.modelRda')
standMethod <- function(object) object@modelRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.modelRda<-'
standGen <- function(x, value) standardGeneric('.settings.modelRda<-')
standMethod <- function(x, value) {x@modelRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.statsRda'
standGen <- function(object) standardGeneric('.settings.statsRda')
standMethod <- function(object) object@statsRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.statsRda<-'
standGen <- function(x, value) standardGeneric('.settings.statsRda<-')
standMethod <- function(x, value) {x@statsRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.regRda'
standGen <- function(object) standardGeneric('.settings.regRda')
standMethod <- function(object) object@regRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.regRda<-'
standGen <- function(x, value) standardGeneric('.settings.regRda<-')
standMethod <- function(x, value) {x@regRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.funcRda'
standGen <- function(object) standardGeneric('.settings.funcRda')
standMethod <- function(object) object@funcRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.funcRda<-'
standGen <- function(x, value) standardGeneric('.settings.funcRda<-')
standMethod <- function(x, value) {x@funcRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.subjectPrefix'
standGen <- function(object) standardGeneric('.settings.subjectPrefix')
standMethod <- function(object) object@subjectPrefix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.subjectPrefix<-'
standGen <- function(x, value) standardGeneric('.settings.subjectPrefix<-')
standMethod <- function(x, value) {x@subjectPrefix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.conditionPrefix'
standGen <- function(object) standardGeneric('.settings.conditionPrefix')
standMethod <- function(object) object@conditionPrefix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.conditionPrefix<-'
standGen <- function(x, value) standardGeneric('.settings.conditionPrefix<-')
standMethod <- function(x, value) {x@conditionPrefix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.modelPrefix'
standGen <- function(object) standardGeneric('.settings.modelPrefix')
standMethod <- function(object) object@modelPrefix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.modelPrefix<-'
standGen <- function(x, value) standardGeneric('.settings.modelPrefix<-')
standMethod <- function(x, value) {x@modelPrefix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.subjectDir'
standGen <- function(object) standardGeneric('.settings.subjectDir')
standMethod <- function(object) object@subjectDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.subjectDir<-'
standGen <- function(x, value) standardGeneric('.settings.subjectDir<-')
standMethod <- function(x, value) {x@subjectDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.conditionDir'
standGen <- function(object) standardGeneric('.settings.conditionDir')
standMethod <- function(object) object@conditionDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.conditionDir<-'
standGen <- function(x, value) standardGeneric('.settings.conditionDir<-')
standMethod <- function(x, value) {x@conditionDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.dataDir'
standGen <- function(object) standardGeneric('.settings.dataDir')
standMethod <- function(object) object@dataDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.dataDir<-'
standGen <- function(x, value) standardGeneric('.settings.dataDir<-')
standMethod <- function(x, value) {x@dataDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.weightsDir'
standGen <- function(object) standardGeneric('.settings.weightsDir')
standMethod <- function(object) object@weightsDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.weightsDir<-'
standGen <- function(x, value) standardGeneric('.settings.weightsDir<-')
standMethod <- function(x, value) {x@weightsDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.avgDir'
standGen <- function(object) standardGeneric('.settings.avgDir')
standMethod <- function(object) object@avgDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.avgDir<-'
standGen <- function(x, value) standardGeneric('.settings.avgDir<-')
standMethod <- function(x, value) {x@avgDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.regDir'
standGen <- function(object) standardGeneric('.settings.regDir')
standMethod <- function(object) object@regDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.regDir<-'
standGen <- function(x, value) standardGeneric('.settings.regDir<-')
standMethod <- function(x, value) {x@regDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.funcDir'
standGen <- function(object) standardGeneric('.settings.funcDir')
standMethod <- function(object) object@funcDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.funcDir<-'
standGen <- function(x, value) standardGeneric('.settings.funcDir<-')
standMethod <- function(x, value) {x@funcDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.betaDir'
standGen <- function(object) standardGeneric('.settings.betaDir')
standMethod <- function(object) object@betaDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.betaDir<-'
standGen <- function(x, value) standardGeneric('.settings.betaDir<-')
standMethod <- function(x, value) {x@betaDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.modelDir'
standGen <- function(object) standardGeneric('.settings.modelDir')
standMethod <- function(object) object@modelDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.modelDir<-'
standGen <- function(x, value) standardGeneric('.settings.modelDir<-')
standMethod <- function(x, value) {x@modelDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.statsDir'
standGen <- function(object) standardGeneric('.settings.statsDir')
standMethod <- function(object) object@statsDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.statsDir<-'
standGen <- function(x, value) standardGeneric('.settings.statsDir<-')
standMethod <- function(x, value) {x@statsDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.modeldatDir'
standGen <- function(object) standardGeneric('.settings.modeldatDir')
standMethod <- function(object) object@modeldatDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.modeldatDir<-'
standGen <- function(x, value) standardGeneric('.settings.modeldatDir<-')
standMethod <- function(x, value) {x@modeldatDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.avgdatFile'
standGen <- function(object) standardGeneric('.settings.avgdatFile')
standMethod <- function(object) object@avgdatFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.avgdatFile<-'
standGen <- function(x, value) standardGeneric('.settings.avgdatFile<-')
standMethod <- function(x, value) {x@avgdatFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.avgWFile'
standGen <- function(object) standardGeneric('.settings.avgWFile')
standMethod <- function(object) object@avgWFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.avgWFile<-'
standGen <- function(x, value) standardGeneric('.settings.avgWFile<-')
standMethod <- function(x, value) {x@avgWFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.avgtstatFile'
standGen <- function(object) standardGeneric('.settings.avgtstatFile')
standMethod <- function(object) object@avgtstatFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.avgtstatFile<-'
standGen <- function(x, value) standardGeneric('.settings.avgtstatFile<-')
standMethod <- function(x, value) {x@avgtstatFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.modelDataFile'
standGen <- function(object) standardGeneric('.settings.modelDataFile')
standMethod <- function(object) object@modelDataFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.modelDataFile<-'
standGen <- function(x, value) standardGeneric('.settings.modelDataFile<-')
standMethod <- function(x, value) {x@modelDataFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.modelnamesRda'
standGen <- function(object) standardGeneric('.settings.modelnamesRda')
standMethod <- function(object) object@modelnamesRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.modelnamesRda<-'
standGen <- function(x, value) standardGeneric('.settings.modelnamesRda<-')
standMethod <- function(x, value) {x@modelnamesRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.residualFile'
standGen <- function(object) standardGeneric('.settings.residualFile')
standMethod <- function(object) object@residualFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.residualFile<-'
standGen <- function(x, value) standardGeneric('.settings.residualFile<-')
standMethod <- function(x, value) {x@residualFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.derivativeFile'
standGen <- function(object) standardGeneric('.settings.derivativeFile')
standMethod <- function(object) object@derivativeFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.derivativeFile<-'
standGen <- function(x, value) standardGeneric('.settings.derivativeFile<-')
standMethod <- function(x, value) {x@derivativeFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.weightFile'
standGen <- function(object) standardGeneric('.settings.weightFile')
standMethod <- function(object) object@weightFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.weightFile<-'
standGen <- function(x, value) standardGeneric('.settings.weightFile<-')
standMethod <- function(x, value) {x@weightFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.lowresFile'
standGen <- function(object) standardGeneric('.settings.lowresFile')
standMethod <- function(object) object@lowresFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.lowresFile<-'
standGen <- function(x, value) standardGeneric('.settings.lowresFile<-')
standMethod <- function(x, value) {x@lowresFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.lowresAvg'
standGen <- function(object) standardGeneric('.settings.lowresAvg')
standMethod <- function(object) object@lowresAvg
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.lowresAvg<-'
standGen <- function(x, value) standardGeneric('.settings.lowresAvg<-')
standMethod <- function(x, value) {x@lowresAvg<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.logFile'
standGen <- function(object) standardGeneric('.settings.logFile')
standMethod <- function(object) object@logFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.logFile<-'
standGen <- function(x, value) standardGeneric('.settings.logFile<-')
standMethod <- function(x, value) {x@logFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'settings'
funcname <-'.settings.version'
standGen <- function(object) standardGeneric('.settings.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.settings.version<-'
standGen <- function(x, value) standardGeneric('.settings.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.nlm.gradtol'
standGen <- function(object) standardGeneric('.options.nlm.gradtol')
standMethod <- function(object) object@nlm.gradtol
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.nlm.gradtol<-'
standGen <- function(x, value) standardGeneric('.options.nlm.gradtol<-')
standMethod <- function(x, value) {x@nlm.gradtol<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.nlm.steptol'
standGen <- function(object) standardGeneric('.options.nlm.steptol')
standMethod <- function(object) object@nlm.steptol
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.nlm.steptol<-'
standGen <- function(x, value) standardGeneric('.options.nlm.steptol<-')
standMethod <- function(x, value) {x@nlm.steptol<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.opt.method'
standGen <- function(object) standardGeneric('.options.opt.method')
standMethod <- function(object) object@opt.method
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.opt.method<-'
standGen <- function(x, value) standardGeneric('.options.opt.method<-')
standMethod <- function(x, value) {x@opt.method<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.opt.lower'
standGen <- function(object) standardGeneric('.options.opt.lower')
standMethod <- function(object) object@opt.lower
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.opt.lower<-'
standGen <- function(x, value) standardGeneric('.options.opt.lower<-')
standMethod <- function(x, value) {x@opt.lower<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.opt.upper'
standGen <- function(object) standardGeneric('.options.opt.upper')
standMethod <- function(object) object@opt.upper
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.opt.upper<-'
standGen <- function(x, value) standardGeneric('.options.opt.upper<-')
standMethod <- function(x, value) {x@opt.upper<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.min.analyticalgrad'
standGen <- function(object) standardGeneric('.options.min.analyticalgrad')
standMethod <- function(object) object@min.analyticalgrad
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.min.analyticalgrad<-'
standGen <- function(x, value) standardGeneric('.options.min.analyticalgrad<-')
standMethod <- function(x, value) {x@min.analyticalgrad<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.min.iterlim'
standGen <- function(object) standardGeneric('.options.min.iterlim')
standMethod <- function(object) object@min.iterlim
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.min.iterlim<-'
standGen <- function(x, value) standardGeneric('.options.min.iterlim<-')
standMethod <- function(x, value) {x@min.iterlim<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.min.routine'
standGen <- function(object) standardGeneric('.options.min.routine')
standMethod <- function(object) object@min.routine
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.min.routine<-'
standGen <- function(x, value) standardGeneric('.options.min.routine<-')
standMethod <- function(x, value) {x@min.routine<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.min.boundlim'
standGen <- function(object) standardGeneric('.options.min.boundlim')
standMethod <- function(object) object@min.boundlim
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.min.boundlim<-'
standGen <- function(x, value) standardGeneric('.options.min.boundlim<-')
standMethod <- function(x, value) {x@min.boundlim<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.start.method'
standGen <- function(object) standardGeneric('.options.start.method')
standMethod <- function(object) object@start.method
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.start.method<-'
standGen <- function(x, value) standardGeneric('.options.start.method<-')
standMethod <- function(x, value) {x@start.method<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.start.maxfac'
standGen <- function(object) standardGeneric('.options.start.maxfac')
standMethod <- function(object) object@start.maxfac
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.start.maxfac<-'
standGen <- function(x, value) standardGeneric('.options.start.maxfac<-')
standMethod <- function(x, value) {x@start.maxfac<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.start.vector'
standGen <- function(object) standardGeneric('.options.start.vector')
standMethod <- function(object) object@start.vector
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.start.vector<-'
standGen <- function(x, value) standardGeneric('.options.start.vector<-')
standMethod <- function(x, value) {x@start.vector<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.chk.method'
standGen <- function(object) standardGeneric('.options.chk.method')
standMethod <- function(object) object@chk.method
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.chk.method<-'
standGen <- function(x, value) standardGeneric('.options.chk.method<-')
standMethod <- function(x, value) {x@chk.method<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.chk.range'
standGen <- function(object) standardGeneric('.options.chk.range')
standMethod <- function(object) object@chk.range
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.chk.range<-'
standGen <- function(x, value) standardGeneric('.options.chk.range<-')
standMethod <- function(x, value) {x@chk.range<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.sw.type'
standGen <- function(object) standardGeneric('.options.sw.type')
standMethod <- function(object) object@sw.type
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.sw.type<-'
standGen <- function(x, value) standardGeneric('.options.sw.type<-')
standMethod <- function(x, value) {x@sw.type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.output.mode'
standGen <- function(object) standardGeneric('.options.output.mode')
standMethod <- function(object) object@output.mode
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.output.mode<-'
standGen <- function(x, value) standardGeneric('.options.output.mode<-')
standMethod <- function(x, value) {x@output.mode<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'options'
funcname <-'.options.version'
standGen <- function(object) standardGeneric('.options.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.options.version<-'
standGen <- function(x, value) standardGeneric('.options.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.fullpath'
standGen <- function(object) standardGeneric('.registration.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.fullpath<-'
standGen <- function(x, value) standardGeneric('.registration.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.filename'
standGen <- function(object) standardGeneric('.registration.filename')
standMethod <- function(object) object@filename
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.filename<-'
standGen <- function(x, value) standardGeneric('.registration.filename<-')
standMethod <- function(x, value) {x@filename<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.linkedfile'
standGen <- function(object) standardGeneric('.registration.linkedfile')
standMethod <- function(object) object@linkedfile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.linkedfile<-'
standGen <- function(x, value) standardGeneric('.registration.linkedfile<-')
standMethod <- function(x, value) {x@linkedfile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.examp2high'
standGen <- function(object) standardGeneric('.registration.examp2high')
standMethod <- function(object) object@examp2high
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.examp2high<-'
standGen <- function(x, value) standardGeneric('.registration.examp2high<-')
standMethod <- function(x, value) {x@examp2high<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.high2stand'
standGen <- function(object) standardGeneric('.registration.high2stand')
standMethod <- function(object) object@high2stand
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.high2stand<-'
standGen <- function(x, value) standardGeneric('.registration.high2stand<-')
standMethod <- function(x, value) {x@high2stand<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.examp2stand'
standGen <- function(object) standardGeneric('.registration.examp2stand')
standMethod <- function(object) object@examp2stand
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.examp2stand<-'
standGen <- function(x, value) standardGeneric('.registration.examp2stand<-')
standMethod <- function(x, value) {x@examp2stand<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.example'
standGen <- function(object) standardGeneric('.registration.example')
standMethod <- function(object) object@example
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.example<-'
standGen <- function(x, value) standardGeneric('.registration.example<-')
standMethod <- function(x, value) {x@example<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.highres'
standGen <- function(object) standardGeneric('.registration.highres')
standMethod <- function(object) object@highres
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.highres<-'
standGen <- function(x, value) standardGeneric('.registration.highres<-')
standMethod <- function(x, value) {x@highres<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.standard'
standGen <- function(object) standardGeneric('.registration.standard')
standMethod <- function(object) object@standard
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.standard<-'
standGen <- function(x, value) standardGeneric('.registration.standard<-')
standMethod <- function(x, value) {x@standard<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.Dex'
standGen <- function(object) standardGeneric('.registration.Dex')
standMethod <- function(object) object@Dex
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.Dex<-'
standGen <- function(x, value) standardGeneric('.registration.Dex<-')
standMethod <- function(x, value) {x@Dex<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.Dhi'
standGen <- function(object) standardGeneric('.registration.Dhi')
standMethod <- function(object) object@Dhi
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.Dhi<-'
standGen <- function(x, value) standardGeneric('.registration.Dhi<-')
standMethod <- function(x, value) {x@Dhi<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.Dst'
standGen <- function(object) standardGeneric('.registration.Dst')
standMethod <- function(object) object@Dst
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.Dst<-'
standGen <- function(x, value) standardGeneric('.registration.Dst<-')
standMethod <- function(x, value) {x@Dst<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.SXhi'
standGen <- function(object) standardGeneric('.registration.SXhi')
standMethod <- function(object) object@SXhi
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.SXhi<-'
standGen <- function(x, value) standardGeneric('.registration.SXhi<-')
standMethod <- function(x, value) {x@SXhi<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.Aex2hi'
standGen <- function(object) standardGeneric('.registration.Aex2hi')
standMethod <- function(object) object@Aex2hi
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.Aex2hi<-'
standGen <- function(x, value) standardGeneric('.registration.Aex2hi<-')
standMethod <- function(x, value) {x@Aex2hi<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.Ahi2st'
standGen <- function(object) standardGeneric('.registration.Ahi2st')
standMethod <- function(object) object@Ahi2st
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.Ahi2st<-'
standGen <- function(x, value) standardGeneric('.registration.Ahi2st<-')
standMethod <- function(x, value) {x@Ahi2st<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.OXst'
standGen <- function(object) standardGeneric('.registration.OXst')
standMethod <- function(object) object@OXst
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.OXst<-'
standGen <- function(x, value) standardGeneric('.registration.OXst<-')
standMethod <- function(x, value) {x@OXst<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'registration'
funcname <-'.registration.version'
standGen <- function(object) standardGeneric('.registration.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.registration.version<-'
standGen <- function(x, value) standardGeneric('.registration.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.fullpath'
standGen <- function(object) standardGeneric('.nifti.fileinfo.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.fullpath<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.filename'
standGen <- function(object) standardGeneric('.nifti.fileinfo.filename')
standMethod <- function(object) object@filename
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.filename<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.filename<-')
standMethod <- function(x, value) {x@filename<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.filetype'
standGen <- function(object) standardGeneric('.nifti.fileinfo.filetype')
standMethod <- function(object) object@filetype
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.filetype<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.filetype<-')
standMethod <- function(x, value) {x@filetype<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.extension'
standGen <- function(object) standardGeneric('.nifti.fileinfo.extension')
standMethod <- function(object) object@extension
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.extension<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.extension<-')
standMethod <- function(x, value) {x@extension<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.gzipped'
standGen <- function(object) standardGeneric('.nifti.fileinfo.gzipped')
standMethod <- function(object) object@gzipped
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.gzipped<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.gzipped<-')
standMethod <- function(x, value) {x@gzipped<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.endian'
standGen <- function(object) standardGeneric('.nifti.fileinfo.endian')
standMethod <- function(object) object@endian
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.endian<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.endian<-')
standMethod <- function(x, value) {x@endian<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.fileinfo'
funcname <-'.nifti.fileinfo.version'
standGen <- function(object) standardGeneric('.nifti.fileinfo.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.fileinfo.version<-'
standGen <- function(x, value) standardGeneric('.nifti.fileinfo.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.sizeof_hdr'
standGen <- function(object) standardGeneric('.nifti.header.sizeof_hdr')
standMethod <- function(object) object@sizeof_hdr
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.sizeof_hdr<-'
standGen <- function(x, value) standardGeneric('.nifti.header.sizeof_hdr<-')
standMethod <- function(x, value) {x@sizeof_hdr<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.data_type'
standGen <- function(object) standardGeneric('.nifti.header.data_type')
standMethod <- function(object) object@data_type
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.data_type<-'
standGen <- function(x, value) standardGeneric('.nifti.header.data_type<-')
standMethod <- function(x, value) {x@data_type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.db_name'
standGen <- function(object) standardGeneric('.nifti.header.db_name')
standMethod <- function(object) object@db_name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.db_name<-'
standGen <- function(x, value) standardGeneric('.nifti.header.db_name<-')
standMethod <- function(x, value) {x@db_name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.extents'
standGen <- function(object) standardGeneric('.nifti.header.extents')
standMethod <- function(object) object@extents
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.extents<-'
standGen <- function(x, value) standardGeneric('.nifti.header.extents<-')
standMethod <- function(x, value) {x@extents<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.session_error'
standGen <- function(object) standardGeneric('.nifti.header.session_error')
standMethod <- function(object) object@session_error
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.session_error<-'
standGen <- function(x, value) standardGeneric('.nifti.header.session_error<-')
standMethod <- function(x, value) {x@session_error<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.regular'
standGen <- function(object) standardGeneric('.nifti.header.regular')
standMethod <- function(object) object@regular
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.regular<-'
standGen <- function(x, value) standardGeneric('.nifti.header.regular<-')
standMethod <- function(x, value) {x@regular<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.dim_info'
standGen <- function(object) standardGeneric('.nifti.header.dim_info')
standMethod <- function(object) object@dim_info
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.dim_info<-'
standGen <- function(x, value) standardGeneric('.nifti.header.dim_info<-')
standMethod <- function(x, value) {x@dim_info<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.dims'
standGen <- function(object) standardGeneric('.nifti.header.dims')
standMethod <- function(object) object@dims
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.dims<-'
standGen <- function(x, value) standardGeneric('.nifti.header.dims<-')
standMethod <- function(x, value) {x@dims<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.intent_p1'
standGen <- function(object) standardGeneric('.nifti.header.intent_p1')
standMethod <- function(object) object@intent_p1
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.intent_p1<-'
standGen <- function(x, value) standardGeneric('.nifti.header.intent_p1<-')
standMethod <- function(x, value) {x@intent_p1<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.intent_p2'
standGen <- function(object) standardGeneric('.nifti.header.intent_p2')
standMethod <- function(object) object@intent_p2
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.intent_p2<-'
standGen <- function(x, value) standardGeneric('.nifti.header.intent_p2<-')
standMethod <- function(x, value) {x@intent_p2<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.intent_p3'
standGen <- function(object) standardGeneric('.nifti.header.intent_p3')
standMethod <- function(object) object@intent_p3
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.intent_p3<-'
standGen <- function(x, value) standardGeneric('.nifti.header.intent_p3<-')
standMethod <- function(x, value) {x@intent_p3<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.intent_code'
standGen <- function(object) standardGeneric('.nifti.header.intent_code')
standMethod <- function(object) object@intent_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.intent_code<-'
standGen <- function(x, value) standardGeneric('.nifti.header.intent_code<-')
standMethod <- function(x, value) {x@intent_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.datatype'
standGen <- function(object) standardGeneric('.nifti.header.datatype')
standMethod <- function(object) object@datatype
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.datatype<-'
standGen <- function(x, value) standardGeneric('.nifti.header.datatype<-')
standMethod <- function(x, value) {x@datatype<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.bitpix'
standGen <- function(object) standardGeneric('.nifti.header.bitpix')
standMethod <- function(object) object@bitpix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.bitpix<-'
standGen <- function(x, value) standardGeneric('.nifti.header.bitpix<-')
standMethod <- function(x, value) {x@bitpix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.slice_start'
standGen <- function(object) standardGeneric('.nifti.header.slice_start')
standMethod <- function(object) object@slice_start
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.slice_start<-'
standGen <- function(x, value) standardGeneric('.nifti.header.slice_start<-')
standMethod <- function(x, value) {x@slice_start<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.pixdim'
standGen <- function(object) standardGeneric('.nifti.header.pixdim')
standMethod <- function(object) object@pixdim
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.pixdim<-'
standGen <- function(x, value) standardGeneric('.nifti.header.pixdim<-')
standMethod <- function(x, value) {x@pixdim<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.vox_offset'
standGen <- function(object) standardGeneric('.nifti.header.vox_offset')
standMethod <- function(object) object@vox_offset
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.vox_offset<-'
standGen <- function(x, value) standardGeneric('.nifti.header.vox_offset<-')
standMethod <- function(x, value) {x@vox_offset<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.scl_slope'
standGen <- function(object) standardGeneric('.nifti.header.scl_slope')
standMethod <- function(object) object@scl_slope
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.scl_slope<-'
standGen <- function(x, value) standardGeneric('.nifti.header.scl_slope<-')
standMethod <- function(x, value) {x@scl_slope<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.scl_inter'
standGen <- function(object) standardGeneric('.nifti.header.scl_inter')
standMethod <- function(object) object@scl_inter
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.scl_inter<-'
standGen <- function(x, value) standardGeneric('.nifti.header.scl_inter<-')
standMethod <- function(x, value) {x@scl_inter<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.slice_end'
standGen <- function(object) standardGeneric('.nifti.header.slice_end')
standMethod <- function(object) object@slice_end
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.slice_end<-'
standGen <- function(x, value) standardGeneric('.nifti.header.slice_end<-')
standMethod <- function(x, value) {x@slice_end<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.slice_code'
standGen <- function(object) standardGeneric('.nifti.header.slice_code')
standMethod <- function(object) object@slice_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.slice_code<-'
standGen <- function(x, value) standardGeneric('.nifti.header.slice_code<-')
standMethod <- function(x, value) {x@slice_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.xyzt_units'
standGen <- function(object) standardGeneric('.nifti.header.xyzt_units')
standMethod <- function(object) object@xyzt_units
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.xyzt_units<-'
standGen <- function(x, value) standardGeneric('.nifti.header.xyzt_units<-')
standMethod <- function(x, value) {x@xyzt_units<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.cal_max'
standGen <- function(object) standardGeneric('.nifti.header.cal_max')
standMethod <- function(object) object@cal_max
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.cal_max<-'
standGen <- function(x, value) standardGeneric('.nifti.header.cal_max<-')
standMethod <- function(x, value) {x@cal_max<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.cal_min'
standGen <- function(object) standardGeneric('.nifti.header.cal_min')
standMethod <- function(object) object@cal_min
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.cal_min<-'
standGen <- function(x, value) standardGeneric('.nifti.header.cal_min<-')
standMethod <- function(x, value) {x@cal_min<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.slice_duration'
standGen <- function(object) standardGeneric('.nifti.header.slice_duration')
standMethod <- function(object) object@slice_duration
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.slice_duration<-'
standGen <- function(x, value) standardGeneric('.nifti.header.slice_duration<-')
standMethod <- function(x, value) {x@slice_duration<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.toffset'
standGen <- function(object) standardGeneric('.nifti.header.toffset')
standMethod <- function(object) object@toffset
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.toffset<-'
standGen <- function(x, value) standardGeneric('.nifti.header.toffset<-')
standMethod <- function(x, value) {x@toffset<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.glmax'
standGen <- function(object) standardGeneric('.nifti.header.glmax')
standMethod <- function(object) object@glmax
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.glmax<-'
standGen <- function(x, value) standardGeneric('.nifti.header.glmax<-')
standMethod <- function(x, value) {x@glmax<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.glmin'
standGen <- function(object) standardGeneric('.nifti.header.glmin')
standMethod <- function(object) object@glmin
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.glmin<-'
standGen <- function(x, value) standardGeneric('.nifti.header.glmin<-')
standMethod <- function(x, value) {x@glmin<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.descrip'
standGen <- function(object) standardGeneric('.nifti.header.descrip')
standMethod <- function(object) object@descrip
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.descrip<-'
standGen <- function(x, value) standardGeneric('.nifti.header.descrip<-')
standMethod <- function(x, value) {x@descrip<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.aux_file'
standGen <- function(object) standardGeneric('.nifti.header.aux_file')
standMethod <- function(object) object@aux_file
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.aux_file<-'
standGen <- function(x, value) standardGeneric('.nifti.header.aux_file<-')
standMethod <- function(x, value) {x@aux_file<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.qform_code'
standGen <- function(object) standardGeneric('.nifti.header.qform_code')
standMethod <- function(object) object@qform_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.qform_code<-'
standGen <- function(x, value) standardGeneric('.nifti.header.qform_code<-')
standMethod <- function(x, value) {x@qform_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.sform_code'
standGen <- function(object) standardGeneric('.nifti.header.sform_code')
standMethod <- function(object) object@sform_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.sform_code<-'
standGen <- function(x, value) standardGeneric('.nifti.header.sform_code<-')
standMethod <- function(x, value) {x@sform_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.quatern_b'
standGen <- function(object) standardGeneric('.nifti.header.quatern_b')
standMethod <- function(object) object@quatern_b
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.quatern_b<-'
standGen <- function(x, value) standardGeneric('.nifti.header.quatern_b<-')
standMethod <- function(x, value) {x@quatern_b<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.quatern_c'
standGen <- function(object) standardGeneric('.nifti.header.quatern_c')
standMethod <- function(object) object@quatern_c
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.quatern_c<-'
standGen <- function(x, value) standardGeneric('.nifti.header.quatern_c<-')
standMethod <- function(x, value) {x@quatern_c<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.quatern_d'
standGen <- function(object) standardGeneric('.nifti.header.quatern_d')
standMethod <- function(object) object@quatern_d
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.quatern_d<-'
standGen <- function(x, value) standardGeneric('.nifti.header.quatern_d<-')
standMethod <- function(x, value) {x@quatern_d<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.qoffset_x'
standGen <- function(object) standardGeneric('.nifti.header.qoffset_x')
standMethod <- function(object) object@qoffset_x
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.qoffset_x<-'
standGen <- function(x, value) standardGeneric('.nifti.header.qoffset_x<-')
standMethod <- function(x, value) {x@qoffset_x<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.qoffset_y'
standGen <- function(object) standardGeneric('.nifti.header.qoffset_y')
standMethod <- function(object) object@qoffset_y
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.qoffset_y<-'
standGen <- function(x, value) standardGeneric('.nifti.header.qoffset_y<-')
standMethod <- function(x, value) {x@qoffset_y<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.qoffset_z'
standGen <- function(object) standardGeneric('.nifti.header.qoffset_z')
standMethod <- function(object) object@qoffset_z
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.qoffset_z<-'
standGen <- function(x, value) standardGeneric('.nifti.header.qoffset_z<-')
standMethod <- function(x, value) {x@qoffset_z<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.srow_x'
standGen <- function(object) standardGeneric('.nifti.header.srow_x')
standMethod <- function(object) object@srow_x
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.srow_x<-'
standGen <- function(x, value) standardGeneric('.nifti.header.srow_x<-')
standMethod <- function(x, value) {x@srow_x<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.srow_y'
standGen <- function(object) standardGeneric('.nifti.header.srow_y')
standMethod <- function(object) object@srow_y
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.srow_y<-'
standGen <- function(x, value) standardGeneric('.nifti.header.srow_y<-')
standMethod <- function(x, value) {x@srow_y<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.srow_z'
standGen <- function(object) standardGeneric('.nifti.header.srow_z')
standMethod <- function(object) object@srow_z
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.srow_z<-'
standGen <- function(x, value) standardGeneric('.nifti.header.srow_z<-')
standMethod <- function(x, value) {x@srow_z<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.intent_name'
standGen <- function(object) standardGeneric('.nifti.header.intent_name')
standMethod <- function(object) object@intent_name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.intent_name<-'
standGen <- function(x, value) standardGeneric('.nifti.header.intent_name<-')
standMethod <- function(x, value) {x@intent_name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.magic'
standGen <- function(object) standardGeneric('.nifti.header.magic')
standMethod <- function(object) object@magic
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.magic<-'
standGen <- function(x, value) standardGeneric('.nifti.header.magic<-')
standMethod <- function(x, value) {x@magic<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.data.type'
standGen <- function(object) standardGeneric('.nifti.header.data.type')
standMethod <- function(object) object@data.type
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.data.type<-'
standGen <- function(x, value) standardGeneric('.nifti.header.data.type<-')
standMethod <- function(x, value) {x@data.type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.data.signed'
standGen <- function(object) standardGeneric('.nifti.header.data.signed')
standMethod <- function(object) object@data.signed
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.data.signed<-'
standGen <- function(x, value) standardGeneric('.nifti.header.data.signed<-')
standMethod <- function(x, value) {x@data.signed<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.fullpath'
standGen <- function(object) standardGeneric('.nifti.header.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.fullpath<-'
standGen <- function(x, value) standardGeneric('.nifti.header.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.filename'
standGen <- function(object) standardGeneric('.nifti.header.filename')
standMethod <- function(object) object@filename
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.filename<-'
standGen <- function(x, value) standardGeneric('.nifti.header.filename<-')
standMethod <- function(x, value) {x@filename<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.filetype'
standGen <- function(object) standardGeneric('.nifti.header.filetype')
standMethod <- function(object) object@filetype
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.filetype<-'
standGen <- function(x, value) standardGeneric('.nifti.header.filetype<-')
standMethod <- function(x, value) {x@filetype<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.extension'
standGen <- function(object) standardGeneric('.nifti.header.extension')
standMethod <- function(object) object@extension
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.extension<-'
standGen <- function(x, value) standardGeneric('.nifti.header.extension<-')
standMethod <- function(x, value) {x@extension<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.gzipped'
standGen <- function(object) standardGeneric('.nifti.header.gzipped')
standMethod <- function(object) object@gzipped
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.gzipped<-'
standGen <- function(x, value) standardGeneric('.nifti.header.gzipped<-')
standMethod <- function(x, value) {x@gzipped<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.endian'
standGen <- function(object) standardGeneric('.nifti.header.endian')
standMethod <- function(object) object@endian
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.endian<-'
standGen <- function(x, value) standardGeneric('.nifti.header.endian<-')
standMethod <- function(x, value) {x@endian<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'nifti.header'
funcname <-'.nifti.header.version'
standGen <- function(object) standardGeneric('.nifti.header.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.nifti.header.version<-'
standGen <- function(x, value) standardGeneric('.nifti.header.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.datavec'
standGen <- function(object) standardGeneric('.fmri.data.datavec')
standMethod <- function(object) object@datavec
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.datavec<-'
standGen <- function(x, value) standardGeneric('.fmri.data.datavec<-')
standMethod <- function(x, value) {x@datavec<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.sizeof_hdr'
standGen <- function(object) standardGeneric('.fmri.data.sizeof_hdr')
standMethod <- function(object) object@sizeof_hdr
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.sizeof_hdr<-'
standGen <- function(x, value) standardGeneric('.fmri.data.sizeof_hdr<-')
standMethod <- function(x, value) {x@sizeof_hdr<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.data_type'
standGen <- function(object) standardGeneric('.fmri.data.data_type')
standMethod <- function(object) object@data_type
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.data_type<-'
standGen <- function(x, value) standardGeneric('.fmri.data.data_type<-')
standMethod <- function(x, value) {x@data_type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.db_name'
standGen <- function(object) standardGeneric('.fmri.data.db_name')
standMethod <- function(object) object@db_name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.db_name<-'
standGen <- function(x, value) standardGeneric('.fmri.data.db_name<-')
standMethod <- function(x, value) {x@db_name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.extents'
standGen <- function(object) standardGeneric('.fmri.data.extents')
standMethod <- function(object) object@extents
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.extents<-'
standGen <- function(x, value) standardGeneric('.fmri.data.extents<-')
standMethod <- function(x, value) {x@extents<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.session_error'
standGen <- function(object) standardGeneric('.fmri.data.session_error')
standMethod <- function(object) object@session_error
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.session_error<-'
standGen <- function(x, value) standardGeneric('.fmri.data.session_error<-')
standMethod <- function(x, value) {x@session_error<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.regular'
standGen <- function(object) standardGeneric('.fmri.data.regular')
standMethod <- function(object) object@regular
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.regular<-'
standGen <- function(x, value) standardGeneric('.fmri.data.regular<-')
standMethod <- function(x, value) {x@regular<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.dim_info'
standGen <- function(object) standardGeneric('.fmri.data.dim_info')
standMethod <- function(object) object@dim_info
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.dim_info<-'
standGen <- function(x, value) standardGeneric('.fmri.data.dim_info<-')
standMethod <- function(x, value) {x@dim_info<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.dims'
standGen <- function(object) standardGeneric('.fmri.data.dims')
standMethod <- function(object) object@dims
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.dims<-'
standGen <- function(x, value) standardGeneric('.fmri.data.dims<-')
standMethod <- function(x, value) {x@dims<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.intent_p1'
standGen <- function(object) standardGeneric('.fmri.data.intent_p1')
standMethod <- function(object) object@intent_p1
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.intent_p1<-'
standGen <- function(x, value) standardGeneric('.fmri.data.intent_p1<-')
standMethod <- function(x, value) {x@intent_p1<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.intent_p2'
standGen <- function(object) standardGeneric('.fmri.data.intent_p2')
standMethod <- function(object) object@intent_p2
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.intent_p2<-'
standGen <- function(x, value) standardGeneric('.fmri.data.intent_p2<-')
standMethod <- function(x, value) {x@intent_p2<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.intent_p3'
standGen <- function(object) standardGeneric('.fmri.data.intent_p3')
standMethod <- function(object) object@intent_p3
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.intent_p3<-'
standGen <- function(x, value) standardGeneric('.fmri.data.intent_p3<-')
standMethod <- function(x, value) {x@intent_p3<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.intent_code'
standGen <- function(object) standardGeneric('.fmri.data.intent_code')
standMethod <- function(object) object@intent_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.intent_code<-'
standGen <- function(x, value) standardGeneric('.fmri.data.intent_code<-')
standMethod <- function(x, value) {x@intent_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.datatype'
standGen <- function(object) standardGeneric('.fmri.data.datatype')
standMethod <- function(object) object@datatype
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.datatype<-'
standGen <- function(x, value) standardGeneric('.fmri.data.datatype<-')
standMethod <- function(x, value) {x@datatype<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.bitpix'
standGen <- function(object) standardGeneric('.fmri.data.bitpix')
standMethod <- function(object) object@bitpix
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.bitpix<-'
standGen <- function(x, value) standardGeneric('.fmri.data.bitpix<-')
standMethod <- function(x, value) {x@bitpix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.slice_start'
standGen <- function(object) standardGeneric('.fmri.data.slice_start')
standMethod <- function(object) object@slice_start
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.slice_start<-'
standGen <- function(x, value) standardGeneric('.fmri.data.slice_start<-')
standMethod <- function(x, value) {x@slice_start<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.pixdim'
standGen <- function(object) standardGeneric('.fmri.data.pixdim')
standMethod <- function(object) object@pixdim
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.pixdim<-'
standGen <- function(x, value) standardGeneric('.fmri.data.pixdim<-')
standMethod <- function(x, value) {x@pixdim<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.vox_offset'
standGen <- function(object) standardGeneric('.fmri.data.vox_offset')
standMethod <- function(object) object@vox_offset
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.vox_offset<-'
standGen <- function(x, value) standardGeneric('.fmri.data.vox_offset<-')
standMethod <- function(x, value) {x@vox_offset<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.scl_slope'
standGen <- function(object) standardGeneric('.fmri.data.scl_slope')
standMethod <- function(object) object@scl_slope
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.scl_slope<-'
standGen <- function(x, value) standardGeneric('.fmri.data.scl_slope<-')
standMethod <- function(x, value) {x@scl_slope<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.scl_inter'
standGen <- function(object) standardGeneric('.fmri.data.scl_inter')
standMethod <- function(object) object@scl_inter
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.scl_inter<-'
standGen <- function(x, value) standardGeneric('.fmri.data.scl_inter<-')
standMethod <- function(x, value) {x@scl_inter<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.slice_end'
standGen <- function(object) standardGeneric('.fmri.data.slice_end')
standMethod <- function(object) object@slice_end
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.slice_end<-'
standGen <- function(x, value) standardGeneric('.fmri.data.slice_end<-')
standMethod <- function(x, value) {x@slice_end<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.slice_code'
standGen <- function(object) standardGeneric('.fmri.data.slice_code')
standMethod <- function(object) object@slice_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.slice_code<-'
standGen <- function(x, value) standardGeneric('.fmri.data.slice_code<-')
standMethod <- function(x, value) {x@slice_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.xyzt_units'
standGen <- function(object) standardGeneric('.fmri.data.xyzt_units')
standMethod <- function(object) object@xyzt_units
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.xyzt_units<-'
standGen <- function(x, value) standardGeneric('.fmri.data.xyzt_units<-')
standMethod <- function(x, value) {x@xyzt_units<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.cal_max'
standGen <- function(object) standardGeneric('.fmri.data.cal_max')
standMethod <- function(object) object@cal_max
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.cal_max<-'
standGen <- function(x, value) standardGeneric('.fmri.data.cal_max<-')
standMethod <- function(x, value) {x@cal_max<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.cal_min'
standGen <- function(object) standardGeneric('.fmri.data.cal_min')
standMethod <- function(object) object@cal_min
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.cal_min<-'
standGen <- function(x, value) standardGeneric('.fmri.data.cal_min<-')
standMethod <- function(x, value) {x@cal_min<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.slice_duration'
standGen <- function(object) standardGeneric('.fmri.data.slice_duration')
standMethod <- function(object) object@slice_duration
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.slice_duration<-'
standGen <- function(x, value) standardGeneric('.fmri.data.slice_duration<-')
standMethod <- function(x, value) {x@slice_duration<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.toffset'
standGen <- function(object) standardGeneric('.fmri.data.toffset')
standMethod <- function(object) object@toffset
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.toffset<-'
standGen <- function(x, value) standardGeneric('.fmri.data.toffset<-')
standMethod <- function(x, value) {x@toffset<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.glmax'
standGen <- function(object) standardGeneric('.fmri.data.glmax')
standMethod <- function(object) object@glmax
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.glmax<-'
standGen <- function(x, value) standardGeneric('.fmri.data.glmax<-')
standMethod <- function(x, value) {x@glmax<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.glmin'
standGen <- function(object) standardGeneric('.fmri.data.glmin')
standMethod <- function(object) object@glmin
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.glmin<-'
standGen <- function(x, value) standardGeneric('.fmri.data.glmin<-')
standMethod <- function(x, value) {x@glmin<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.descrip'
standGen <- function(object) standardGeneric('.fmri.data.descrip')
standMethod <- function(object) object@descrip
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.descrip<-'
standGen <- function(x, value) standardGeneric('.fmri.data.descrip<-')
standMethod <- function(x, value) {x@descrip<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.aux_file'
standGen <- function(object) standardGeneric('.fmri.data.aux_file')
standMethod <- function(object) object@aux_file
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.aux_file<-'
standGen <- function(x, value) standardGeneric('.fmri.data.aux_file<-')
standMethod <- function(x, value) {x@aux_file<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.qform_code'
standGen <- function(object) standardGeneric('.fmri.data.qform_code')
standMethod <- function(object) object@qform_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.qform_code<-'
standGen <- function(x, value) standardGeneric('.fmri.data.qform_code<-')
standMethod <- function(x, value) {x@qform_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.sform_code'
standGen <- function(object) standardGeneric('.fmri.data.sform_code')
standMethod <- function(object) object@sform_code
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.sform_code<-'
standGen <- function(x, value) standardGeneric('.fmri.data.sform_code<-')
standMethod <- function(x, value) {x@sform_code<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.quatern_b'
standGen <- function(object) standardGeneric('.fmri.data.quatern_b')
standMethod <- function(object) object@quatern_b
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.quatern_b<-'
standGen <- function(x, value) standardGeneric('.fmri.data.quatern_b<-')
standMethod <- function(x, value) {x@quatern_b<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.quatern_c'
standGen <- function(object) standardGeneric('.fmri.data.quatern_c')
standMethod <- function(object) object@quatern_c
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.quatern_c<-'
standGen <- function(x, value) standardGeneric('.fmri.data.quatern_c<-')
standMethod <- function(x, value) {x@quatern_c<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.quatern_d'
standGen <- function(object) standardGeneric('.fmri.data.quatern_d')
standMethod <- function(object) object@quatern_d
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.quatern_d<-'
standGen <- function(x, value) standardGeneric('.fmri.data.quatern_d<-')
standMethod <- function(x, value) {x@quatern_d<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.qoffset_x'
standGen <- function(object) standardGeneric('.fmri.data.qoffset_x')
standMethod <- function(object) object@qoffset_x
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.qoffset_x<-'
standGen <- function(x, value) standardGeneric('.fmri.data.qoffset_x<-')
standMethod <- function(x, value) {x@qoffset_x<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.qoffset_y'
standGen <- function(object) standardGeneric('.fmri.data.qoffset_y')
standMethod <- function(object) object@qoffset_y
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.qoffset_y<-'
standGen <- function(x, value) standardGeneric('.fmri.data.qoffset_y<-')
standMethod <- function(x, value) {x@qoffset_y<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.qoffset_z'
standGen <- function(object) standardGeneric('.fmri.data.qoffset_z')
standMethod <- function(object) object@qoffset_z
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.qoffset_z<-'
standGen <- function(x, value) standardGeneric('.fmri.data.qoffset_z<-')
standMethod <- function(x, value) {x@qoffset_z<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.srow_x'
standGen <- function(object) standardGeneric('.fmri.data.srow_x')
standMethod <- function(object) object@srow_x
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.srow_x<-'
standGen <- function(x, value) standardGeneric('.fmri.data.srow_x<-')
standMethod <- function(x, value) {x@srow_x<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.srow_y'
standGen <- function(object) standardGeneric('.fmri.data.srow_y')
standMethod <- function(object) object@srow_y
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.srow_y<-'
standGen <- function(x, value) standardGeneric('.fmri.data.srow_y<-')
standMethod <- function(x, value) {x@srow_y<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.srow_z'
standGen <- function(object) standardGeneric('.fmri.data.srow_z')
standMethod <- function(object) object@srow_z
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.srow_z<-'
standGen <- function(x, value) standardGeneric('.fmri.data.srow_z<-')
standMethod <- function(x, value) {x@srow_z<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.intent_name'
standGen <- function(object) standardGeneric('.fmri.data.intent_name')
standMethod <- function(object) object@intent_name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.intent_name<-'
standGen <- function(x, value) standardGeneric('.fmri.data.intent_name<-')
standMethod <- function(x, value) {x@intent_name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.magic'
standGen <- function(object) standardGeneric('.fmri.data.magic')
standMethod <- function(object) object@magic
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.magic<-'
standGen <- function(x, value) standardGeneric('.fmri.data.magic<-')
standMethod <- function(x, value) {x@magic<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.data.type'
standGen <- function(object) standardGeneric('.fmri.data.data.type')
standMethod <- function(object) object@data.type
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.data.type<-'
standGen <- function(x, value) standardGeneric('.fmri.data.data.type<-')
standMethod <- function(x, value) {x@data.type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.data.signed'
standGen <- function(object) standardGeneric('.fmri.data.data.signed')
standMethod <- function(object) object@data.signed
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.data.signed<-'
standGen <- function(x, value) standardGeneric('.fmri.data.data.signed<-')
standMethod <- function(x, value) {x@data.signed<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.fullpath'
standGen <- function(object) standardGeneric('.fmri.data.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.fullpath<-'
standGen <- function(x, value) standardGeneric('.fmri.data.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.filename'
standGen <- function(object) standardGeneric('.fmri.data.filename')
standMethod <- function(object) object@filename
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.filename<-'
standGen <- function(x, value) standardGeneric('.fmri.data.filename<-')
standMethod <- function(x, value) {x@filename<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.filetype'
standGen <- function(object) standardGeneric('.fmri.data.filetype')
standMethod <- function(object) object@filetype
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.filetype<-'
standGen <- function(x, value) standardGeneric('.fmri.data.filetype<-')
standMethod <- function(x, value) {x@filetype<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.extension'
standGen <- function(object) standardGeneric('.fmri.data.extension')
standMethod <- function(object) object@extension
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.extension<-'
standGen <- function(x, value) standardGeneric('.fmri.data.extension<-')
standMethod <- function(x, value) {x@extension<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.gzipped'
standGen <- function(object) standardGeneric('.fmri.data.gzipped')
standMethod <- function(object) object@gzipped
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.gzipped<-'
standGen <- function(x, value) standardGeneric('.fmri.data.gzipped<-')
standMethod <- function(x, value) {x@gzipped<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.endian'
standGen <- function(object) standardGeneric('.fmri.data.endian')
standMethod <- function(object) object@endian
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.endian<-'
standGen <- function(x, value) standardGeneric('.fmri.data.endian<-')
standMethod <- function(x, value) {x@endian<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'fmri.data'
funcname <-'.fmri.data.version'
standGen <- function(object) standardGeneric('.fmri.data.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.fmri.data.version<-'
standGen <- function(x, value) standardGeneric('.fmri.data.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'wald'
funcname <-'.wald.consts'
standGen <- function(object) standardGeneric('.wald.consts')
standMethod <- function(object) object@consts
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.wald.consts<-'
standGen <- function(x, value) standardGeneric('.wald.consts<-')
standMethod <- function(x, value) {x@consts<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'wald'
funcname <-'.wald.stats'
standGen <- function(object) standardGeneric('.wald.stats')
standMethod <- function(object) object@stats
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.wald.stats<-'
standGen <- function(x, value) standardGeneric('.wald.stats<-')
standMethod <- function(x, value) {x@stats<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'wald'
funcname <-'.wald.df1'
standGen <- function(object) standardGeneric('.wald.df1')
standMethod <- function(object) object@df1
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.wald.df1<-'
standGen <- function(x, value) standardGeneric('.wald.df1<-')
standMethod <- function(x, value) {x@df1<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'wald'
funcname <-'.wald.df2'
standGen <- function(object) standardGeneric('.wald.df2')
standMethod <- function(object) object@df2
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.wald.df2<-'
standGen <- function(x, value) standardGeneric('.wald.df2<-')
standMethod <- function(x, value) {x@df2<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'wald'
funcname <-'.wald.pvalues'
standGen <- function(object) standardGeneric('.wald.pvalues')
standMethod <- function(object) object@pvalues
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.wald.pvalues<-'
standGen <- function(x, value) standardGeneric('.wald.pvalues<-')
standMethod <- function(x, value) {x@pvalues<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'functional'
funcname <-'.functional.fullpath'
standGen <- function(object) standardGeneric('.functional.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.functional.fullpath<-'
standGen <- function(x, value) standardGeneric('.functional.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'functional'
funcname <-'.functional.functionaldata'
standGen <- function(object) standardGeneric('.functional.functionaldata')
standMethod <- function(object) object@functionaldata
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.functional.functionaldata<-'
standGen <- function(x, value) standardGeneric('.functional.functionaldata<-')
standMethod <- function(x, value) {x@functionaldata<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'functional'
funcname <-'.functional.filename'
standGen <- function(object) standardGeneric('.functional.filename')
standMethod <- function(object) object@filename
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.functional.filename<-'
standGen <- function(x, value) standardGeneric('.functional.filename<-')
standMethod <- function(x, value) {x@filename<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'functional'
funcname <-'.functional.linkedfiles'
standGen <- function(object) standardGeneric('.functional.linkedfiles')
standMethod <- function(object) object@linkedfiles
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.functional.linkedfiles<-'
standGen <- function(x, value) standardGeneric('.functional.linkedfiles<-')
standMethod <- function(x, value) {x@linkedfiles<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'functional'
funcname <-'.functional.timings'
standGen <- function(object) standardGeneric('.functional.timings')
standMethod <- function(object) object@timings
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.functional.timings<-'
standGen <- function(x, value) standardGeneric('.functional.timings<-')
standMethod <- function(x, value) {x@timings<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'functional'
funcname <-'.functional.version'
standGen <- function(object) standardGeneric('.functional.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.functional.version<-'
standGen <- function(x, value) standardGeneric('.functional.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.modelname'
standGen <- function(object) standardGeneric('.model.modelname')
standMethod <- function(object) object@modelname
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.modelname<-'
standGen <- function(x, value) standardGeneric('.model.modelname<-')
standMethod <- function(x, value) {x@modelname<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.modelpath'
standGen <- function(object) standardGeneric('.model.modelpath')
standMethod <- function(object) object@modelpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.modelpath<-'
standGen <- function(x, value) standardGeneric('.model.modelpath<-')
standMethod <- function(x, value) {x@modelpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.modeldatapath'
standGen <- function(object) standardGeneric('.model.modeldatapath')
standMethod <- function(object) object@modeldatapath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.modeldatapath<-'
standGen <- function(x, value) standardGeneric('.model.modeldatapath<-')
standMethod <- function(x, value) {x@modeldatapath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.residualFile'
standGen <- function(object) standardGeneric('.model.residualFile')
standMethod <- function(object) object@residualFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.residualFile<-'
standGen <- function(x, value) standardGeneric('.model.residualFile<-')
standMethod <- function(x, value) {x@residualFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.derivativeFile'
standGen <- function(object) standardGeneric('.model.derivativeFile')
standMethod <- function(object) object@derivativeFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.derivativeFile<-'
standGen <- function(x, value) standardGeneric('.model.derivativeFile<-')
standMethod <- function(x, value) {x@derivativeFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.weightFile'
standGen <- function(object) standardGeneric('.model.weightFile')
standMethod <- function(object) object@weightFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.weightFile<-'
standGen <- function(x, value) standardGeneric('.model.weightFile<-')
standMethod <- function(x, value) {x@weightFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.modelDataFile'
standGen <- function(object) standardGeneric('.model.modelDataFile')
standMethod <- function(object) object@modelDataFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.modelDataFile<-'
standGen <- function(x, value) standardGeneric('.model.modelDataFile<-')
standMethod <- function(x, value) {x@modelDataFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.fullmodelDataFile'
standGen <- function(object) standardGeneric('.model.fullmodelDataFile')
standMethod <- function(object) object@fullmodelDataFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.fullmodelDataFile<-'
standGen <- function(x, value) standardGeneric('.model.fullmodelDataFile<-')
standMethod <- function(x, value) {x@fullmodelDataFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.modelFile'
standGen <- function(object) standardGeneric('.model.modelFile')
standMethod <- function(object) object@modelFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.modelFile<-'
standGen <- function(x, value) standardGeneric('.model.modelFile<-')
standMethod <- function(x, value) {x@modelFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.optionsFile'
standGen <- function(object) standardGeneric('.model.optionsFile')
standMethod <- function(object) object@optionsFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.optionsFile<-'
standGen <- function(x, value) standardGeneric('.model.optionsFile<-')
standMethod <- function(x, value) {x@optionsFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.startFile'
standGen <- function(object) standardGeneric('.model.startFile')
standMethod <- function(object) object@startFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.startFile<-'
standGen <- function(x, value) standardGeneric('.model.startFile<-')
standMethod <- function(x, value) {x@startFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.logFile'
standGen <- function(object) standardGeneric('.model.logFile')
standMethod <- function(object) object@logFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.logFile<-'
standGen <- function(x, value) standardGeneric('.model.logFile<-')
standMethod <- function(x, value) {x@logFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.convergence'
standGen <- function(object) standardGeneric('.model.convergence')
standMethod <- function(object) object@convergence
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.convergence<-'
standGen <- function(x, value) standardGeneric('.model.convergence<-')
standMethod <- function(x, value) {x@convergence<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.iterates'
standGen <- function(object) standardGeneric('.model.iterates')
standMethod <- function(object) object@iterates
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.iterates<-'
standGen <- function(x, value) standardGeneric('.model.iterates<-')
standMethod <- function(x, value) {x@iterates<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.minimum'
standGen <- function(object) standardGeneric('.model.minimum')
standMethod <- function(object) object@minimum
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.minimum<-'
standGen <- function(x, value) standardGeneric('.model.minimum<-')
standMethod <- function(x, value) {x@minimum<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.estimates'
standGen <- function(object) standardGeneric('.model.estimates')
standMethod <- function(object) object@estimates
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.estimates<-'
standGen <- function(x, value) standardGeneric('.model.estimates<-')
standMethod <- function(x, value) {x@estimates<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.gradient'
standGen <- function(object) standardGeneric('.model.gradient')
standMethod <- function(object) object@gradient
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.gradient<-'
standGen <- function(x, value) standardGeneric('.model.gradient<-')
standMethod <- function(x, value) {x@gradient<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.hessian'
standGen <- function(object) standardGeneric('.model.hessian')
standMethod <- function(object) object@hessian
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.hessian<-'
standGen <- function(x, value) standardGeneric('.model.hessian<-')
standMethod <- function(x, value) {x@hessian<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.params'
standGen <- function(object) standardGeneric('.model.params')
standMethod <- function(object) object@params
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.params<-'
standGen <- function(x, value) standardGeneric('.model.params<-')
standMethod <- function(x, value) {x@params<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.modeltype'
standGen <- function(object) standardGeneric('.model.modeltype')
standMethod <- function(object) object@modeltype
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.modeltype<-'
standGen <- function(x, value) standardGeneric('.model.modeltype<-')
standMethod <- function(x, value) {x@modeltype<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.sandwichmethod'
standGen <- function(object) standardGeneric('.model.sandwichmethod')
standMethod <- function(object) object@sandwichmethod
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.sandwichmethod<-'
standGen <- function(x, value) standardGeneric('.model.sandwichmethod<-')
standMethod <- function(x, value) {x@sandwichmethod<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.varcov'
standGen <- function(object) standardGeneric('.model.varcov')
standMethod <- function(object) object@varcov
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.varcov<-'
standGen <- function(x, value) standardGeneric('.model.varcov<-')
standMethod <- function(x, value) {x@varcov<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.warnings'
standGen <- function(object) standardGeneric('.model.warnings')
standMethod <- function(object) object@warnings
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.warnings<-'
standGen <- function(x, value) standardGeneric('.model.warnings<-')
standMethod <- function(x, value) {x@warnings<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.fit'
standGen <- function(object) standardGeneric('.model.fit')
standMethod <- function(object) object@fit
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.fit<-'
standGen <- function(x, value) standardGeneric('.model.fit<-')
standMethod <- function(x, value) {x@fit<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.wald'
standGen <- function(object) standardGeneric('.model.wald')
standMethod <- function(object) object@wald
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.wald<-'
standGen <- function(x, value) standardGeneric('.model.wald<-')
standMethod <- function(x, value) {x@wald<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.regions'
standGen <- function(object) standardGeneric('.model.regions')
standMethod <- function(object) object@regions
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.regions<-'
standGen <- function(x, value) standardGeneric('.model.regions<-')
standMethod <- function(x, value) {x@regions<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.startval'
standGen <- function(object) standardGeneric('.model.startval')
standMethod <- function(object) object@startval
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.startval<-'
standGen <- function(x, value) standardGeneric('.model.startval<-')
standMethod <- function(x, value) {x@startval<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.proctime'
standGen <- function(object) standardGeneric('.model.proctime')
standMethod <- function(object) object@proctime
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.proctime<-'
standGen <- function(x, value) standardGeneric('.model.proctime<-')
standMethod <- function(x, value) {x@proctime<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.valid'
standGen <- function(object) standardGeneric('.model.valid')
standMethod <- function(object) object@valid
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.valid<-'
standGen <- function(x, value) standardGeneric('.model.valid<-')
standMethod <- function(x, value) {x@valid<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.name'
standGen <- function(object) standardGeneric('.model.name')
standMethod <- function(object) object@name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.name<-'
standGen <- function(x, value) standardGeneric('.model.name<-')
standMethod <- function(x, value) {x@name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.fullpath'
standGen <- function(object) standardGeneric('.model.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.fullpath<-'
standGen <- function(x, value) standardGeneric('.model.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.betafiles'
standGen <- function(object) standardGeneric('.model.betafiles')
standMethod <- function(object) object@betafiles
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.betafiles<-'
standGen <- function(x, value) standardGeneric('.model.betafiles<-')
standMethod <- function(x, value) {x@betafiles<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.weightfiles'
standGen <- function(object) standardGeneric('.model.weightfiles')
standMethod <- function(object) object@weightfiles
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.weightfiles<-'
standGen <- function(x, value) standardGeneric('.model.weightfiles<-')
standMethod <- function(x, value) {x@weightfiles<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.avgdatfile'
standGen <- function(object) standardGeneric('.model.avgdatfile')
standMethod <- function(object) object@avgdatfile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.avgdatfile<-'
standGen <- function(x, value) standardGeneric('.model.avgdatfile<-')
standMethod <- function(x, value) {x@avgdatfile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.avgWfile'
standGen <- function(object) standardGeneric('.model.avgWfile')
standMethod <- function(object) object@avgWfile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.avgWfile<-'
standGen <- function(x, value) standardGeneric('.model.avgWfile<-')
standMethod <- function(x, value) {x@avgWfile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.avgtstatFile'
standGen <- function(object) standardGeneric('.model.avgtstatFile')
standMethod <- function(object) object@avgtstatFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.avgtstatFile<-'
standGen <- function(x, value) standardGeneric('.model.avgtstatFile<-')
standMethod <- function(x, value) {x@avgtstatFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.n'
standGen <- function(object) standardGeneric('.model.n')
standMethod <- function(object) object@n
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.n<-'
standGen <- function(x, value) standardGeneric('.model.n<-')
standMethod <- function(x, value) {x@n<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.mask'
standGen <- function(object) standardGeneric('.model.mask')
standMethod <- function(object) object@mask
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.mask<-'
standGen <- function(x, value) standardGeneric('.model.mask<-')
standMethod <- function(x, value) {x@mask<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.ss'
standGen <- function(object) standardGeneric('.model.ss')
standMethod <- function(object) object@ss
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.ss<-'
standGen <- function(x, value) standardGeneric('.model.ss<-')
standMethod <- function(x, value) {x@ss<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.regDir'
standGen <- function(object) standardGeneric('.model.regDir')
standMethod <- function(object) object@regDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.regDir<-'
standGen <- function(x, value) standardGeneric('.model.regDir<-')
standMethod <- function(x, value) {x@regDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.regRda'
standGen <- function(object) standardGeneric('.model.regRda')
standMethod <- function(object) object@regRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.regRda<-'
standGen <- function(x, value) standardGeneric('.model.regRda<-')
standMethod <- function(x, value) {x@regRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.funcDir'
standGen <- function(object) standardGeneric('.model.funcDir')
standMethod <- function(object) object@funcDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.funcDir<-'
standGen <- function(x, value) standardGeneric('.model.funcDir<-')
standMethod <- function(x, value) {x@funcDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.funcRda'
standGen <- function(object) standardGeneric('.model.funcRda')
standMethod <- function(object) object@funcRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.funcRda<-'
standGen <- function(x, value) standardGeneric('.model.funcRda<-')
standMethod <- function(x, value) {x@funcRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.runs'
standGen <- function(object) standardGeneric('.model.runs')
standMethod <- function(object) object@runs
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.runs<-'
standGen <- function(x, value) standardGeneric('.model.runs<-')
standMethod <- function(x, value) {x@runs<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.dataHeader'
standGen <- function(object) standardGeneric('.model.dataHeader')
standMethod <- function(object) object@dataHeader
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.dataHeader<-'
standGen <- function(x, value) standardGeneric('.model.dataHeader<-')
standMethod <- function(x, value) {x@dataHeader<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'model'
funcname <-'.model.version'
standGen <- function(object) standardGeneric('.model.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.model.version<-'
standGen <- function(x, value) standardGeneric('.model.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.name'
standGen <- function(object) standardGeneric('.data.name')
standMethod <- function(object) object@name
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.name<-'
standGen <- function(x, value) standardGeneric('.data.name<-')
standMethod <- function(x, value) {x@name<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.fullpath'
standGen <- function(object) standardGeneric('.data.fullpath')
standMethod <- function(object) object@fullpath
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.fullpath<-'
standGen <- function(x, value) standardGeneric('.data.fullpath<-')
standMethod <- function(x, value) {x@fullpath<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.betafiles'
standGen <- function(object) standardGeneric('.data.betafiles')
standMethod <- function(object) object@betafiles
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.betafiles<-'
standGen <- function(x, value) standardGeneric('.data.betafiles<-')
standMethod <- function(x, value) {x@betafiles<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.weightfiles'
standGen <- function(object) standardGeneric('.data.weightfiles')
standMethod <- function(object) object@weightfiles
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.weightfiles<-'
standGen <- function(x, value) standardGeneric('.data.weightfiles<-')
standMethod <- function(x, value) {x@weightfiles<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.avgdatfile'
standGen <- function(object) standardGeneric('.data.avgdatfile')
standMethod <- function(object) object@avgdatfile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.avgdatfile<-'
standGen <- function(x, value) standardGeneric('.data.avgdatfile<-')
standMethod <- function(x, value) {x@avgdatfile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.avgWfile'
standGen <- function(object) standardGeneric('.data.avgWfile')
standMethod <- function(object) object@avgWfile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.avgWfile<-'
standGen <- function(x, value) standardGeneric('.data.avgWfile<-')
standMethod <- function(x, value) {x@avgWfile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.avgtstatFile'
standGen <- function(object) standardGeneric('.data.avgtstatFile')
standMethod <- function(object) object@avgtstatFile
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.avgtstatFile<-'
standGen <- function(x, value) standardGeneric('.data.avgtstatFile<-')
standMethod <- function(x, value) {x@avgtstatFile<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.n'
standGen <- function(object) standardGeneric('.data.n')
standMethod <- function(object) object@n
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.n<-'
standGen <- function(x, value) standardGeneric('.data.n<-')
standMethod <- function(x, value) {x@n<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.mask'
standGen <- function(object) standardGeneric('.data.mask')
standMethod <- function(object) object@mask
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.mask<-'
standGen <- function(x, value) standardGeneric('.data.mask<-')
standMethod <- function(x, value) {x@mask<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.ss'
standGen <- function(object) standardGeneric('.data.ss')
standMethod <- function(object) object@ss
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.ss<-'
standGen <- function(x, value) standardGeneric('.data.ss<-')
standMethod <- function(x, value) {x@ss<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.regDir'
standGen <- function(object) standardGeneric('.data.regDir')
standMethod <- function(object) object@regDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.regDir<-'
standGen <- function(x, value) standardGeneric('.data.regDir<-')
standMethod <- function(x, value) {x@regDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.regRda'
standGen <- function(object) standardGeneric('.data.regRda')
standMethod <- function(object) object@regRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.regRda<-'
standGen <- function(x, value) standardGeneric('.data.regRda<-')
standMethod <- function(x, value) {x@regRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.funcDir'
standGen <- function(object) standardGeneric('.data.funcDir')
standMethod <- function(object) object@funcDir
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.funcDir<-'
standGen <- function(x, value) standardGeneric('.data.funcDir<-')
standMethod <- function(x, value) {x@funcDir<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.funcRda'
standGen <- function(object) standardGeneric('.data.funcRda')
standMethod <- function(object) object@funcRda
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.funcRda<-'
standGen <- function(x, value) standardGeneric('.data.funcRda<-')
standMethod <- function(x, value) {x@funcRda<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.runs'
standGen <- function(object) standardGeneric('.data.runs')
standMethod <- function(object) object@runs
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.runs<-'
standGen <- function(x, value) standardGeneric('.data.runs<-')
standMethod <- function(x, value) {x@runs<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.dataHeader'
standGen <- function(object) standardGeneric('.data.dataHeader')
standMethod <- function(object) object@dataHeader
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.dataHeader<-'
standGen <- function(x, value) standardGeneric('.data.dataHeader<-')
standMethod <- function(x, value) {x@dataHeader<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'data'
funcname <-'.data.version'
standGen <- function(object) standardGeneric('.data.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.data.version<-'
standGen <- function(x, value) standardGeneric('.data.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.current'
standGen <- function(object) standardGeneric('.sequence.current')
standMethod <- function(object) object@current
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.current<-'
standGen <- function(x, value) standardGeneric('.sequence.current<-')
standMethod <- function(x, value) {x@current<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.regions'
standGen <- function(object) standardGeneric('.sequence.regions')
standMethod <- function(object) object@regions
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.regions<-'
standGen <- function(x, value) standardGeneric('.sequence.regions<-')
standMethod <- function(x, value) {x@regions<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.mnames'
standGen <- function(object) standardGeneric('.sequence.mnames')
standMethod <- function(object) object@mnames
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.mnames<-'
standGen <- function(x, value) standardGeneric('.sequence.mnames<-')
standMethod <- function(x, value) {x@mnames<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.fit'
standGen <- function(object) standardGeneric('.sequence.fit')
standMethod <- function(object) object@fit
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.fit<-'
standGen <- function(x, value) standardGeneric('.sequence.fit<-')
standMethod <- function(x, value) {x@fit<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.minimum'
standGen <- function(object) standardGeneric('.sequence.minimum')
standMethod <- function(object) object@minimum
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.minimum<-'
standGen <- function(x, value) standardGeneric('.sequence.minimum<-')
standMethod <- function(x, value) {x@minimum<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.best'
standGen <- function(object) standardGeneric('.sequence.best')
standMethod <- function(object) object@best
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.best<-'
standGen <- function(x, value) standardGeneric('.sequence.best<-')
standMethod <- function(x, value) {x@best<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'sequence'
funcname <-'.sequence.valid'
standGen <- function(object) standardGeneric('.sequence.valid')
standMethod <- function(object) object@valid
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.sequence.valid<-'
standGen <- function(x, value) standardGeneric('.sequence.valid<-')
standMethod <- function(x, value) {x@valid<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'mnames'
funcname <-'.mnames.experiment'
standGen <- function(object) standardGeneric('.mnames.experiment')
standMethod <- function(object) object@experiment
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.mnames.experiment<-'
standGen <- function(x, value) standardGeneric('.mnames.experiment<-')
standMethod <- function(x, value) {x@experiment<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'mnames'
funcname <-'.mnames.subject'
standGen <- function(object) standardGeneric('.mnames.subject')
standMethod <- function(object) object@subject
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.mnames.subject<-'
standGen <- function(x, value) standardGeneric('.mnames.subject<-')
standMethod <- function(x, value) {x@subject<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'mnames'
funcname <-'.mnames.condition'
standGen <- function(object) standardGeneric('.mnames.condition')
standMethod <- function(object) object@condition
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.mnames.condition<-'
standGen <- function(x, value) standardGeneric('.mnames.condition<-')
standMethod <- function(x, value) {x@condition<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'mnames'
funcname <-'.mnames.mnames'
standGen <- function(object) standardGeneric('.mnames.mnames')
standMethod <- function(object) object@mnames
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.mnames.mnames<-'
standGen <- function(x, value) standardGeneric('.mnames.mnames<-')
standMethod <- function(x, value) {x@mnames<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'arfcorrelation'
funcname <-'.arfcorrelation.timebyreg'
standGen <- function(object) standardGeneric('.arfcorrelation.timebyreg')
standMethod <- function(object) object@timebyreg
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.arfcorrelation.timebyreg<-'
standGen <- function(x, value) standardGeneric('.arfcorrelation.timebyreg<-')
standMethod <- function(x, value) {x@timebyreg<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'arfcorrelation'
funcname <-'.arfcorrelation.corr'
standGen <- function(object) standardGeneric('.arfcorrelation.corr')
standMethod <- function(object) object@corr
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.arfcorrelation.corr<-'
standGen <- function(x, value) standardGeneric('.arfcorrelation.corr<-')
standMethod <- function(x, value) {x@corr<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'arfcorrelation'
funcname <-'.arfcorrelation.corr.pval'
standGen <- function(object) standardGeneric('.arfcorrelation.corr.pval')
standMethod <- function(object) object@corr.pval
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.arfcorrelation.corr.pval<-'
standGen <- function(x, value) standardGeneric('.arfcorrelation.corr.pval<-')
standMethod <- function(x, value) {x@corr.pval<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'arfcorrelation'
funcname <-'.arfcorrelation.pacorr'
standGen <- function(object) standardGeneric('.arfcorrelation.pacorr')
standMethod <- function(object) object@pacorr
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.arfcorrelation.pacorr<-'
standGen <- function(x, value) standardGeneric('.arfcorrelation.pacorr<-')
standMethod <- function(x, value) {x@pacorr<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'arfcorrelation'
funcname <-'.arfcorrelation.pacorr.pval'
standGen <- function(object) standardGeneric('.arfcorrelation.pacorr.pval')
standMethod <- function(object) object@pacorr.pval
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.arfcorrelation.pacorr.pval<-'
standGen <- function(x, value) standardGeneric('.arfcorrelation.pacorr.pval<-')
standMethod <- function(x, value) {x@pacorr.pval<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'arfcorrelation'
funcname <-'.arfcorrelation.num.corr'
standGen <- function(object) standardGeneric('.arfcorrelation.num.corr')
standMethod <- function(object) object@num.corr
setGeneric(funcname,standGen,package='arf3DS4')
setMethod(funcname,classname,standMethod)
slotreplace <-'.arfcorrelation.num.corr<-'
standGen <- function(x, value) standardGeneric('.arfcorrelation.num.corr<-')
standMethod <- function(x, value) {x@num.corr<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
