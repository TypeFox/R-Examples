##################################################################################
# Set Class: rasclassRaster
##################################################################################
setClass('rasclassRaster',

representation(
# Header
ncols     = 'numeric',
nrows     = 'numeric',
xllcorner = 'numeric',
yllcorner = 'numeric',
cellsize  = 'numeric',
NAvalue   = 'numeric',

# Data
grid = 'numeric'
),

prototype = list(ncols = NULL, nrows = NULL, xllcorner = NULL, yllcorner = NULL, cellsize = NULL, NAvalue = NULL, grid = NULL)
)


##################################################################################
# Set Class: rasclass
##################################################################################
setOldClass('multinom')
setOldClass('mlp')
setOldClass('svm.formula')
setOldClass('randomForest.formula')

setClass('rasclass',

representation(
# Input Variables
path = 'character',
data = 'data.frame',
samplename = 'character',
formula = 'character',
call = 'call',

# Memory saving skeleton
gridSkeleton = 'rasclassRaster',

# Training vector for splitfraction
training = 'logical',

# Reconstruction type specifics
maximumLikelihood = 'list',
randomForest = 'randomForest.formula',
logit = 'multinom',
neuralNetwork = 'mlp',
supportVector = 'svm.formula',

# Predicted Variables
predictedGrid = 'rasclassRaster',

# Accuracy and reconstruction statistics
overallAccuracy = 'numeric',
accuracyMatrix = 'matrix',
kappa = 'numeric'),

prototype = list(
path = NULL,
data = NULL,
samplename = NULL,
formula = NULL,
call = NULL,
gridSkeleton = new('rasclassRaster'),
training = NULL,
maximumLikelihood = NULL,
randomForest = NULL,
logit = NULL,
neuralNetwork = NULL,
supportVector = NULL,
predictedGrid = new('rasclassRaster'),
overallAccuracy = NULL,
accuracyMatrix = NULL,
kappa = NULL))
