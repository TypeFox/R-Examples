scalevalue<-function(key.scale, z){
# scalevalue - computes scale of detection function for a given 
#              set of scale covariates and parameters, using log link
# Args:
#  key.scale      scale parameters
#  z              design matrix for scale covariates
    exp(z %*% key.scale)
}
