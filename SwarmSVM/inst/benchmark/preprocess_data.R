source('utils.R')
require(Matrix)

########################
# For Clustered SVM
########################

# svmguide1
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/binary/svmguide1',
              'svmguide1')
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/binary/svmguide1.t',
              'svmguide1.t')
svmguide1 = read.libsvm('svmguide1')
svmguide1.t = read.libsvm('svmguide1.t')
svmguide1[,-1] = rowl2norm(svmguide1[,-1])
svmguide1.t[,-1] = rowl2norm(svmguide1.t[,-1])
svmguide1 = list(svmguide1,svmguide1.t)
save(svmguide1, file='svmguide1.RData',compress = 'xz')

# ijcnn1
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/binary/ijcnn1.bz2',
              'ijcnn1.bz2')
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/binary/ijcnn1.t.bz2',
              'ijcnn1.t.bz2')
ijcnn1 = read.libsvm(bzfile('ijcnn1.bz2'))
ijcnn1.t = read.libsvm(bzfile('ijcnn1.t.bz2'))
ijcnn1[,-1] = rowl2norm(ijcnn1[,-1])
ijcnn1.t[,-1] = rowl2norm(ijcnn1.t[,-1])
ijcnn1 = list(ijcnn1,ijcnn1.t)
save(ijcnn1,file='ijcnn1.RData',compress = 'xz')

# usps
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/multiclass/usps.bz2',
              'usps.bz2')
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/multiclass/usps.t.bz2',
              'usps.t.bz2')
usps = read.libsvm(bzfile('usps.bz2'))
usps.t = read.libsvm(bzfile('usps.t.bz2'))
usps[,-1] = rowl2norm(usps[,-1])
usps.t[,-1] = rowl2norm(usps.t[,-1])
usps[,1] = as.numeric(usps[,1]%%2==0)
usps.t[,1] = as.numeric(usps.t[,1]%%2==0)
usps = list(usps,usps.t)
save(usps,file='usps.RData',compress = 'xz')

# mnist
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/multiclass/mnist.bz2',
              'mnist.bz2')
download.file('http://www.csie.ntu.edu.tw/%7Ecjlin/libsvmtools/datasets/multiclass/mnist.t.bz2',
              'mnist.t.bz2')
mnist = read.libsvm(bzfile('mnist.bz2'), dims=c(60000,782))
mnist.t = read.libsvm(bzfile('mnist.t.bz2'),dim=c(10000,782))
mnist[,-1] = rowl2norm(mnist[,-1])
mnist.t[,-1] = rowl2norm(mnist.t[,-1])

mnist38 = mnist[which(mnist[,1]==3 | mnist[,1]==8),]
mnist38.t = mnist.t[which(mnist.t[,1]==3 | mnist.t[,1]==8),]
mnist38[,1] = as.numeric(mnist38[,1]==3)
mnist38.t[,1] = as.numeric(mnist38.t[,1]==3)

mnist49 = mnist[which(mnist[,1]==4 | mnist[,1]==9),]
mnist49.t = mnist.t[which(mnist.t[,1]==4 | mnist.t[,1]==9),]
mnist49[,1] = as.numeric(mnist49[,1]==4)
mnist49.t[,1] = as.numeric(mnist49.t[,1]==4)

mnistoe = mnist
mnistoe[,1] = as.numeric(as.vector(mnist[,1])%%2==0)
mnistoe.t = mnist.t
mnistoe.t[,1] = as.numeric(as.vector(mnist.t[,1])%%2==0)

mnist = list(mnist38,mnist38.t,mnist49,mnist49.t,mnistoe,mnistoe.t)
save(mnist, file='mnist.RData',compress = 'xz')


########################
# For DC SVM
########################

# ijcnn1
download.file('http://www.sfu.ca/~hetongh/data/ijcnn1.train',
              'ijcnn1.train')
download.file('http://www.sfu.ca/~hetongh/data/ijcnn1.t',
              'ijcnn1.t')
ijcnn1 = read.libsvm('ijcnn1.train')
ijcnn1.t = read.libsvm('ijcnn1.t')
# for (i in 2:ncol(ijcnn1)) {
#   if (sd(ijcnn1[,i])>0)
#     ijcnn1[,i] = (ijcnn1[,i]-min(ijcnn1[,i]))/(max(ijcnn1[,i])-min(ijcnn1[,i]))
# }
# for (i in 2:ncol(ijcnn1.t)) {
#   if (sd(ijcnn1.t[,i])>0)
#     ijcnn1.t[,i] = (ijcnn1.t[,i]-min(ijcnn1.t[,i]))/(max(ijcnn1.t[,i])-min(ijcnn1.t[,i]))
# }
ijcnn1 = list(ijcnn1,ijcnn1.t)
save(ijcnn1,file='ijcnn1.dcsvm.RData',compress = 'xz')

# covtype.binary
download.file('http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/covtype.libsvm.binary.scale.bz2',
              'covtype.libsvm.binary.scale.bz2')
covtype = read.libsvm(bzfile('covtype.libsvm.binary.scale.bz2'))
set.seed(1024)
n = nrow(covtype)
ind = sample(n,n*0.2)
covtype.t = covtype[ind,]
covtype = covtype[-ind,]
covtype = list(covtype,covtype.t)
save(covtype,file = 'covtype.RData',compress = 'xz')

########################
# For Gater SVM
########################

# Toy example 
set.seed(1024)
train.1 = cbind(runif(333,-1.7,-0.7),
                runif(333,0.7,1.7))
train.2 = cbind(runif(333,-0.5,0.5),
                runif(333,-0.5,0.5))
train.3 = cbind(runif(334,0.7,1.7),
                runif(334,-1.7,-0.7))
y = c(rep(1,333),rep(0,333),rep(1,334))
train = rbind(train.1, train.2, train.3)
train = cbind(train, y)

test.1 = cbind(runif(3333,-1.7,-0.7),
               runif(3333,0.7,1.7))
test.2 = cbind(runif(3333,-0.5,0.5),
               runif(3333,-0.5,0.5))
test.3 = cbind(runif(3334,0.7,1.7),
               runif(3334,-1.7,-0.7))
y = c(rep(1,3333), rep(0,3333), rep(1,3334))
test = rbind(test.1, test.2, test.3)
test = cbind(test, y)

toydata = list(train, test)
save(toydata, file='toydata.RData', compress='xz')

# covtype
download.file('http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass/covtype.bz2',
              'covtype.bz2')
covtype = read.libsvm(bzfile('covtype.bz2'))
set.seed(1024)
n = nrow(covtype)
covtype[,1] = as.numeric(covtype[,1]==2)
covtype = covtype[,-2]
max.train = apply(covtype,2,max)
for (i in which(max.train!=1)) {
  covtype[,i] = covtype[,i]/max.train[i]
  cat(i,'\r')
}
save(covtype,file = 'covtype.mult.RData',compress = 'xz')




