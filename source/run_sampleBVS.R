require('BVS')
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/fitBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS2.R")
data=read.table("bvs.datamatrix.simulation1.txt")
cov=as.matrix(read.table("bvs.cov.simulation1.txt"))
region=as.matrix(read.table("bvs.region.simulation1.txt"))
#fit=glm(data[,1]~data[,9991]+data[,9992]+data[,9993]+data[,9994]+data[,9995],family="binomial")
result=sampleBVS(cov=cov,inform=TRUE,mult.regions=TRUE,regions=region,iter=bvs.iter,data=data,rare=TRUE,status.file="/dev/null")
summary.result = summaryBVS(result)
pdf("bvs.simulation1.modelsnp.pdf")
plotBVS(summary.result,num.models = 100, num.snps = 100, type = "s")
dev.off()
pdf("bvs.simulation1.modelregion.pdf")
plotBVS(summary.result,num.models = 100, num.regions = 10,type = "r")
dev.off()