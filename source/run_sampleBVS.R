require('BVS')
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/fitBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS2.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/summaryBVS.R")
data=read.table("bvs.datamatrix.simulation1.txt")
cov=as.matrix(read.table("bvs.cov.simulation1.txt"))
region=as.matrix(read.table("bvs.region.simulation1.txt"))
#fit=glm(data[,1]~data[,9991]+data[,9992]+data[,9993]+data[,9994]+data[,9995],family="binomial")
result=sampleBVS(cov=cov,inform=TRUE,mult.regions=TRUE,regions=region,iter=BVS_ITERATION,data=data,rare=TRUE,status.file="/dev/null")
summary.result = summaryBVS(result,data,rare=TRUE,mult.regions=TRUE,cov=cov,regions=region,inform=TRUE,burnin=BVS_BURNIN)
#sapply(c("Global","Marg.RBF"),function(i){
#       write(paste(i,summary.result[[which(names(summary.result)==i)]]),
#	     "bvs.result.simulation1.txt",append=TRUE)
#})
write(paste(summary.result$Global,paste(summary.result$Marg.RBF,collapse=" "),sep=" "),
      "bvs.result.simulation1.txt")
