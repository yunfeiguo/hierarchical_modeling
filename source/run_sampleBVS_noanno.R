require('BVS')
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/fitBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/summaryBVS.R")
data=read.table("bvs.datamatrix.simulation1.txt")
cov=as.matrix(read.table("bvs.cov.simulation1.txt"))
region=as.matrix(read.table("bvs.region.simulation1.txt"))
result=sampleBVS(mult.regions=TRUE,regions=region,iter=BVS_ITERATION,data=data,rare=TRUE,status.file="/dev/null")
summary.result = summaryBVS(result,data,rare=TRUE,mult.regions=TRUE,regions=region,burnin=BVS_BURNIN)
write(summary.result$Global,
      "bvs_noanno.sampleBVS.result.simulation1.txt")
write(paste(summary.result$Marg.RBF,collapse=" "),
      append=TRUE,"bvs_noanno.sampleBVS.result.simulation1.txt")
write(paste(summary.result$MargBF,collapse=" "),
      append=TRUE,"bvs_noanno.sampleBVS.result.simulation1.txt")
