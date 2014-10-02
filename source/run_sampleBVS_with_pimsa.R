args <- commandArgs(TRUE)
pimsa = args[1]
require('BVS')
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/fitBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/sampleBVS2.R")
source("/home/yunfeiguo/projects/Hierarchical_modeling/ASHG2014/source/summaryBVS.R")
system(paste(pimsa,"example_pimsa_settings.simulation1.xml"),wait=TRUE,ignore.stdout=TRUE)
data=read.table("bvs.datamatrix.simulation1.txt")
cov=as.matrix(read.table("bvs.cov.simulation1.txt"))
region=as.matrix(read.table("bvs.region.simulation1.txt"))
start_snp = read.table("models.out",as.is=TRUE,skip=PIMSA_ITERATION)
start_model = rep(0,dim(data)[2]-1)
if(start_snp[1,2] != '-1')
{
    start_snp = as.numeric(unlist(strsplit(start_snp[1,2],';')))
    #first index is -1 (for intercept), also the index is 0-based
    start_model[start_snp[-1]+1] = 1
}
result=sampleBVS2(start=start_model,cov=cov,inform=TRUE,mult.regions=TRUE,regions=region,iter=BVS_ITERATION,data=data,rare=TRUE,status.file="/dev/null")
summary.result = summaryBVS(result,data,rare=TRUE,mult.regions=TRUE,cov=cov,regions=region,inform=TRUE,burnin=BVS_BURNIN)
write(paste(summary.result$Global,paste(summary.result$Marg.RBF,collapse=" "),sep=" "),
      "bvs.result.simulation1.txt")
