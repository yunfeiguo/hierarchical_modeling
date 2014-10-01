args <- commandArgs(TRUE)
file = args[1]
#format:
#N	TIME	GLOBAL_BF	BF1	BF2...
data = read.table(file,fill=TRUE,header=FALSE)
total = length(data[,2])
success =sum(data[,2]!="ERROR") 
fail =sum(data[,2]=="ERROR") 

cat("Total simulations:",total,"\n")
if(total != 0)
{
    cat("Successful simulations:",success,"(",success/total,")\n")
    cat("Failed simulations:",fail,"(",fail/total,")\n")
    if(dim(data)[2]>=4)
    {
	n_hiBF = sum(apply(data[data[,2]!="ERROR",],1,function(i){i=as.numeric(i);max(i[4:length(i)])==i[4]}))
	cat("Simulations with 1st region having highest BF:", n_hiBF,"(",n_hiBF/total,")\n")
	cat("Global BF:\n")
	print(summary(data[,3]))
	cat("BF1:\n")
	print(summary(data[,4]))
    }
}
