require('pROC')
##########################SUBROUTINES###################################
make_roc <- function(file="roc.pdf",names,data,outcome,...)
{
    pdf(file)
    auc.val = NULL
    for(i in 1:length(args))
    {
	roc.obj = roc(outcome,data[i,],na.rm=TRUE)
	if(i == 1)
	{
	    plot(roc.obj,col=col[i])
	} else
	{
	    lines.roc(roc.obj,col=col[i])
	}
	auc.val = c(auc.val,roc.obj$auc)
    }
    legend("bottomright", legend=paste(args,"AUC:",format(auc.val,digits=2)),
	   col=col, lwd=2)
    dev.off()
}
########################MAIN#############################################
args <- commandArgs(TRUE)
globalBF=NULL
RBF=NULL
margBF=NULL
if(length(args)>0)
{
    for (i in 1:length(args))
    {
	#format:
	#N	TIME	GLOBAL_BF	10	REGIONAL_BF	BF1	BF2...	MARGINAL_BF	BF1	BF2...
	file = args[i]
	data = read.table(file,fill=TRUE,header=FALSE)
	#make sure we got equal risk-containing and risk-free data sets
	#a lot code depends on this spec!!!!
        stopifnot(dim(data)[1] %% 2 == 0)
	#rearrange based on simulation number
	data = data[order(data[,1]),]
	success.idx = data[,2]!="ERROR"
	total = length(data[,2])
	success =sum(success.idx) 
	fail =sum(!success.idx) 
	#For lines with ERROR, 3rd column and beyond are filled with blanks,
	#we can convert the blanks into NAs to suppress warnings in as.numeric
	data[!success.idx,3:dim(data)[2]]=NA

	globalBF.idx = which(data[success.idx,][1,] == c('GLOBAL_BF'))
	RBF.idx = which(data[success.idx,][1,] == c('REGIONAL_BF'))
	margBF.idx = which(data[success.idx,][1,] == c('MARGINAL_BF'))
	stopifnot(margBF.idx>RBF.idx)

	globalBF.idx = globalBF.idx + 1
	RBF.idx = (RBF.idx+1):(margBF.idx-1)
	margBF.idx = (margBF.idx+1):(dim(data)[2])
	#choose the model with median globalBF for evaluating RBF and margBF
	idx.median_globalBF_risk = 0
	idx.median_globalBF_norisk = 0
	#by default, all NAs will appear at end in order()'s result
	#so we don't have to worry about them if we only choose the midpoint
	idx.median_globalBF_risk = order(data[1:(dim(data)[1]/2),globalBF.idx])[floor(dim(data)[1]/4)]
	idx.median_globalBF_norisk = dim(data)[1]/2+order(data[(dim(data)[1]/2+1):dim(data)[1],globalBF.idx])[floor(dim(data)[1]/4)]

	globalBF = rbind(globalBF,as.numeric(data[,globalBF.idx]))
	#same situation as margBF
	RBF = rbind(RBF,as.numeric(c(data[idx.median_globalBF_risk,RBF.idx],data[idx.median_globalBF_norisk,RBF.idx])))
	#data has been ordered, the first half contains no risk variants, the second half does
	#therefore we would have to construct a vector with length equal to twice
	#the number of variants, with first half from data sets containing risk variants
	#2nd half not
	margBF = rbind(margBF,as.numeric(c(data[idx.median_globalBF_risk,margBF.idx],data[idx.median_globalBF_norisk,margBF.idx])))


	cat("File:",file,"\n")
	cat("Total simulations:",total,"\n")
	if(total != 0)
	{
	    cat("Successful simulations:",success,"(",100*success/total,"%)\n")
	    cat("Failed simulations:",fail,"(",100*fail/total,"%)\n")
	    if(dim(data)[2]>=6)
	    {
		n_hiBF = sum(apply(data[success.idx,RBF.idx],1,function(i){i=as.numeric(i);max(i)==i[1]}))
		cat("Successful simulations with 1st region having highest BF:", n_hiBF,"(",100*n_hiBF/success,"%)\n")
		cat("Global BF:\n")
		print(summary(as.numeric(data[success.idx,globalBF.idx])))
		cat("BF1:\n")
		print(summary(as.numeric(data[success.idx,RBF.idx][,1])))
	    }
	}
    }

    #plot ROC reports
    col = rainbow(length(args))
    n_risk_region = 5
    n_risk_var = 20
    #globalBF, assume first half is associated, second half is not
    stopifnot(dim(globalBF)[2] %% 2 == 0)
    make_roc(file="globalBF_ROC.pdf",names=args,data=globalBF,outcome=c(rep(1,dim(globalBF)[2]/2),rep(0,dim(globalBF)[2]/2)),col)
    #RBF(assume 1st region is always associated risk, no other region is associated with risk)
    make_roc(file="RBF_ROC.pdf",names=args,data=RBF,outcome=c(rep(1,n_risk_region),rep(0,dim(RBF)[2]/2-n_risk_region),rep(0,dim(RBF)[2]/2)),col)
    #marginal BF (each variant), assume only last 10 variants are risk variants
    stopifnot(dim(margBF)[2]>=10)
    make_roc(file="margBF_ROC.pdf",names=args,data=margBF,outcome=c(rep(0,dim(margBF)[2]/2-n_risk_var),rep(1,n_risk_var),rep(0,dim(margBF)[2]/2)),col)
}
