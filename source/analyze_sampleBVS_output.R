#CONSTANTS
    n_risk_region = 5
    n_region = 20
    n_risk_var = 20
    region_order = paste("pathway",1:n_region,sep="")
require('pROC')
##########################SUBROUTINES###################################
convert_region_output = function(x)
{
    return (
	    t(apply(x,1,function(i)
		    {
			i.split = unlist(strsplit(i,":"))
			stopifnot(length(i.split) %% 2 == 0)
			name = i.split[-seq(2,length(i.split),by=2)]
			val = i.split[seq(2,length(i.split),by=2)]
			return(val[match(name,region_order)])
		    })
	    )
	    )
}
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
make_boxplot = function(file,names,data,...)
{
    pdf(file)
    op <- par(no.readonly = TRUE)
    transform.data = data.frame(group=rep(args,(lapply(data,length))),outcome=as.numeric(unlist(data)))
    for(i in 1:length(args))
    {
	boxplot(outcome~group,transform.data,xlab='Group',ylab='Outcome')
    }
    par(op)
    dev.off()
}
########################MAIN#############################################
args <- commandArgs(TRUE)
globalBF=NULL
RBF=NULL
margBF=NULL
data.time=list()
if(length(args)>0)
{
    for (i in 1:length(args))
    {
	#format:
	#N	TIME	GLOBAL_BF	10	REGIONAL_BF	BF1	BF2...	MARGINAL_BF	BF1	BF2...
	file = args[i]
	data = read.table(file,fill=TRUE,header=FALSE,as.is=TRUE)
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

	#consider situation where we have unequal lengths
	data.time[[args[i]]] = data[success.idx,][,2]

	globalBF.idx = which(data[success.idx,][1,] == c('GLOBAL_BF'))
	RBF.idx = which(data[success.idx,][1,] == c('REGIONAL_BF'))
	margBF.idx = which(data[success.idx,][1,] == c('MARGINAL_BF'))
	stopifnot(margBF.idx>RBF.idx)
	globalPosterior.idx = which(data[success.idx,][1,] == c('GLOBAL_POSTERIOR'))
	regionalPosterior.idx = which(data[success.idx,][1,] == c('REGIONAL_POSTERIOR'))
	margPosterior.idx = which(data[success.idx,][1,] == c('MARGINAL_POSTERIOR'))
	stopifnot(margPosterior.idx>regionalPosterior.idx)

	globalBF.idx = globalBF.idx + 1
	RBF.idx = (RBF.idx+1):(margBF.idx-1)
	margBF.idx = (margBF.idx+1):(globalPosterior.idx-1)
	globalPosterior.idx = globalPosterior.idx + 1
	regionalPosterior.idx = (regionalPosterior.idx+1):(margPosterior.idx-1)
	margPosterior.idx = (margPosterior.idx+1):(dim(data)[2])

	#tranform RBF and regionalPosterior
	#old format: pathway1:0.1	pathway2:0.2
	#new format: 0.1	0.2
	data[success.idx,RBF.idx] = convert_region_output(data[success.idx,RBF.idx])
	data[success.idx,regionalPosterior.idx] = convert_region_output(data[success.idx,regionalPosterior.idx])
	##################MEDIAN GLOBAL BF##############################
	##choose the model with median globalBF for evaluating RBF and margBF
	#idx.median_globalBF_risk = 0
	#idx.median_globalBF_norisk = 0
	##by default, all NAs will appear at end in order()'s result
	##so we don't have to worry about them if we only choose the midpoint
	#idx.median_globalBF_risk = order(data[1:(dim(data)[1]/2),globalBF.idx])[floor(dim(data)[1]/4)]
	#idx.median_globalBF_norisk = dim(data)[1]/2+order(data[(dim(data)[1]/2+1):dim(data)[1],globalBF.idx])[floor(dim(data)[1]/4)]

	#globalBF = rbind(globalBF,as.numeric(data[,globalBF.idx]))
	##same situation as margBF
	#RBF = rbind(RBF,as.numeric(c(data[idx.median_globalBF_risk,RBF.idx],data[idx.median_globalBF_norisk,RBF.idx])))
	##data has been ordered, the first half contains no risk variants, the second half does
	##therefore we would have to construct a vector with length equal to twice
	##the number of variants, with first half from data sets containing risk variants
	##2nd half not
	#margBF = rbind(margBF,as.numeric(c(data[idx.median_globalBF_risk,margBF.idx],data[idx.median_globalBF_norisk,margBF.idx])))
	##################MEDIAN GLOBAL BF##############################

	##################MAX GLOBAL BF##############################
	#choose the model with max globalBF for evaluating RBF and margBF
	idx.max_globalBF_risk = 0
	idx.max_globalBF_norisk = 0
	#by default, all NAs will appear at end in order()'s result
	#so we don't have to worry about them if we only choose the max
	idx.max_globalBF_risk = order(data[1:(dim(data)[1]/2),globalBF.idx],decreasing=TRUE)[1]
	idx.max_globalBF_norisk = dim(data)[1]/2+order(data[(dim(data)[1]/2+1):dim(data)[1],globalBF.idx],decreasing=TRUE)[1]

	globalBF = rbind(globalBF,as.numeric(data[,globalBF.idx]))
	#same situation as margBF
	RBF = rbind(RBF,as.numeric(c(data[idx.max_globalBF_risk,RBF.idx],data[idx.max_globalBF_norisk,RBF.idx])))
	#data has been ordered, the first half contains no risk variants, the second half does
	#therefore we would have to construct a vector with length equal to twice
	#the number of variants, with first half from data sets containing risk variants
	#2nd half not
	margBF = rbind(margBF,as.numeric(c(data[idx.max_globalBF_risk,margBF.idx],data[idx.max_globalBF_norisk,margBF.idx])))
	##################MAX GLOBAL BF##############################

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
	    #the following is not really top models
	    #they are data sets that return highest globalBF
	    #for ( j in 1:10)
	    #{
	    #    topmodel.idx = order(data[,globalBF.idx],decreasing=TRUE)[j]
	    ##output top 10 regions in top 10 models
	    #    top10.region = order(data[topmodel.idx,RBF.idx], decreasing=TRUE)[1:10]
	    #    #cat("Top 10 regions in top",j,"model: ",paste(top10.region,format(data[topmodel.idx,top10.region],digits=1),sep=":"),"\n")
	    #    #output top top50 variants in top 10 models
	    #    top20.variant = order(data[topmodel.idx,margBF.idx], decreasing=TRUE)[1:20]
	    #    #cat("Top 50 variants in top",j,"model: ",paste(top50.variant,format(data[topmodel.idx,top50.variant],digits=1),sep=":"),"\n")
	    #    #cat("Model",j,": percentage of TRUE positives for regions: ",sum(top10.region %in% (1:5))/5,"\n")
	    #    #cat("Model",j,": percentage of TRUE positives for variants: ",sum(top50.variant %in% (9980:10000))/20,"\n")
	    #    cat(sum(top20.variant %in% (9980:10000))/20,"\n")
	    #}
	    
	    #output median posterior probabilities for top 20 variants
	    median_marg_posterior = apply(data[success.idx,margPosterior.idx],2,median)
	    top20.variant = order(median_marg_posterior,decreasing=TRUE)[1:20]
	    top20.variant.posterior = median_marg_posterior[top20.variant]
	    cat("Top 20 variants (median posterior):",paste(top20.variant,format(top20.variant.posterior,digits=1),sep=":"),"\n")
	        #cat("Model",j,": percentage of TRUE positives for regions: ",sum(top10.region %in% (1:5))/5,"\n")
	        #cat("Model",j,": percentage of TRUE positives for variants: ",sum(top50.variant %in% (9980:10000))/20,"\n")
	        #cat(sum(top20.variant %in% (9980:10000))/20,"\n")
	}
    }

    #plot ROC reports
    col = rainbow(length(args))
    #globalBF, assume first half is associated, second half is not
    stopifnot(dim(globalBF)[2] %% 2 == 0)
    make_roc(file="globalBF_ROC.pdf",names=args,data=globalBF,outcome=c(rep(1,dim(globalBF)[2]/2),rep(0,dim(globalBF)[2]/2)),col)
    #RBF(assume 1st region is always associated risk, no other region is associated with risk)
    make_roc(file="RBF_ROC.pdf",names=args,data=RBF,outcome=c(rep(1,n_risk_region),rep(0,dim(RBF)[2]/2-n_risk_region),rep(0,dim(RBF)[2]/2)),col)
    #marginal BF (each variant), assume only last 10 variants are risk variants
    stopifnot(dim(margBF)[2]>=10)
    make_roc(file="margBF_ROC.pdf",names=args,data=margBF,outcome=c(rep(0,dim(margBF)[2]/2-n_risk_var),rep(1,n_risk_var),rep(0,dim(margBF)[2]/2)),col)
    #draw time boxplot for comparison
    make_boxplot(file="boxplot.pdf",names=args,data=data.time)
}
