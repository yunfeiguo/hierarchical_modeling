##########################PARAMETERS#############################
#based on the study results from 202 genes sequencing in 14,002 samples
#most variants are rare or private
n_replicate = 200 #number of simulation data sets
n_case = 500 #number of cases
n_ctrl = 500 #number of controls
#n_var = 10000
n_var = 10000 #total number of variants
prop_very_common = 0 #proportion of very common variants
prop_common = 0 #proportion of common variants
#proportion of rare variants in all variants
prop_rare = 1 #proportion of rare variants
maf_very_common = 0.1 #MAF of very common variants
maf_common = 0.05 #MAF of common variants
maf_rare = 0.005 #MAF of rare variants
n_annot = 10 #number of annotations (predictor-level covariates)
#annot_cover percentage of all variants will have annotation
annot_cover = 0.05
#annot_cover_risk percentage of risk variants covered by annotation
annot_cover_risk = 0.8
n_region = 10 #number of regions(genes/pathways)
n_risk_region = 1 #number of regions associated with disease
n_risk_very_common_var = 0 #number of very common risk variants
n_risk_common_var = 0 #number of common risk variants
n_risk_rare_var = 10 #number of rare risk variants
or_very_common_var = 1.05 #OR of very common risk variants
or_common_var = 1.2 #OR of common risk variants
or_rare_var = 10 #OR of rare risk variants
#assume subjects are independent, variants are independent (no LD)
#we only simulate rare variants, so ignore any code related to common variants
################################SUBROUTINES###################################
#expand genotype counts to actual genotypes (0,1,2 minor alleles)
expand.geno = function(x)
{
    return (
	    sapply(1:nrow(x),function(i)
		   {
		       result = rep(NA,sum(x[i,]))
		       for(j in 1:3)
		       {
			   which = which(is.na(result))
			   if(length(which) <= 1)
			   {
			       if(x[i,j] > 0)
				   result[which] = 3-j
			   } else
			   {
			       result[sample(which,x[i,j])] = 3-j
			   }
		       }
		       return(result)
		   }
	    )
	    )
}
#generate genotypes following Hardy-Weinberg equilibrium
#modified from HWData function in HardyWeinberg package
gen_hw_genotype  = function (nm = 100, n = 100, f = 0, p = runif(1), 
			     pfixed = FALSE, exactequilibrium = FALSE, pdist = "runif", ...) 
    #refer to manual of HardyWeinberg package for help
{
    Xt <- NULL
    if (length(p) == 1) 
	p <- rep(p, nm)
    stopifnot(length(n) == 1)
    if (length(n) == 1) 
	n <- rep(n, nm)
    if (length(f) == 1) 
	f <- rep(f, nm)
    ps <- p
    if (!pfixed) {
	for (i in 1:nm) {
	    if (is.null(p[i])) {
		if (pdist == "runif") 
		    ps[i] <- runif(1, ...)
		else {
		    if (pdist == "rbeta") 
			ps[i] <- rbeta(1, ...)
		    else stop("unknown value for pdist")
		}
	    }
	    else ps <- p
	    if (!exactequilibrium) 
		X <- t(rmultinom(1, size = n[i], prob = c(ps[i]^2 + 
							  f[i] * ps[i] * (1 - ps[i]), (1 - f[i]) * 2 * 
							  ps[i] * (1 - ps[i]), (1 - ps[i])^2 + f[i] * 
							  ps[i] * (1 - ps[i]))))
	    else X <- c(n[i] * ps[i]^2 + n[i] * f[i] * ps[i] * 
			(1 - ps[i]), n[i] * (1 - f[i]) * 2 * ps[i] * 
			(1 - ps[i]), n[i] * (1 - ps[i])^2 + n[i] * f[i] * 
			ps[i] * (1 - ps[i]))
	    Xt <- rbind(Xt, X)
	}
    }
    else {
	if (is.null(p)) 
	    p <- rep(0.5, nm)
	nAll <- 2 * n
	nA <- round(p * nAll, digits = 0)
	nB <- nAll - nA
	Xt <- NULL
	for (i in 1:nm) {
	    Pop <- c(rep(1, nA[i]), rep(0, nB[i]))
	    sam <- matrix(sample(Pop, nAll), ncol = 2, byrow = TRUE)
	    status <- apply(sam, 1, sum)
	    nAA <- sum(status == 2)
	    nAB <- sum(status == 1)
	    nBB <- sum(status == 0)
	    if ((nAA + nAB + nBB) != n[i]) 
		stop("HWData: error")
	    Xt <- rbind(Xt, c(nAA, nAB, nBB))
	}
    }
    #A is our allele of interest
    #B is the other allele
    colnames(Xt) <- c("AA", "AB", "BB")
    #expand Xt to the actual genotype (0,1,2 minor alleles)
    return(expand.geno(Xt))
}
#inverse of logit
invlogit = function (x)
{
    return(exp(x)/(1+exp(x)))
}
#generate actual genotypes for each subject at risk loci
#assume logistic regression model, binary outcome
#nm number of markers
#n number of subjects
#case: vector indicating case or not, recycled
#or: odds ratio vector, recycled
#maf: minor allele frequency, recycled
#beta0: constant in regression model
gen_riskvar_genotype = function (nm,n,case=FALSE,or,maf,beta0=0)
{
    #convert OR to logOR
    lor = rep(log(or),nm)[1:nm]
    beta0 = rep(beta0,nm)[1:nm]
    maf = rep(maf,nm)[1:nm]
    if(!case)
    {
	lor = -lor
	beta0 = -beta0
    }

    #each variant of each subject is independent from another
    return(
	   t(sapply(1:n,function(j)
		    {
			sapply(1:nm,function(i)
			       {
				   #P(variant|case) = P(case|variant)P(variant)/sum(P(case|variant)P(variant))
				   #sum(P(case|variant)P(variant))
				   p0 = (invlogit(sum(c(beta0[i],lor[i])*c(1,0)))*((1-maf[i])^2))
				   p1 = (invlogit(sum(c(beta0[i],lor[i])*c(1,1)))*(2*maf[i]*(1-maf[i])))
				   p2 = (invlogit(sum(c(beta0[i],lor[i])*c(1,2)))*(maf[i]^2))
				   p.sum = sum(p0,p1,p2)
				   variant_assign = t(rmultinom(1,1,prob=c(p0,p1,p2)/p.sum))
				   return( which(variant_assign == 1)-1)
			       }
			)
		    })
	   )
	   )
}

################################MAIN#########################################
for (i in 1:n_replicate)
{
    #we need prepare input for both pimsa and BVS
    #files for pimsa
    #0=missing,1=homozygous wildtype,2=heterozygous,3=homozygous mutant
    pimsafile.trait = paste("pimsa.trait.simulation",i,".txt",sep="")
    pimsafile.snplist = paste("pimsa.snplist.simulation",i,".txt",sep="")
    pimsafile.study_geno = paste("pimsa.study_geno.simulation",i,".txt",sep="")
    #pimsafile.study_env_var = paste("pimsa.study_env_var.simulation",i,".txt",sep="")
    pimsafile.initmodel = paste("pimsa.initmodel.simulation",i,".txt",sep="")
    pimsafile.z_matrix = paste("pimsa.z_matrix.simulation",i,".txt",sep="")
    pimsafile.a_matrix = paste("pimsa.a_matrix.simulation",i,".txt",sep="")
    #files for BVS
    #number of copies of minor alleles (0|1|2)
    bvsfile.datamatrix = paste("bvs.datamatrix.simulation",i,".txt",sep="")
    bvsfile.cov = paste("bvs.cov.simulation",i,".txt",sep="")
    bvsfile.region = paste("bvs.region.simulation",i,".txt",sep="")

    #prepare data
    #put all risk variants at the end
    n_risk_var = n_risk_rare_var + n_risk_common_var + n_risk_very_common_var
    idx.risk = (n_var-n_risk_var+1):(n_var)
    case.data = matrix(NA,nrow=n_case,ncol=n_var+1)
    ctrl.data = matrix(NA,nrow=n_ctrl,ncol=n_var+1)
    case.data[,1] = 1
    ctrl.data[,1] = 0
    #outcome and genotype (coded as 0,1,2 minor alleles)
    stopifnot(n_risk_rare_var<=n_var*prop_rare)
    hw_genotype = gen_hw_genotype(nm=n_var*prop_rare-n_risk_rare_var,n=n_case+n_ctrl,p=maf_rare)
    case.data[,-1] = cbind(hw_genotype[1:n_case,], 
			   gen_riskvar_genotype (nm=n_risk_rare_var,n=n_case,case=TRUE,or=or_rare_var,maf=maf_rare,beta0=0))
    ctrl.data[,-1] = cbind(hw_genotype[(n_case+1):(n_case+n_ctrl),],
			   gen_riskvar_genotype (nm=n_risk_rare_var,n=n_ctrl,case=FALSE,or=or_rare_var,maf=maf_rare,beta0=0))
    phenogeno.data = rbind(case.data,ctrl.data)
    rm(hw_genotype) #remove unnecessary memory
    #genotype and annotation
    #make sure to put all risk variants in annotation besides other variants in annotation
    #cov.data is pxq matrix, p for number of variants, q for number of annotations
    cov.data = matrix(0,nrow=n_var,ncol=n_annot)
    #assume sample will work (i.e. do not consider conditions that make sample fail)
    for (annot_i in 1:n_annot)
    {
	cov.data[sample(1:nrow(cov.data),nrow(cov.data)*annot_cover),annot_i] = 1
	cov.data[sample(idx.risk,n_risk_var*annot_cover_risk),annot_i] = 1
    }
    #genotype and region
    #make sure to put all risk variants into one region to maximize power
    #p-dimesional vector identifying region for p variants
    region.level = paste("pathway",1:n_region,sep="")
    region.data = rep(NA,n_var)
    region.data[-idx.risk] = region.level
    #put all risk variants into region 1
    region.data[idx.risk] = region.level[1]


    ###############output data for pimsa######################
    pimsadata.trait = phenogeno.data[,1]
    write.table(pimsadata.trait,pimsafile.trait,quote=FALSE,row.names=FALSE,col.names=FALSE)
    pimsadata.snplist = paste("snp",1:ncol(phenogeno.data[,-1]),sep="")
    write.table(pimsadata.snplist,pimsafile.snplist,quote=FALSE,row.names=FALSE,col.names=FALSE)
    #0=missing,1=homozygous wildtype,2=heterozygous,3=homozygous mutant
    pimsadata.study_geno = apply(phenogeno.data[,-1]+1,1,function(x){paste(x,collapse="")})
    write.table(pimsadata.study_geno,pimsafile.study_geno,quote=FALSE,row.names=FALSE,col.names=FALSE)
    ##set all environmental variables to 0
    #pimsadata.study_env_var = rep(0,nrow(phenogeno.data))
    #write.table(pimsadata.study_env_var,pimsafile.study_env_var,quote=FALSE,row.names=FALSE,col.names=FALSE)
    #just randomly choose 10 SNPs to include, SNP index is 0-based
    pimsadata.initmodel = sample(1:ncol(phenogeno.data),10)-1
    write.table(pimsadata.initmodel,pimsafile.initmodel,quote=FALSE,row.names=FALSE,col.names=FALSE)
    #z_matrix is fixed effect matrix, ie annotation
    pimsadata.z_matrix = cov.data
    write.table(pimsadata.z_matrix,pimsafile.z_matrix,quote=FALSE,row.names=FALSE,col.names=FALSE)
    #a_matrix is random effect matrix ie region
    pimsadata.a_matrix = model.matrix(~region.data-1)
    write.table(pimsadata.a_matrix,pimsafile.a_matrix,quote=FALSE,row.names=FALSE,col.names=FALSE)
    rm(pimsadata.trait,pimsadata.snplist,pimsadata.study_geno,pimsadata.initmodel,pimsadata.z_matrix,pimsadata.a_matrix)


    ####################output data for BVS#########################
    write.table(phenogeno.data,bvsfile.datamatrix ,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(cov.data,bvsfile.cov,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(region.data,bvsfile.region,quote=FALSE,row.names=FALSE,col.names=FALSE)
    cat("Simulation data set",i,"done.\n")
}
