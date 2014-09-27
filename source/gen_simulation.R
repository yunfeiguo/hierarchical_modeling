##########################PARAMETERS#############################
#based on the study results from 202 genes sequencing in 14,002 samples
#most variants are rare or private
n_replicate = 1000
n_case = 500
n_ctrl = 500
n_var = 10000
prop_very_common = 0
prop_common = 0
prop_rare = 0
maf_very_common = 0.1
maf_common = 0.05
maf_rare = 0.005
n_annot = 2
#annot_cover percentage of all variants will have annotation
annot_cover = 0.05
n_region = 10
n_risk_region = 1
n_risk_very_common_var = 0
n_risk_common_var = 0
n_risk_rare_var = 10
lod_very_common_var = 1.05
lod_common_var = 1.2
lod_rare_var = 20
#assume subjects are independent, variants are independent (no LD)
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
			   cat(which,"\n")
			   if(length(which) <= 1)
			   {
			       if(x[i,j] > 0)
				   result[which] = 3-j
			   } else
			   {
			       cat(sample(which,x[i,j]),"\n")
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
gen_hw_genotype  = function (nm = 100, n = rep(100, nm), f = rep(0, nm), p = runif(nm), 
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

################################MAIN#########################################
for (i in 1:n_replicate)
{
    #we need prepare input for both pimsa and BVS
    #files for pimsa
    #0=missing,1=homozygous wildtype,2=heterozygous,3=homozygous mutant
    pimsafile.trait = paste("pimsa.trait.simulation",i,".txt",sep="")
    pimsafile.snplist = paste("pimsa.snplist.simulation",i,".txt",sep="")
    pimsafile.study_geno = paste("pimsa.study_geno.simulation",i,".txt",sep="")
    pimsafile.study_env_var = paste("pimsa.study_env_var.simulation",i,".txt",sep="")
    pimsafile.initmodel = paste("pimsa.initmodel.simulation",i,".txt",sep="")
    pimsafile.z_matrix = paste("pimsa.z_matrix.simulation",i,".txt",sep="")
    pimsafile.a_matrix = paste("pimsa.a_matrix.simulation",i,".txt",sep="")
    #files for BVS
    #number of copies of minor alleles (0|1|2)
    bvsfile.datamatrix = paste("bvs.datamatrix.simulation",i,".txt",sep="")
    bvsfile.cov = paste("bvs.cov.simulation",i,".txt",sep="")
    bvsfile.region = paste("bvs.region.simulation",i,".txt",sep="")

    #prepare data
    case.data = matrix(NA,nrow=n_case,ncol=n_var+1)
    ctrl.data = matrix(NA,nrow=n_ctrl,ncol=n_var+1)
    case.data[,1] = 1
    ctrl.data[,1] = 0
    #outcome and genotype (coded as 0,1,2 minor alleles)
    stopifnot(n_risk_rare_var>n_var*prop_rare)
    n_risk_rare_var
    lod_rare_var
    rpois(n_var*prop_rare-n_risk_rare_var,maf_rare)
    #genotype and annotation
    #genotype and region

    data = rbind(case.data,ctrl.data)
    #output data for pimsa
    #output data for BVS
}
