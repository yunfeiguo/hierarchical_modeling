fitBVS <- function(Z,data,forced=NULL,cov=NULL,a1=NULL,rare=FALSE,mult.regions=FALSE,regions=NULL,hap=FALSE,inform=FALSE,which=NULL,which.char=NULL){
	    case = names(data)[1]
	    #number of predictors
	    p = dim(data)[2] - 1
	    #number of observations
	    n = dim(data)[1]
	    #number of predictors in the model
	    pgamma = sum(Z!=0)
	    coef = rep(0,p)
	    #number of predicotr-level covariates
	    q = dim(cov)[2]
	    #number of observations of predictor-level covariates
	    r = dim(cov)[1]
	    if(length(cov)==0){
	    	q=1
	    	r=p
	    }
	    #BMU (assuming no effect of predictor-level covariate) model baseline
	    a0 = qnorm((1-2^(-1/r)))
	    #which.ind represents the number of coeffecients or independent variables in the glm model
	    which.ind = 1:p
	    if(rare==TRUE&mult.regions==FALSE)
	    {
	    	coef=0
	    	which.ind=1
	    	region.ind = as.matrix(rep(1,length(Z)))
	    }
	    if(rare==TRUE&mult.regions==TRUE)
	    {
	        coef=rep(0,length(unique(regions)))
		#when rare==TRUE and mult.regions==TRUE, number of independent variables in the model
		#is determined by number of unique regions
	        which.ind=1:length(unique(regions))
	        regions = as.factor(regions)
		#in '~regions-1', '-1' means removal of intercept term
	        region.ind = model.matrix(~regions-1)	
	        colnames(region.ind) = paste("G",c(1:length(levels(region.ind))),sep="")
	    }
	    #get names of covariates and predictors in the model
	    predictors = c(names(forced),colnames(data[,-1])[Z>0])	
	    data.i = data
	    
	    
	    	
	    ##Check if sampled new model already exists in which matrix
	    new.model.ind = c(1:length(which.char))[which.char==paste(Z,collapse="")]
	    
	    ##If sampled model already exists get the results from which
	    if(length(new.model.ind)>0)
	    {
	    	results.new = which[new.model.ind[1],]
	    	coef = as.numeric(results.new[which.ind])
	    	a1.old = as.numeric(results.new[-which.ind][1:q])
	    	fitness.old = as.numeric(results.new[-which.ind][q+1])
	    	logPrM.old = as.numeric(results.new[-which.ind][q+2]) 
	        ll = fitness.old + logPrM.old
	    }
	    
	    ##If sampled model does not already exist	
	    if(length(new.model.ind)==0)
	    {	
	      ##If we are doing an analysis with rare variants make risk index and new predict vector
	      if(rare == TRUE & pgamma>0)
	      {
		  #region.ind*as.numeric(z) will result in a new matrix with same dimension as region.ind
		  #this new matrix has rows filled with Zeros where the corresponding predictor is
		  #excluded from current model
		  #risk.index is a vector with same length as which.ind
		  #its row can be interpreted as the vector sum of alleles in each region for 
		  #each observation/subject
	    	  risk.index = as.matrix(data[,-1])%*%(region.ind*as.numeric(Z))
	    	  risk.index = as.matrix(risk.index)
	    	  
	    	  ##remove null indices 
	    	  which.ind = which.ind[apply(risk.index,2,sum)>0]
	    	  col.names =  paste("RI",colnames(region.ind)[apply(risk.index,2,sum)>0],sep=".")
	    	  risk.index = as.matrix(risk.index[,apply(risk.index,2,sum)>0])
	    	  colnames(risk.index) = col.names
	    	  
	    	  data.i = cbind(data,risk.index)
	    	  predictors=c(names(forced),colnames(risk.index))
	      }
	    	  
	      ##If we are doing a haplotype analysis estimate haplotypes and new predict vector
	      #do not consider haplotype for now--09/21/2014
	      if(hap==TRUE & pgamma>1)
	      {
	          G <- data[, colnames(data[,-1])[Z>0]]
			  X.haps = hapBVS(G,min.Hap.freq=.02)
			  data.i <- cbind(data, X.haps)
			  predictors = c(names(forced),names(X.haps))
	        }
	    
	      ##Attach forced variables to data matrix  	
	      if(length(forced)>0)
	      {
	          data.i = cbind(data.i,forced)
	      }   
	      
	        
	      ##Make fmla formula for glm    
	      if(length(predictors)==0)
	      {
	    	  fmla = as.formula(paste(case,"~1",sep=""))
	      }
	      if(length(predictors)>0)
	      {
	          fmla = as.formula(paste(case,"~",paste(predictors,collapse="+")))
	      }
	        
	      ##Fit glm
	      fit = suppressWarnings(glm(fmla,data = data.i,family="binomial"))
	      
	      
	      #!!!don't understand why we need penalty here
	      ##Calculate Log Likelihood
	      penalty = length(predictors) - length(forced)
	      ll = (deviance(fit) + penalty*2)/2
	    
	      ##Calculate coef. vector
	      if((rare==FALSE & hap==FALSE))
	      {
	          coef[Z!=0] = fit$coef[(length(fit$coef)-pgamma+1):length(fit$coef)]
	      }
	      if(rare==TRUE & pgamma>0)
	      {
	    	  coef[which.ind] = fit$coef[(2+length(forced)):length(fit$coef)]
	      }
	    }
	     
 
	    ##Calculate prior on the model
	    ##If inform ==TRUE use probit specification and if mult.regions==TRUE use region level probit specification
	    if(inform==TRUE)
	    {
	      eta = a0+as.matrix(cov)%*%as.matrix(a1)
          lprob.inc = pnorm(0,mean=eta,lower.tail=FALSE,log.p=TRUE)
	      lprob.ninc = pnorm(0,mean=eta,lower.tail=TRUE,log.p=TRUE)
	      logPrM = t(Z)%*%lprob.inc + t(1-Z)%*%lprob.ninc
	    }
	      
	      	   
	    ##BetaBinomial Prior
	    BetaBinomial = function(p,pgamma)
	    {
		#the following statement is actually trying to calculate factorial
		#gmma function is faster and less memory-demanding than direct factorial calculation
		#why p+p instead of p???
	    	logPrM = log(p) - lgamma(p+ p + 1) + lgamma(1 + pgamma) + lgamma(p + p - pgamma)
	    	return(logPrM)
	    }
	    	
	    if(inform==FALSE)
	    {
	      logPrM = BetaBinomial(p=p,pgamma=pgamma)
	      }
	                  
           
	    
       	fitness = ll - logPrM
       	results = c(coef,fitness,logPrM)
       	names(results) = c(paste("coef.",c(1:length(coef)),sep=""),"fitness","logPrM")
       	return(results)
    	}
