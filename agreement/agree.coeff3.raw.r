#								  AGREE.COEFF3.RAW.R
#						 		 (September 2, 2016)
#Description: This script file contains a series of R functions for computing various agreement coefficients
#	      for multiple raters (2 or more) when the input data file is in the form of nxr matrix or data frame showing 
#             the actual ratings each rater (column) assigned to each subject (in row). That is n = number of subjects, and r = number of raters.
#             A typical table entry (i,g) represents the rating associated with subject i and rater g. 
#Author: Kilem L. Gwet, Ph.D. (gwet@agreestat.com)
#-----------------------------------------------------------------
#     EXAMPLES OF SIMPLE CALLS OF THE MAIN FUNCTIONS:
# > gwet.ac1.raw(YourRatings)       # to obtain gwet's AC1 coefficient
# > fleiss.kappa.raw(YourRatings)   # to obtain fleiss' unweighted generalized kappa coefficient
# > krippen.alpha.raw(YourRatings)  # to obtain krippendorff's unweighted alpha coefficient
# > conger.kappa.raw(YourRatings)   # to obtain conger's unweighted generalized kappa coefficient
# > bp.coeff.raw(YourRatings)       # to obtain Brennan-Prediger unweighted coefficient
#

#===========================================================================================
#gwet.ac1.raw: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in the presence of high
#		agreement." British Journal of Mathematical and Statistical Psychology, 61, 29-48.
#============================================================================================
gwet.ac1.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings)
  if (is.character(ratings.mat)){ratings.mat <- toupper(ratings.mat)}
  
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(na.omit(ratings.mat)))
  if (is.numeric(categ.init))
    categ <- sort(categ.init)
  else {
    ratings.mat<-trim(ratings.mat)
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- sort(categ.init[nchar(categ.init)>0])
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
        k.mis <-(ratings.mat==categ[k])
        in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
        agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  gwet.ac1 <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-gwet.ac1) * (pe.ivec-pe)/(1-pe)
  
  var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - gwet.ac1)^2)
  stderr <- sqrt(var.ac1)# ac1's standard error
  p.value <- 2*(1-pt(abs(gwet.ac1/stderr),n-1))
  
  lcb <- gwet.ac1 - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,gwet.ac1 + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE) {
    if (weights=="unweighted") {
       cat("Gwet's AC1 Coefficient\n")
       cat('======================\n')	
       cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
       cat('AC1 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
       cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
       cat('P-value: ',p.value,'\n')
    }
    else {
       cat("Gwet's AC2 Coefficient\n")
       cat('==========================\n')	
       cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
	 cat('AC2 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
       cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
       cat('P-value: ',p.value,'\n')
	 cat('\n')
       if (!is.numeric(weights)) {
	    cat('Weights: ', weights,'\n')
          cat('---------------------------\n') 	
	 }
       else{
	    cat('Weights: Custom Weights\n')
          cat('---------------------------\n')
	 }
	 print(weights.mat)
    }
  }
  invisible(c(pa,pe,gwet.ac1,stderr,p.value))
}

#=====================================================================================
#fleiss.kappa.raw: This function computes Fleiss' generalized kappa coefficient (see Fleiss(1971)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John Wiley & Sons.
#======================================================================================
fleiss.kappa.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(ratings.mat)}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(na.omit(ratings.mat)))
  if (is.numeric(categ.init))
    categ <- sort(categ.init)
  else {
    ratings.mat<-trim(ratings.mat)
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- sort(categ.init[nchar(categ.init)>0])
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating fleiss's generalized kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  fleiss.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  kappa.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.vec
  kappa.ivec.x <- kappa.ivec - 2*(1-fleiss.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.fleiss <- ((1-f)/(n*(n-1))) * sum((kappa.ivec.x - fleiss.kappa)^2)
  stderr <- sqrt(var.fleiss)# kappa's standard error
  p.value <- 2*(1-pt(abs(fleiss.kappa/stderr),n-1))
  
  lcb <- fleiss.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,fleiss.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE){
    cat("Fleiss' Kappa Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Fleiss kappa coefficient:',fleiss.kappa,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
    if (weights!="unweighted") {
	 cat('\n')
       if (!is.numeric(weights)) {
	    cat('Weights: ', weights,'\n')
          cat('---------------------------\n')
	 }
       else{
	    cat('Weights: Custom Weights\n')
          cat('---------------------------\n')
	 }
	 print(weights.mat)	
    }
  }
  invisible(c(pa,pe,fleiss.kappa,stderr,p.value))
}


#=====================================================================================
#krippen.alpha.raw: This function computes Krippendorff's alpha coefficient (see Krippendorff(1970, 1980)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The algorithm used to compute krippendorff's alpha is very different from anything that was published on this topic. Instead,
#it follows the equations presented by K. Gwet (2012)
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. (2012). Handbook of Inter-Rater Reliability: The Definitive Guide to Measuring the Extent of Agreement Among 
#	Multiple Raters, 3rd Edition.  Advanced Analytics, LLC; 3rd edition (March 2, 2012)
#Krippendorff (1970). "Bivariate agreement coefficients for reliability of data." Sociological Methodology,2,139-150
#Krippendorff (1980). Content analysis: An introduction to its methodology (2nd ed.), New-bury Park, CA: Sage.
#======================================================================================
krippen.alpha.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(ratings.mat)}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(na.omit(ratings.mat)))
  if (is.numeric(categ.init))
    categ <- sort(categ.init)
  else {
    ratings.mat<-trim(ratings.mat)
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- sort(categ.init[nchar(categ.init)>0])
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating krippendorff's alpha coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  agree.mat<-agree.mat[(ri.vec>=2),]
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(agree.mat)
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
print(n)
paprime <-sum(sum.q/(ri.mean*(ri.vec-1)))/n
print(paprime)

  pa <- (1-epsi)*sum(sum.q/(ri.mean*(ri.vec-1)))/n + epsi

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/ri.mean))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  krippen.alpha <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.mean*(ri.vec-1)
  pa.ivec <- sum.q/den.ivec
  pa.v <- mean(pa.ivec)
  pa.ivec <- pa.ivec-pa.v*(ri.vec-ri.mean)/ri.mean

  krippen.ivec <- (pa.ivec-pe)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec

  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.mean - sum(pi.vec) * (ri.vec-ri.mean)/ri.mean
  krippen.ivec.x <- krippen.ivec - 2*(1-krippen.alpha) * (pe.ivec-pe)/(1-pe)
  
  var.krippen <- ((1-f)/(n*(n-1))) * sum((krippen.ivec.x - krippen.alpha)^2)
  stderr <- sqrt(var.krippen)# alpha's standard error
  p.value <- 2*(1-pt(abs(krippen.alpha/stderr),n-1))
  
  lcb <- krippen.alpha - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,krippen.alpha + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE){
    cat("Krippendorff's Alpha Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Krippendorff alpha coefficient:',krippen.alpha,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
    if (weights!="unweighted") {
	 cat('\n')
       if (!is.numeric(weights)) {
	    cat('Weights: ', weights,'\n')
          cat('---------------------------\n')
	 }
       else{
	    cat('Weights: Custom Weights\n')
          cat('---------------------------\n')
	 }
	print(weights.mat)
    }
  }
  invisible(c(pa,pe,krippen.alpha,stderr,p.value))
}


#===========================================================================================
#conger.kappa.raw: Conger's kappa coefficient (see Conger(1980)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Conger, A. J. (1980), ``Integration and Generalization of Kappas for Multiple Raters,"
#		Psychological Bulletin, 88, 322-328.
#======================================================================================
conger.kappa.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(ratings.mat)}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(na.omit(ratings.mat)))
  if (is.numeric(categ.init))
	  categ <- sort(categ.init)
  else {
  	ratings.mat<-trim(ratings.mat)
  	categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
   	categ <- sort(categ.init[nchar(categ.init)>0])
  }

  q <- length(categ)
  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
  	k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	agree.mat[,k] <- in.categ.k%*%rep(1,r) 
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
	
  # creating the rxq rater-category matrix representing the distribution of subjects by rater and category

  classif.mat <- matrix(0,nrow=r,ncol=q)
  for(k in 1:q){
	with.mis <-(t(ratings.mat)==categ[k])
	without.mis <- replace(with.mis,is.na(with.mis),FALSE)
	classif.mat[,k] <- without.mis%*%rep(1,n)
  }
  # calculating conger's kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  ng.vec <- classif.mat%*%rep(1,q)
  pgk.mat <- classif.mat/(ng.vec%*%rep(1,q))
  p.mean.k <- (t(pgk.mat)%*%rep(1,r))/r 
  s2kl.mat <- (t(pgk.mat)%*%pgk.mat - r * p.mean.k%*%t(p.mean.k))/(r-1)
  pe <- sum(weights.mat * (p.mean.k%*%t(p.mean.k) -  s2kl.mat/r))
  conger.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of conger's kappa coefficient
  
  bkl.mat <- (weights.mat+t(weights.mat))/2
  pe.ivec1 <- r*(agree.mat%*%t(t(p.mean.k)%*%bkl.mat))
  pe.ivec2 = rep(0,n)
  
  lamda.ig.mat=matrix(0,n,r)
  if (is.numeric(ratings.mat)){
	epsi.ig.mat <-1-is.na(ratings.mat)
	epsi.ig.mat <- replace(epsi.ig.mat,is.na(epsi.ig.mat),FALSE)
  }else{
	epsi.ig.mat <- 1-(ratings.mat=="")
	epsi.ig.mat <- replace(epsi.ig.mat,is.na(epsi.ig.mat),FALSE)
  }
  for(k in 1:q){
	lamda.ig.kmat=matrix(0,n,r)
	for(l in 1:q){
		delta.ig.mat <- (ratings.mat==categ[l])
		delta.ig.mat <- replace(delta.ig.mat,is.na(delta.ig.mat),FALSE)
		lamda.ig.kmat <- lamda.ig.kmat + weights.mat[k,l] * (delta.ig.mat - (epsi.ig.mat - rep(1,n)%*%t(ng.vec/n)) * (rep(1,n)%*%t(pgk.mat[,l])))
	}
	lamda.ig.kmat = lamda.ig.kmat*(rep(1,n)%*%t(n/ng.vec))
	lamda.ig.mat = lamda.ig.mat+ lamda.ig.kmat*(r*mean(pgk.mat[,k]) - rep(1,n)%*%t(pgk.mat[,k]))
  }
  pe.ivec <- (lamda.ig.mat%*%rep(1,r)) / (r*(r-1))
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  pe.r2 <- pe*(ri.vec>=2)
  conger.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe) 
  conger.ivec.x <- conger.ivec - 2*(1-conger.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.conger <- ((1-f)/(n*(n-1))) * sum((conger.ivec.x - conger.kappa)^2)
  stderr <- sqrt(var.conger)# conger's kappa standard error
  p.value <- 2*(1-pt(abs(conger.kappa/stderr),n-1))
  
  lcb <- conger.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,conger.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE) {
    cat("Conger's Kappa Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement: ',pa,'Percent chance agreement: ',pe,'\n')
    cat("Conger's kappa coefficient: ",conger.kappa,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
    if (weights!="unweighted") {
	 cat('\n')
       if (!is.numeric(weights)) {
          cat('Weights: ', weights,'\n')
          cat('---------------------------\n')
       }
       else{
	    cat('Weights: Custom Weights\n')
          cat('---------------------------\n')
       }
	 print(weights.mat)
    }
  }
  invisible(c(pa,pe,conger.kappa,stderr,p.value))
}

#===========================================================================================
#bp.coeff.raw: Brennan-Prediger coefficient (see Brennan & Prediger(1981)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, misuses, and alternatives."
#           Educational and Psychological Measurement, 41, 687-699.
#======================================================================================
bp.coeff.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(ratings.mat)}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(na.omit(ratings.mat)))
  if (is.numeric(categ.init))
    categ <- sort(categ.init)
  else {
    ratings.mat<-trim(ratings.mat)
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- sort(categ.init[nchar(categ.init)>0])
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) / (q^2)
  bp.coeff <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  bp.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  var.bp <- ((1-f)/(n*(n-1))) * sum((bp.ivec - bp.coeff)^2)
  stderr <- sqrt(var.bp)# BP's standard error
  p.value <- 2*(1-pt(abs(bp.coeff/stderr),n-1))
  
  lcb <- bp.coeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,bp.coeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE) {
    cat("Brennan-Prediger Coefficient\n")
    cat('============================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('B-P coefficient:',bp.coeff,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
    if (weights!="unweighted") {
       if (!is.numeric(weights)) {
	    cat('\n')
	    cat('Weights: ', weights,'\n')
	              cat('---------------------------\n')	
	    print(weights.mat)
	 }
       else{
	    cat('Weights: Custom Weights\n')
          cat('---------------------------\n')
	 }
	 print(weights.mat)
    }
  }
  invisible(c(pa,pe,bp.coeff,stderr,p.value))
}

#
#-----  Additional functions needed to run the main functions. If the main functions must be included in another R script, then
# the user will need to add these additional functions to the new script file.
#
# ==============================================================
# trim(x): This is an r function for trimming leading and trealing blanks
# ==============================================================
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) 
}
# ==============================================================
# The following functions generate various weight matrices used 
# in the weighted or unweighted analyses.
# ==============================================================
identity.weights<-function(categ){
  weights<-diag(length(categ))
  return (weights)
}

quadratic.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) { 
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-(categ.vec[k]-categ.vec[l])^2/(xmax-xmin)^2 
    }
  }
  return (weights)
}

linear.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) { 
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-abs(categ.vec[k]-categ.vec[l])/abs(xmax-xmin)
    }
  }
  return (weights)
}
#--------------------------------
radical.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) { 
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-sqrt(abs(categ.vec[k]-categ.vec[l]))/sqrt(abs(xmax-xmin))
    }
  }
  return (weights)
}

#--------------------------------
ratio.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) { 
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-((categ.vec[k]-categ.vec[l])/(categ.vec[k]+categ.vec[l]))^2 / ((xmax-xmin)/(xmax+xmin))^2
    }
  }
  return (weights)
}

#--------------------------------
circular.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) { 
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  U = xmax-xmin+1
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- (sin(pi*(categ.vec[k]-categ.vec[l])/U))^2
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}

#--------------------------------
bipolar.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) { 
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      if (k!=l)
        weights[k,l] <- (categ.vec[k]-categ.vec[l])^2 / (((categ.vec[k]+categ.vec[l])-2*xmin)*(2*xmax-(categ.vec[k]+categ.vec[l])))
      else weights[k,l] <- 0
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}
#--------------------------------
ordinal.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  categ.vec<-1:length(categ)
  for(k in 1:q){
    for(l in 1:q){
      nkl <- max(k,l)-min(k,l)+1
      weights[k,l] <- nkl * (nkl-1)/2
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}
