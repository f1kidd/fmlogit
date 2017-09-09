#' Estimate Fractional Multinomial Logit Models
#' 
#' Used to estimate fractional multinomial logit models using quasi-maximum
#' likelihood estimations following Papke and Wooldridge(1996).
#' @param y the dependent variable (N*J). Can be a matrix or a named data frame.
#'   The first column of the matrix is automatically treated as the baseline.
#' @param X independent variable (N*K). Can be a matrix or a named data frame.
#'   If there is no intercept term in the X, then an intercept term is
#'   automatically added.
#' @param beta0 Initial value for beta used in optimization. Uses a 1*K(J-1)
#'   vector. Default to a vector of zeros.
#' @param MLEmethod Method of optimization. Goes into
#'   \code{maxLik(method=MLEmethod))}. Choose from "NR","BFGS","CG","BHHH","SANN",or "NM".  
#'   Default to "CG", the conjugate gradients method. See Details. 
#' @param maxit Maximum number of iteration.
#' @param abstol Tolerence.
#' @param ... additional parameters that goes into \code{optim()}
#' @return The function returns an object of class "fmlogit". Use \code{effects}, \code{predict}, 
#'  \code{residual}, \code{fitted} to extract various useful features of the value returned by 
#' \code{fmlogit}. 
#' @return An object of class "fmlogit" contains the following components: 
#' @return \code{estimates}   A list of matrices containing parameter estimates,
#'   standard errors, and hypothesis testing results.
#' @return \code{baseline}    The baseline choice
#' @return \code{likelihood}  The likelihood value
#' @return \code{conv_code}   Convergence diagnostics code. 
#' @return \code{convergence} Convergence messages. 
#' @return \code{count}       Provides dataset information
#' @return \code{y}           The dependent variable data frame.
#' @return \code{X}           The independent variable data frame. Augmented by factor dummy transformation
#' , constant term added. 
#' @return \code{coefficient} Matrix of estimated coefficients. Augmented with the baseline coefficient
#' (which is a vector of zeros). 
#' @return \code{vcov}        A list of matrices containing the robust variance covariance matrix for each choice
#' variable. 
#' @details The fractional multinomial model is the expansion of the multinomial
#' logit to fractional responses. Unlike standard multinomial logit models,
#' which only considers 0-1 respones, fractional multinomial model considers the
#' case where the response variable is fractions that sums up to one. Examples
#' of these type of data are, percentages of budget spent in education, defense,
#' public health; fractions of a population that have middle school, high
#' school, college, or post college education, etc.
#' 
#' This function follows Papke and Wooldridge(1996)'s paper, in which they
#' proposed a quasi-maximum likelihood estimator for fractional response data.
#' The likelihood function used here is a standard multinomial likelihood
#' function, see \url{http://maartenbuis.nl/software/likelihoodFmlogit.pdf} for
#' the likelihood used here. Robust standard errors are provided following Papke
#' and Wooldridge(1996), in which they proposed an asymptotically consistent
#' estimator of variance.
#' 
#' Maximization is done by calling \code{\link{maxLik}}. maxLik is a wrapper function 
#' for different maximization methods in R. This include most methods provided by \code{\link{maxLik}},  
#' but also other methods such as BHHH(Berndt-Hall-Hall-Hausman). 
#' 
#' MLE convergence can be a problem in R, especially if dataset is large with many explanatory variables. 
#' It is recommended to call CG(Conjugate Gradients) or BHHH(Berndt-Hall-Hall-Hausman).
#' Conjugate gradients method is usually faster, but could lead to non-convergence under 
#' certain scenarios. BHHH is slower, but has better convergence performance.
#' 
#' 
#' @examples 
#' data = spending
#' X = data[,2:5]
#' y = data[,6:11]
#' results1 = fmlogit(y,X)
#' @references Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods
#'   for fractional response variables with an application to 401(k) plan
#'   participation rates. J. Appl. Econ., 11: 619-632.
#' @export fmlogit




fmlogit <- function(y,X,beta0=NULL,MLEmethod="CG",maxit=500000,abstol=0.00001,...){

start.time=proc.time()  

# y has to be numerical. X can be numerical or factor/categorical
Xclass = sapply(X,class)

# Get factor and character columns
Xfac = which(Xclass %in% c("factor","character"))

# Check for categorical variables and generate dummies
if( length(Xfac > 0)){
  Xfacnames = colnames(X)[Xfac]
  strformFac = paste(Xfacnames,collapse="+")
  # This creates a formula ~dummies, which goes into model.matrix to generate dummies
  Xdum = model.matrix(as.formula(paste("~",strformFac,sep="")),data=X)
  X = cbind(X,Xdum); X = X[,-Xfac] 
}

Xnames = colnames(X); ynames = colnames(y)
X = as.matrix(X); y = as.matrix(y)
n = dim(X)[1]; j=dim(y)[2]; k=dim(X)[2]

# handles missing value
xy = cbind(X,y)
xy = na.omit(xy)
X = xy[,1:k]; y=xy[,(k+1):(k+j)]
n = dim(X)[1]
remove(xy)

# remove pre-existing constant variables
if( length(Xfac > 0 )) X = X[,sapply(X,function(x) length(unique(x))!=1)]

# add constant term if necessary
X = cbind(X,rep(1,n))
# here k is No. of explanatories, without the constant term. 
knew = qr(X)$rank # rank of X matrix
k=knew-1
if(knew<dim(X)[2]){ # X has constant column
  X = X[,1:knew]
} else{ # Augment constant column
  if(length(Xnames)==k) Xnames = c(Xnames,"constant"); colnames(X) = Xnames
}

QMLE<-function(betas){
  # This is the overall(sum) log likelihood. 
  # betas here is as.vector(beta), which is (k+1)*(j-1)*1
  betas = matrix(betas,nrow=j-1,byrow=T)
  betamat = rbind(rep(0,k+1),betas) # augment with baseline betas of zero.
  llf = 0
  for(i in 1:j){
    L = y[,i] * ((X %*% betamat[i,] ) - log(rowSums(exp(X %*% t(betamat)))))
    llf = llf + sum(L)
  }
  return(llf)
}

QMLE_Obs<-function(betas){
  # This is the individual specific log likelihood. Required if using BHHH to maximize. 
  # betas here is as.vector(beta), which is (k+1)*(j-1)*1
  betas = matrix(betas,nrow=j-1,byrow=T)
  betamat = rbind(rep(0,k+1),betas) # augment with baseline betas of zero.
  
  llf = rep(0,n)
  for(i in 1:j){
    L = y[,i] * ((X %*% betamat[i,] ) - log(rowSums(exp(X %*% t(betamat)))))
    llf = llf + L
  }
  return(llf)
}

if(length(beta0)==0) beta0 = rep(0,(k+1)*(j-1))
if(length(beta0) != (k+1)*(j-1)) { # user gives a wrong beta0 vector
  beta0 = rep(0,(k+1)*(j-1)) # default to 0
  warning("Wrong length of beta0 given. Use default setting instead.")
}
# Here starts the MLE loop. 
# obviously, using different method of optimization can have big impact on parameter estimates. 
# Using "CG" (conjugate gradients) method yields same estimates as Stata does. 
opt<-maxLik(QMLE_Obs,start=beta0,method=MLEmethod, control=list(iterlim=maxit,tol=abstol),...)

# reformat beta matrix
betamat = matrix(opt$estimate,ncol=k+1,byrow=T) #(j-1)*k
betamat_aug = rbind(rep(0,k+1),betamat) # augment beta matrix with baseline
colnames(betamat_aug) = Xnames; rownames(betamat_aug)=ynames

# robust standard error, following PW(1996)
sigmat = matrix(nrow=j-1,ncol=k+1)

vcov = list()
for(i in 1:j){
# start calculation  
sum_expxb = rowSums(exp(X %*% t(betamat_aug))) # sum of the exp(x'b)s
expxb = exp(X %*% betamat_aug[i,]) # individual exp(x'b)
G = expxb / sum_expxb # exp(X'bj) / sum^J(exp(X'bj))
g = (expxb * sum_expxb - expxb^2) / sum_expxb^2 # derivative of the logit function

# Here the diagonal of A is the 'standard' standard error
# hat(A) = sum hat(gi)^2 * xi'xi / hat(Gi)(1-hat(Gi))
# or, Xtilde = X * sqrt(g^2/G(1-G)), A = Xtilde'Xtilde
X_a = X * as.vector(sqrt(g^2/(G*(1-G))))
A = t(X_a) %*% X_a

# robust standard error, again following PW(1996)
mu = y[,i] - G
X_b = X * as.vector(mu * g / G / (1-G))
B = t(X_b) %*% X_b
Var_b = solve(A) %*% B %*% solve(A)
std_b = sqrt(diag(Var_b))
# std_b= sqrt(diag(solve(A))) is the "unrobust" standard error. 
vcov[[i]] = Var_b
if(i>1) sigmat[i-1,] = std_b
}

# generating hypothesis testing tables.
listmat = list()
for(i in 1:(j-1)){
  tabout = matrix(ncol=4,nrow=k+1)
  tabout[,1:2] = t(rbind(betamat[i,],sigmat[i,]))
  tabout[,3] = tabout[,1] / tabout[,2]
  tabout[,4] = 2*(1-pnorm(abs(tabout[,3])))
  colnames(tabout) = c("estimate","std","z","p-value")
  if(length(Xnames)>0) rownames(tabout) = Xnames
  listmat[[i]] = tabout
}
if(length(ynames)>0) names(listmat) = ynames[2:j]

outlist = list()
outlist$estimates = listmat
outlist$baseline = ynames[1]
outlist$likelihood = opt$maximum
outlist$conv_code = opt$code
outlist$convergence = paste(opt$type, paste(as.character(opt$iterations),"iterations"), opt$message,sep=",")
outlist$count = c(Obs=n,Explanatories=k,Choices=j)
outlist$y = y
outlist$X = X
outlist$coefficient = betamat_aug
names(vcov) = ynames
outlist$vcov = vcov

print(paste("Fractional logit model estimation completed. Time:",round(proc.time()[3]-start.time[3],1),"seconds"))

return(structure(outlist,class="fmlogit"))
}





