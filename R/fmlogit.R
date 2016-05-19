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
#'   \code{optim(method=MLEmethod))}. Default to "CG", the conjugate gradients
#'   method.
#' @param maxit Maximum number of iteration.
#' @param abstol Tolerence.
#' @param ... additional parameters that goes into \code{optim()}
#' @return The function returns a list of elements below:
#' @return \code{$estimates} A list of matrices containing parameter estimates,
#'   standard errors, and hypothesis testing results.
#' @return \code{$likelihood} The likelihood value
#' @return \code{$convergence} Convergence diagnostics code
#' @return \code{$count} Provides dataset information
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
#' Maximization is done by calling \code{\link{optim}}. It seems that difference
#' choice of maximization method may yield drastically different parameter
#' estimates. Using the conjugate gradients method yields the same results as
#' using Stata's fmlogit package by Maarten Buis.
#' 
#' @examples 
#' require(foreign)
#' data = read.dta("http://fmwww.bc.edu/repec/bocode/c/citybudget.dta")
#' X = data[,2:5]
#' y = data[,6:11]
#' results1 = fmlogit(y,X)
#' @references Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods
#'   for fractional response variables with an application to 401(k) plan
#'   participation rates. J. Appl. Econ., 11: 619-632.
#' @export fmlogit




fmlogit <- function(y,X,beta0=NULL,MLEmethod="CG",maxit=50000,abstol=0.0001,...){

Xnames = colnames(X); ynames = colnames(y)
X = as.matrix(X); y = as.matrix(y)
n = dim(X)[1]; j=dim(y)[2]; k=dim(X)[2]

# handles missing value
xy = cbind(X,y)
xy = na.omit(xy)
X = xy[,1:k]; y=xy[,(k+1):(k+j)]
n = dim(X)[1]
remove(xy)

# add constant term if necessary
X = cbind(X,rep(1,n))
# here k is No. of explanatories, without the constant term. 
knew = qr(X)$rank # rank of X matrix
k=knew-1
if(knew<dim(X)[2]){ # X has constant column
  X = X[,1:knew]
} else{ # Augment constant column
  if(length(Xnames)==k) Xnames = c(Xnames,"constant")
}

QMLE<-function(betas,y,X){
  # betas here is as.vector(beta), which is (k+1)*(j-1)*1
  betas = matrix(betas,nrow=j-1,byrow=T)
  betamat = rbind(rep(0,k+1),betas) # augment with baseline betas of zero.
  
  llf = 0
  for(i in 1:j){
    L = y[,i] * ((X %*% betamat[i,] ) - log(rowSums(exp(X %*% t(betamat)))))
    llf = llf + sum(L)
  }
  return(-llf)
}

if(length(beta0)==0) beta0 = rep(0,(k+1)*(j-1))
if(length(beta0) != (k+1)*(j-1)) { # user gives a wrong beta0 vector
  beta0 = rep(0,(k+1)*(j-1)) # default to 0
  warning("Wrong length of beta0 given. Use default setting instead.")
}
# Here starts the MLE loop. 
# obviously, using different method of optimization can have big impact on parameter estimates. 
# Using "CG" (conjugate gradients) method yields same estimates as Stata does. 
opt<-optim(beta0,QMLE,gr=NULL,y,X,method=MLEmethod,
           control=list(maxit=maxit,abstol=abstol))

# reformat beta matrix
betamat = matrix(opt$par,ncol=k+1,byrow=T) #(j-1)*k
betamat_aug = rbind(rep(0,k+1),betamat) # augment beta matrix with baseline

# robust standard error, following PW(1996)
sigmat = matrix(nrow=j-1,ncol=k+1)
for(i in 2:j){
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
sigmat[i-1,] = std_b
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
outlist$likelihood = -opt$value
outlist$convergence = opt$convergence
outlist$count = c(Obs=n,Explanatories=k,Choices=j)

return(outlist)
}

