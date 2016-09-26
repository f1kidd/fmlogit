#' Extract fitted values, residuals, and residuals
#' 
#' @name fitted.fmlogit
#' @aliases residuals.fmlogit
#' @aliases predicted.fmlogit
#' Extract fitted dependent variable from a fractional multinomial logit model. 
#' @param object A "fmlogit" object.
#' @param newdata A new X matrix to perform model prediction. If Null, default to the original dataset. 
#' X can be a vector with length k, or a matrix with k columns, where k is the number of explanatory 
#' variables in the original model. 
#' @param newbeta A new augmented matrix of coefficients that can be used to predict outcome variables. 
#' Feeds into object$coefficient, which contains the baseline coefficient. Useful for constructing
#' confidence intervals via simulation or bootstrapping. 
#' @examples 
#' #results1 = fmlogit(y,X)
#' fitted(results1)
#' residuals(results1)
#' predict(results1)
#' # predict using the first observation from the original dataset.
#' predict(results1,X[1,])
#' @rdname fitted.fmlogit
#' @export fitted.fmlogit
#' 


fitted.fmlogit <-function(object){
  j=length(object$estimates)+1; k=dim(object$estimates[[1]])[1]; N=dim(object$y)[1]
  betamat_aug = object$coefficient; X=object$X; y=object$y
  sum_expxb = rowSums(exp(X %*% t(betamat_aug))) # sum of the exp(x'b)s
  yhat = y
  for(i in 1:j){
    expxb = exp(X %*% betamat_aug[i,]) # individual exp(x'b)
    yhat[,i] = expxb / sum_expxb
  }
  return(as.data.frame(yhat))
}

#' @rdname fitted.fmlogit
#' @export residuals.fmlogit
#' 
residuals.fmlogit <- function(object){
  yhat = fitted(object)
  return(as.data.frame(object$y-yhat))
}

#' @rdname fitted.fmlogit
#' @export predict.fmlogit
#' 
predict.fmlogit <- function(object,newdata=NULL,newbeta = NULL){
  if(length(newdata)==0) return(fitted(object))
  if(length(newbeta)>0) object$coefficient = newbeta
  j=length(object$estimates)+1; k=dim(object$estimates[[1]])[1]; N=dim(object$y)[1]
  betamat_aug = object$coefficient;
  newdata = as.matrix(newdata)
  if(length(newdata) == dim(newdata)[1]) newdata = t(newdata) # vector
  if(k != dim(newdata)[2]+1) stop(paste("Dimension of newdata is wrong. Should be",k-1,"instead of",dim(newdata)[2]))
  X = cbind(newdata,1); N = dim(X)[1]
  yhat = matrix(ncol=j,nrow=N); colnames(yhat) = colnames(object$y)
  sum_expxb = rowSums(exp(X %*% t(betamat_aug))) # sum of the exp(x'b)s
  for(i in 1:j){
    expxb = exp(X %*% betamat_aug[i,]) # individual exp(x'b)
    yhat[,i] = expxb / sum_expxb
  }
  return(as.data.frame(yhat))
}