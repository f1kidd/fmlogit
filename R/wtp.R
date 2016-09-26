#' "Willingness to Pay" for fmlogit models
#' 
#' Calculates the willingness to pay for fractional multinomial logit models.  
#' 
#' @param object An "fmlogit.margins" object.
#' @param wtp.vec A 1*J vector that contains the willingness to pay for each choice j. 
#' @param varlist A string vector which provides the name of variables to calculate 
#' the wtp. If missing, all variables in object will be calculated. 
#' @return A matrix containing the estimates, standard error, z-stats, and p-value. 
#' @details This function calculates the aggregate effect of a variable on the 
#' "willingness to pay" by linearly multiplying coefficients with ex-ante (arbitary) 
#' willingness to pay numbers associated with each choice. 
#' 
#' Suppose there are three choices A,B,C, each with a willingness to pay (or cost, profit, budget),
#' of 100, 200, and 300. The discrete effect of variable X on A,B and C are 0.5, 0.5, and -1, with 
#' standard error 0.2, 0.3 and 0.5. The aggregated discrete effect of X on the total willingness 
#' to pay (or cost), is thus 100*0.5 + 200*0.5 + 300*(-1) = -150. And the standard error can be also
#' calculated to be 162.8, assuming that the standard error is independent. 
#' A simple z-test is provided to test whether the aggregate effect is different from zero. 
#' 
#' Note that if the input fmlogit.margins object has no standard error computation, then no standard error
#' @examples
#' #results1 = fmlogit(y,X)
#' #effects1 = effects(results1,effect="marginal",se=T)
#' # assume that the WTP = 1,2,3,...J for each choice j. 
#' wtp(effects1,seq(1:nrow(effects1$effects))) 
#' 
#'

wtp = function(object,wtp.vec,varlist=NULL){
  j=nrow(object$effects); K=ncol(object$effects)
  Xnames = colnames(object$effects); ynames = rownames(object$effects)
  if(length(varlist)==0){
    varlist=Xnames
    var_colNo = c(1:K)
    k = length(var_colNo)
  }else{
    var_colNo = which(varlist %in% Xnames)
    k = length(var_colNo)
  }
  if(length(wtp.vec)!=j) stop("Wrong length of wtp.vec. Please check specification again.")
  # wtp calcs
  betamat = object$effects[,varlist]; semat = object$se[,varlist]
  wtp_mean = wtp.vec %*% betamat
  if(length(object$R)>0){
  wtp_se = sqrt(wtp.vec^2 %*% semat^2)
  # output tables
  tabout = matrix(ncol=4,nrow=k)
  tabout[,1] = wtp_mean
  tabout[,2] = wtp_se
  tabout[,3] = tabout[,1] / tabout[,2]
  tabout[,4] = 2*(1-pnorm(abs(tabout[,3])))
  colnames(tabout) = c("estimate","std","z","p-value")
  rownames(tabout) = varlist
  return(tabout)
  }else return(wtp_mean)
}
  