#' Plot marginal or discrete effects of willingness to pay
#' 
#' Plot marginal or discrete effects of willingness to pay, potentially against another variable
#' 
#' @param object An "fmlogit" object.
#' @param varlist A string vector which provides the name of variables to plot the effect.
#'  If missing, all variables in object will be plotted.
#' @param X The covariates matrix. Recommend to use element X from the fmlogit object. 
#' @param y The covariates matrix. Recommend to use element y from the fmlogit object. 
#' @param against A vector with the same length as the number of observations in the model. 
#' Serve as the x-axis in the plots.
#' @param against.x A character string, Supply the column name in the X matrix to be plot against.
#' @param against.y A character string, Supply the column name in the y matrix to be plot against.
#' @param group.x A character string. Supply the column name in the X matrix to be grouped upon. 
#' @param group.by A character string. Supply additional algebra emposed on the group variable. 
#' @param mfrow A numeric vector with two elements. Specify the number of rows and columns in a panel.
#' Similar to par(mfrow=c()). Default to Null, and the program will choose a square panel. 
#' @param plot.show If true, the plot will be created. Otherwise the function returns raw data that can be
#' used to create user-specified (fancier) plots. 
#' @return Panel plots of effects vs. chosen variables
#' @details 
#' This function provides a visualization tool for potentially heterogeneous marginal and discrete effects.
#' The function lets the user to plot marginal effects to detect any patterns in the effects, in itself
#' and against other variables. The plot also allows visualization of sub-groups in data, which can be
#' very useful to visualize categorical and dummy variables. 
#' 
#' The functions takes an fmlogit.margins object, created by the effects(fmlogit) function. Note that since 
#' the plotting requires marginal effects for all observations, the object should be created by choosing 
#' \code{marg.type="aveacr"}, the average across method for effects calculation. 
#' 
#' Additional parameters including \code{varlist}, a vector of string variable names to be plotted. \code{X}
#'  and \code{y}, the dependent and independent variable matrix in the original regression model. 
#'  
#'  \code{against}, \code{against.x}, and \code{against.y} allows different variables to be chosen
#'  as the x-axis. \code{against} directly supplies the vector to be plotted against, whereas \code{against.x}
#'  and \code{against.y} supplies variable names in the original dataset. Note that the user has to provide
#'  \code{X} and \code{y} in order to use the column name option, respectively. 
#'  
#'  \code{group.x} supplies the column name in the X matrix to be grouped by. The plot will be able to 
#'  differentiate different groups by colors. Additionally, the user can supply a string to \code{group.by},
#'  which provides a algebra method that will be evaluated on the group vector. For example, choose 
#'  \code{group.x = "a"} and \code{group.by= ">0"} will create two groups, one with X$a>0, and one with X$a
#'  <=0
#' @examples  
#' # Not running
#' # results1 = fmlogit(y,X)
#' # effect1 = effects(results1,effect="marginal",marg.type="aveacr")
#' 
#' # Plot only takes effects with marg.type="aveacr". 
#' plot(effect1,X=results1$X,against.x = "popdens", group = "tot", groupby = ">3")
#' @export plot.fmlogit


plot.fmlogit = function(object,wtp.vec,varlist, against=NULL,mfrow=NULL,t=500,effect=c("discrete","marginal"),
                        type="l",plot.show=T,...){
  K = ncol(object$X); j = ncol(object$y); N = nrow(object$X); 
  Xnames = colnames(object$X) ; ynames = colnames(object$y)
  X = object$X; y=object$y
  
  # determine variable list
  var_colNo = which(Xnames %in% varlist)
  k = length(var_colNo)
  
  if(is.null(mfrow)){
    js = ceiling(sqrt(k))
    jr = ifelse(js*(js-1)>=k,js-1,js)
  }else{
    jr = mfrow[1]; js = mfrow[2]
  }
  
  if(!is.null(against)) {
    ag_No = which(Xnames == against)
    if(length(ag_No)==0) stop(paste("The against vector specified,",against,
                                    "is not in the list of explanatory variables. Please check again."))
    ag_min = min(X[,ag_No]); ag_max = max(X[,ag_No])
    ag_vec = seq(ag_min,ag_max,length.out = t)
    wtp_mat = matrix(nrow=t,ncol=k)
    colnames(wtp_mat) = varlist
    for(i in 1:t){
      newdata = colMeans(X[,-K])
      newdata[ag_No] = ag_vec[i]
      wtp_mat[i,] = wtp(effects(object,effect=effect,se=F,varlist=varlist,at=newdata),wtp.vec)[[1]]
    }
  }else{
    against="ObsNo"
    ag_vec=1:N
    wtp_mat = matrix(nrow=N,ncol=k)
    colnames(wtp_mat) = varlist
    for(i in 1:N){
      newdata = X[i,-K]
      wtp_mat[i,] = wtp(effects(object,effect=effect,se=F,varlist=varlist,at=newdata),wtp.vec)[[1]]
    }
  }
  # plotting
  if(plot.show){
    par(mfrow=c(jr,js))
    if(is.null(type)){type="l"} # default to line plot. 
    for(i in 1:k){
      plot(ag_vec,wtp_mat[,i],xlab=against,ylab=paste(effect,"effect of", varlist[i]),...)
    }}
  return(list(ag_vec,wtp_mat))
}
