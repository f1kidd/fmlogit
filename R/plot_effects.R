#' "Willingness to Pay" for fmlogit models
#' 
#' Calculates the willingness to pay for fractional multinomial logit models.  
#' 
#' @param object An "fmlogit.margins" object.
#' @param varlist A string vector which provides the name of variables to calculate 
#' the wtp. If missing, all variables in object will be calculated. 
#' @param against A character string of an index to plot against. Default to Null, and uses
#' "natural" order of the dataset.
#' @param X The covariates matrix. Recommend to use element X from the fmlogit object. 
#' @param y The covariates matrix. Recommend to use element y from the fmlogit object. 
#' @param mfrow A two-element numeric vector indicating plot panel row and column. 
#' Feeds into par(mfrow). Please make sure that the size of the panel is larger than the
#' number of choices. If NULL, use default method to generate square panels. 
#' @return Plots of effects vs. chosen variables
#' @details 
#' 



plot.fmlogit.margins = function(object,varlist=NULL,X=NULL,y=NULL, 
                                against=NULL,against.x=NULL,against.y=NULL,
                                group=NULL, group.by=NULL,
                                mfrow=NULL){
  if(is.null(object[["marg.list"]])) stop("Please choose marg.type=aveacr when calculating effects")
  k = ncol(object$effects); j = nrow(object$effects) ;N = nrow(object$marg.list[[1]]); 
  Xnames = colnames(object$effects) ; ynames = rownames(object$effects)
  
  # determine variable list
  if(length(varlist)==0){
    varlist=Xnames
    var_colNo = 1:k
  }else{
    var_colNo = which(Xnames %in% varlist)
    k = length(var_colNo)
  }
  if(k==0) stop("Variable list not matched. Please check your varlist input.")
  
  # determine panel size
  if(is.null(mfrow)){
    js = ceiling(sqrt(j))
    jr = ifelse(js*(js-1)>j,js-1,js)
    mfrow=c(jr,js)
  }
  
  # determine plotting x axis. 
  if(is.null(against) & is.null(against.x) & is.null(against.y)) against=1:N
  if(is.null(agaist.x)==F) against = X[,against.x]
  if(is.null(agaist.y)==F) against = y[,against.y]
  
  # determine group variables
  if(is.null(group.x) & is.null(group.by)) group=NULL
  if(is.null(group.x)==F) group = X[,group.x]
  if(is.null(group.by)==F) group = eval(parse(text=group.by))
  
  for(c in var_colNo){
    for(i in 1:j){
      temp.data = cbind(object$marg.list[[c]][,i],against)
      colnames(temp.data)[1] = ynames[i]
      if(is.null(group)==F) temp.data = cbind(temp.data,group)
    }
  }
}