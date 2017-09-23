#' Plot marginal or discrete effects, at each observation & for each choice
#' 
#' Plot the desired effect at each observed value for each choice
#' 
#' @param object An "fmlogit.margins" object.
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
#' @export plot.fmlogit.margins



plot.fmlogit.margins = function(object,varlist=NULL,X=NULL,y=NULL, 
                                against=NULL,against.x=NULL,against.y=NULL,
                                group.x=NULL, group.algebra=NULL,
                                mfrow=NULL){
  require(ggplot2)
  require(grid)
  
  if(is.null(object[["marg.list"]])) stop("Please choose marg.type=aveacr when calculating effects")
  k = ncol(object$effects); j = nrow(object$effects); N = nrow(object$marg.list[[1]]); 
  Xnames = colnames(object$effects) ; ynames = rownames(object$effects)
  X = object$X; y=object$y
  
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
  }else{
    jr = mfrow[1]; js = mfrow[2]
  }
  
  # determine plotting x axis. 
  if(is.null(against) & is.null(against.x) & is.null(against.y)) {M.against=1:N; ag.name = "ObsNo"}
  if(is.null(against.x)==F) {M.against = X[,against.x]; ag.name = against.x}
  if(is.null(against.y)==F) {M.against = y[,against.y]; ag.name = against.y}
  
  
  # determine group variables
  if(is.null(group.x) & is.null(group.algebra)) {M.group=NULL; g.name=NULL}
  if(is.null(group.x)==F) {M.group = X[,group.x]; g.name.display <- g.name <- group.x;}
  if(is.null(group.algebra)==F) {
    M.group = eval(parse(text=paste("X[,",'"',group.x,'"',"]",group.algebra,sep="")))
    M.group = ifelse(M.group,"Yes","No")
    g.name = group.x
    g.name.display = paste(group.x,group.algebra,sep="")
    }
  
  for(c in var_colNo){
    ggplot()
    pushViewport(viewport(layout = grid.layout(jr, js)))
    temp.data = cbind(object$marg.list[[c]],M.against)
    temp.data = as.data.frame(temp.data)
    colnames(temp.data) = c(colnames(object$marg.list[[c]]),ag.name)
    if(is.null(M.group)==F){
      temp.data = cbind(temp.data,as.factor(M.group))
      colnames(temp.data)[-1] = g.name}
    for(i in 1:j){
      g <- ggplot(temp.data,aes_string(ag.name,ynames[i],color=g.name)) + geom_point() 
      g <- g + geom_hline(yintercept = 0) + theme_classic() + ggtitle(paste("Effects on", Xnames[c]))
      if(is.null(M.group)==F) g <- g + theme(legend.title = element_text(colour="black"))+
        scale_color_discrete(name=g.name.display)
      print(g,vp = viewport(layout.pos.row = ifelse(i%%jr==0,jr,i%%jr), layout.pos.col = (i-1) %/%js + 1) )
  }
}}