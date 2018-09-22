#' Generate summary tables for fmlogit objects
#' 
#' Generate tables of coefficient estimates, partial effects, and willingness to pay from
#' fmlogit-type objects. 
#' 
#' @name summary.fmlogit
#' @aliases summary.fmlogit.margins
#' @aliases summary.fmlogit.wtp
#' 
#' @param object an object with class "fmlogit", "fmlogit.margins", or "fmlogit.wtp". 
#' @param varlist select a subset of variable names to be processed. Default to NULL, of which all variables will
#' be processed.
#' @param sepline whether the output table uses separate lines for coefficients and standard errors. 
#' @param digits number of digits to be signifed. Default to show 3 digits. 
#' @param add.info whether to add additional descriptive information to the output. 
#' @param list whether to output a list object, or a single data frame. 
#' @param sigcode the significance code to be used. Has to be a three-component vector. 
#' @return Either a list (for display purposes) or a data.frame (for csv output purposes). If list return (which is
#' the default) is selected, then the list will contain 4 components: $estimates the estimate; $N number of 
#' observations, $llf value of the log-likelihood function; and $baseline the name of the baseline choice. 
#' 
#' @details This module provides summary methods for three fmlogit objects: \code{fmlogit}, \code{fmlogit.margins}
#' , and \code{fmlogit.wtp}. 
#' 
#' The summary method offers several options to the users. The user can choose for a list output \code{list=T}, which is
#'  good for display and quoting purposes, or a data frame output \code{list=F}, which is good for table outputs. The user
#' can also specify whether to provide additional information other than the parameter estimates, whether to use 
#' seperate lines for the estimates and the standard errors (which mimics the output style in Stata),
#'  as well as the significance code. 
#' 
#' @examples 
#' # generate fmlogit summary
#' #results1 = fmlogit(y,X)
#' 
#' # generate marginal effects summary
#' #effects1 = effects(results1,effect="marginal")
#' summary(effects1)
#' 
#' # generate latex style output
#' # require(xtable)
#' xtable(summary(effects1,list=F,sepline=T))
#' @rdname summary.fmlogit
#' @export summary.fmlogit

############
# generate fmlogit style table
###########

summary.fmlogit = function(object,varlist=NULL,sepline=F,digits=3,add.info=T,list=T,sigcode=c(0.05,0.01,0.001),
                           print=F){
  # define significance code first. 
  asterisk = function(x,k=sigcode){
    if(x>k[1]) return("")
    if(x>k[2]) return("*")
    if(x>k[3]){return("**")}else
    {return("***")}
  }
  # main text  
  # pre matters
  if(!class(object)=="fmlogit") stop("Expect an fmlogit object. Wrong object type given.")
  ynames = names(object[[1]]); Xnames = rownames(object[[1]][[1]])
  if(length(varlist)==0) varlist=Xnames
  var_colNo = which(Xnames %in% varlist)
  j = object$count[3]; K = length(var_colNo)
  if(K < length(varlist)) warning("Some variables requested are not in the variable list. Those variables are omitted.")
  varlist = Xnames[var_colNo]
  # generating tables
  if(!sepline){
    store_mat = matrix(ncol=j-1,nrow=K)
    colnames(store_mat)=ynames
    rownames(store_mat)=Xnames[var_colNo]
    for(i in 1:(j-1)){
      temp_data = signif(object$estimates[[i]][var_colNo,],digits=digits)
      store_mat[,i]=apply(temp_data, 1, function(x) paste(x[1],"(",x[2],")",asterisk(x[4]),sep=""))    
    }}else{
      store_beta = store_se = matrix(ncol=j-1,nrow=K)   
      colnames(store_beta)=ynames
      rownames(store_beta)=varlist
      for(i in 1:(j-1)){
        temp_data = signif(object$estimates[[i]][var_colNo,],digits=digits)
        store_beta[,i]=apply(temp_data,1, function(x) paste(x[1],asterisk(x[4]),sep=""))
        store_se[,i]=apply(temp_data, 1, function(x) paste("(",x[2],")",sep=""))
      }
      for(i in 1:K){
        if(i==1) store_mat=matrix(ncol=j-1)
        store_mat = rbind(store_mat,store_beta[i,],store_se[i,])
      }
      store_mat=store_mat[-1,]
      rownames(store_mat) = rep(" ",length=nrow(store_mat))
      rownames(store_mat)[seq(1,K*2,2)] = varlist
    }
  # output matters
  sig.print = paste("Significance code: 0", "'***'", sigcode[3], "'**'", sigcode[2], "'*'", sigcode[1], "' ", 1)
  if(add.info){
    nc = paste("N=",object$count[1],sep="")
    llf = paste("log pseudo-likelihood=",round(object$likelihood,digits=2),sep="")
    bl = paste("Baseline choice:", object$baseline)
  }
  if(list){
    outlist = list(estimates=store_mat)
    if(add.info){
      outlist$N = nc
      outlist$llf = llf
      outlist$baseline = bl
      outlist$sigcode = sig.print
    }
    if(print){print(outlist)}
    return(outlist)
  }else{
    if(add.info){
      info = matrix(ncol=j-1,nrow=4)
      info[,1] = c(nc,llf,bl,sig.print)
      store_mat = rbind(store_mat,info)
    }
    if(print){print(store_mat)}
    return(as.data.frame(store_mat))
  }
}

##########
# summary for fmlogit.margins
##########

#' @rdname summary.fmlogit
#' @export summary.fmlogit.margins

summary.fmlogit.margins = function(object,varlist=NULL,sepline=F,digits=3,add.info=T,list=T,sigcode=c(0.05,0.01,0.001),
                                   print=F){
  # define significance code first. 
  asterisk = function(x,k=sigcode){
    if(x>k[1]) return("")
    if(x>k[2]) return("*")
    if(x>k[3]){return("**")}else
    {return("***")}
  }
  # main text  
  if(!class(object)=="fmlogit.margins") stop("Expect an fmlogit.margins object. Wrong object type given.")
  ynames = rownames(object[[1]]); Xnames = colnames(object[[1]])
  if(length(varlist)==0) varlist=Xnames
  var_colNo = which(Xnames %in% varlist)
  j = length(ynames); K = length(var_colNo)
  if(K < length(varlist)) warning("Some variables requested are not in the variable list. Those variables are omitted.")
  varlist = Xnames[var_colNo]
  
  # table process
  if(object$R==0) sepline=FALSE
  if(!sepline){
    store_mat = matrix(ncol=j,nrow=K)
    colnames(store_mat)=ynames
    rownames(store_mat)=Xnames
    if(object$R>0){
      for(i in var_colNo){
        temp_data = signif(object$ztable[[i]],digits=digits)
        store_mat[i,]=apply(temp_data, 1, function(x) paste(x[1],"(",x[2],")",asterisk(x[4]),sep=""))    
      }}else{
        store_mat = signif(t(object$effects),digits=digits)
      }
  }else{
    store_beta = store_se = matrix(ncol=j,nrow=K)   
    colnames(store_beta)=ynames
    rownames(store_beta)=Xnames
    for(i in var_colNo){
      temp_data = signif(object$ztable[[i]],digits=digits)
      store_beta[i,]=apply(temp_data,1, function(x) paste(x[1],asterisk(x[4]),sep=""))
      store_se[i,]=apply(temp_data, 1, function(x) paste("(",x[2],")",sep=""))
    }
    for(i in 1:K){
      if(i==1) store_mat=matrix(ncol=j)
      store_mat = rbind(store_mat,store_beta[i,],store_se[i,])
    }
    store_mat=store_mat[-1,]
    rownames(store_mat) = rep("",length=nrow(store_mat))
    rownames(store_mat)[seq(1,K*2,2)] = varlist
  }
  # output matters
  sig.print = paste("Significance code: 0", "'***'", sigcode[3], "'**'", sigcode[2], "'*'", sigcode[1], "' ", 1)
  if(add.info){
    expl = object$expl
  }
  if(list){
    outlist = list(estimates=store_mat)
    if(add.info){
      outlist$expl = expl
      outlist$sigcode = sig.print
    }
    if(print){print(outlist)}
    return(outlist)
  }else{
    if(add.info){
      info = matrix(ncol=j,nrow=2)
      info[,1] = c(expl,sig.print)
      store_mat = rbind(store_mat,info)
    }
    if(print){print(store_mat)}
    return(as.data.frame(store_mat))
  }
}

############
# generate willingness to pay tables
############

#' @rdname summary.fmlogit
#' @export summary.fmlogit.wtp

summary.fmlogit.wtp = function(object,varlist=NULL,sepline=F,digits=3,sigcode=c(0.05,0.01,0.001),
                               print=F){
  # define significance code first. 
  asterisk = function(x,k=sigcode){
    if(x>k[1]) return("")
    if(x>k[2]) return("*")
    if(x>k[3]){return("**")}else
    {return("***")}
  }
  # main text  
  if(!class(object)=="fmlogit.wtp") stop("Expect an fmlogit.wtp object. Wrong object type given.")
  if(colnames(object$wtp)[1]!="estimate") return(object$wtp) # no need to summary. 
  Xnames = rownames(object$wtp)
  if(length(varlist)==0) varlist=Xnames
  var_colNo = which(Xnames %in% varlist)
  K = length(var_colNo)
  if(K < length(varlist)) warning("Some variables requested are not in the variable list. Those variables are omitted.")
  varlist = Xnames[var_colNo]
  sig.print = paste("Significance code: 0", "'***'", sigcode[3], "'**'", sigcode[2], "'*'", sigcode[1], "' ", 1)
  if(!sepline){
  # table process
  store_mat = apply(signif(object$wtp[var_colNo,],digits=digits), 1, function(x) paste(x[1],"(",x[2],")",asterisk(x[4]),sep=""))
  store_mat = as.data.frame(store_mat)
  colnames(store_mat)=NULL
  # output matters
  }else{          
    store_beta=apply(signif(object$wtp[var_colNo,],digits=digits),1, function(x) paste(x[1],asterisk(x[4]),sep=""))
    store_se=apply(signif(object$wtp[var_colNo,],digits=digits), 1, function(x) paste("(",x[2],")",sep=""))
    for(i in 1:K){
      if(i==1) store_mat=vector()
      store_mat = c(store_mat,store_beta[i],store_se[i])
    }
    names(store_mat) = rep("",length=length(store_mat))
    names(store_mat)[seq(1,K*2,2)] = varlist
  }
  if(print){print(store_mat);print(sig.print)}    
  return(store_mat)
}
