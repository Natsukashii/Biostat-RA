library(gam)
my_plot.gam <-
  function(data, x, deg, scope, which=1:p,rugplot=TRUE, ylims, factor = 1.0, type="binary", ...){
    p = ncol(data$X)
    gam_fit=x
    x=data$X
    gam_scope = scope
    terms = attr(gam_fit$terms, "term.labels")
    linear = vector(length = p)
    lin_pattern = paste0(names(gam_scope),'$')
    for (each in 1:length(lin_pattern)){
      if (length(grep(lin_pattern[each], terms))>0) {linear[each]=T}
      else {linear[each]=F}
    }
    linear = which(linear==TRUE)
    
    non_linear = vector(length = p)
    non_lin_pattern = paste0('s\\(',names(gam_scope), '\\,')
    for (each in 1:length(non_lin_pattern)){
      if (length(grep(non_lin_pattern[each], terms))>0) {non_linear[each]=T}
      else {non_linear[each]=F}
    }
    non_linear = which(non_linear==TRUE)
    zero = setdiff(c(1:p), union(linear, non_linear))
    pred = predict(gam_fit, type = "terms")
    colvals=rep("blue",p)
    colvals[linear]="green"
    colvals[non_linear]="red"
    
    
    y = matrix(NA, nrow = nrow(data$X), ncol = p)
    for (j in 1:p){
      y[,j] = data$U[,((j-1)*deg+1):(j*deg)] %*% data$beta[((j-1)*deg+1):(j*deg)] + data$X[,j] * data$alpha[j]
    }  
    if(missing(ylims))ylims=range(y)
    if (type == "binary"){
      y_bin = predict(gam_fit, type = "response") > 0.5
      y_bin = as.numeric(y_bin)
    }
    for (each in 1:p){
      if (each %in% linear){
        o=order(data$X[,each])
        index = grep(lin_pattern[each], colnames(pred))
        plot(data$X[o,each], pred[o,index], type='l', ylab=paste("f(v",each,")",sep=""),xlab=paste("v",each,sep=""),ylim = ylims,col=colvals[each], lwd=2)
        lines(data$X[o,each], y[o,each])
        if (type == "binary"){
          par(new = T)
          plot(x[o,each], jitter(y_bin[o], factor = factor), pch='.', type = "p", xaxt = "n", yaxt = "n", ylab = "", xlab = "", lty = 3)
          axis(side = 4)
        }
        if(rugplot)rug(x[,each])
      } 
      if (each %in% non_linear){
        o=order(data$X[,each])
        index = grep(non_lin_pattern[each], colnames(pred))
        plot(data$X[o,each], pred[o,index], type='l', ylab=paste("f(v",each,")",sep=""),xlab=paste("v",each,sep=""),ylim = ylims,col=colvals[each], lwd=2)
        lines(data$X[o,each], y[o,each])
        if (type == "binary"){
          par(new = T)
          plot(x[o,each], jitter(y_bin[o], factor = factor), pch='.', type = "p", xaxt = "n", yaxt = "n", ylab = "", xlab = "", lty = 3)
          axis(side = 4)
        }
        if(rugplot)rug(x[,each])
      } 
      if (each %in% zero){
        o=order(data$X[,each])
        plot(data$X[o,each],rep(0, nrow(data$X)), type = "l", ylab=paste("f(v",each,")",sep=""),xlab=paste("v",each,sep=""),ylim = ylims,col=colvals[each], lwd=2)
        lines(data$X[o,each], y[o,each])
        if (type == "binary"){
          par(new = T)
          plot(x[o,each], jitter(y_bin[o], factor = factor), pch='.', type = "p", xaxt = "n", yaxt = "n", ylab = "", xlab = "", lty = 3)
          axis(side = 4)
        }
        if(rugplot)rug(x[,each])
      } 
    }
  }