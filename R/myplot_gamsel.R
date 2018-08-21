library(gamsel)
my_plot.gamsel <-
  function(data, deg, x, newx, index,which=1:p,rugplot=TRUE,ylims, factor = 1.0, type="binary", ...){
    gamsel.out=x
    x=newx
    if(missing(index))index=length(gamsel.out$lambdas)
    pJ=dim(gamsel.out$alpha)
    p=pJ[1]
    maxJ=pJ[2]
    nonlin <- getActive(gamsel.out, index,type="nonlinear")[[1]]
    active_linear <- getActive(gamsel.out, index,type="linear")[[1]]
    lin=setdiff(active_linear,nonlin)
    zero=setdiff(1:p,union(nonlin,lin))
    if(which[1]=="nonzero")which=seq(p)[-zero]
    if(which[1]=="nonlinear")which=nonlin
    if(is.null(which)){warning("Nothing to plot with this choice of argument which");return()}
    colvals=rep("blue",p)
    colvals[lin]="green"
    colvals[nonlin]="red"
    termmat=drop(predict.gamsel(gamsel.out,x,index=index,type="terms"))
    lambda=gamsel.out$lambda[index]
    
    # y = matrix(NA, nrow = nrow(data$X), ncol = p)
    # for (j in 1:p){
    #   y[,j] = data$U[,((j-1)*deg+1):(j*deg)] %*% data$beta[((j-1)*deg+1):(j*deg)] + data$X[,j] * data$alpha[j]
    # }
    y = data$y_mat
    if(missing(ylims))ylims=range(y)
    if (type == "binary"){
      y_bin = predict.gamsel(gamsel.out,x,index=index,type = "response") > 0.5
      y_bin = as.numeric(y_bin)
    }
    # if (type == "counts"){
    #   y_count = predict.gamsel(gamsel.out,x,index=index,type = "response") > 0.5
    #   y_count = as.numeric(y_bin)
    # }
    for(ell in which){
      o=order(x[,ell])
      plot(x[o,ell], termmat[o,ell], type='l', ylab=paste("f(v",ell,")",sep=""),xlab=paste("v",ell,sep=""),ylim=ylims,col=colvals[ell],lwd=2,...)
      o2 = order(data$X[,ell])
      lines(data$X[o2,ell], y[o2,ell])
      if (type == "binary"){
        par(new = T)
        plot(x[o,ell], jitter(y_bin[o], factor = factor), pch='.', type = "p", xaxt = "n", yaxt = "n", ylab = "", xlab = "", lty = 3)
        axis(side = 4)
      }
      if(rugplot)rug(x[,ell])
    }
  }