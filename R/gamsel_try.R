library(tictoc)
library(gamsel)
library(ggplot2)
library(dplyr)
source("control_gendata.R")
source("my_gendata.R")
source("myplot_gamsel.R")

# plot_true <-
#   function(data, p, deg, ylims){
#     y = matrix(NA, nrow = nrow(data$XX), ncol = p)
#     for (j in 1:p){
#       y[,j] = data$U[,((j-1)*deg+1):(j*deg)] %*% data$beta[((j-1)*deg+1):(j*deg)] + data$X[,j] * data$alpha[j]
#     }  
#     if(missing(ylims))ylims=range(y)
#     for (j in 1:p){
#       o = order(data$X[,j])
#       plot(data$X[o,j],y[o,j], type = 'l', ylab=paste("true_f(v",j,")",sep=""),xlab=paste("v",j,sep=""),lwd=2, ylim = ylims)
#     }
#     return(list(ylims = ylims, y=y))
#   }
beta_6 = c(4, -4, 0, 0, 0)
beta_7 = c(8, -24, 16, 0, 0)
beta_8 = c(-3/2, 3, 0, 0, 0)
beta_9 = c(-200/9, 1100/9, -200, 100, 0)
beta_10 = c(3/2, -3, 0, 0, 0)
beta_nonlinear = c(beta_6, beta_7, beta_8, beta_9, beta_10)


gamma0 = 0.4
degree = 5
sample_size =  10000
no_var = 25
fixed_beta = c(rep(0, degree*5), beta_nonlinear, rep(0, degree*(no_var-10)))
data = control_gendata(n=sample_size, p=no_var,k.lin=5,k.nonlin=5,deg=degree,sigma=0.5,
                  fixed_X = matrix(runif(sample_size*no_var), sample_size, no_var), fixed_beta = fixed_beta)
bases = pseudo.bases(data$X, degree=10, df=8)


### Fit the logistic model
# tic("GAMSEL model fitting")
# gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
# gamsel.bincv=cv.gamsel(data$X, data$yb, bases=bases, family = "binomial", gamma = gamma0)
# toc()
# 
# par(mfrow=c(3,5), mars(1,1,1,1))
# my_plot.gamsel(data=data, deg = degree,
#                gamsel.binout, newx=data$X,index=gamsel.bincv$index.min,
#                which = 1:15, rugplot=F, factor = 0.2, type = "binary")

## Define lambda myself
tic("GAMSEL model fitting")
gamsel.binout = gamsel(data$X, data$yb, num_lambda = 50, lambda = seq(5,0.1, length.out = 50), bases = bases, family = "binomial", gamma = gamma0)
toc()
gamsel.bincv=cv.gamsel(data$X, data$yb, lambda = seq(5,0.1, length.out = 50), bases=bases, family = "binomial", gamma = gamma0)
# toc()
par(mfrow=c(3,5), mars(1,1,1,1))
my_plot.gamsel(data=data, deg = degree,
               gamsel.binout, newx=data$X,index=50,
               which = 1:15, rugplot=F, factor = 0.2, type = "binary")

# par(mfrow=c(5,5), mars(1,1,1,1))
# for (num in c(10,20,30,40,50)){
#   my_plot.gamsel(data=data, deg = degree,
#                  gamsel.binout, newx=data$X,index=num,
#                  which = 6:10, rugplot=F, factor = 0.1, type = "binary")
# }




### Fit the linear model

# gamsel.out = gamsel(data$X, data$y, bases = bases, gamma = gamma0)
# gamsel.cv=cv.gamsel(data$X,data$y,bases=bases, gamma = gamma0)
# par(mfrow=c(3,5), mars(1,1,1,1))
# my_plot.gamsel(data=data, deg = degree,
#                gamsel.out,newx=data$X,index=gamsel.cv$index.1se, which = 1:15, type = "notbinary")