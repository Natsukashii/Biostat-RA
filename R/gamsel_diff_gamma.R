library(gamsel)
library(ggplot2)
library(dplyr)
source("control_gendata.R")
source("my_gendata.R")
source("myplot_gamsel.R")

beta_6 = c(4, -4, 0, 0, 0)
beta_7 = c(8, -24, 16, 0, 0)
beta_8 = c(-3/2, 3, 0, 0, 0)
beta_9 = c(-200/9, 1100/9, -200, 100, 0)
beta_10 = c(3/2, -3, 0, 0, 0)
beta_nonlinear = c(beta_6, beta_7, beta_8, beta_9, beta_10)

gamma0 = 0.4
degree = 5
sample_size =  2000
no_var = 25
fixed_beta = c(rep(0, degree*5), beta_nonlinear, rep(0, degree*(no_var-10)))
data = control_gendata(n=2000, p=25,k.lin=5,k.nonlin=5,deg=degree,sigma=0.5,
                  fixed_X = matrix(runif(sample_size*no_var), sample_size, no_var), fixed_beta = fixed_beta)
bases = pseudo.bases(data$X, degree=10, df=5)
s
### Fit the logistic model

# gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
# gamsel.bincv=cv.gamsel(data$X, data$yb, bases=bases, family = "binomial", gamma = gamma0)
# par(mfrow=c(3,5), mars(1,1,1,1))
# my_plot.gamsel(data=data, deg = degree,
#                gamsel.binout, newx=data$X,index=gamsel.bincv$index.1se,
#                which = 1:15, rugplot=F, factor = 0.2, type = "binary")

### Fit the linear model

# gamsel.out = gamsel(data$X, data$y, bases = bases, gamma = gamma0)
# gamsel.cv=cv.gamsel(data$X,data$y,bases=bases, gamma = gamma0)
# par(mfrow=c(3,5), mars(1,1,1,1))
# my_plot.gamsel(data=data, deg = degree,
#                gamsel.out,newx=data$X,index=gamsel.cv$index.1se, which = 1:15, type = "notbinary")

 
### Different gamma
gamma_list = c(0.1,0.4,0.5,0.7,0.9)
par(mfrow=c(3,5), mars(1,1,1,1))
for (each in gamma_list){
  gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = each)
  gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = each)
  my_plot.gamsel(data=data, deg = degree,
                 gamsel.binout,newx=data$X,index=gamsel.bincv$index.1se,
                 which = 6:10, rugplot=F, factor = 0.2, type = "binary")
}
