library(tictoc)
library(dplyr)
library(gam)
source("my_gendata.R")
source("control_gendata.R")
source("myplot_gamsel.R")
source("myplot_gam.R")
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

tbl_data = as_tibble(data.frame(cbind(data$X,data$yb)))


### Fit the logistic model using GAM
gam_scope = gam.scope(tbl_data, response = 26, arg = "df=8", form = F)
var_c = c()
for (i in 1:25) {
  var_c = c(var_c, gam_scope[[i]][3])
}
tic("model fitting")
full_fit = gam(as.formula(paste("X26~", paste(var_c, collapse = "+"))), family="binomial", tbl_data)
toc()
par(mfrow=c(3,5), mars(1,1,1,1))
my_plot.gam(data, full_fit, deg = degree, gam_scope, which = 1:15, factor = 0.1)


 ### Fit the logistic model using step.Gam
# gam_scope = gam.scope(tbl_data, response = 26, arg = c("df=3", "df=4", "df=5"))
# tic("step.Gam model fitting")
# null_fit = gam(X26 ~ 1, family="binomial", tbl_data)
# step_fit = step.Gam(null_fit, gam_scope)
# toc()
# par(mfrow=c(3,5))
# my_plot.gam(data, step_fit, deg = 5, gam_scope, which = 1:15, factor = 0.1)





# ### Compared to the model using GAMSEL
# bases = pseudo.bases(data$X, degree=10, df=5)
# gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
# gamsel.bincv=cv.gamsel(data$X, data$yb, bases=bases, family = "binomial", gamma = gamma0)
# par(mfrow=c(5,5), mars(1,1,1,1))
# my_plot.gamsel(data=data, deg = degree,
#                gamsel.binout, newx=data$X,index=gamsel.bincv$index.1se, 
#                which = 1:25, rugplot=F, factor = 0.1, type = "binary")