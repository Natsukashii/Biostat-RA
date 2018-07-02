library(dplyr)
source("my_gendata.R")
source("myplot_gamsel.R")
source("myplot_gam.R")
beta_nonlinear = c(1, 20, 0.2, -15, 0.1)
gamma0 = 0.4
degree = 5
sample_size =  2000
no_var = 25
data = my_gendata(n=2000, p=25,k.lin=5,k.nonlin=5,deg=degree,sigma=0.5,
                  fixed_X = matrix(runif(sample_size*no_var), sample_size, no_var),fixed = F)
tbl_data = as_tibble(data.frame(cbind(data$X,data$yb)))


### Fit the logistic model using step.Gam

null_fit = gam(X26 ~ 1, family="binomial", tbl_data)
gam_scope = gam.scope(tbl_data, response = 26, arg = "df=5")
gam_fit = step.Gam(null_fit, gam_scope)

par(mfrow=c(5,5))
my_plot.gam(data, gam_fit, deg = 5, gam_scope, factor = 0.1)

### Compared to the model using GAMSEL
bases = pseudo.bases(data$X, degree=10, df=5)
gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
gamsel.bincv=cv.gamsel(data$X, data$yb, bases=bases, family = "binomial", gamma = gamma0)
par(mfrow=c(5,5), mars(1,1,1,1))
my_plot.gamsel(data=data, deg = degree,
               gamsel.binout, newx=data$X,index=gamsel.bincv$index.1se, 
               which = 1:25, rugplot=F, factor = 0.1, type = "binary")