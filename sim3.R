library(gamsel)
library(doParallel)
source("control_gendata.R")
### Set the paramters
beta_6 = c(4, -4, 0, 0, 0)
beta_7 = c(8, -24, 16, 0, 0)
beta_8 = c(-3/2, 3, 0, 0, 0)
beta_9 = c(-200/9, 1100/9, -200, 100, 0)
beta_10 = c(3/2, -3, 0, 0, 0)
beta_nonlinear = c(beta_6, beta_7, beta_8, beta_9, beta_10)
beta_nonlinear = c(beta_nonlinear, -beta_nonlinear)

degree = 5
sample_size =  2000
no_var = 30
fixed_beta = c(rep(0, degree*10), beta_nonlinear, rep(0, degree*(no_var-10)))

gamma0 = 0.1
nrep = 500
res = c()

## Parallel Computing
no_cores = detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
res = foreach(i=1:nrep, .combine = rbind, .packages = 'gamsel') %dopar% {
  data = control_gendata(n=sample_size, p=30,k.lin=10,k.nonlin=10,deg=degree,sigma=0.5,
                         fixed_X = matrix(runif(sample_size*no_var), sample_size, no_var), fixed_beta = fixed_beta)
  bases = pseudo.bases(data$X, degree=10, df=5)
  
  ### Fit the logistic model
  
  gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
  gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = gamma0)
  
  lin_list = getActive(gamsel.binout, index= gamsel.bincv$index.1se , type="linear")[[1]]
  nonlin_list = getActive(gamsel.binout, index= gamsel.bincv$index.1se , type="nonlinear")[[1]]
  lin_list = setdiff(lin_list, nonlin_list)
  zero_list = setdiff(c(1:30), c(lin_list,nonlin_list))
  
  
  lin_ppv = sum(c(1:10) %in% lin_list)/length(lin_list)
  non_ppv = sum(c(11:20) %in% nonlin_list)/length(nonlin_list)
  zero_ppv = sum(c(21:30) %in% zero_list)/length(zero_list)
  
  lin_sen = sum(lin_list %in% c(1:10))/10
  non_sen = sum(nonlin_list %in% c(11:20))/10
  zero_sen = sum(zero_list%in% c(21:30))/10
  
  lin_in_non = sum(c(11:20) %in% lin_list)/10
  zero_in_non = sum(c(21:30) %in% zero_list)/10
  lin_in_zero = sum(c(21:30) %in% lin_list)/10
  non_in_zero = sum(c(21:30) %in% nonlin_list)/10 
  zero_in_lin = sum(c(1:10) %in% zero_list)/10
  non_in_lin = sum(c(1:10) %in% nonlin_list)/10
  a_res = data.frame(cbind(lin_ppv,non_ppv,zero_ppv,lin_sen,non_sen,zero_sen,
                           lin_in_non,zero_in_non, lin_in_zero,non_in_zero,zero_in_lin, non_in_lin))
  a_res
}
stopCluster(cl)
colMeans(res,na.rm = TRUE)
