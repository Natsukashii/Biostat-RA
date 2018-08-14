library(gamsel)
library(doParallel)
gamma0 = 0.4
nrep = 500
res = c()
## Parallel Computing
no_cores = detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
res = foreach(i=1:nrep, .combine = rbind, .packages = 'gamsel') %dopar% {
  data = gendata(n=500, p=30,k.lin=10,k.nonlin=10,deg=5,sigma=0.5)
  bases = pseudo.bases(data$X, degree=10, df=5)
  
  ### Fit the logistic model
  
  gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
  gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = gamma0)
  
  lin_list = getActive(gamsel.binout, index= gamsel.bincv$index.1se , type="linear")[[1]]
  nonlin_list = getActive(gamsel.binout, index= gamsel.bincv$index.1se , type="nonlinear")[[1]]
  lin_list = setdiff(lin_list, nonlin_list)
  zero_list = setdiff(c(1:30), c(lin_list,nonlin_list))
  
  lin_corr = sum(c(1:10) %in% lin_list)/10
  non_corr = sum(c(11:20) %in% nonlin_list)/10
  zero_corr = sum(c(21:30) %in% zero_list)/10
  
  lin_in_non = sum(c(11:20) %in% lin_list)/10
  zero_in_non = sum(c(21:30) %in% zero_list)/10
  lin_in_zero = sum(c(21:30) %in% lin_list)/10
  non_in_zero = sum(c(21:30) %in% nonlin_list)/10 
  zero_in_lin = sum(c(1:10) %in% zero_list)/10
  non_in_lin = sum(c(1:10) %in% nonlin_list)/10
  a_res = data.frame(cbind(lin_corr,non_corr,zero_corr,lin_in_non,zero_in_non, lin_in_zero,non_in_zero,zero_in_lin, non_in_lin))
  a_res
}
stopCluster(cl)
colMeans(res)
