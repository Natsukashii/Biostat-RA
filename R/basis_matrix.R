library(splines)
library(gamsel)
library(assist)
### Generate truncated DR basis
n = 50
x = rnorm(n)
sort.x = sort(x)
pU = basis.gen(sort.x, df=6, degree = 10)
par(mfrow = c(3,4))
plot(sort.x,pU[,1],type = 'l')
plot(sort.x,pU[,2],type = 'l')
plot(sort.x,pU[,3],type = 'l')
plot(sort.x,pU[,4],type = 'l')
plot(sort.x,pU[,5],type = 'l')
plot(sort.x,pU[,6],type = 'l')
plot(sort.x,pU[,7],type = 'l')
plot(sort.x,pU[,8],type = 'l')
plot(sort.x,pU[,9],type = 'l')
plot(sort.x,pU[,10],type = 'l')

### Generate DR basis
smooth.matrix = function(x, df){
 n = length(x);
 A = matrix(0, n, n);
 for(i in 1:n){
     y = rep(0, n); y[i]=1;
     yi = smooth.spline(x, y, df = df)$y;
     A[,i]= yi;
  }
  lambda = smooth.spline(x, y, df = df)$lambda 
  return(list("mat" = (A+t(A))/2, "lam"=lambda))
}
sm = smooth.matrix(sort.x, df=6)
sm_mat = sm$mat
lambda = sm$lam
# sum(diag(sm_mat))

inv.sm_mat = solve(sm_mat)
I = diag(rep(1, n))
K = (inv.sm_mat-I)/lambda
eig.K = eigen(K)
eig.K$values
U = eig.K$vectors[,49:1]
par(mfrow = c(3,4))
plot(sort.x,U[,1],type = 'l')
plot(sort.x,U[,2],type = 'l')
plot(sort.x,U[,3],type = 'l')
plot(sort.x,U[,4],type = 'l')
plot(sort.x,U[,5],type = 'l')
plot(sort.x,U[,6],type = 'l')
plot(sort.x,U[,7],type = 'l')
plot(sort.x,U[,8],type = 'l')
plot(sort.x,U[,9],type = 'l')
plot(sort.x,U[,10],type = 'l')


