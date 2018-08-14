library(gamsel)
gamma0 = 0.4
### Simulate data

data = gendata(n=200, p=25,k.lin=6,k.nonlin=4,deg=5,sigma=0.5)
bases = pseudo.bases(data$X, degree=10, df=5)

### Fit the model

gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
par(mfrow=c(1,2),mar=c(5,4,3,1))
summary(gamsel.binout)

gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = gamma0)
par(mfrow=c(1,1))
plot(gamsel.bincv)

par(mfrow=c(5,5))
plot(gamsel.binout,newx=data$X,index=gamsel.bincv$index.1se)


### Try different n
data = gendata(n=500, p=25,k.lin=6,k.nonlin=4,deg=5,sigma=0.5)
bases = pseudo.bases(data$X, degree=10, df=5)
gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = gamma0)
par(mfrow=c(5,5))
plot(gamsel.binout,newx=data$X,index=gamsel.bincv$index.1se)

data = gendata(n=1000, p=25,k.lin=6,k.nonlin=4,deg=5,sigma=0.5)
bases = pseudo.bases(data$X, degree=10, df=5)
gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial",gamma = gamma0)
gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = gamma0)
par(mfrow=c(5,5))
plot(gamsel.binout,newx=data$X,index=gamsel.bincv$index.1se)

data = gendata(n=10000, p=25,k.lin=6,k.nonlin=4,deg=5,sigma=0.1)
bases = pseudo.bases(data$X, degree=10, df=5)
gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial", gamma = gamma0)
gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial", gamma = gamma0)
par(mfrow=c(5,5))
par(mar=c(1,1,1,1))
plot(gamsel.binout,newx=data$X,index=gamsel.bincv$index.1se)

y_hat = predict.gamsel(gamsel.binout, newdata = data$X, index = gamsel.bincv$index.1se, type="response")>0.5
mean(data$yb !=y_hat) 
