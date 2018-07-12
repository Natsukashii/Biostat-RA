library(gamsel)
### Simulate data

data = gendata(n=2000, p=30,k.lin=6,k.nonlin=4,deg=5,sigma=0.5)

bases = pseudo.bases(data$X, degree=10, df=5)
U=do.call("cbind",bases)

gamsel.out = gamsel(data$X, data$y, bases = bases)
par(mfrow=c(1,2),mar=c(5,4,3,1))
summary(gamsel.out)
gamsel.cv=cv.gamsel(data$X,data$y,bases=bases)
par(mfrow=c(1,1))
plot(gamsel.cv)
par(mfrow=c(4,6))

plot(gamsel.out,newx=data$X,index=5, which = c(1,3,9,10,15,16))
plot(gamsel.out,newx=data$X,index=15, which = c(1,3,9,10,15,16))
plot(gamsel.out,newx=data$X,index=25, which = c(1,3,9,10,15,16))
plot(gamsel.out,newx=data$X,index=40, which = c(1,3,9,10,15,16))

### Simulate binary data

gamsel.binout = gamsel(data$X, data$yb, bases = bases, family = "binomial")
summary(gamsel.binout)

par(mfrow=c(4,6))
plot(gamsel.binout,newx=data$X,index=5, which = c(1,3,9,10,15,16))
plot(gamsel.binout,newx=data$X,index=15, which = c(1,3,9,10,15,16))
plot(gamsel.binout,newx=data$X,index=25, which = c(1,3,9,10,15,16))
plot(gamsel.binout,newx=data$X,index=40, which = c(1,3,9,10,15,16))



gamsel.bincv=cv.gamsel(data$X,data$yb,bases=bases, family = "binomial")
par(mfrow=c(1,1))
plot(gamsel.bincv)

par(mfrow=c(6,5))
plot(gamsel.binout,newx=data$X,index=42)
