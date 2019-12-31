# I. Introduction
library(boot)
data("urine")
# Cleaning the data
urine1 <- urine[complete.cases(urine[ , c(1:7)]),] # removing all rows with "NA" in any category
#dim(urine1) # 7 variables
#urine1$r <- as.factor(urine1$r)
#head(urine1)

# II. Exploratory Data Analysis
pairs(urine1[-c(1)], pch = 19)
cor(urine1)
cor(urine1$osmo, urine1$urea)

par(mfrow=c(2,4))
hist(urine1$r)
hist(urine1$gravity)
hist(urine1$ph)
hist(urine1$osmo)
hist(urine1$cond)
hist(urine1$urea)
hist(urine1$calc)

# III. Methods

## Frequentist Regression Model
fit <- glm(r ~ ., family = binomial, data = urine1)
fit.sum <- summary(fit)
coef(fit.sum)

fit1 <- glm(r ~ gravity + ph + cond + urea + calc, family = binomial, data = urine1)
fit1.sum <- summary(fit1)
coef(fit1.sum)

fit2 <- glm(r ~ gravity + cond + urea + calc, family = binomial, data = urine1)
fit2.sum <- summary(fit2)
coef(fit2.sum)

## Generalized Linear Model in the Bayesian Framework
   
# A function to evaluate the log of the posterior density
    logP=function(y,X,b,b0,varB){
      Xb=X%*%b
      theta=exp(Xb)/(1+exp(Xb))
      logLik=sum( dbinom(x=y,p=theta,size=1,log=T)  )
      logPrior=sum(  dnorm(x=b,sd=sqrt(varB),mean=b0,log=T))
      return(logLik+logPrior)
    }
    
    logisticRegressionBayes=function(y,X,nIter=20000,V=.02,varB=rep(10000,ncol(X)),b0=rep(0,ncol(X))){
      ####### Arguments #######################
      # y  a vector with 0/1 values
      # X  incidence matrix of effects
      # b0,varB, the prior mean and prior variance bj~N(b0[j],varB[j])
      # V the variance of the normal distribution used to generate candidates~N(b[i-1],V)
      # nIter: number of iterations of the sampler
      # Details: generates samples from the posterior distribution of a logistic regression
      # using a Metropolis algorithm
      #########################################
      
      # A matrix to store samples
      p=ncol(X)
      B=matrix(nrow=nIter,ncol=p)
      colnames(B)=colnames(X)
      
      # A vector to trace acceptance
      accept=matrix(nrow=nIter,ncol=p,NA)
      accept[1,]=TRUE 
      
      # Initialize
      B[1,]=0
      B[1,1]=log(mean(y)/(1-mean(y)))
      b=B[1,]
      for(i in 2:nIter){
        for(j in 1:p){
          candidate=b
          candidate[j]=rnorm(mean=b[j],sd=sqrt(V),n=1)
          logP_current=logP(y,X,b0=b0,varB=varB,b=b)
          logP_candidate=logP(y,X,b0=b0,varB=varB,b=candidate)
          r=min(1,exp(logP_candidate-logP_current))
          delta=rbinom(n=1,size=1,p=r)
          accept[i,j]=delta
          if(delta==1){ b[j]=candidate[j] }
        }
        B[i,]=b
        #if(i%%1000==0){
        #  message(" Iteration ",i)
        #}
      }
      return(list(B=B,accept=accept))
    }
    
    Z = as.matrix(model.matrix(~ osmo + cond + urea + calc, data = urine1)) #[,-1]
    # collecting 40,000 samples
    samples = logisticRegressionBayes(y=urine1$r, X=cbind(Z), nIter=20000)
    # discarding 20,000 iterations as burn-in
    #cbind(fit1$coef, colMeans(samples$B[-c(1:5000),])) 
    library(coda)
    samples$B = as.mcmc(samples$B[-c(1:10000),])

summary(samples$B)

## trace plots of a few parameters (plotting 10,000 iterations after burning in 10,000 iterations)
par(mfrow=c(2,2))
# osmo
plot(samples$B[1:10000, 2],type='o',main="osmo") # betas
# cond
plot(samples$B[1:10000, 3],type='o',main="cond") # betas
# urea
plot(samples$B[1:10000, 4],type='o',main="urea") # betas
# calc
plot(samples$B[1:10000, 5],type='o',main="calc") # betas

effectiveSize(samples$B)
sqrt(.02/effectiveSize(samples$B))

contrastA <- c(0,mean(urine1$osmo),mean(urine1$cond),mean(urine1$urea),mean(urine1$calc)) # 0 for intercept, average values of osmo, cond, urea, and calc
tmpA <- samples$B%*%contrastA
# plot histogram on top of density function
hist(tmpA, freq = FALSE, main = "Posterior Density")
lines(density(tmpA))
# 95% posterior credibility region
abline(v=0)
abline(v = HPDinterval(as.mcmc(tmpA), p=.95), col = 2, lty = 2) 