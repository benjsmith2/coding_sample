# Simulation to help determine the sample size needed for a DNAm ~ PTSD 
# association study. 

library(ggplot2)

set.seed(19700101)

N = seq(50,300,25) # number of individuals 
R2 = seq(0.05,.4,0.05) # R^2 value be used

Pars = expand.grid(N=N,R2=R2,Power = NA) # create matrix to store output
nRep = 20000 # number of repetitions to run at each sample size and R^2

i = 1
j = 1
for (j in 1:nrow(Pars)) {
  pval = rep(NA,nRep)
  for (i in 1:nRep) {
    x = rnorm(n = Pars$N[j])
    b=1
    signal = x*b
    vSignal = var(signal)
    #simulate the amount of noise (error) to add to achieve a given R^2 
    error = rnorm(n = Pars$N[j], mean = 0, 
                  sd = sqrt(vSignal / Pars$R2[j]*(1-Pars$R2[j])))
    y = signal + error # y = error for t1 error
    fit = lm(y~x)
    pval[[i]] = summary(fit)$coef[2,4]
  }
  Pars$Power[j] = mean(pval<(0.05/nRep))
}

Pars$N <- as.character(Pars$N)

ggplot(data = Pars, aes(x = R2, y = Power, group = N))+ 
  geom_point(aes(color = N))+
  geom_line(aes(color = N))

Pars$N <- as.numeric(Pars$N)
Pars$R2 <- as.character(Pars$R2)

ggplot(data = Pars, aes(x = N, y = Power, group = R2))+ 
  geom_point(aes(color = R2))+
  geom_line(aes(color = R2))


mean(pval<(0.05/1e5))

