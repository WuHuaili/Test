set.seed(1000)

r = 0.05
q = 0
rf = 0.04
sigma1 = 0.1
sigma2 = 0.15
gamma = -0.04
lambda = 1.5
s0 = 6000
e0 = 0.0096
K = 60
T = 1
N = 100000
delta = 0.01
rho = -0.25

St = matrix(s0, nrow = N, ncol = 101)
Et = matrix(e0, nrow = N, ncol = 101)

randNum = matrix(rnorm(N*100*2,0,1),nrow = N)
Wt = sqrt(delta)*randNum[,1:100]
Bt = rho*Wt + sqrt(delta)*sqrt(1-rho^2)*randNum[,101:200]
Jt = matrix(rpois(N*100,lambda*delta),nrow = N)

for(i in 1:100){
  St[,i+1] = St[,i]*(1+(r-q)*delta+sigma1*Wt[,i]+gamma*Jt[,i])
  Et[,i+1] = Et[,i]+(r-rf)*Et[,i]*delta+sigma2*Et[,i]*Bt[,i]
}

payoff = rep(0,N)
for(i in 1:N){
  payoff[i] = max(St[i,101]*Et[i,101]-K,0)*exp(-r*T)
}
mean(payoff)