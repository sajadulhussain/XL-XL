# Normalization condition
norm_cond<-function(theta,beta) {
  fun <- function(x) dXLXL(x, theta,beta)
  return (as.numeric(integrate(Vectorize(fun), lower = 0, upper = Inf)[1]))}

#CDF of XL distribution
pXL <- function(x,theta){
  return(1-(1+theta*x/(theta+1)^2)*exp(-theta*x))}

#PDF of XL distribution
dXL <- function(x,theta){
  return(theta*theta*(2+theta+x)*exp(-theta*x)/(1+theta)/(1+theta))}

#CDF of XLXL distribution
pXLXL <- function(x,theta,beta){
  return(pXL(x,theta)/(pXL(x,theta)+1-pXL(x,beta)))}

#PDF of XLXL distribution
dXLXL <- function(x,theta,beta){
  L=dXL(x,theta)*(1-pXL(x,beta))+pXL(x,theta)*dXL(x,beta)
  M=(pXL(x,theta)+1-pXL(x,beta))^2
  return(L/M)}

# Survival function of XLXL distribution
sXLXL <- function(x,theta,beta){
  return(1-pXLXL(x,theta,beta))}

# Hazard rate function of XLXL distribution
hXLXL <- function(x, theta,beta){
    return(dXLXL(x,theta,beta)/sXLXL(x,theta,beta))}

# Quantile function of XLXL distribution
qXLXL=function(p,theta,beta){
  u1 = function(x,theta,beta) (pXLXL(x,theta,beta)-p)
  return(uniroot(u1, c(0,1e06), tol = 0.0000000001, theta=theta, beta=beta)$root)}

# Pseudo-random number generator
rXLXL =function(n,theta,beta) {
  x=numeric(n)
  for (i in 1:n) x[i]=qXLXL(runif(1,0,1),theta,beta)
  return(sort(x))}

# Ordinary moments
mXLXL=function(k,theta,beta) {
  return(integral(function(x) x^k*dXLXL(x,theta,beta), 0, Inf, reltol = 1e-12, method = "Simpson"))}

# Incomplete moments
imXLXL=function(k,t,theta,beta) {
  return(integral(function(x) x^k*dXLXL(x,theta,beta), 0, t, reltol = 1e-12, method = "Simpson"))}

# Characteristics
chXLXL=function(theta,beta){
  x=numeric(5)
  # mean
  x[1]=mXLXL(1,theta,beta)
  #variance
  x[2]=mXLXL(2,theta,beta)-x[1]^2
  # coefficient of variation
  x[3]=sqrt(x[2])/x[1]
  w1=mXLXL(3,theta,beta)-3*mXLXL(1,theta,beta)*mXLXL(2,theta,beta)+2* x[1]^3
  #skewness
  x[4]=w1/(x[2])^(1.5) 
  w2=mXLXL(4,theta,beta)-4*mXLXL(1,theta,beta)*mXLXL(3,theta,beta)+6*mXLXL(1,theta,beta)^2*
    mXLXL(2,theta,beta)-3*x[1]^4
  #kurtosis
  x[5]=w2/(x[2])^2 
  return(x)}

# Moments generating function
mgfXLXL=function(t,theta,beta) {
  return(integral(function(x) exp(t*x)*dXLXL(x,theta,beta), 0, Inf, reltol = 1e-12, method = "Simpson"))}

#Pdf of order statistics
dosXLXL=function(x,r,n,theta,beta) {
  pdf=dXLXL(x,theta,beta)
  cdf=pXLXL(x,theta,beta)
  return(fact(n)/fact(r-1)/fact(n-r)*pdf*cdf^(r-1)*(1-cdf)^(n-r))}

# Moments of order statistics
mosXLXL=function(k,r,n,theta,beta) {
  return(integral(function(x) x^k*dosAPTZLD(x,r,n,theta,beta), 0, Inf, reltol = 1e-12, method = "Simpson"))
}

