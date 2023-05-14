library(MASS)
petrol
data=petrol[,2:6]
data
n=length(data[,1])
##a
eigen(cov(data))
lamda1.hat=eigen(cov(data))$values[1]
lamda1.hat
sum(eigen(cov(data))$values)
lambda1=function(data) (lamda1.hat=eigen(cov(data))$values[1]/sum(eigen(cov(data))$values))
lambda1(data)

jknife(data,lambda1,fxn_dim=2)

##b
bootstrap_CIs(data,fxn=lambda1,fxn_dim=2,B=1000,alfa=0.05)

##c
head(data)
detach(data)
attach(data)
lr=lm(Y~ SG+VP+V10+EP)
summary(lr)
anova(lr)

y=data[,5]
X=as.matrix(cbind(1,data[,1:4]))
B = 2000;index.star = matrix(sample(1:n,B*n,replace=T),nrow=B);

Err.hat.632.etc = function(X,y,B=2000){
   n=length(y)
  X=as.matrix(X)
  y=as.vector(y)
  Err.hat.cv1 = 0
   for ( i in 1:n ){
     beta.hat.jack = solve(t(X[-i,])%*%X[-i,])%*%t(X[-i,])%*%y[-i]
     Err.hat.cv1 = Err.hat.cv1+(y[i]-X[i,]%*%beta.hat.jack)^2
     }
   Err.hat.cv1 = Err.hat.cv1/n
   beta.hat = solve(t(X)%*%X)%*%t(X)%*%y
   err.bar = mean((y-X%*%beta.hat)^2)
   numerators = denominators = rep(0,n)
   for ( b in 1:B ){
     ind.star = index.star[b,]
     X.star = X[ind.star,]; y.star = y[ind.star]
     beta.hat.star = solve(t(X.star)%*%X.star)%*%t(X.star)%*%y.star
     for ( i in 1:n ) if (!(i%in%ind.star)){
       numerators[i] = numerators[i]+(y[i]-X[i,]%*%beta.hat.star)^2
       denominators[i] = denominators[i]+1 }
     }
   Err.hat.1 = mean(numerators/denominators)
   Err.hat.632 = 0.368*err.bar+0.632*Err.hat.1
   gamma.hat = mean(outer(1:n,1:n,function(i,j) (y[i]-X[j,]%*%beta.hat)^2))
   R.hat = (Err.hat.1-err.bar)/(gamma.hat-err.bar)
   if (R.hat<0||R.hat>1) R.hat = 0
   w.hat = 0.632/(1-0.368*R.hat)
   Err.hat.632plus = (1-w.hat)*err.bar+w.hat*Err.hat.1
   c(err.bar,Err.hat.632,Err.hat.632plus,Err.hat.cv1)
}
?combn
combn(16)
Err.hat.632.etc((X[,1]),y)
Err.hat.632.etc((X[,c(1,2)]),y)
Err.hat.632.etc((X[,c(1,3)]),y)
Err.hat.632.etc((X[,c(1,4)]),y)
Err.hat.632.etc((X[,c(1,5)]),y)
Err.hat.632.etc((X[,c(1,2,3)]),y)
Err.hat.632.etc((X[,c(1,2,4)]),y)
Err.hat.632.etc((X[,c(1,2,5)]),y)

Err.hat.632.etc((X[,c(1,3,4)]),y)
Err.hat.632.etc((X[,c(1,2,3,4)]),y)

for (i in 0:4){
  a=combn(2:5,i)
  for (j in 1:ncol(a)) {
  asd=as.vector(a[,j])
  cat(c(1,asd),':', Err.hat.632.etc(X[,c(1,asd)],y),'\n')
  }
}

head(data)
combn(letters[1:4], 3)
as.vector(combn(2:5,2)[,2])

oti=combn(2:5,1)
oti2=as.vector(combn[,4])

ncol(combn(2:5,1))
as.vector(combn[,j])
X
c(1,c(2,2))


##d
mean_ratio=function(data) mean(data[,1])/mean(data[,2])
mean_diff=function(x) mean(x[,1])-mean(x[,2])
petrol
y_a=petrol$Y[petrol$No=='A']
y_a
y_d=petrol$Y[petrol$No=='D']

######μχ/μy προσεγγιση##########
data=cbind(y_a,y_d)
var(y_d)
bootstrap_CIs(data,fxn=mean_ratio,fxn_dim = 2)
test_bootstrap_CIs(data,fxn=mean_ratio,fxn_dim = 2,mu0=3/2)
1.4231+c(-1,1)*qnorm(0.975)*0.0631
#### hist(bootstrap_CIs(data,fxn=mean_ratio,fxn_dim = 2)$r.star,breaks=seq(1,2,0.1)) ###
t = (mean_ratio(data)-mu0)/ sqrt(var(data)/n)
mean(y_a/y_d)
mean(y_a)/mean(y_d)

########χ=3/2y προσσέγιση
mu0=0


data=cbind(y_a,3/2*y_d)
var(data[,2])
####Welch test
welch_test=function(data){
  t = (mean(data[,1])-mean(data[,2]))/sqrt(var(data[,1])/length(data[,1])+var(data[,2])/length(data[,2]))
  nu1=((var(data[,1])/length(data[,1]))/(var(data[,1])/length(data[,1])+var(data[,2])/length(data[,2])))^2/(length(data[,2])-1)
  nu2 = ((var(data[,2])/length(data[,1]))/(var(data[,1])/length(data[,1])+var(data[,2])/length(data[,2])))^2/(length(data[,2])-1)
  df = 1/(nu1+nu2)
  return(2*pt(-abs(t),df))
}
welch_test(data)

mu_test=function(data,B=2000,alfa=0.05,mu0=0){
  data=as.matrix(data)
  x=data[,1];y=data[,2];
  n=length(x)
  m=length(y)
  t = (mean(x)-mean(y)-mu0)/(sqrt(var(x)/n +var(y)/m))
  xy.bar=mean(c(x,y))
  x.tilde = x-mean(x)+xy.bar
  y.tilde = y-mean(y)+xy.bar
  set.seed(1)
  x.star=replicate(B,x.tilde[sample(1:n,n,replace = TRUE)])
  set.seed(50900)
  y.star=replicate(B,y.tilde[sample(1:m,m,replace = TRUE)])
  
  #euresh tou bootstrap CI
  bootstrap_cint=NULL
  
  t.star = NULL
  for (b in 1:B){
    t.star[b] = (mean(x.star[,b])-mean(y.star[,b]))/(sqrt(var(x.star[,b])/n+ var(y.star[,b])/m))
    
  }
  #λόγω μικρού n κάποια bootstrap δείγματα μπορεί να περιέχουν 4 φορές την ίδια παρατήρηση επομένως να δίνουν διασπορά 0 και το t να γίνεται άπειρο#
  #τα αφαιρούμε από το δείγμα#
  #αναμένονται να είναι (1/4)^4*Β=0.004*Β τέτοια δείγματα
  #t.star=t.star[!is.na(t.star) & !is.infinite(t.star)]
  bootstrap_cint[1]=t-mu0-quantile(t.star,c(1-alfa/2,alfa/2),type=1,names=F)[1]*sd(t.star)
  bootstrap_cint[2]=t-mu0-quantile(t.star,c(1-alfa/2,alfa/2),type=1,names=F)[2]*sd(t.star)
  list(est=mean(t.star),p_value=2*min(mean(t.star<=t+mu0),mean(t.star>=t+mu0)),CIs=bootstrap_cint)
}
mu_test(data)

(1/4)^4
2000*0.004